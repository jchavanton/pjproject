/* $Id$ */
/* 
 * Copyright (C) 2008-2011 Teluu Inc. (http://www.teluu.com)
 * Copyright (C) 2003-2008 Benny Prijono <benny@prijono.org>
 * Copyright (C) 2024      Julien Chavanton <jchavanton@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
 */
#include <pjmedia/wav_port.h>
#include <pjmedia/alaw_ulaw.h>
#include <pjmedia/errno.h>
#include <pjmedia/wave.h>
#include <pj/assert.h>
#include <pj/file_access.h>
#include <pj/file_io.h>
#include <pj/log.h>
#include <pj/pool.h>
#include <pj/string.h>


#define THIS_FILE	    "tone_detector.c"
#define SIGNATURE	    PJMEDIA_SIG_PORT_WAV_WRITER

// begin addition
#define PI 3.14159265358979323846
#define MAX_CORRELATION_SAMPLES 10000
#define SOUND_CARD_MIN_LATENCY 64
#define REF_TONE 0
#define READ_TONE 1
#define TONE_440HZ 0
#define TONE_480HZ 1

typedef struct correlation_variable {
	int val[MAX_CORRELATION_SAMPLES];
	float sum;
	float std;
	float avg;
} correlation_variable_t;

typedef struct correlation {
	correlation_variable_t var[2];
	int count;
	float sum_cross_deviation;
	float avg_sum_cross_deviation;
} correlation_t;

static float get_corr_coeff(correlation_t *c) {
	int count = c->count;
	int i=0;
	float corr_coef;
	if (count == 0)
		return 0.0f;
	c->var[REF_TONE].avg = c->var[REF_TONE].sum / count;
	c->var[READ_TONE].avg = c->var[READ_TONE].sum / count;
	c->sum_cross_deviation = 0;
	c->var[REF_TONE].std = 0;
	c->var[READ_TONE].std = 0;
	for (i=0;i<count;i++) {
		c->var[REF_TONE].std += pow(c->var[REF_TONE].val[i] - c->var[REF_TONE].avg, 2);
		c->var[READ_TONE].std += pow(c->var[READ_TONE].val[i] - c->var[READ_TONE].avg, 2);
		c->sum_cross_deviation += (c->var[REF_TONE].val[i]-c->var[REF_TONE].avg) * (c->var[READ_TONE].val[i]-c->var[READ_TONE].avg);
	}
	c->var[REF_TONE].std = sqrt(c->var[REF_TONE].std/count);
	c->var[READ_TONE].std = sqrt(c->var[READ_TONE].std/count);
	/*
    	PJ_LOG(4,(THIS_FILE, 
	      "File writer '%.*s' created: samp.rate=%d, bufsize=%uKB",
	      (int)fport->base.info.name.slen,
	      fport->base.info.name.ptr,
	      PJMEDIA_PIA_SRATE(&fport->base.info),
	      fport->bufsize / 1000));
	     */
	PJ_LOG(4,(THIS_FILE, "sum[%.3f|%.3f]avg[%.3f|%.3f]std[%.3f|%.3f] [%.2f/%d]",
              c->var[REF_TONE].sum, c->var[READ_TONE].sum, c->var[REF_TONE].avg, c->var[READ_TONE].avg, c->var[REF_TONE].std, c->var[READ_TONE].std, c->sum_cross_deviation, count);
	c->avg_sum_cross_deviation = c->sum_cross_deviation / count);

	if (c->var[REF_TONE].std == 0 || c->var[READ_TONE].std == 0)
		return 0.0f;
	corr_coef = c->avg_sum_cross_deviation / (c->var[REF_TONE].std * c->var[READ_TONE].std);
	return corr_coef;
}

typedef struct goertzel_state {
	float Q2;
	float Q1;
	float sine;
	float cosine;
	float coeff;
	float sampling_rate;
	int block_size;
	int block_processed_samples;
	float target_hz;
	float last_block_real_part;
	float last_block_imag_part;
	int block_processed;
} goertzel_state_t;

/* Call this routine before every "block" (size=N) of samples. */
static void goertzel_reset(goertzel_state_t* g) {
	g->Q2 = 0;
	g->Q1 = 0;
	g->block_processed_samples = 0;
}

/* Call this once, to precompute the constants. */
static void goertzel_init(goertzel_state_t* g, float target_hz, int block_size, float sampling_rate) {
	int k;
	float float_n;
	float omega;
	g->sampling_rate = sampling_rate;
	g->block_size = block_size;
	g->target_hz = target_hz;
	float_n = (float) block_size;
	k = (int) (0.5 + ((float_n * g->target_hz) / g->sampling_rate));
	omega = (2.0 * PI * k) / float_n;
	g->sine = sin(omega);
	g->cosine = cos(omega);
	g->coeff = 2.0 * g->cosine;
	g->block_processed=0;
	goertzel_reset(g);
}

/* Call this routine after every block to get the complex result. */
static void goertzel_get_complex(goertzel_state_t* g) {
	g->last_block_real_part = (g->Q1 - g->Q2 * g->cosine);
	g->last_block_imag_part = (g->Q2 * g->sine);
}

static void update_corr_reference(goertzel_state_t* g, correlation_t* c, int delay_offset) {
	int block_processed = 0;
	int sum=0;
	int i=0;
	for (i = 0 ; i < g->block_processed; i++) {
		int ms, ms_delayed;
		int reference = 1000;
		block_processed++;
		ms = block_processed*g->block_size/(int)(g->sampling_rate/1000);
		ms_delayed = ms - delay_offset;

		if (ms_delayed%5000 > 2000) {
			reference = 0;
		}
		sum += reference;
		c[TONE_440HZ].var[REF_TONE].val[i] = reference;
		c[TONE_480HZ].var[REF_TONE].val[i] = reference;
	}
	c[TONE_440HZ].var[REF_TONE].sum = sum;
	c[TONE_480HZ].var[REF_TONE].sum = sum;
}

static int goertzel_block_check(goertzel_state_t* g, correlation_t* c) {
	if (g->block_processed_samples == g->block_size) {
		float magnitude_squared, magnitude, db;

		goertzel_get_complex(g);
		goertzel_reset(g);
		g->block_processed++;
		magnitude_squared = g->last_block_real_part*g->last_block_real_part + g->last_block_imag_part*g->last_block_imag_part;
		magnitude = sqrt(magnitude_squared);

		db = 10*log10f(magnitude/32768);
		PJ_LOG(4,(THIS_FILE,"goertzel_block_check[%d] [%fHz] mag[%0.0f] corr:dB[%0.2f]", g->block_processed, g->target_hz, magnitude, db));
		if (g->block_processed <= MAX_CORRELATION_SAMPLES) {
			c->var[READ_TONE].val[c->count] = (int)(db*100);
			c->var[READ_TONE].sum += (int)(db*100);
			c->count++;
			return 1;
		}
	} else {
		g->block_processed_samples++;
	}
	return 0;
}

/* Call this routine for every sample. */
static void goertzel_process_sample(goertzel_state_t* g, uint16_t sample) {
	float Q0;
	Q0 = g->coeff * g->Q1 - g->Q2 + (float) sample;
	g->Q2 = g->Q1;
	g->Q1 = Q0;
}

typedef struct detector_state {
	float target_fz[2]; // the two target Fz will that we will search in the signal
	correlation_t c[2];  // the correlation data between each frequency and the reference signal
	goertzel_state_t g[2]; // the goertzel analysis for each frequency, currently 440Hz and 480Hz
	// MSBufferizer *buf;
	int rate;
	int framesize;
	int frame_ms;
#if SIG_DUMP == 1
	FILE *recfile;
#endif
} detector_state_t;

// static void detector_init(MSFilter *f) {
// static void detector_init(pj_pool_t *pool, pjmedia_port *p) {
static void detector_init(detector_state_t *s) {
	int goertzel_block_size = 256; // 32 ms //  8000/1000 8*40 = 320
    	// detector_state_t *s=PJ_POOL_ZALLOC_T(pool, detector_state_t);
	// p->s.buf=ms_bufferizer_new();
	s->rate=8000;
	s->target_fz[0] = 480.0f; // US 480Hz
	s->target_fz[1] = 440.0f; // US 440Hz
	memset(&s->c[TONE_440HZ],0,sizeof(correlation_t));
	memset(&s->c[TONE_480HZ],0,sizeof(correlation_t));
	goertzel_init(&s->g[0], s->target_fz[0], goertzel_block_size, s->rate);
	goertzel_init(&s->g[1], s->target_fz[1], goertzel_block_size, s->rate);
	s->frame_ms=20;
	s->framesize=2*(s->frame_ms*s->rate)/1000;
	//f->data=s;
//#if SIG_DUMP == 1
//{
//	char *fname = ms_strdup_printf("%s/ringbacktone-detector-%p.raw", SIG_DUMP_PREFIX, f);
//	s->recfile = fopen(fname, "w");
//	if(s->recfile) {
//		ms_message("tone detector recording[%s][%p]", fname, s->recfile);
//	} else {
//		ms_message("tone detector error opening file recording[%s][%p]", fname, s->recfile);
//	}
//	ms_free(fname);
//}
//#endif
}

static void detector_uninit(pjmedia_port *p) {
	// detector_state_t *s=(detector_state_t *)f->data;
#if SIG_DUMP == 1
	if(s->recfile)
		fclose(s->recfile);
//	ms_message("tone detector recording file closed");
#endif
	// ms_bufferizer_destroy (s->buf);
	// ms_free(f->data);
}

// static int detector_set_rate(MSFilter *f, void *arg){
// 	detector_state_t *s=(detector_state_t *)f->data;
// 	s->rate = *((int*) arg);
// 	return 0;
// }

struct file_port {
    pjmedia_port     base;
    pjmedia_wave_fmt_tag fmt_tag;
    pj_uint16_t	     bytes_per_sample;

    pj_size_t	     bufsize;
    char	    *buf;
    char	    *writepos;
    pj_size_t	     total;

    pj_oshandle_t    fd;

    pj_size_t	     cb_size;
    pj_status_t	   (*cb)(pjmedia_port*, void*);
    pj_bool_t	     subscribed;
    pj_bool_t	     cb_called;
    void	   (*cb2)(pjmedia_port*, void*);
    detector_state_t  state;
};

static pj_status_t file_put_frame(pjmedia_port *this_port, pjmedia_frame *frame);
static pj_status_t file_get_frame(pjmedia_port *this_port, pjmedia_frame *frame);
static pj_status_t file_on_destroy(pjmedia_port *this_port);


/*
 * Create file writer port.
 */
PJ_DEF(pj_status_t) pjmedia_tone_detector_port_create( pj_pool_t *pool,
						     const char *filename,
						     unsigned sampling_rate,
						     unsigned channel_count,
						     unsigned samples_per_frame,
						     unsigned bits_per_sample,
						     unsigned flags,
						     pj_ssize_t buff_size,
						     pjmedia_port **p_port )
{
    struct file_port *fport;
    pjmedia_wave_hdr wave_hdr;
    pj_ssize_t size;
    pj_str_t name;
    pj_status_t status;

    /* Check arguments. */
    PJ_ASSERT_RETURN(pool && filename && p_port, PJ_EINVAL);

    /* Only supports 16bits per sample for now.
     * See flush_buffer().
     */
    PJ_ASSERT_RETURN(bits_per_sample == 16, PJ_EINVAL);


    /* Create file port instance. */
    fport = PJ_POOL_ZALLOC_T(pool, struct file_port);
    PJ_ASSERT_RETURN(fport != NULL, PJ_ENOMEM);

    /* Initialize port info. */
    pj_strdup2(pool, &name, filename);
    pjmedia_port_info_init(&fport->base.info, &name, SIGNATURE,
			   sampling_rate, channel_count, bits_per_sample,
			   samples_per_frame);

    detector_init(&fport->state);
//		goertzel_process_sample(&fport->state.g[0], *src);

    // init 
    fport->base.get_frame = &file_get_frame;
    fport->base.put_frame = &file_put_frame;
    fport->base.on_destroy = &file_on_destroy;

    if (flags == PJMEDIA_FILE_WRITE_ALAW) {
	fport->fmt_tag = PJMEDIA_WAVE_FMT_TAG_ALAW;
	fport->bytes_per_sample = 1;
    } else if (flags == PJMEDIA_FILE_WRITE_ULAW) {
	fport->fmt_tag = PJMEDIA_WAVE_FMT_TAG_ULAW;
	fport->bytes_per_sample = 1;
    } else {
	fport->fmt_tag = PJMEDIA_WAVE_FMT_TAG_PCM;
	fport->bytes_per_sample = 2;
    }

    /* Open file in write and read mode.
     * We need the read mode because we'll modify the WAVE header once
     * the recording has completed.
     */
    status = pj_file_open(pool, filename, PJ_O_WRONLY, &fport->fd);
    if (status != PJ_SUCCESS)
	return status;

    /* Initialize WAVE header */
    pj_bzero(&wave_hdr, sizeof(pjmedia_wave_hdr));
    wave_hdr.riff_hdr.riff = PJMEDIA_RIFF_TAG;
    wave_hdr.riff_hdr.file_len = 0; /* will be filled later */
    wave_hdr.riff_hdr.wave = PJMEDIA_WAVE_TAG;

    wave_hdr.fmt_hdr.fmt = PJMEDIA_FMT_TAG;
    wave_hdr.fmt_hdr.len = 16;
    wave_hdr.fmt_hdr.fmt_tag = (pj_uint16_t)fport->fmt_tag;
    wave_hdr.fmt_hdr.nchan = (pj_int16_t)channel_count;
    wave_hdr.fmt_hdr.sample_rate = sampling_rate;
    wave_hdr.fmt_hdr.bytes_per_sec = sampling_rate * channel_count * 
				     fport->bytes_per_sample;
    wave_hdr.fmt_hdr.block_align = (pj_uint16_t)
				   (fport->bytes_per_sample * channel_count);
    wave_hdr.fmt_hdr.bits_per_sample = (pj_uint16_t)
				       (fport->bytes_per_sample * 8);

    wave_hdr.data_hdr.data = PJMEDIA_DATA_TAG;
    wave_hdr.data_hdr.len = 0;	    /* will be filled later */


    /* Convert WAVE header from host byte order to little endian
     * before writing the header.
     */
    pjmedia_wave_hdr_host_to_file(&wave_hdr);


    /* Write WAVE header */
    if (fport->fmt_tag != PJMEDIA_WAVE_FMT_TAG_PCM) {
	pjmedia_wave_subchunk fact_chunk;
	pj_uint32_t tmp = 0;

	fact_chunk.id = PJMEDIA_FACT_TAG;
	fact_chunk.len = 4;

	PJMEDIA_WAVE_NORMALIZE_SUBCHUNK(&fact_chunk);

	/* Write WAVE header without DATA chunk header */
	size = sizeof(pjmedia_wave_hdr) - sizeof(wave_hdr.data_hdr);
	status = pj_file_write(fport->fd, &wave_hdr, &size);
	if (status != PJ_SUCCESS) {
	    pj_file_close(fport->fd);
	    return status;
	}

	/* Write FACT chunk if it stores compressed data */
	size = sizeof(fact_chunk);
	status = pj_file_write(fport->fd, &fact_chunk, &size);
	if (status != PJ_SUCCESS) {
	    pj_file_close(fport->fd);
	    return status;
	}
	size = 4;
	status = pj_file_write(fport->fd, &tmp, &size);
	if (status != PJ_SUCCESS) {
	    pj_file_close(fport->fd);
	    return status;
	}

	/* Write DATA chunk header */
	size = sizeof(wave_hdr.data_hdr);
	status = pj_file_write(fport->fd, &wave_hdr.data_hdr, &size);
	if (status != PJ_SUCCESS) {
	    pj_file_close(fport->fd);
	    return status;
	}
    } else {
	size = sizeof(pjmedia_wave_hdr);
	status = pj_file_write(fport->fd, &wave_hdr, &size);
	if (status != PJ_SUCCESS) {
	    pj_file_close(fport->fd);
	    return status;
	}
    }

    /* Set buffer size. */
    if (buff_size < 1) buff_size = PJMEDIA_FILE_PORT_BUFSIZE;
    fport->bufsize = buff_size;

    /* Check that buffer size is greater than bytes per frame */
    pj_assert(fport->bufsize >= PJMEDIA_PIA_AVG_FSZ(&fport->base.info));


    /* Allocate buffer and set initial write position */
    fport->buf = (char*) pj_pool_alloc(pool, fport->bufsize);
    if (fport->buf == NULL) {
	pj_file_close(fport->fd);
	return PJ_ENOMEM;
    }
    fport->writepos = fport->buf;

    /* Done. */
    *p_port = &fport->base;

    PJ_LOG(4,(THIS_FILE, 
	      "File writer '%.*s' created: samp.rate=%d, bufsize=%uKB",
	      (int)fport->base.info.name.slen,
	      fport->base.info.name.ptr,
	      PJMEDIA_PIA_SRATE(&fport->base.info),
	      fport->bufsize / 1000));

    return PJ_SUCCESS;
}


/*
 * Get current writing position. 
 */
PJ_DEF(pj_ssize_t) pjmedia_tone_detector_port_get_pos( pjmedia_port *port )
{
    struct file_port *fport;

    /* Sanity check */
    PJ_ASSERT_RETURN(port, -PJ_EINVAL);

    /* Check that this is really a writer port */
    PJ_ASSERT_RETURN(port->info.signature == SIGNATURE, -PJ_EINVALIDOP);

    fport = (struct file_port*) port;

    return fport->total;
}


#if !DEPRECATED_FOR_TICKET_2251
/*
 * Register callback.
 */
PJ_DEF(pj_status_t) pjmedia_tone_detector_port_set_cb( pjmedia_port *port,
				pj_size_t pos,
				void *user_data,
			        pj_status_t (*cb)(pjmedia_port *port,
						  void *usr_data))
{
    struct file_port *fport;

    /* Sanity check */
    PJ_ASSERT_RETURN(port && cb, PJ_EINVAL);

    /* Check that this is really a writer port */
    PJ_ASSERT_RETURN(port->info.signature == SIGNATURE, PJ_EINVALIDOP);

    PJ_LOG(1, (THIS_FILE, "pjmedia_tone_detector_port_set_cb() is deprecated. "
    	       "Use pjmedia_tone_detector_port_set_cb2() instead."));

    fport = (struct file_port*) port;

    fport->cb_size = pos;
    fport->base.port_data.pdata = user_data;
    fport->cb = cb;

    return PJ_SUCCESS;
}
#endif


/*
 * Register callback.
 */
PJ_DEF(pj_status_t) pjmedia_tone_detector_port_set_cb2(pjmedia_port *port,
				pj_size_t pos,
				void *user_data,
			        void (*cb)(pjmedia_port *port,
					   void *usr_data))
{
    struct file_port *fport;

    /* Sanity check */
    PJ_ASSERT_RETURN(port && cb, PJ_EINVAL);

    /* Check that this is really a writer port */
    PJ_ASSERT_RETURN(port->info.signature == SIGNATURE, PJ_EINVALIDOP);

    fport = (struct file_port*) port;

    fport->cb_size = pos;
    fport->base.port_data.pdata = user_data;
    fport->cb2 = cb;
    fport->cb_called = PJ_FALSE;

    return PJ_SUCCESS;
}


#if defined(PJ_IS_BIG_ENDIAN) && PJ_IS_BIG_ENDIAN!=0
    static void swap_samples(pj_int16_t *samples, unsigned count)
    {
	unsigned i;
	for (i=0; i<count; ++i) {
	    samples[i] = pj_swap16(samples[i]);
	}
    }
#else
#   define swap_samples(samples,count)
#endif

/*
 * Flush the contents of the buffer to the file.
 */
static pj_status_t flush_buffer(struct file_port *fport)
{
    pj_ssize_t bytes = fport->writepos - fport->buf;
    pj_status_t status;

    /* Convert samples to little endian */
    swap_samples((pj_int16_t*)fport->buf, bytes/fport->bytes_per_sample);

    /* Write to file. */
    status = pj_file_write(fport->fd, fport->buf, &bytes);

    /* Reset writepos */
    fport->writepos = fport->buf;

    return status;
}

static pj_status_t file_on_event(pjmedia_event *event,
                                 void *user_data)
{
    struct file_port *fport = (struct file_port*)user_data;

    if (event->type == PJMEDIA_EVENT_CALLBACK) {
	if (fport->cb2)
	    (*fport->cb2)(&fport->base, fport->base.port_data.pdata);
    }
    
    return PJ_SUCCESS;
}

#define MIN(a,b) ((a)>(b) ? (b) : (a))
static int get_result(detector_state_t *s){
//	detector_state_t *s=(detector_state_t *)f->data;
	float best_corr_coef=-1.0f;
	float corr_coef_1=0.0f;
	float corr_coef_2=0.0f;
	int delay_adj = SOUND_CARD_MIN_LATENCY; // we start with 64ms, this is considered the lowest delay echo delay (between speaker and mic)
	int delay_inc = s->g->block_size / (s->g->sampling_rate/1000); // delay increment in ms is equivalent to the goertzel block size
//	rbt_rbt_echo_detection_data_t *data = (rbt_rbt_echo_detection_data_t *) arg;
	
	// int echo_delay_ms;
	for(delay_adj=SOUND_CARD_MIN_LATENCY; delay_adj< 500; delay_adj += delay_inc) {
		update_corr_reference(s->g, s->c, delay_adj);
		corr_coef_1 = get_corr_coeff(&s->c[TONE_440HZ]);
		corr_coef_2 = get_corr_coeff(&s->c[TONE_480HZ]);
		PJ_LOG(4,(THIS_FILE,"get_corr_coef %dms[+%dms] best[%.3f] [%.3f|%.3f]", delay_adj, delay_inc, best_corr_coef, corr_coef_1, corr_coef_2));
		if ( best_corr_coef < MIN(corr_coef_1, corr_coef_2)) {
			best_corr_coef = MIN(corr_coef_1,corr_coef_2);
	//		echo_delay_ms =delay_adj;
		} 
	}

		PJ_LOG(4,(THIS_FILE, "get_corr_coef result[%.3f] echo_delay[%dms]", best_corr_coef, delay_adj));
//	data->duration_ms = s->g->block_size / (s->g->sampling_rate/1000) * s->g->block_processed;
//	data->echo_delay_ms = echo_delay_ms;
//	data->corr_coef = best_corr_coef; // We return the worst one to minize false positif echo detection
	return 0;
}

/*
 * Put a frame into the buffer. When the buffer is full, flush the buffer
 * to the file.
 */
static pj_status_t file_put_frame(pjmedia_port *this_port, 
				  pjmedia_frame *frame)
{
    struct file_port *fport = (struct file_port *)this_port;
    pj_size_t frame_size;

    if (fport->fmt_tag == PJMEDIA_WAVE_FMT_TAG_PCM)
	frame_size = frame->size;
    else
	frame_size = frame->size >> 1;

    /* Flush buffer if we don't have enough room for the frame. */
    if (fport->writepos + frame_size > fport->buf + fport->bufsize) {
	pj_status_t status;
	status = flush_buffer(fport);
	if (status != PJ_SUCCESS)
	    return status;
    }

    /* Check if frame is not too large. */
    PJ_ASSERT_RETURN(fport->writepos+frame_size <= fport->buf+fport->bufsize,
		     PJMEDIA_EFRMFILETOOBIG);

    /* Copy frame to buffer. */
    if (fport->fmt_tag == PJMEDIA_WAVE_FMT_TAG_PCM) {
	unsigned i;
	pj_memcpy(fport->writepos, frame->buf, frame->size);
	pj_int16_t *src = (pj_int16_t*)frame->buf;
        for (i = 0; i < frame_size/2; ++i) {
        	goertzel_process_sample(&fport->state.g[0], *src);
		goertzel_process_sample(&fport->state.g[1], *src);
		goertzel_block_check(&fport->state.g[0], &fport->state.c[TONE_440HZ]);
		goertzel_block_check(&fport->state.g[1], &fport->state.c[TONE_480HZ]);
		*src++;
	    }
    } else {
	unsigned i;
	pj_int16_t *src = (pj_int16_t*)frame->buf;
	pj_uint8_t *dst = (pj_uint8_t*)fport->writepos;

	if (fport->fmt_tag == PJMEDIA_WAVE_FMT_TAG_ULAW) {
	    for (i = 0; i < frame_size; ++i) {
		*dst++ = pjmedia_linear2ulaw(*src++);
	    }
	} else {
	    for (i = 0; i < frame_size; ++i) {
		*dst++ = pjmedia_linear2alaw(*src++);
	    }
	}

    }
    fport->writepos += frame_size;

    /* Increment total written, and check if we need to call callback */
    fport->total += frame_size;
    if (fport->total >= fport->cb_size) {
	if (fport->cb2) {
	    if (!fport->subscribed) {
	        pj_status_t status;

	    	status = pjmedia_event_subscribe(NULL, &file_on_event,
	    				         fport, fport);
	    	fport->subscribed = (status == PJ_SUCCESS)? PJ_TRUE:
	    			    PJ_FALSE;
	    }

	    if (fport->subscribed && !fport->cb_called) {
	    	pjmedia_event event;

	    	/* To prevent the callback from being called more than once. */
	    	fport->cb_called = PJ_TRUE;

	    	pjmedia_event_init(&event, PJMEDIA_EVENT_CALLBACK,
	                      	   NULL, fport);
	    	pjmedia_event_publish(NULL, fport, &event,
	                              PJMEDIA_EVENT_PUBLISH_POST_EVENT);
	    }
	} else if (fport->cb) {
	    pj_status_t (*cb)(pjmedia_port*, void*);
	    pj_status_t status;

	    cb = fport->cb;
	    fport->cb = NULL;

	    status = (*cb)(this_port, this_port->port_data.pdata);
	    return status;
	}
    }

    return PJ_SUCCESS;
}

/*
 * Get frame, basicy is a no-op operation.
 */
static pj_status_t file_get_frame(pjmedia_port *this_port, 
				  pjmedia_frame *frame)
{
    PJ_UNUSED_ARG(this_port);
    PJ_UNUSED_ARG(frame);
    return PJ_EINVALIDOP;
}



/*
 * Close the port, modify file header with updated file length.
 */
static pj_status_t file_on_destroy(pjmedia_port *this_port)
{
    enum { FILE_LEN_POS = 4, DATA_LEN_POS = 40 };
    struct file_port *fport = (struct file_port *)this_port;
    pj_off_t file_size;
    pj_ssize_t bytes;
    pj_uint32_t wave_file_len;
    pj_uint32_t wave_data_len;
    pj_status_t status;
    pj_uint32_t data_len_pos = DATA_LEN_POS;

    get_result(&fport->state);

    if (fport->subscribed) {
    	pjmedia_event_unsubscribe(NULL, &file_on_event, fport, fport);
    	fport->subscribed = PJ_FALSE;
    }

    /* Flush remaining buffers. */
    if (fport->writepos != fport->buf) 
	flush_buffer(fport);

    /* Get file size. */
    status = pj_file_getpos(fport->fd, &file_size);
    if (status != PJ_SUCCESS) {
        pj_file_close(fport->fd);
	return status;
    }

    /* Calculate wave fields */
    wave_file_len = (pj_uint32_t)(file_size - 8);
    wave_data_len = (pj_uint32_t)(file_size - sizeof(pjmedia_wave_hdr));

#if defined(PJ_IS_BIG_ENDIAN) && PJ_IS_BIG_ENDIAN!=0
    wave_file_len = pj_swap32(wave_file_len);
    wave_data_len = pj_swap32(wave_data_len);
#endif

    /* Seek to the file_len field. */
    status = pj_file_setpos(fport->fd, FILE_LEN_POS, PJ_SEEK_SET);
    if (status != PJ_SUCCESS) {
        pj_file_close(fport->fd);
	return status;
    }

    /* Write file_len */
    bytes = sizeof(wave_file_len);
    status = pj_file_write(fport->fd, &wave_file_len, &bytes);
    if (status != PJ_SUCCESS) {
        pj_file_close(fport->fd);
	return status;
    }

    /* Write samples_len in FACT chunk */
    if (fport->fmt_tag != PJMEDIA_WAVE_FMT_TAG_PCM) {
	enum { SAMPLES_LEN_POS = 44};
	pj_uint32_t wav_samples_len;

	/* Adjust wave_data_len & data_len_pos since there is FACT chunk */
	wave_data_len -= 12;
	data_len_pos += 12;
	wav_samples_len = wave_data_len;

	/* Seek to samples_len field. */
	status = pj_file_setpos(fport->fd, SAMPLES_LEN_POS, PJ_SEEK_SET);
        if (status != PJ_SUCCESS) {
            pj_file_close(fport->fd);
	    return status;
        }

	/* Write samples_len */
	bytes = sizeof(wav_samples_len);
	status = pj_file_write(fport->fd, &wav_samples_len, &bytes);
	if (status != PJ_SUCCESS) {
            pj_file_close(fport->fd);
	    return status;
        }
    }

    /* Seek to data_len field. */
    status = pj_file_setpos(fport->fd, data_len_pos, PJ_SEEK_SET);
    if (status != PJ_SUCCESS) {
        pj_file_close(fport->fd);
	return status;
    }

    /* Write file_len */
    bytes = sizeof(wave_data_len);
    status = pj_file_write(fport->fd, &wave_data_len, &bytes);
    if (status != PJ_SUCCESS) {
        pj_file_close(fport->fd);
	return status;
    }

    /* Close file */
    status = pj_file_close(fport->fd);
    if (status != PJ_SUCCESS)
	return status;

    /* Done. */
    return PJ_SUCCESS;
}

