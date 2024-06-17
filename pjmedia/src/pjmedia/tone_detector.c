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

#define PI 3.14159265358979323846
#define TONE_440HZ 0
#define TONE_480HZ 1

static const float energy_min_threshold=0.01f;

typedef struct goertzel_state{
        float coef;
} goertzel_state_t;

static void goertzel_state_init(goertzel_state_t *gs, int frequency, int sampling_frequency){
        gs->coef=(float)2*(float)cos(2*M_PI*((float)frequency/(float)sampling_frequency));
}

static float goertzel_state_run(goertzel_state_t *gs,int16_t  *samples, int nsamples, float total_energy){
        int i;
        float tmp;
        float q1=0;
        float q2=0;
        float freq_en;

        for(i=0;i<nsamples;++i){
                tmp=q1;
                q1=(gs->coef*q1) - q2 + (float)samples[i];
                q2=tmp;
        }

        freq_en= (q1*q1) + (q2*q2) - (q1*q2*gs->coef);
        /* return a relative frequency energy compared over the total signal energy */
        return freq_en/(total_energy*(float)nsamples*0.5f);
}

static float compute_energy(int16_t *samples, int nsamples){
        float en=0;
        int i;
        for(i=0;i<nsamples;++i){
                float s=(float)samples[i];
                en+=s*s;
        }
        return en;
}

struct tone_detector_port {
    pjmedia_port     base;
    pj_size_t	     cb_size;
    pj_bool_t	     subscribed;
    pj_bool_t	     cb_called;
    goertzel_state_t    state[2];
    void	     (*cb)(pjmedia_port*, void*);
};

static pj_status_t process_frame(pjmedia_port *this_port, pjmedia_frame *frame);
// static pj_status_t file_get_frame(pjmedia_port *this_port, pjmedia_frame *frame);
static pj_status_t tone_detector_on_destroy(pjmedia_port *this_port);

/*
 * Register callback.
 */
static pj_status_t pjmedia_tone_detector_port_set_cb(pjmedia_port *port,
				void *user_data,
			        void (*cb)(pjmedia_port *port, void *usr_data))
{
    return PJ_SUCCESS;
}

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
						     pjmedia_port **p_port,
						     pj_status_t (*cb)(pjmedia_port *port, void *usr_data),
						     void *cb_user_data)
{
    struct tone_detector_port *td_port;
    pjmedia_wave_hdr wave_hdr;
    pj_ssize_t size;
    pj_str_t name;
    pj_status_t status;

    /* Check arguments. */
    PJ_ASSERT_RETURN(pool && filename && p_port, PJ_EINVAL);

    /* Only supports 16bits per sample for now. */
    PJ_ASSERT_RETURN(bits_per_sample == 16, PJ_EINVAL);

    /* Create file port instance. */
    td_port = PJ_POOL_ZALLOC_T(pool, struct tone_detector_port);
    PJ_ASSERT_RETURN(td_port != NULL, PJ_ENOMEM);

    /* Initialize port info. */
    pj_strdup2(pool, &name, filename);
    pjmedia_port_info_init(&td_port->base.info, &name, SIGNATURE,
			   sampling_rate, channel_count, bits_per_sample,
			   samples_per_frame);

    goertzel_state_init(&td_port->state[0], 480.0f, sampling_rate);
    goertzel_state_init(&td_port->state[1], 440.0f, sampling_rate);
    
    td_port->base.put_frame = &process_frame;
    td_port->base.on_destroy = &tone_detector_on_destroy;

    *p_port = &td_port->base;

    td_port->base.port_data.pdata = cb_user_data;
    td_port->cb = cb;
    td_port->cb_called = PJ_FALSE;
    return PJ_SUCCESS;
}

static pj_status_t process_frame(pjmedia_port *this_port, pjmedia_frame *frame) 
{
	struct tone_detector_port *td_port = (struct tone_detector_port *)this_port;
	if (!this_port || !frame || td_port->cb_called) {
		return PJ_FALSE;
	}
	pj_int16_t *src = (pj_int16_t*)frame->buf;
	pj_bool_t t1 = PJ_FALSE;
	pj_bool_t t2 = PJ_FALSE;
	float en=compute_energy((int16_t*)frame->buf,frame->size/2);
	if (en>energy_min_threshold*(32767.0*32767.0*0.7)){
		goertzel_state_t *gs=&td_port->state[0];
		float freq_en1 = goertzel_state_run(gs,(int16_t*) src,frame->size/2,en);
		if (freq_en1 >= 0.4f) {
			t1 = PJ_TRUE;
		}
		gs=&td_port->state[1];
		float freq_en2 = goertzel_state_run(gs,(int16_t*) src,frame->size/2,en);
		if (freq_en2 >= 0.4f) {
			t2 = PJ_TRUE;
		}
		PJ_LOG(4,(THIS_FILE, "process_frame energy[%f] Hz1[%f] Hz2[%f]", en, freq_en1, freq_en2));
	} else {
		PJ_LOG(4,(THIS_FILE, "process_frame energy: %f", en));
	}
	if (t1 && t2 && td_port->cb) {
		(*td_port->cb)(&td_port->base, td_port->base.port_data.pdata);
		td_port->cb_called = PJ_TRUE;
	}
	return PJ_SUCCESS;
}

/*
 * Close the port, modify file header with updated file length.
 */
static pj_status_t tone_detector_on_destroy(pjmedia_port *this_port)
{
	PJ_LOG(4,(THIS_FILE, "tone detector on destroy"));
	return PJ_SUCCESS;
}

