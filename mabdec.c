#include "mabdec.h"
#include "matrix.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>

//
//  PRIVATE VARIABLES
//

static struct mabdec dec;

//
//  PRIVATE FUNCTIONS
//

static inline void mabdec_phaseDifference(void)
{
    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        float ph = (*dec.ph[i]);
        for (unsigned int j = i + 1; j < MAB_PORTS; j++) {
            dec.pha[i][j] = ph - (*dec.ph[j]);
        }
    }
}

static inline void mabdec_s_vector(struct matrix* restrict s)
{
    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        s->data[i][0] = dec.i[i];
    }

    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        for (unsigned int j = i + 1; j < MAB_PORTS; j++) {
            float res = dec.K[i][j] * dec.pha[i][j] * (1.0F - 2.0F * fabs(dec.pha[i][j]));
            s->data[i][0] -= (*dec.v[j]) * res;
            s->data[j][0] += (*dec.v[i]) * res;
        }
    }

    s->rows = MAB_PORTS;
    s->cols = 1;
}

static inline void mabdec_g_matrix(struct matrix* restrict w)
{
    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        w->data[i][i] = 0.0F;
    }

    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        for (unsigned int j = i + 1; j < MAB_PORTS; j++) {
            float res = dec.K[i][j] * (4.0F * fabs(dec.pha[i][j]) - 1.0F);
            w->data[i][j]  = (*dec.v[j]) * res;
            w->data[j][i]  = (*dec.v[i]) * res;
            w->data[i][i] -= w->data[i][j];
            w->data[j][j] -= w->data[j][i];
        }
    }

    w->rows = MAB_PORTS;
    w->cols = MAB_PORTS;
}

static inline void mabdec_addIncrements(struct matrix* restrict inc) 
{
    for (unsigned int i = 0; i < MAB_PORTS; i++) {
		*dec.ph[i] += inc->data[i][0];
	}
}

//
//  PUBLIC FUNCTIONS
//

void mabdec_init(void)
{
	memset(&dec, 0, sizeof(struct mabdec));
}

void mabdec_constUpdate(void)
{
    float Leq = 0.0F;
    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        float n = dec.n[i] * dec.n[i];
        Leq += (n / dec.Ls[i]);
        Leq += (n / dec.Lm[i]);
    }
    Leq = 1.0F / Leq;

    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        for (unsigned int j = i + 1; j < MAB_PORTS; j++) {
            float n = dec.n[i]  * dec.n[j];
            float L = dec.Ls[i] * dec.Ls[j];
            dec.K[i][j] = n / (dec.fs * (L / Leq));
        }
    }
}

void mabdec_phaseUpdate(void)
{
    static struct matrix tmp[3];

    mabdec_phaseDifference();
	mabdec_g_matrix(&tmp[0]);
	matrix_geninv  (&tmp[1], &tmp[0]);
	mabdec_s_vector(&tmp[0]);
	matrix_multiply(&tmp[2], &tmp[1], &tmp[0]);
    mabdec_addIncrements(&tmp[2]);
}

void mabdec_setFs(float freq)
{
    dec.fs = freq;
}

void mabdec_setLs(unsigned int id, float value)
{
    if (id < MAB_PORTS) {
        dec.Ls[id] = value;
    }
}

void mabdec_setLm(unsigned int id, float value)
{
    if (id < MAB_PORTS) {
        dec.Lm[id] = value;
    }
}

void mabdec_setRatio(unsigned int id, float value)
{
    if (id < MAB_PORTS) {
        dec.n[id] = value;
    }
}

void mabdec_setVoltagePtr(unsigned int id, float* ptr)
{
    if (id < MAB_PORTS) {
        dec.v[id] = ptr;
    }
}
void mabdec_setPhasePtr(unsigned int id, float* ptr)
{
    if (id < MAB_PORTS) {
	    dec.ph[id] = ptr;
    }
}

void mabdec_setRequest(unsigned int id, float req)
{
    if (id < MAB_PORTS) {
	    dec.i[id] = req;
    }
}
