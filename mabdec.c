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

static void mabdec_s_matrix(struct matrix* s)
{
    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        float result = 0.0F;
        for (unsigned int j = 0; j < MAB_PORTS; j++) {
            if (i != j) {
                float pha = (*dec.ph[i]) - (*dec.ph[j]);
                result += (*dec.v[j]) * dec.K[i][j] * pha * (1.0F - 2.0F * fabs(pha));
            }
        }
        s->data[i][0] = -1.0F * (result - dec.i[i]);
    }

    s->rows = MAB_PORTS;
    s->cols = 1;
}

static void mabdec_w_matrix(struct matrix* w)
{
    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        for (unsigned int j = 0; j < MAB_PORTS; j++) {
            float result = 0.0F;
            if (i == j) {
                for (unsigned int k = 0; k < MAB_PORTS; k++) {
                    if (i != k) {
                        result += (*dec.v[k]) * dec.K[i][k] * (1.0F - 4.0F * fabs((*dec.ph[i]) - (*dec.ph[k])));
                    }
                }
            }
            else {
                result = (*dec.v[j]) * dec.K[i][j] * (4.0F * fabs((*dec.ph[i]) - (*dec.ph[j])) - 1.0F);
            }
            w->data[i][j] = result;
        }
    }
    w->rows = MAB_PORTS;
    w->cols = MAB_PORTS;
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
        float n2 = dec.n[i] * dec.n[i];
        Leq += (n2 / dec.Ls[i]) + (n2 / dec.Lm[i]);
    }
    Leq = 1.0F / Leq;

    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        for (unsigned int j = 0; j < MAB_PORTS; j++) {
            if (i != j) {
                dec.K[i][j] = (1.0F / dec.fs) * ((dec.n[i] * dec.n[j]) / ((dec.Ls[i] * dec.Ls[j]) / Leq));
            }
            else {
                dec.K[i][j] = 0.0F;
            }
        }
    }
}

void mabdec_phaseUpdate(void)
{
    static struct matrix tmp[3];

	mabdec_w_matrix(&tmp[0]);
	matrix_geninv  (&tmp[1], &tmp[0]);
	mabdec_s_matrix(&tmp[0]);
	matrix_multiple(&tmp[2], &tmp[1], &tmp[0]);

	for (unsigned int i = 0; i < MAB_PORTS; i++) {
		*dec.ph[i] += tmp[2].data[i][0];
	}
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
