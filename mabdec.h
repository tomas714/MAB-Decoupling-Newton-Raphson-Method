#ifndef MABDEC_H
#define MABDEC_H

#include "config.h"

struct mabdec 
{
    float  Ls[MAB_PORTS];
    float  Lm[MAB_PORTS];
    float   n[MAB_PORTS];
    float   K[MAB_PORTS][MAB_PORTS];
    float   i[MAB_PORTS];
    float*  v[MAB_PORTS];
    float* ph[MAB_PORTS];
    float  fs;
};

void mabdec_init(void);
void mabdec_constUpdate(void);
void mabdec_phaseUpdate(void);

void mabdec_setFs(float freq);
void mabdec_setLs(unsigned int id, float value);
void mabdec_setLm(unsigned int id, float value);
void mabdec_setRatio(unsigned int id, float value);
void mabdec_setVoltagePtr(unsigned int id, float* ptr);
void mabdec_setPhasePtr(unsigned int id, float* ptr);
void mabdec_setRequest(unsigned int id, float req);

#endif /* MABDEC_H */
