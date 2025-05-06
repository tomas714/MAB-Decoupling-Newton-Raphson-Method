#include <stdio.h>
#include <stdlib.h>

#include "mabdec.h"

int main(void)
{
    float voltage[MAB_PORTS] = {24.0F, 24.0F, 24.0F};
    float   phase[MAB_PORTS] = { 0.0F,  0.0F,  0.0F};

    //
    //  SETTINGS
    //

    mabdec_init();
    mabdec_setFs(100e3);
    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        mabdec_setVoltagePtr(i, (float*)&voltage[i]);
        mabdec_setPhasePtr(i, (float*)&phase[i]);
        mabdec_setRatio(i, 6.0F/3.0F);
        mabdec_setLm(i, 600e-6F);
        mabdec_setLs(i, 1e-6F);
    }
    mabdec_constUpdate();

    //
    //  UPDATE PHASE
    //

    mabdec_setRequest(0, 10.0F);
    mabdec_setRequest(1, -5.0F);
    mabdec_setRequest(2, -5.0F);
    
    for (unsigned int i = 0; i < 3; i++) {
        mabdec_phaseUpdate();
    }

    //
    //  PRINTING RESULT
    //

    for (unsigned int i = 0; i < MAB_PORTS; i++) {
        printf("ph%u = %+.5f\n", i+1, phase[i]);
    }

    return 0;
}
