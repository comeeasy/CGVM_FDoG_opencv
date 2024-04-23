#pragma once

#ifndef __CONV_CLR_H__
#define __CONV_CLR_H__

#include "mtoolkits.h"

#include <math.h>
#include "main.h"

void conv_RGBtoLAB(V3DF** buf, int width, int height);

void conv_LABtoRGB(V3DF** buf, int width, int height);

void RGBtoLAB(float iR, float iG, float iB, float *oL, float *oa, float *ob);
void LABtoRGB(float iL, float ia, float ib, float *oR, float *oG, float *oB);

#endif