#pragma once


#ifndef __CONV_CLR_H__
#define __CONV_CLR_H__

#include "mtoolkits.h"

#include <math.h>

#define R_			0		//color
#define G_			1
#define B_			2
#define X			0
#define Y			1
#define Z			2


//변환행렬
//M3DF m_RGBtoXYZ = {0.490f, 0.310f, 0.200f, 0.177f, 0.812f, 0.011f, 0.000f, 0.010f, 0.990f};		//RcGcBc
// Matrix3DF m_RGBtoXYZ = {0.607f, 0.174f, 0.200f, 0.299f, 0.587f, 0.114f, 0.000f, 0.066f, 1.116f};		//RnGnBn
//M3DF m_RGBtoXYZ = {0.393f, 0.365f, 0.192f, 0.212f, 0.701f, 0.087f, 0.019f, 0.112f, 0.958f};		//RsGsBs

//M3DF m_XYZtoRGB = {2.36486f, -0.897081f, -0.467783f, -0.515564f, 1.42727f, 0.0882959f, 0.00520772f, -0.0144169f, 1.00921f};		//RcGcBc
// Matrix3DF m_XYZtoRGB = {1.91046f, -0.53394f, -0.287834f, -0.984436f, 1.9985f, -0.027726f, 0.0582193f, -0.118191f, 0.897697f};		//RnGnBn
//M3DF m_XYZtoRGB = {3.50969f, -1.74031f, -0.545358f, -1.06828f, 1.97725f, 0.0345393f, 0.0552852f, -0.196645f, 1.05062f};		//RsGsBs

//reference white
//V3DF rWhite = {1.0f, 1.0f, 1.0f};
//V3DF rWhite = {0.46f, 0.46f, 0.46f};
// V3DF rWhite = {0.964221f, 1.0f, 0.825211f};


void conv_RGBtoLAB(V3DF** buf, int width, int height);

void conv_LABtoRGB(V3DF** buf, int width, int height);

void RGBtoLAB(float iR, float iG, float iB, float *oL, float *oa, float *ob);
void LABtoRGB(float iL, float ia, float ib, float *oR, float *oG, float *oB);
void RGBtoXYZ(float iR, float iG, float iB, float *oX, float *oY, float *oZ);
void XYZtoRGB(float iX, float iY, float iZ, float *oR, float *oG, float *oB);
void XYZtoLUV(float iX, float iY, float iZ, float *oL, float *oU, float *oV);
void LUVtoXYZ(float iL, float iU, float iV, float *oX, float *oY, float *oZ);
void LUVtoRGB(float iL, float iU, float iV, float *oR, float *oG, float *oB);
void RGBtoLUV(float iR, float iG, float iB, float *oL, float *oU, float *oV);

void conv_RGBtoLUV(V3DF** buf, int width, int height);		//RGB -> LUV
void conv_LUVtoRGB(V3DF** buf, int width, int height);		//LUV -> RGB
void conv_RGBtoXYZ(V3DF** buf, int width, int height);		//RGB -> XYZ
void conv_XYZtoRGB(V3DF** buf, int width, int height);		//XYZ -> RGB
void conv_XYZtoLUV(V3DF** buf, int width, int height);		//XYZ -> LUV
void conv_LUVtoXYZ(V3DF** buf, int width, int height);		//LUV -> XYZ

#endif