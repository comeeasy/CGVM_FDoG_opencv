#pragma once

#ifndef __COMMON_FUNC_H__
#define __COMMON_FUNC_H__

#include "mtoolkits.h"


void init_buf( V3DF **OutImg, V3DF **InImg, int width, int height);
void init_buf( V3DF **OutImg, cv::Mat InImg, int width, int height);

void GetValAtPoint_V3DF(V3DF **ival, int w, int h, float px, float py, V3DF oval);
void GetValAtPoint(float **ival, int w, int h, float px, float py, float *oval);

void BuildGaussianMask(int n, float s, float *m);
void GaussianFilter(float **InImg, float **OutImg, int width, int height);

void ScalarBuf(float *OutBuf, float scalar, int length);
void DifferenceBuf(float *OutBuf, float *InBuf1, float *InBuf2, int length);

#endif