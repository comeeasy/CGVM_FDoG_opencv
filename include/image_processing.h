#pragma once

#ifndef __IMAGE_PROCESSING_H__
#define __IMAGE_PROCESSING_H__

#include <time.h>

#include "mtoolkits.h"
#include "common_func.h"
#include "conv_clr.h"

#include "FlowPath.h"

#define NEGATIVE 0
#define POSITIVE 1




// Flow-Based Bilateral Filter
void apply_FBL_filter( 
    V3DF** src_cim, 
	FlowPath** src_fpath,
    V3DF** &cimFBL,
    int w, int h,
    float sigma_e, float gamma_e, float sigma_g, float gamma_g,	int iteration,
    int threshold_T
);

void get_gradient(
	float** src_gray_im, 
	V3DF** dst_grad, 
	int w, int h, 
	float grad_thr
);
void get_gradient(
	cv::Mat &src_gray_img,
	cv::Mat &dst_grad,
	int w, int h, 
	float grad_thr
);


void get_tangent(
	V3DF** src_grad, 
	V3DF** dst_tangent, 
	int w, int h
);
void get_tangent(
	cv::Mat &src_grad, 
	cv::Mat &dst_tangent, 
	int w, int h
);

void get_ETF(
	V3DF** src_grad, V3DF** &dst_etf,
	int nbhd, int w, int h
);
void get_ETF(
	cv::Mat &src_grad, cv::Mat &dst_etf,
	int nbhd, int w, int h
);

void get_flow_path(
	V3DF** src_etf, 
	V3DF** src_grad, 
	FlowPath** dst_fpath, 
	int w, int h, int threshold_S, int threshold_T
);
void get_flow_path(
	cv::Mat &src_etf, 
	cv::Mat &src_grad, 
	FlowPath** dst_fpath, 
	int w, int h, int threshold_S, int threshold_T
);

void get_coherent_line(
	float** src_gray_im, V3DF** src_etf, FlowPath** src_fpath, float** dst_imCL,
 	int w, int h, float threshold_T, float CL_tanh_he_thr,
	float sigmaC, float sigmaM, float P, int iterations
);
void get_coherent_line(
	cv::Mat &src_gray_im, cv::Mat &src_etf, FlowPath** src_fpath, cv::Mat &dst_imCL,
 	int w, int h, float threshold_T, float CL_tanh_he_thr,
	float sigmaC, float sigmaM, float P, int iterations
);

void cl_set_flow_at_point_S(
	V3DF** src_etf, FlowPath* dst_fpath,
	int px, int py, int w, int h, int threshold_S /*, int* nPoints, V3DF *Alpha*/
);
void cl_set_flow_at_point_S(
	cv::Mat &src_etf, FlowPath* dst_fpath,
	int px, int py, int w, int h, int threshold_S /*, int* nPoints, V3DF *Alpha*/
);

void cl_set_flow_at_point_T(
	V3DF** src_grad, FlowPath* dst_fpath,
	int px, int py, int w, int h, int threshold_T /*int *nPoints, V3DF *Beta*/
);
void cl_set_flow_at_point_T(
	cv::Mat &src_grad, FlowPath* dst_fpath,
	int px, int py, int w, int h, int threshold_T /*int *nPoints, V3DF *Beta*/
);

void V3DF_interpolate(
	V3DF** src, V3DF dst, 
	float px, float py, int w, int h
);
void V3DF_interpolate(
	cv::Mat &src, V3DF dst, 
	float px, float py, int w, int h
);




// void get_vector_path(
// 	V3DF** src_etf, float** src_imCL, 
// 	FlowPath** dst_vpath, 
// 	int w, int h, int threshold_S
// );

// void vpath_set_flow_at_point_S(
// 	V3DF** src_etf, float** src_imCL,
// 	FlowPath* dst_vpath,
// 	int px, int py, int w, int h, int threshold_S
// );

// void vectorize_vector_path(
	
// 	int** dst_skltn
// );

// void _vectorize_build_skltn(
	
// );

//	vector length square
// float vlsquare ( V3DF v )
// {
// 	return vdot ( v, v );
// }

// //	vector length
// float vlength ( V3DF v )
// {
// 	return sqrt ( vlsquare ( v ) );
// }

// //	vector scalar multiplication
// void vscalar ( V3DF v, float sca )
// {
// 	v[0] *= sca;
// 	v[1] *= sca;
// 	v[2] *= sca;
// }
// void vscalar ( V3DI v, int sca )
// {
// 	v[0] *= sca;
// 	v[1] *= sca;
// 	v[2] *= sca;
// }

// //	vector add and assignment
// //	dest = dest + src
// void vadd ( V3DF dest, V3DF src ) 
// {
// 	dest[0] += src[0];
// 	dest[1] += src[1];
// 	dest[2] += src[2];
// }


// void RGBtoLAB(float iR, float iG, float iB, float *oL, float *oa, float *ob);
// void conv_RGBtoLAB(V3DF** buf, int width, int height);

#endif