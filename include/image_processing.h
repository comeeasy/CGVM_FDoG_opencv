#pragma once

#ifndef __IMAGE_PROCESSING_H__
#define __IMAGE_PROCESSING_H__

#include <time.h>
#include <thread>
#include <vector>

#include "mtoolkits.h"
#include "common_func.h"
#include "conv_clr.h"

#include "FlowPath.h"
#include "Label.h"

#define NEGATIVE 0
#define POSITIVE 1


// Gradients
//=====================================
//	Legacy code
//=====================================
void get_gradient(float** src_gray_im, V3DF** dst_grad, int w, int h, float grad_thr);
void get_gradient(cv::Mat &src_gray_img, cv::Mat &dst_grad, int w, int h, float grad_thr);
void get_gradient(const cv::Mat& src, cv::Mat& dst, float grad_thr);

//====================================
// Usage Example
//
// cv::Mat get_gradient(cv::Mat &src_gray_img, float grad_thr);
//====================================
//
// cv::Mat img_gray = cv::imread( input_img_path, cv::IMREAD_GRAYSCALE );
// if ( !img_gray.data )
// {
//     printf("No image data \n");
//     return -1;
// }
// img_gray.convertTo(img_gray, CV_32F);
// cv::normalize(img_gray, img_gray, 0.0f, 1.0f, cv::NORM_MINMAX);
//
// cv::Mat grad = get_gradient(img_gray, grad_thr);
//
//====================================
cv::Mat get_gradient(cv::Mat &src_gray_img, float grad_thr);


// Tangemt
//=====================================
//	Legacy code
//=====================================
void get_tangent(V3DF** src_grad, V3DF** dst_tangent, int w, int h);
void get_tangent(cv::Mat &src_grad, cv::Mat &dst_tangent, int w, int h);

//====================================
// Usage Example
//
// cv::Mat get_tangent(cv::Mat &src_grad);
//====================================
//
// cv::Mat grad = get_gradient(img_gray, grad_thr);
// 
// **Check get_gradient function usage (line 31)
//
//====================================
cv::Mat get_tangent(cv::Mat &src_grad);

// ETF (Edge Tangent Flow)
void get_ETF(V3DF** src_grad, V3DF** &dst_etf, int nbhd, int w, int h);
void get_ETF(cv::Mat &src_grad, cv::Mat &dst_etf, int nbhd, int w, int h);
void get_ETF(cv::Mat &src_grad, cv::Mat &dst_etf,int nbhd, int w, int h, int num_workers);

cv::Mat get_ETF(cv::Mat &src_grad, cv::Mat &src_tangent, int nbhd, int iteration);
cv::Mat get_ETF(cv::Mat &src_grad, cv::Mat &src_tangent, int nbhd, int iteration, size_t num_workers);

// Flow Path
void get_flow_path(V3DF** src_etf, V3DF** src_grad, FlowPath** dst_fpath, int w, int h, int threshold_S, int threshold_T);

void get_flow_path(cv::Mat &src_etf, cv::Mat &src_grad, FlowPath** dst_fpath, int w, int h, int threshold_S, int threshold_T);
void get_flow_path(cv::Mat &src_etf, cv::Mat &src_grad, FlowPath** dst_fpath, int w, int h, int threshold_S, int threshold_T, int num_workers);

// Coherent Line
void get_coherent_line(float** src_gray_im, V3DF** src_etf, FlowPath** src_fpath, float** dst_imCL, int w, int h, float threshold_T, float CL_tanh_he_thr, float sigmaC, float sigmaM, float P, int iterations);
void get_coherent_line(cv::Mat &src_gray_im, cv::Mat &src_etf, FlowPath** src_fpath, cv::Mat &dst_imCL, int w, int h, float threshold_T, float CL_tanh_he_thr, float sigmaC, float sigmaM, float P, int iterations);
void get_coherent_line(cv::Mat &src_gray_im, cv::Mat &src_etf, FlowPath** src_fpath, cv::Mat &dst_imCL, int w, int h, float threshold_T, float CL_tanh_he_thr, float sigmaC, float sigmaM, float P, int iterations, int num_workers);

cv::Mat get_coherent_line(cv::Mat &src_gray_im, cv::Mat &src_etf, FlowPath** src_fpath, float threshold_T, float CL_tanh_he_thr, float sigmaC, float sigmaM, float P, int iterations);
cv::Mat get_coherent_line(cv::Mat &src_gray_im, cv::Mat &src_etf, FlowPath** src_fpath, float threshold_T, float CL_tanh_he_thr, float sigmaC, float sigmaM, float P, int iterations, int num_workers);

void cl_set_flow_at_point_S(V3DF** src_etf, FlowPath* dst_fpath, int px, int py, int w, int h, int threshold_S);
void cl_set_flow_at_point_S(cv::Mat &src_etf, FlowPath* dst_fpath, int px, int py, int w, int h, int threshold_S);

void cl_set_flow_at_point_T(V3DF** src_grad, FlowPath* dst_fpath, int px, int py, int w, int h, int threshold_T);
void cl_set_flow_at_point_T(cv::Mat &src_grad, FlowPath* dst_fpath, int px, int py, int w, int h, int threshold_T);

void V3DF_interpolate(V3DF** src, V3DF dst, float px, float py, int w, int h);
void V3DF_interpolate(cv::Mat &src, V3DF dst, float px, float py, int w, int h);

// Flow-Based Bilateral Filter
//=====================================
//	Legacy code
//=====================================
void apply_FBL_filter(V3DF** src_cim, FlowPath** src_fpath, V3DF** cimFBL, int w, int h, float sigma_e, float gamma_e, float sigma_g, float gamma_g, int threshold_T, int iteration);
void apply_FBL_filter(cv::Mat &src_cim, FlowPath** src_fpath, cv::Mat &cimFBL, int w, int h, float sigma_e, float gamma_e, float sigma_g, float gamma_g, int threshold_T, int iteration);
void apply_FBL_filter_parellel(cv::Mat &src_cim, FlowPath** src_fpath, cv::Mat &cimFBL, int w, int h, float sigma_e, float gamma_e, float sigma_g, float gamma_g, int threshold_T, int iteration, size_t num_workers);
//====================================

//====================================
// Usage Example
//====================================
//
//
//
//
//
//
//====================================
cv::Mat apply_FBL_filter(cv::Mat &src_cim, FlowPath** src_fpath, float sigma_e, float gamma_e, float sigma_g, float gamma_g, int threshold_T, int iteration);
cv::Mat apply_FBL_filter(cv::Mat &src_cim, FlowPath** src_fpath, float sigma_e, float gamma_e, float sigma_g, float gamma_g, int threshold_T, int iteration, size_t num_workers);

// Segmentation
void get_segmentation (cv::Mat &InImg, int** dst_sgid, cv::Mat &dst_cimSG, int hs, int hr, int M, int w, int h);
void get_segmentation (V3DF** InImg, int** dst_sgid, V3DF** dst_cimSG, int hs, int hr, int M, int w, int h);

// Merging FBl result + segmentation
void sum_FDoG_FBL(cv::Mat &src_imCL, cv::Mat &src_cimSG, cv::Mat &dst_cimsum, int w, int h);

#endif