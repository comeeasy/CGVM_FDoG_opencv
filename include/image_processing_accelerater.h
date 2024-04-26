#pragma once

#ifndef __IMAGE_PROCESSING_ACCELERATER_H__
#define __IMAGE_PROCESSING_ACCELERATER_H__

#include <opencv2/opencv.hpp>

#include "common_func.h"


class GradientOperator {
public:
    GradientOperator(const cv::Mat& src, cv::Mat& dst, float grad_thr)
    : src_(src), dst_(dst), grad_thr_(grad_thr) {
        this->w = this->src_.cols;
        this->h = this->src_.rows;

        this->imG                       = new float *[this->w];
        this->src_gray_img_arr          = new float *[w];
        for (int i=0; i<this->w; i++){
            this->imG[i] 			    = new float[h];
            this->src_gray_img_arr[i]   = new float[h];
        }

        for(int i=0; i<this->w; ++i) {
            for(int j=0; j<this->h; ++j) {
                this->src_gray_img_arr[i][j] = this->src_.at<float>(j, i);
            }
        }
        
        /* Declare Sobel masks */
        GxMask[0][2] = -1.0f; GxMask[1][2] = 0.0f; GxMask[2][2] = 1.0f;
        GxMask[0][1] = -2.0f; GxMask[1][1] = 0.0f; GxMask[2][1] = 2.0f;
        GxMask[0][0] = -1.0f; GxMask[1][0] = 0.0f; GxMask[2][0] = 1.0f;

        GyMask[0][2] =  1.0f; GyMask[1][2] =  2.0f; GyMask[2][2] =  1.0f;
        GyMask[0][1] =  0.0f; GyMask[1][1] =  0.0f; GyMask[2][1] =  0.0f;
        GyMask[0][0] = -1.0f; GyMask[1][0] = -2.0f; GyMask[2][0] = -1.0f;

        GaussianFilter(this->src_gray_img_arr, this->imG, this->w, this->h);

    }

    void operator()(const cv::Range& range) const {
        int row,col;
        int colOffset,rowOffset;
        int colTotal,rowTotal;
        float Gx, Gy;

        for (row = range.start; row < range.end; row++) {
            for (col = 0; col < src_.cols; col++) {
                // Compute Gx, Gy using your existing logic here
                // Example:
                Gx = 0.0f; Gy = 0.0f;
                // Apply Sobel filter, check bounds etc.
                /* Calculate the sum of the Sobel mask times the nine surrounding pixels in the x and y direction */
                for (colOffset=-1; colOffset<=1; colOffset++) {
                    for (rowOffset=-1; rowOffset<=1; rowOffset++) {
                        colTotal = col+colOffset;
                        rowTotal = row+rowOffset;
                        //���ó��(��Ī)
                        if (rowTotal<0) rowTotal *= (-1.0f);
                        else if (rowTotal>=h)	rowTotal = h-(rowTotal-(h-1));
                        if (colTotal<0) colTotal *= (-1.0f);
                        else if (colTotal>=w)	colTotal = w-(colTotal-(w-1));

                        Gx = Gx + (imG[colTotal][rowTotal] * GxMask[colOffset+1][rowOffset+1]);//??
                        Gy = Gy + (imG[colTotal][rowTotal] * GyMask[colOffset+1][rowOffset+1]);//??
                    }
                }

                // Store results
                cv::Vec3f& pixel = dst_.at<cv::Vec3f>(row, col);
                float magnitude = sqrt(Gx * Gx + Gy * Gy);
                if (magnitude < grad_thr_) {
                    pixel = cv::Vec3f(0, 0, 0);
                } else {
                    pixel[0] = Gx;
                    pixel[1] = Gy;
                    pixel[2] = magnitude;
                }
            }
        }
    }

private:
    const cv::Mat& src_;
    cv::Mat& dst_;
    float grad_thr_;

    float **imG;
	float **src_gray_img_arr;

    int w, h;
    float GxMask[3][3] = {0.0f};
    float GyMask[3][3] = {0.0f};;
};

#endif