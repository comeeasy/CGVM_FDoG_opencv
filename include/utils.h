#pragma once

#ifndef __UTILS_H__
#define __UTILS_H__

#include <string>
#include <filesystem>
#include <opencv2/opencv.hpp>

#include "FlowPath.h"


namespace utils
{
    void save_fpath_to_npy(std::string output_dir, std::string input_img_path, cv::Mat &img, FlowPath** fpath, int w, int h, int threshold_S);

    cv::Mat read_RGB_normalized_image(std::string path);
    cv::Mat read_Gray_normalized_image(std::string path);

    void save_image(std::string output_dir, std::string input_img_path, cv::Mat norm_img, std::string postfix);

} // namespace utils


#endif