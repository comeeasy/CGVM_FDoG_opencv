#include "utils.h"



void utils::save_RGBpixel_of_fpath_to_npy(std::string output_dir, std::string input_img_path, cv::Mat &img, FlowPath** fpath, int w, int h, int threshold_S) {

    float px_old, py_old, px_cur, py_cur;

    std::vector<float> data;
    for(int x=0; x<w; ++x) {
        for(int y=0; y<h; ++y) {
            // Threshold_S is defined out of this space.
            // Please be careful to use.
            int sn = fpath[x][y].sn;
            for(int r=0; r<threshold_S; ++r) {
                if (r < sn) { // insert fpath pixel values if there is fpath 
                    px_cur = fpath[x][y].Alpha[r][0];	py_cur = fpath[x][y].Alpha[r][1];
                
                    cv::Vec3f &pixel = img.at<cv::Vec3f>(py_cur, px_cur);
                    for(int c=0; c<3; ++c) { // 3 for RGB
                        data.push_back(pixel[c]);                    
                    }
                } else { // if fpath is shorter than threshold_S, then fill it with 0s.
                    for(int c=0; c<3; ++c) { // 3 for RGB
                        data.push_back(1.0f);                    
                    }
                }
            }
        }
    }

    npy::npy_data_ptr<float> d;
	d.data_ptr = data.data();
	d.shape = { (unsigned long)w, (unsigned long)h, (unsigned long)threshold_S, 3 };
	d.fortran_order = false; // optional

    std::filesystem::path input_file_path(input_img_path);
	std::filesystem::path base_output_path(output_dir);
	std::string base_name = input_file_path.stem().string();

	std::filesystem::path final_output_dir = base_output_path / "fpath_npzs";
	std::filesystem::create_directories(final_output_dir);

    std::string final_path;
	final_path = final_output_dir.string() + "/" + base_name + "_fpath_of_infodraw.npy";
	npy::write_npy(final_path, d);
}

void utils::save_Graypixel_of_fpath_to_npy(std::string output_dir, std::string input_img_path, cv::Mat &img, FlowPath** fpath, int w, int h, int threshold_S) {
    float px_old, py_old, px_cur, py_cur;
    float pixel;

    std::vector<float> data;
    for(int x=0; x<w; ++x) {
        for(int y=0; y<h; ++y) {
            // Threshold_S is defined out of this space.
            // Please be careful to use.
            int sn = fpath[x][y].sn;
            for(int r=0; r<threshold_S; ++r) {
                if (r < sn) { // insert fpath pixel values if there is fpath 
                    px_cur = fpath[x][y].Alpha[r][0];	py_cur = fpath[x][y].Alpha[r][1];
                    
                    // Must save interpolated value
                    pixel = float_interpolate(img, px_cur, py_cur, w, h);

                    // ??
                    // pixel = img.at<float>(py_cur, px_cur);

                    data.push_back(pixel);
                } else { // if fpath is shorter than threshold_S, then fill it with 0s.
                    data.push_back(1.0f);                    
                }
            }
        }
    }

    npy::npy_data_ptr<float> d;
	d.data_ptr = data.data();
	d.shape = { (unsigned long)w, (unsigned long)h, (unsigned long)threshold_S};
	d.fortran_order = false; // optional

    std::filesystem::path input_file_path(input_img_path);
	std::filesystem::path base_output_path(output_dir);
	std::string base_name = input_file_path.stem().string();

	std::filesystem::path final_output_dir = base_output_path / "fpath_npzs";
	std::filesystem::create_directories(final_output_dir);

    std::string final_path;
	final_path = final_output_dir.string() + "/" + base_name + "_fpath_of_infodraw.npy";
	npy::write_npy(final_path, d);
}


void utils::save_fpath(std::string output_dir, std::string input_img_path, FlowPath** fpath, int w, int h, int threshold_S) {
    float px_old, py_old, px_cur, py_cur;

    std::vector<float> data;
    for(int x=0; x<w; ++x) {
        for(int y=0; y<h; ++y) {
            // Threshold_S is defined out of this space.
            // Please be careful to use.
            int sn = fpath[x][y].sn;
            for(int r=0; r<threshold_S; ++r) {
                if (r < sn) { // insert fpath pixel values if there is fpath 
                    px_cur = fpath[x][y].Alpha[r][0];	py_cur = fpath[x][y].Alpha[r][1];
                    data.push_back(px_cur);
                    data.push_back(py_cur);
                } else { // if fpath is shorter than threshold_S, then fill it with 0s.
                    data.push_back(-1.0f);
                    data.push_back(-1.0f);
                }
            }
        }
    }

    npy::npy_data_ptr<float> d;
	d.data_ptr = data.data();
	d.shape = { (unsigned long)w, (unsigned long)h, (unsigned long)threshold_S, 2};
	d.fortran_order = false; // optional

    std::filesystem::path input_file_path(input_img_path);
	std::filesystem::path base_output_path(output_dir);
	std::string base_name = input_file_path.stem().string();

	std::filesystem::path final_output_dir = base_output_path / "fpath_npzs";
	std::filesystem::create_directories(final_output_dir);

    std::string final_path;
	final_path = final_output_dir.string() + "/" + base_name + "_fpath.npy";
	npy::write_npy(final_path, d);
}

cv::Mat utils::read_RGB_normalized_image(std::string path) {
    cv::Mat image = cv::imread( path, cv::IMREAD_COLOR );
    if ( !image.data )
    {
        printf("No image data \n");
        exit(1);
    }
	cv::cvtColor(image, image, cv::COLOR_BGR2RGB);
	image.convertTo(image, CV_32FC3);
	cv::normalize(image, image, 0.0f, 1.0f, cv::NORM_MINMAX);

    return image;
}

cv::Mat utils::read_Gray_normalized_image(std::string path) {
    cv::Mat img_gray = cv::imread( path, cv::IMREAD_GRAYSCALE );
    if ( !img_gray.data )
    {
        printf("No image data \n");
        exit(1);
    }
	img_gray.convertTo(img_gray, CV_32F);
	cv::normalize(img_gray, img_gray, 0.0f, 1.0f, cv::NORM_MINMAX);

    return img_gray;
}

void utils::save_image(std::string output_dir, std::string input_img_path, cv::Mat norm_img, std::string postfix) {
    // pixel values of image must have range 0 to 1 
    // normalized image must be passed
    // condition check logic must be added later.

    std::filesystem::path input_file_path(input_img_path);
	std::filesystem::path base_output_path(output_dir);
	std::string base_name = input_file_path.stem().string();

	std::filesystem::path final_output_dir = base_output_path / base_name;
	std::filesystem::create_directories(final_output_dir);

	std::string final_path;
	final_path = final_output_dir.string() + "/" + base_name + postfix + ".png";
	cv::imwrite(final_path, norm_img * 255);
}
