#include "main.h"


int main1(int argc, char** argv) {
	std::string input_img_path(argv[1]); 
	std::string seed_img_path(argv[2]);
	std::string output_dir(argv[3]);

	size_t num_workers = std::thread::hardware_concurrency();

	// ========== Adjust below parameters manually ==========
	// Set Direction parameters 
    float CL_sigma_c_line_width = 3.0f;
    float CL_sigma_m_line_coherence = 4.0f;
    float CL_tanh_he_thr = 0.998f;
	float CL_P = 0.99;
	int CL_iterations = 2;

	// grad parameters
    float GRAD_grad_thr = 0.0001f;

	// flow path parameters
	int FPath_threshold_S = 21;
	int FPath_threshold_T = 11;

	// etf parameters
	int ETF_nbhd = 5;
	int ETF_iteration = 2;
	// =======================================================

	// Read input image and infodraw images
	cv::Mat image = utils::read_RGB_normalized_image(input_img_path);
	cv::Mat img_gray = utils::read_Gray_normalized_image(input_img_path);

	cv::Mat seed_img = utils::read_RGB_normalized_image(seed_img_path);
	cv::Mat seed_img_gray = utils::read_Gray_normalized_image(seed_img_path);
	if (image.cols != seed_img.cols || image.rows != seed_img.rows) {
		printf("Size of image and pred_sketch must be the same.\n");
		exit(1);
	}
	int w = image.cols, h = image.rows;

	// ======================================
	// 1. Get gradient map
	// 2. Get tangent map
	// 3. Get ETF map
	// 4. Get flow pathes for input image
	// 5. Get Coherent Line image from input_gray image
	// 6. Get FBL filtered infodraw image
	// ++++++++++++++++++++++++++++++++++++++

	cv::Mat grad = get_gradient(img_gray, GRAD_grad_thr);
	cv::Mat tangent = get_tangent(grad);
	cv::Mat etf = get_ETF(grad, tangent, ETF_nbhd, ETF_iteration, num_workers);

	FlowPath** fpath = new FlowPath*[w];
	for(int i=0; i<w; ++i) {
		fpath[i] = new FlowPath[h];
		for (int j=0; j<h; j++) {
			fpath[i][j].Init();
		}
	}
	get_flow_path(etf, grad, fpath, w, h, FPath_threshold_S, FPath_threshold_T, num_workers);

	cv::Mat imCL = get_coherent_line(seed_img_gray, etf, fpath, FPath_threshold_T, CL_tanh_he_thr, CL_sigma_c_line_width, CL_sigma_m_line_coherence, CL_P, CL_iterations, num_workers);
	// std::string prefix = "_imCL_" + std::to_string(i);
	std::string prefix = "";
	utils::save_image(output_dir, input_img_path, imCL, prefix);

	// free fpath memory
	for(int i=0; i<w; ++i) {
		delete [] fpath[i];
	}
	delete [] fpath;

    return 0;
}

int main2(int argc, char** argv) {
	std::string input_img_path(argv[1]); 
	std::string output_dir(argv[2]);

	size_t num_workers = std::thread::hardware_concurrency();

	// ========== Adjust below parameters manually ==========
	// Set Direction parameters 
    float CL_sigma_c_line_width = 3.0f;
    float CL_sigma_m_line_coherence = 4.0f;
    float CL_tanh_he_thr = 0.998f;
	float CL_P = 0.99;
	int CL_iterations = 2;

	// grad parameters
    float GRAD_grad_thr = 0.0001f;

	// flow path parameters
	int FPath_threshold_S = 21;
	int FPath_threshold_T = 11;

	// etf parameters
	int ETF_nbhd = 5;
	int ETF_iteration = 2;
	// =======================================================

	// Read input image and infodraw images
	cv::Mat image = utils::read_RGB_normalized_image(input_img_path);
	cv::Mat img_gray = utils::read_Gray_normalized_image(input_img_path);
	int w = image.cols, h = image.rows;
	// ======================================
	// 1. Get gradient map
	// 2. Get tangent map
	// 3. Get ETF map
	// 4. Get flow pathes for input image
	// 5. Get Coherent Line image from input_gray image
	// 6. Get FBL filtered infodraw image
	// ++++++++++++++++++++++++++++++++++++++

	cv::Mat grad = get_gradient(img_gray, GRAD_grad_thr);
	cv::Mat tangent = get_tangent(grad);
	cv::Mat etf = get_ETF(grad, tangent, ETF_nbhd, ETF_iteration, num_workers);

	FlowPath** fpath = new FlowPath*[w];
	for(int i=0; i<w; ++i) {
		fpath[i] = new FlowPath[h];
		for (int j=0; j<h; j++) {
			fpath[i][j].Init();
		}
	}
	get_flow_path(etf, grad, fpath, w, h, FPath_threshold_S, FPath_threshold_T, num_workers);

	cv::Mat imCL = get_coherent_line(img_gray, etf, fpath, FPath_threshold_T, CL_tanh_he_thr, CL_sigma_c_line_width, CL_sigma_m_line_coherence, CL_P, CL_iterations, num_workers);
	// std::string prefix = "_imCL_" + std::to_string(i);
	std::string prefix = "";
	utils::save_image(output_dir, input_img_path, imCL, prefix);

	// free fpath memory
	for(int i=0; i<w; ++i) {
		delete [] fpath[i];
	}
	delete [] fpath;

    return 0;
}

int main(int argc, char** argv ) {
	if (argc < 1) {
		printf("Usage1: ./getCL <origin_img> <pred_sketch> <output_dir>\n");
		printf("Usage2: ./getCL <origin_img> <output_dir>\n");
		exit(1);
	}

	if (argc == 4) {
		printf("main1\n");
		return main1(argc, argv);
	} else if (argc == 3) {
		printf("main2\n");
		return main2(argc, argv);
	}

	return -1;
}