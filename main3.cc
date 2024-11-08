#include "main.h"


int main(int argc, char** argv ) {
	if (argc < 2) {
		printf("Usage: ./FDoG_GUI <input_img> <output_dir>");
		exit(1);
	}

	// set inputs
	std::string input_img_path(argv[1]); 
	std::string output_dir = argv[2];

	size_t num_workers = std::thread::hardware_concurrency();

	// ========== Adjust below parameters manually ==========
	// Set Direction parameters 
    float CL_sigma_c_line_width = 1.0f;
    float CL_sigma_m_line_coherence = 5.0f;
    float CL_tanh_he_thr = 0.998f;
	float CL_P = 0.99;
	int CL_iterations = 3;

	// grad parameters
    float GRAD_grad_thr = 0.0001f;

	// flow path parameters
	int FPath_threshold_S = 21;
	int FPath_threshold_T = 11;

	// etf parameters
	int ETF_nbhd = 5;
	int ETF_iteration = 2;

	// FBL image parameters
	float FBL_sigma_e = 2.0f;
	float FBL_gamma_e = 10.0f;
	float FBL_sigma_g = 0.3f;
	float FBL_gamma_g = 10.0f;
	int FBL_threshold_T = 3;
	int FBL_iteration = 1;
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
	cv::Mat cimFBL = apply_FBL_filter(image, fpath, FBL_sigma_e, FBL_gamma_e, FBL_sigma_g, FBL_gamma_g, FBL_threshold_T, FBL_iteration, num_workers);

	std::string prefix = "";
	utils::save_image(output_dir, input_img_path, cimFBL, prefix);

	// free fpath memory
	for(int i=0; i<w; ++i) {
		delete [] fpath[i];
	}
	delete [] fpath;

    return 0;
}