#include "main.h"




int main(int argc, char** argv ) {
	if (argc < 3) {
		printf("Usage: ./FDoG_GUI <src_img>\nExample: ./FDoG_GUI ./inputs/chainsaw.bmp");
		return 1;
	}
 
	std::string input_img_path(argv[1]); 
	std::string infodraw_img_path(argv[2]);

	// Set Direction parameters 
    // int CL_is_rgb_or_lab = 1;
    // float CL_l_weight = 0.9f;
    float CL_sigma_c_line_width = 1.0f;
    float CL_sigma_m_line_coherence = 5.0f;
    float CL_tanh_he_thr = 0.998f;
    float grad_thr = 0.0001f;
	int threshold_S = 21;
	int threshold_T = 11;
	float P = 0.99;
	int CL_iterations = 3;
	int nbhd = 5;

	cv::Mat image = cv::imread( input_img_path, cv::IMREAD_COLOR );
    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }
	cv::cvtColor(image, image, cv::COLOR_BGR2RGB);
	image.convertTo(image, CV_32FC3);
	cv::normalize(image, image, 0.0f, 1.0f, cv::NORM_MINMAX);

    cv::Mat img_gray = cv::imread( input_img_path, cv::IMREAD_GRAYSCALE );
    if ( !img_gray.data )
    {
        printf("No image data \n");
        return -1;
    }
	img_gray.convertTo(img_gray, CV_32F);
	cv::normalize(img_gray, img_gray, 0.0f, 1.0f, cv::NORM_MINMAX);

    cv::Size size = image.size();
    int w = size.width;
    int h = size.height;
	int cols = image.cols;
	int rows = image.rows;

	printf("w: %d, h: %d, cols: %d, rows: %d\n", w, h, cols, rows);

	cv::Mat grad = get_gradient(img_gray, grad_thr);
	std::cout << "gradient doen." << std::endl;

	cv::Mat tangent = get_tangent(grad);
	std::cout << "tangent doen." << std::endl;

	int iter_etf = 2;
	size_t num_workers = std::thread::hardware_concurrency();
	// size_t num_workers = 32;
	std::cout << "num_workers: " << num_workers << std::endl;

	cv::Mat etf = get_ETF(grad, tangent, nbhd, iter_etf, num_workers);
	// cv::Mat etf = get_ETF(grad, tangent, nbhd, iter_etf);
	std::cout << "etf done." << std::endl;

	FlowPath** fpath = new FlowPath*[w];
	FlowPath** vpath = new FlowPath*[w];
	for(int i=0; i<w; ++i) {
		fpath[i] = new FlowPath[h];
		vpath[i] = new FlowPath[h];
		for (int j=0; j<h; j++) {
			fpath[i][j].Init();
			vpath[i][j].Init();
		}
	}
	get_flow_path(etf, grad, fpath, w, h, threshold_S, threshold_T, num_workers);
	// get_flow_path(etf, grad, fpath, w, h, threshold_S, threshold_T);
	std::cout << "flow path done." << std::endl;

	cv::Mat imCL = get_coherent_line(img_gray, etf, fpath, threshold_T, CL_tanh_he_thr, CL_sigma_c_line_width, CL_sigma_m_line_coherence, P, CL_iterations, num_workers);
	// cv::Mat imCL = get_coherent_line(img_gray, etf, fpath, threshold_T, CL_tanh_he_thr, CL_sigma_c_line_width, CL_sigma_m_line_coherence, P, CL_iterations);
	std::cout << "coherent line done." << std::endl;
	
	cv::Mat infodraw_img = cv::imread(infodraw_img_path, cv::IMREAD_COLOR);
	if ( !infodraw_img.data ) {
        printf("No image data \n");
        exit(1);
    }
	cv::Size info_size = infodraw_img.size();
	if (size != info_size) {
		printf("input image and info draw image must be matched!\n");
		exit(1);
	}

	cv::cvtColor(infodraw_img, infodraw_img, cv::COLOR_BGR2RGB);
	infodraw_img.convertTo(infodraw_img, CV_32FC3);
	cv::normalize(infodraw_img, infodraw_img, 0.0f, 1.0f, cv::NORM_MINMAX);
	

	// Set FBL params
	float FBL_sigma_e = 2.0f;
	float FBL_gamma_e = 50.0f;
	float FBL_sigma_g = 2.0f;
	float FBL_gamma_g = 10.0f;
	int FBL_threshold_T = 3;
	int FBL_iteration = 5;
	cv::Mat cimFBL = apply_FBL_filter(infodraw_img, fpath, FBL_sigma_e, FBL_gamma_e, FBL_sigma_g, FBL_gamma_g, FBL_threshold_T, FBL_iteration, num_workers);
	// cv::Mat cimFBL = apply_FBL_filter(infodraw_img, fpath, FBL_sigma_e, FBL_gamma_e, FBL_sigma_g, FBL_gamma_g, FBL_threshold_T, FBL_iteration);

	int** sgid = new int*[w];
	for(int i=0; i<w; ++i) {
		sgid[i] = new int[h];
	}
	// cv::Mat cimSG = cv::Mat::zeros(h, w, CV_32FC3);
	// int hs=10, hr=20, M=300;
	// get_segmentation(cimFBL, sgid, cimSG, hs, hr, M, w, h);

	// cv::Mat cimSum = cv::Mat::zeros(h, w, CV_32FC3);
	// sum_FDoG_FBL(imCL, cimSG, cimSum, w, h); // 없어도 괜찮음.

	//== draw images==
	// origin image
	// cv::namedWindow("Original image", cv::WINDOW_NORMAL );
	// cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
    // cv::imshow("Original image", image);

	// // gray image
	// cv::namedWindow("Gray image", cv::WINDOW_GUI_NORMAL);
    // cv::imshow("Gray image", img_gray);

	// // gradient
	// cv::namedWindow("Gradient Image", cv::WINDOW_NORMAL );
	// cv::normalize(grad, grad, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32FC3);
	// cv::imshow("Gradient Image", grad);

	// // tangent 
	// cv::namedWindow("Tangent Image", cv::WINDOW_NORMAL );
	// cv::normalize(tangent, tangent, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32FC3);
	// cv::imshow("Tangent Image", tangent); 

	// // ETF
	// cv::namedWindow("ETF Image", cv::WINDOW_NORMAL );
	// cv::normalize(etf, etf, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32FC3);
	// cv::imshow("ETF Image", etf); 

	// // CL
	// cv::namedWindow("CL Image", cv::WINDOW_NORMAL );
	// cv::normalize(imCL, imCL, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32F);
	// cv::imshow("CL Image", imCL); 

	// // infodraw
	// cv::namedWindow("InfoDraw Image", cv::WINDOW_NORMAL );
	// cv::normalize(infodraw_img, infodraw_img, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32F);
	// cv::imshow("InfoDraw Image", infodraw_img); 

	// // FBL
	// cv::namedWindow("FBL Image", cv::WINDOW_NORMAL );
	// cv::normalize(cimFBL, cimFBL, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32F);
	// cv::imshow("FBL Image", cimFBL);

	// image save
	std::filesystem::path input_file_path(input_img_path);
	std::filesystem::path base_output_path("results");
	std::string base_name = input_file_path.stem().string();

	std::filesystem::path output_dir = base_output_path / base_name;
	std::filesystem::create_directories(output_dir);

	std::string final_path;
	final_path = output_dir.string() + "/" + base_name + "_imCL.png";
	cv::imwrite(final_path, imCL * 255);
	std::cout << final_path << " is saved." << std::endl;
	final_path = output_dir.string() + "/" + base_name + "_FBL.png";
	cv::imwrite(final_path, cimFBL * 255);
	std::cout << final_path << " is saved." << std::endl;
	final_path = output_dir.string() + "/" + base_name + "_grad.png";
	cv::imwrite(final_path, grad * 255);
	std::cout << final_path << " is saved." << std::endl;
	final_path = output_dir.string() + "/" + base_name + "_tangent.png";
	cv::imwrite(final_path, tangent * 255);
	std::cout << final_path << " is saved." << std::endl;
	final_path = output_dir.string() + "/" + base_name + "_etf.png";
	cv::imwrite(final_path, etf * 255);
	std::cout << final_path << " is saved." << std::endl;
	final_path = output_dir.string() + "/" + base_name + "_origin.png";
	cv::imwrite(final_path, image * 255);
	std::cout << final_path << " is saved." << std::endl;
	final_path = output_dir.string() + "/" + base_name + "_infodraw.png";
	cv::imwrite(final_path, infodraw_img * 255);
	std::cout << final_path << " is saved." << std::endl;

	// Interactive
	

	cv::Mat imgs[7];
	imgs[utils::ORIGIN] = image;
	imgs[utils::INFODRAW] = infodraw_img;
	imgs[utils::GRAD] = grad;
	imgs[utils::TANGENT] = tangent;
	imgs[utils::ETF] = etf;
	imgs[utils::CL] = imCL;
	imgs[utils::FBL] = cimFBL;
	utils::interactive_monitor(imgs, fpath);

	for(int i=0; i<w; ++i) {
		delete [] fpath[i];
		delete [] vpath[i];
		delete [] sgid[i];
	}
	delete [] fpath;
	delete [] vpath;
	delete [] sgid;

    return 0;
}