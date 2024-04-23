#include "main.h"
 


int main(int argc, char** argv ) {
	if (argc < 2) {
		printf("Usage: ./FDoG_GUI <src_img>\nExample: ./FDoG_GUI ./inputs/chainsaw.bmp");
		return 1;
	}
 
	// Set Direction parameters 
    int CL_is_rgb_or_lab = 1;
    int CL_itration = 3;
    float CL_l_weight = 0.9f;
    float CL_sigma_c_line_width = 1.0f;
    float CL_sigma_m_line_coherence = 5.0f;
    float CL_tanh_he_thr = 0.998f;
    float grad_thr = 0.0001f;
	int threshold_S = 21;
	int threshold_T = 11;
	float P = 0.99;
	int CL_iterations = 3;
	int nbhd = 5;

	// Set FBL params
	float FBL_sigma_e = 2.0f;
	float FBL_gamma_e = 50.0f;
	float FBL_sigma_g = 2.0f;
	float FBL_gamma_g = 10.0f;
	int FBL_threshold_T = threshold_T;
	int FBL_iteration = 3;


    cv::Mat __image, _image, image;
    __image = cv::imread( argv[1], cv::IMREAD_COLOR );
    if ( !__image.data )
    {
        printf("No image data \n");
        return -1;
    }
	__image.convertTo(_image, CV_32FC3);
	cv::normalize(_image, image, 0.0f, 1.0f, cv::NORM_MINMAX);

    cv::Mat __img_gray, _img_gray, img_gray;
    __img_gray = cv::imread( argv[1], cv::IMREAD_GRAYSCALE );
    if ( !__img_gray.data )
    {
        printf("No image data \n");
        return -1;
    }
	__img_gray.convertTo(_img_gray, CV_32F);
	cv::normalize(_img_gray, img_gray, 0.0f, 1.0f, cv::NORM_MINMAX);

    cv::Size size = image.size();
    int w = size.width;
    int h = size.height;
	int cols = image.cols;
	int rows = image.rows;

	printf("w: %d, h: %d, cols: %d, rows: %d\n", w, h, cols, rows);

	// cv::Mat3f grad; // origin
	cv::Mat grad = cv::Mat::zeros(h, w, CV_32FC3);

	// get_gradient(img_gray, grad, w, h, grad_thr); // origin
	get_gradient(img_gray, grad, w, h, grad_thr);
	std::cout << "gradient doen." << std::endl;

	cv::Mat tangent = cv::Mat::zeros(h, w, CV_32FC3);
	get_tangent(grad, tangent, w, h);
	std::cout << "tangent doen." << std::endl;

	cv::Mat etf = tangent.clone();
	get_ETF(grad, etf, nbhd, w, h);
	get_ETF(grad, etf, nbhd, w, h);
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
	std::cout << "fpath initialization done." << std::endl;

	get_flow_path(etf, grad, fpath, w, h, threshold_S, threshold_T);
	std::cout << "flow path done." << std::endl;

	cv::Mat imCL = cv::Mat::zeros(h, w, CV_32F);
	get_coherent_line(img_gray, etf, fpath, imCL, w, h, threshold_T, CL_tanh_he_thr, CL_sigma_c_line_width, CL_sigma_m_line_coherence, P, CL_iterations);
	std::cout << "coherent line done." << std::endl;
	
	cv::Mat infodraw_img = cv::imread(argv[2], cv::IMREAD_COLOR);
	infodraw_img.convertTo(infodraw_img, CV_32FC3);
	cv::normalize(infodraw_img, infodraw_img, 0.0f, 1.0f, cv::NORM_MINMAX);
	cv::Size info_size = infodraw_img.size();
	if (size != info_size) {
		printf("input image and info draw image must be matched!\n");
		exit(1);
	}
	
	cv::Mat FBL_res = cv::Mat::zeros(info_size, CV_32FC3);
	apply_FBL_filter(infodraw_img, fpath, FBL_res, w, h, FBL_sigma_e, FBL_gamma_e, FBL_sigma_g, FBL_gamma_g, FBL_iteration, FBL_threshold_T);


	//== draw images==
	// origin image
	cv::namedWindow("Original image", cv::WINDOW_AUTOSIZE );
    cv::imshow("Original image", image);

	// gray image
	cv::namedWindow("Gray image", cv::WINDOW_AUTOSIZE );
    cv::imshow("Gray image", img_gray);

	// gradient
	cv::namedWindow("Gradient Image", cv::WINDOW_AUTOSIZE );
	cv::normalize(grad, grad, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32FC3);
	cv::imshow("Gradient Image", grad);

	// tangent 
	cv::namedWindow("Tangent Image", cv::WINDOW_AUTOSIZE );
	cv::normalize(tangent, tangent, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32FC3);
	cv::imshow("Tangent Image", tangent); 

	// ETF
	cv::namedWindow("ETF Image", cv::WINDOW_AUTOSIZE );
	cv::normalize(etf, etf, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32FC3);
	cv::imshow("ETF Image", etf); 

	// CL
	cv::namedWindow("CL Image", cv::WINDOW_AUTOSIZE );
	cv::normalize(imCL, imCL, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32F);
	cv::imshow("CL Image", imCL); 

	// FBL
	cv::namedWindow("FBL Image", cv::WINDOW_AUTOSIZE );
	cv::normalize(FBL_res, FBL_res, 0.0f, 1.0f, cv::NORM_MINMAX, CV_32F);
	cv::imshow("FBL Image", FBL_res);



    cv::waitKey(0);

	for(int i=0; i<w; ++i) {
		delete [] fpath[i];
		delete [] vpath[i];
	}
	delete [] fpath;
	delete [] vpath;

    return 0;
}