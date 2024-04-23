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

    cv::Mat __image, _image, image;
    __image = cv::imread( argv[1], cv::IMREAD_COLOR );
    if ( !__image.data )
    {
        printf("No image data \n");
        return -1;
    }
	__image.convertTo(_image, CV_32FC3);
	cv::normalize(_image, image, 0.0f, 1.0f, cv::NORM_MINMAX);
    cv::namedWindow("Display Image", cv::WINDOW_AUTOSIZE );
    cv::imshow("Display Image", image);
    cv::waitKey(0);

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
	
	cv::Mat normalizedImage;
	cv::Mat displayImage;

	cv::normalize(grad, normalizedImage, 0, 255, cv::NORM_MINMAX, CV_32FC3);
	normalizedImage.convertTo(displayImage, CV_8UC3);

	cv::imshow("Display Image", img_gray);
    cv::waitKey(0);
	cv::imshow("Display Image", displayImage);
    cv::waitKey(0);

	std::cout << "gradient doen." << std::endl;

	cv::Mat tangent = cv::Mat::zeros(h, w, CV_32FC3);
	get_tangent(grad, tangent, w, h);
	std::cout << "tangent doen." << std::endl;

	cv::normalize(tangent, normalizedImage, 0, 255, cv::NORM_MINMAX, CV_32FC3);
	normalizedImage.convertTo(displayImage, CV_8UC3);
	cv::imshow("Display Image", displayImage); 
    cv::waitKey(0);

	// init_buf(etf, tangent, w, h);
	cv::Mat etf = tangent.clone();

	get_ETF(grad, etf, nbhd, w, h);
	get_ETF(grad, etf, nbhd, w, h);
	std::cout << "etf done." << std::endl;

	cv::normalize(etf, normalizedImage, 0, 255, cv::NORM_MINMAX, CV_32FC3);
	normalizedImage.convertTo(displayImage, CV_8UC3);
	cv::imshow("Display Image", displayImage); 
    cv::waitKey(0);

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
	
	cv::Mat normalizedImage_gray;
	cv::normalize(imCL, normalizedImage_gray, 0, 255, cv::NORM_MINMAX, CV_32F);
	normalizedImage_gray.convertTo(displayImage, CV_8UC1);
	cv::imshow("Display Image", displayImage); 
    cv::waitKey(0);





    return 0;
}