#include "main.h"
 


int main(int argc, char** argv ) {
    if ( argc != 2 )
    {
        printf("usage: DisplayImage.out <Image_Path>\n");
        return -1;
    }
 
    cv::Mat image;
    image = cv::imread( argv[1], cv::IMREAD_COLOR );
 
    if ( !image.data )
    {
        printf("No image data \n");
        return -1;
    }
    cv::namedWindow("Display Image", cv::WINDOW_AUTOSIZE );
    cv::imshow("Display Image", image);
    
    cv::waitKey(0);
 
    if (argc < 2) {
		printf("Usage: ./FDoG_GUI <src_img>\nExample: ./FDoG_GUI ./inputs/chainsaw.bmp");
		return 1;
	}

	std::string input_img_name = argv[1];
	std::cout << input_img_name << std::endl;

    cv::Mat img_c, img_gray;
    img_c = cv::imread( argv[1], cv::IMREAD_COLOR );
    if ( !img_c.data )
    {
        printf("No image data \n");
        return -1;
    }

    img_gray = cv::imread( argv[1], cv::IMREAD_GRAYSCALE );
    if ( !img_gray.data )
    {
        printf("No image data \n");
        return -1;
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

    cv::Size size = img_c.size();
    int w = size.width;
    int h = size.height;
	int cols = img_c.cols;
	int rows = img_c.rows;

	printf("w: %d, h: %d, cols: %d, rows: %d\n", w, h, cols, rows);

	// cv::Mat3f grad; // origin
	cv::Mat grad = cv::Mat::zeros(h, w, CV_32FC3);

	// get_gradient(img_gray, grad, w, h, grad_thr); // origin
	get_gradient(img_gray, grad, w, h, grad_thr);
	
	cv::Mat normalizedImage;
	cv::normalize(grad, normalizedImage, 0, 255, cv::NORM_MINMAX, CV_32FC3);
	cv::Mat displayImage;
	normalizedImage.convertTo(displayImage, CV_8UC3);

	cv::imshow("Display Image", img_gray);
    cv::waitKey(0);
	cv::imshow("Display Image", displayImage);
    cv::waitKey(0);

	std::cout << "gradient doen." << std::endl;

	cv::Mat tangent = cv::Mat::zeros(w, h, CV_32FC3);
	get_tangent(grad, tangent, w, h);
	std::cout << "tangent doen." << std::endl;

	// int nbhd = 5;
	// init_buf(etf, tangent, w, h);
	// get_ETF(grad, etf, nbhd, w, h);
	// get_ETF(grad, etf, nbhd, w, h);
	// std::cout << "etf done." << std::endl;


	// FlowPath** fpath = new FlowPath*[w];
	// FlowPath** vpath = new FlowPath*[w];
	// for(int i=0; i<w; ++i) {
	// 	fpath[i] = new FlowPath[h];
	// 	vpath[i] = new FlowPath[h];
	// 	for (int j=0; j<h; j++) {
	// 		fpath[i][j].Init();
	// 		vpath[i][j].Init();
	// 	}
	// }
	// std::cout << "fpath initialization done." << std::endl;

	// get_flow_path(etf, grad, fpath, w, h, threshold_S, threshold_T);
	// std::cout << "flow path done." << std::endl;

	// float** imCL = new float*[w];
	// for(int i=0; i<w; ++i) {
	// 	imCL[i] = new float[h];
	// }
	
	// get_coherent_line(img.im, etf, fpath, imCL, w, h, threshold_T, CL_tanh_he_thr, CL_sigma_c_line_width, CL_sigma_m_line_coherence, P, CL_iterations);
	// std::cout << "coherent line done." << std::endl;
	
	// for(int i=0; i<w; ++i) {
	// 	delete [] grad[i];
	// 	delete [] tangent[i];
	// 	delete [] etf[i];
	// }

	// delete grad;
	// delete tangent;
	// delete etf; 


    return 0;
}