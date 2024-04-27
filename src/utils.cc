#include "utils.h"


void draw_RGB_image(cv::Mat RGB_image, std::string title) {
    cv::namedWindow(title, cv::WINDOW_NORMAL );
	cv::cvtColor(RGB_image, RGB_image, cv::COLOR_RGB2BGR);
    cv::imshow(title, RGB_image);
}

void utils::draw_fpath(int x, int y) {
	float px_old, py_old, px_cur, py_cur;
	int pixel_size = 10;

	px_old = x, py_old = y;

	for(int r=0; r<utils::g_fpath[x][y].sn; ++r) {
		px_cur = utils::g_fpath[x][y].Alpha[r][0];	py_cur = utils::g_fpath[x][y].Alpha[r][1];
		cv::rectangle(utils::g_monitor_image, cv::Point(pixel_size*r, 0), cv::Point(pixel_size*(r+1), pixel_size), utils::g_images_for_monitor[utils::g_image_idx].at<cv::Vec3f>(py_old, px_old), -1);
		cv::line(utils::g_monitor_image, cv::Point(px_old, py_old), cv::Point(px_cur, py_cur), cv::Scalar(0, 0, 255), 1);
		cv::imshow("img", utils::g_monitor_image);
		px_old = px_cur; py_old = py_cur;
	}
}

void utils::on_mouse(int event, int x, int y, int flags, void*) {
	// 디버깅을 위해 매개변수로 들어온 flags 출력
	switch (event) {
    // 좌클릭이 되는 순간 ptOld에 좌표 저장
	case cv::EVENT_LBUTTONDOWN:
		switch (utils::g_image_idx) {
			case utils::FBL:
			case utils::INFODRAW:
				draw_fpath(x, y);
				break;
			default:
				utils::g_ptOld = cv::Point(x, y);
		}
		
	case cv::EVENT_LBUTTONUP:
		//cout << "EVENT_LBUTTONDOWN: " << x << ", " << y << endl;
		if (utils::g_ptNew.x != 0 && utils::g_ptNew.y != 0) {
			cv::rectangle(utils::g_monitor_image, utils::g_ptOld, utils::g_ptNew, cv::Scalar(255, 0, 0), utils::g_toggleFill);
			cv::imshow("img", utils::g_monitor_image);
			utils::g_ptNew = cv::Point(0, 0);
		}
		break;
	default:
		break;
	}	
}

void set_image(cv::Mat img) {
	utils::g_monitor_image = img.clone();
	cv::imshow("img", utils::g_monitor_image);
}

void utils::interactive_monitor(cv::Mat* imgs, FlowPath** fpath) {
	utils::g_images_for_monitor = imgs;
	utils::g_fpath = fpath;
	
	cv::namedWindow("img", cv::WINDOW_NORMAL);

	set_image(imgs[0]);

	// 마우스 콜백 함수를 등록하는 함수
	//void setMouseCallback(const String & winname, MouseCallback onMouse, void* userdata = 0);
	//winname : 이벤트 처리를 할 창의 이름, onMouse : 마우스콜백 함수 이름, userdata : 콜백함수에 전달한 데이터의 포인터로 전달할 값이 없으면 비움
	cv::setMouseCallback("img", utils::on_mouse);
	cv::imshow("img", utils::g_monitor_image);

	while (true) {
		int keycode = cv::waitKey();

		switch (keycode) {
			// 숫자 입력시 이미지 변환
			case '1': utils::g_image_idx = utils::ORIGIN;
				      set_image(imgs[utils::g_image_idx]); 	 break;
			case '2': utils::g_image_idx = utils::INFODRAW;
					  set_image(imgs[utils::g_image_idx]); 	 break;
			case '3': utils::g_image_idx = utils::GRAD;
					  set_image(imgs[utils::g_image_idx]); 	 break;
			case '4': utils::g_image_idx = utils::TANGENT;
					  set_image(imgs[utils::g_image_idx]); 	 break;
			case '5': utils::g_image_idx = utils::ETF;
					  set_image(imgs[utils::g_image_idx]); 	 break;
			case '6': utils::g_image_idx = utils::CL;
					  set_image(imgs[utils::g_image_idx]); 	 break;
			case '7': utils::g_image_idx = utils::FBL;
					  set_image(imgs[utils::g_image_idx]); 	 break;

			case 'c':
			case 'C':
				utils::g_monitor_image = imgs[utils::g_image_idx].clone();
				cv::imshow("img", utils::g_monitor_image);
				break;

			case 'q':
			case 'Q':
			case 27:
				return;
		}
	}
}