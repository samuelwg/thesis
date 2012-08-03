
#include <pthread.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <libfreenect/libfreenect.h>
#include <cstdio>
#include <iostream>
#include <GL/gl.h>


pthread_t freenect_thread;
volatile int die = 0;

int g_argc;
char **g_argv;

int window;

// back: owned by libfreenect (implicit for depth)
// mid: owned by callbacks, "latest frame ready"
// front: owned by GL, "currently being drawn"
uint8_t *depth_mid, *depth_front;
uint8_t *rgb_back, *rgb_mid, *rgb_front;

GLuint gl_depth_tex;
GLuint gl_rgb_tex;

int freenect_angle = 0;
int freenect_led;

freenect_video_format requested_format = FREENECT_VIDEO_RGB;
freenect_video_format current_format = FREENECT_VIDEO_RGB;

pthread_cond_t gl_frame_cond = PTHREAD_COND_INITIALIZER;
int got_rgb = 0;
int got_depth = 0;
unsigned got_frames = 0;


cv::Mat depthMat(cv::Size(640,480),CV_16UC1);
cv::Mat rgbMat(cv::Size(640,480),CV_8UC3,cv::Scalar(0));
pthread_t fnkt_thread;
freenect_device *f_dev;
pthread_mutex_t buf_mutex = PTHREAD_MUTEX_INITIALIZER;
freenect_context *f_ctx;
pthread_cond_t frame_cond = PTHREAD_COND_INITIALIZER;

// Review: Those is maybe due to a rename in libfreenect
#define FREENECT_DEPTH_SIZE FREENECT_DEPTH_11BIT_SIZE
#define FREENECT_RGB_SIZE FREENECT_VIDEO_RGB_SIZE
#define FREENECT_FORMAT_RGB FREENECT_VIDEO_RGB
#define FREENECT_FORMAT_11_BIT FREENECT_DEPTH_11BIT

void *freenect_threadfunc(void* arg) {
	std::cout << "freenect thread"<<std::endl;
	while(!die && freenect_process_events(f_ctx) >= 0 ) {}
	std::cout << "freenect die"<<std::endl;
	return NULL;
}


void depth_cb(freenect_device *dev, void * depth, uint32_t timestamp)
{
	pthread_mutex_lock(&buf_mutex);
 
	//copy to ocv buf...
	memcpy(depthMat.data, depth, FREENECT_DEPTH_SIZE);
 
	got_frames++;
	pthread_cond_signal(&frame_cond);
	pthread_mutex_unlock(&buf_mutex);
}
 
void rgb_cb(freenect_device *dev, void *rgb, uint32_t timestamp)
{
	pthread_mutex_lock(&buf_mutex);
	got_frames++;
	//copy to ocv_buf..
	memcpy(rgbMat.data, rgb, FREENECT_RGB_SIZE);
 
	pthread_cond_signal(&frame_cond);
	pthread_mutex_unlock(&buf_mutex);
}





int main(int argc, char **argv)
{
	int res;
 
	g_argc = argc;
	g_argv = argv;
 
	if (freenect_init(&f_ctx, NULL) < 0) {
		printf("freenect_init() failed\n");
		return 1;
	}
 
	freenect_set_log_level(f_ctx, FREENECT_LOG_INFO);
 
	int nr_devices = freenect_num_devices (f_ctx);
	printf ("Number of devices found: %d\n", nr_devices);
 
	int user_device_number = 0;
	if (argc > 1)
		user_device_number = atoi(argv[1]);
 
	if (nr_devices < 1)
		return 1;
 
	if (freenect_open_device(f_ctx, &f_dev, user_device_number) < 0) {
		printf("Could not open device\n");
		return 1;
	}
 
	freenect_set_tilt_degs(f_dev,freenect_angle);
	freenect_set_led(f_dev,LED_RED);
	freenect_set_depth_callback(f_dev, depth_cb);
	freenect_set_video_callback(f_dev, rgb_cb);
	freenect_set_video_format(f_dev, FREENECT_FORMAT_RGB);
	freenect_set_depth_format(f_dev, FREENECT_FORMAT_11_BIT);
 
	freenect_start_depth(f_dev);
	freenect_start_video(f_dev);
	res = pthread_create(&fnkt_thread, NULL, freenect_threadfunc, NULL);
	if (res) {
		printf("pthread_create failed\n");
		return 1;
	}
 	unsigned fr=0;
	cv::Mat depthf;
	while (!die) {
		fr++;
 
		cv::imshow("rgb", rgbMat);
		depthMat.convertTo(depthf, CV_8UC1, 255.0/2048.0);
		cv::imshow("depth",depthf);
 
		char k = cvWaitKey(5);
		if( k == 27 ) break;
	}
	die = 1;
 
	printf("-- done!\n");
 
	cvDestroyWindow("rgb");
	cvDestroyWindow("depth");
 
	pthread_join(fnkt_thread, NULL);
	pthread_exit(NULL);
}






























