#include "Application.h"
#include "qt_opengl_framework.h"
#include <cmath>
#include <map>
#include <ctime>
#include <cstdlib>

Application::Application()
{
	srand(time(NULL));
}
Application::~Application()
{

}

//****************************************************************************
//
// * ��l�e���A�����Ntust.png����
// 
//============================================================================
void Application::createScene( void )
{
	
	ui_instance = Qt_Opengl_Framework::getInstance();
	
}

//****************************************************************************
//
// * ���}���w����
// 
//============================================================================
void Application::openImage( QString filePath )
{
	mImageSrc.load(filePath);
	mImageDst.load(filePath);

	renew();

	img_data = mImageSrc.bits();
	img_width = mImageSrc.width();
	img_height = mImageSrc.height();

	ui_instance->ui.label->setFixedHeight(img_height);
	ui_instance->ui.label->setFixedWidth(img_width);
}
//****************************************************************************
//
// * ��s�e��
// 
//============================================================================
void Application::renew()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageDst));

	std::cout << "Renew" << std::endl;
}

//****************************************************************************
//
// * �e����l��
// 
//============================================================================
void Application::reload()
{
	ui_instance = Qt_Opengl_Framework::getInstance();

	ui_instance->ui.label->clear();
	ui_instance->ui.label->setPixmap(QPixmap::fromImage(mImageSrc));
}

//****************************************************************************
//
// * �x�s����
// 
//============================================================================
void Application::saveImage(QString filePath )
{
	mImageDst.save(filePath);
}

//****************************************************************************
//
// * �N���ɸ���ഫ��RGB��m���
// 
//============================================================================
unsigned char* Application::To_RGB( void )
{
	unsigned char *rgb = new unsigned char[img_width * img_height * 3];
	int i, j;

	if (! img_data )
		return NULL;

	// Divide out the alpha
	for (i = 0; i < img_height; i++)
	{
		int in_offset = i * img_width * 4;
		int out_offset = i * img_width * 3;

		for (j = 0 ; j < img_width ; j++)
		{
			RGBA_To_RGB(img_data + (in_offset + j*4), rgb + (out_offset + j*3));
		}
	}

	return rgb;
}

void Application::RGBA_To_RGB( unsigned char *rgba, unsigned char *rgb )
{
	const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

	unsigned char  alpha = rgba[3];

	if (alpha == 0)
	{
		rgb[0] = BACKGROUND[0];
		rgb[1] = BACKGROUND[1];
		rgb[2] = BACKGROUND[2];
	}
	else
	{
		float	alpha_scale = (float)255 / (float)alpha;
		int	val;
		int	i;

		for (i = 0 ; i < 3 ; i++)
		{
			val = (int)floor(rgba[i] * alpha_scale);
			if (val < 0)
				rgb[i] = 0;
			else if (val > 255)
				rgb[i] = 255;
			else
				rgb[i] = val;
		}
	}
}
//------------------------Color------------------------

unsigned char Application::toGray(unsigned char* rgb, int offset) {
	return 0.299 * rgb[offset + rr] + 0.587 * rgb[offset + gg] + 0.114 * rgb[offset + bb];
}

///////////////////////////////////////////////////////////////////////////////
//
//  Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Gray()
{
	unsigned char *rgb = To_RGB();

	for (int i=0; i<img_height; i++)
	{
		for (int j=0; j<img_width; j++)
		{
			int offset_rgb = i*img_width*3+j*3;
			int offset_rgba = i*img_width*4+j*4;
			unsigned char gray = toGray(rgb, offset_rgb);

			for (int k=0; k<3; k++)
				img_data[offset_rgba+k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}
	
	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Uniform()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			int level[] = { 4, 8, 8 };

			for (int k = 0; k < 3; k++) {
				int q = 256 / level[k];
				img_data[offset_rgba + k] = (rgb[offset_rgb + k] / q) * q;
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Quant_Populosity()
{
	unsigned char *rgb = this->To_RGB();

	vector<pair<int, int>> list = vector<pair<int, int>>(65536, pair<int, int>(0, 0));

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;

			int color = 0;

			for (int k = 0; k < 3; k++) {
				color <<= 5;
				color += ((rgb[offset_rgb + k] & 0b11111000) >> 3);
			}

			list[color].first = color;
			list[color].second++;
		}
	}

	sort(list.begin(), list.begin() + 32768, [](const pair<int, int>& a, const pair<int, int>& b) {
		return a.second > b.second;
	});

	int selected[256];

	for (int i = 0; i < 256; i++) {
		uchar b = (list[i].first >> 10) << 3;
		uchar g = ((list[i].first >> 2) & 0b11111000);
		uchar r = (list[i].first & 0b11111) << 3;

		selected[i] = (b << 16) + (g << 8) + r;
	}

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;

			int mDist = INT_MAX;

			int c = ((rgb[offset_rgb + bb] & 0b11111000) << 8) + ((rgb[offset_rgb + gg] & 0b11111100) << 3) + ((rgb[offset_rgb + rr] & 0b11111000) >> 3);

			if (c >= list.size()) {
				continue;
			}

			if (list[c].second != -1) {

				int mIdx = 0;

				for (int k = 0; k < 256; k++) {
					uchar b = selected[k] >> 16;
					uchar g = (selected[k] >> 8) & 255;
					uchar r = selected[k] & 255;

					int dist = (pow(b - rgb[offset_rgb + bb], 2) + pow(g - rgb[offset_rgb + gg], 2) + pow(r - rgb[offset_rgb + rr], 2));
					
					if (dist < mDist) {
						mDist = dist;
						mIdx = k;
					}
				}

				list[c].first = selected[mIdx];

				list[c].second = -1;
			}

			img_data[offset_rgba] = list[c].first >> 16;
			img_data[offset_rgba + 1] = (list[c].first >> 8) & 255;
			img_data[offset_rgba + 2] = list[c].first & 255;

			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Dithering------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Threshold()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			unsigned char gray = toGray(rgb, offset_rgb) >= 128 ? WHITE : BLACK;
			
			for (int k = 0; k < 3; k++)
				img_data[offset_rgba + k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Random()
{
	unsigned char *rgb = this->To_RGB();

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			unsigned char gray = (toGray(rgb, offset_rgb) + rand() % 102 - 51) >= 128 ? WHITE : BLACK;
			for (int k = 0; k < 3; k++)
				img_data[offset_rgba + k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_FS()
{
	unsigned char *rgb = this->To_RGB();

	vector<double> gray;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			gray.push_back(toGray(rgb, offset_rgb));
		}
	}
	double distribution[] = { 7.0 / 16.0, 3.0 / 16.0, 5.0 / 16.0, 1.0 / 16.0 };
	for (int i = 0; i < img_height; i++)
	{
		int x[] = { 1, -1, 0, 1 };
		int y[] = { 0, 1, 1, 1 };
		if (i % 2) {
			for (int j = 0; j < 4; j++) {
				x[j] *= -1;
			}
		}
		for (int j = ((i % 2) ? img_width - 1 : 0); j >= 0 && j < img_width; j += ((i % 2) ? 1 : 1))
		{
			int offset_gray = i * img_width + j;
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			unsigned char value = gray[offset_gray] >= 128 ? WHITE : BLACK;
			double error = gray[offset_gray] - value;
			for (int k = 0; k < 3; k++)
				img_data[offset_rgba + k] = value;
			img_data[offset_rgba + aa] = WHITE;
			for (int k = 0; k < 4; k++) {
				if (i + y[k] >= 0 && i + y[k] < img_height && j + x[k] >= 0 && j + x[k] < img_width) {
					gray[(i + y[k]) * img_width + j + x[k]] += error * distribution[k];
				}
			}
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Bright()
{
	unsigned char *rgb = this->To_RGB();

	double avg = 0;
	vector<int> brightness;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			unsigned char gray = toGray(rgb, offset_rgb);
			brightness.push_back(gray);
			avg += gray;
		}
	}
	sort(brightness.begin(), brightness.end());
	avg /= img_height * img_width;
	int threshold = brightness[(int)((1 - avg / 255) * (img_height * img_width - 1))];
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			unsigned char gray = toGray(rgb, offset_rgb) >= threshold ? WHITE : BLACK;
			
			for (int k = 0; k < 3; k++)
				img_data[offset_rgba + k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Cluster()
{
	unsigned char *rgb = this->To_RGB();

	double mat[4][4] = {
		{0.7059, 0.3529, 0.5882, 0.2353},
		{0.0588, 0.9412, 0.8235, 0.4118},
		{0.4706, 0.7647, 0.8824, 0.1176},
		{0.1765, 0.5294, 0.2941, 0.6471}
	};

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			unsigned char gray = ((double)toGray(rgb, offset_rgb) / 255) < mat[i % 4][j % 4] ? BLACK : WHITE;
			
			for (int k = 0; k < 3; k++)
				img_data[offset_rgba + k] = gray;
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Dither_Color()
{
	unsigned char *rgb = this->To_RGB();

	vector<double> RGB;
	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			for (int k = 0; k < 3; k++) {
				RGB.push_back(rgb[offset_rgb + k]);
			}
		}
	}

	double distribution[] = { 7.0 / 16.0, 3.0 / 16.0, 5.0 / 16.0, 1.0 / 16.0 };
	int level[] = { 4, 8, 8 };
	for (int i = 0; i < img_height; i++)
	{
		int x[] = { 1, -1, 0, 1 };
		int y[] = { 0, 1, 1, 1 };
		if (i % 2) {
			for (int j = 0; j < 4; j++) {
				x[j] *= -1;
			}
		}
		for (int j = ((i % 2) ? img_width - 1 : 0); j >= 0 && j < img_width; j += ((i % 2) ? 1 : 1))
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			unsigned char color[3];
			double error[3];
			for (int k = 0; k < 3; k++) {
				int q = 256 / level[k];
				color[k] = ((unsigned char)RGB[offset_rgb + k] / q) * q;
				error[k] = RGB[offset_rgb + k] - color[k];
				img_data[offset_rgba + k] = color[k];
			}
			img_data[offset_rgba + aa] = WHITE;
			for (int k = 0; k < 4; k++) {
				if (i + y[k] >= 0 && i + y[k] < img_height && j + x[k] >= 0 && j + x[k] < img_width) {
					for (int l = 0; l < 3; l++) {
						RGB[(i + y[k]) * img_width * 3 + (j + x[k]) * 3 + l] += error[l] * distribution[k];
					}
				}
			}
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Filter------------------------

///////////////////////////////////////////////////////////////////////////////
//
//     Filtering the img_data array by the filter from the parameters
//
///////////////////////////////////////////////////////////////////////////////
void Application::filtering( double filter[][5] )
{
	unsigned char *rgb = this->To_RGB();

	double weight = 0;
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 5; j++) {
			weight += filter[i][j];
		}
	}

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			for (int k = 0; k < 3; k++) {
				double sum = 0;
				for (int m = -2; m <= 2; m++) {
					for (int n = -2; n <= 2; n++) {
						if (i + m >= 0 && i + m < img_height && j + n >= 0 && j + n < img_width) {
							sum += rgb[(i + m) * img_width * 3 + (j + n) * 3 + k] * filter[m + 2][n + 2];
						}
					}
				}
				img_data[offset_rgba + k] = sum / weight;
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::filtering( double **filter, int N )
{
	unsigned char *rgb = this->To_RGB();

	double weight = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			weight += filter[i][j];
		}
	}

	for (int i = 0; i < img_height; i++)
	{
		for (int j = 0; j < img_width; j++)
		{
			int offset_rgb = i * img_width * 3 + j * 3;
			int offset_rgba = i * img_width * 4 + j * 4;
			for (int k = 0; k < 3; k++) {
				double sum = 0;
				for (int m = -N / 2; m <= N / 2; m++) {
					for (int n = -N / 2; n <= N / 2; n++) {
						if (i + m >= 0 && i + m < img_height && j + n >= 0 && j + n < img_width) {
							sum += rgb[(i + m) * img_width * 3 + (j + n) * 3 + k] * filter[m + N / 2][n + N / 2];
						}
					}
				}
				img_data[offset_rgba + k] = sum / weight;
			}
			img_data[offset_rgba + aa] = WHITE;
		}
	}

	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Box()
{
	double filter[5][5] = {
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1}
	};
	filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Bartlett()
{
	double filter[5][5] = {
		{1, 2, 3, 2, 1},
		{2, 4, 6, 4, 2},
		{3, 6, 9, 6, 3},
		{2, 4, 6, 4, 2},
		{1, 2, 3, 2, 1},
	};
	filtering(filter);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian()
{
	Filter_Gaussian_N(5);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Gaussian_N( unsigned int N )
{
	double **filter = new double*[N];
	for (int i = 0; i < N; i++) {
		filter[i] = new double[N];
	}
	for (int i = 0; i < N / 2; i++) {
		for (int j = 0; j < N / 2; j++) {
			filter[i][j] = filter[i][N - j - 1] = filter[N - i - 1][j] = filter[N - i - 1][N - j - 1] = exp(-(pow(i - N / 2, 2) + pow(j - N / 2, 2)));
		}
	}
	filtering(filter, N);
}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Edge()
{

}
///////////////////////////////////////////////////////////////////////////////
//
//  Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Filter_Enhance()
{
	unsigned char *rgb = this->To_RGB();



	delete[] rgb;
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Size------------------------

///////////////////////////////////////////////////////////////////////////////
//
//  Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Half_Size()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Double_Size()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}
///////////////////////////////////////////////////////////////////////////////
//
//  resample_src for resize and rotate
//
///////////////////////////////////////////////////////////////////////////////
void Application::resample_src(int u, int v, float ww, unsigned char* rgba)
{

}

///////////////////////////////////////////////////////////////////////////////
//
//  Scale the image dimensions by the given factor.  The given factor is 
//	assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Resize( float scale )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//////////////////////////////////////////////////////////////////////////////
//
//  Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Rotate( float angleDegrees )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

//------------------------Composing------------------------


void Application::loadSecondaryImge( QString filePath )
{
	mImageSrcSecond.load(filePath);

	renew();

	img_data2 = mImageSrcSecond.bits();
	img_width2 = mImageSrcSecond.width();
	img_height2 = mImageSrcSecond.height();
}

//////////////////////////////////////////////////////////////////////////
//
//	Composite the image A and image B by Over, In, Out, Xor and Atom. 
//
//////////////////////////////////////////////////////////////////////////
void Application::Comp_image( int tMethod )
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Over()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_In()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Out()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Atop()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::Comp_Xor()
{
	if (img_height == img_height2 && img_width == img_width2)
	{

	}
	else
	{
		std::cout << "Images not the same size" << std::endl;
	}
}

//------------------------NPR------------------------

///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void Application::NPR_Paint()
{
	mImageDst = QImage(img_data, img_width, img_height, QImage::Format_ARGB32 );
	renew();
}

void Application::NPR_Paint_Layer( unsigned char *tCanvas, unsigned char *tReferenceImage, int tBrushSize )
{

}

///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void Application::Paint_Stroke( const Stroke& s )
{
	int radius_squared = (int)s.radius * (int)s.radius;
	for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) 
	{
		for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) 
		{
			int x_loc = (int)s.x + x_off;
			int y_loc = (int)s.y + y_off;

			// are we inside the circle, and inside the image?
			if ((x_loc >= 0 && x_loc < img_width && y_loc >= 0 && y_loc < img_height)) 
			{
				int dist_squared = x_off * x_off + y_off * y_off;
				int offset_rgba = (y_loc * img_width + x_loc) * 4;

				if (dist_squared <= radius_squared) 
				{
					img_data[offset_rgba + rr] = s.r;
					img_data[offset_rgba + gg] = s.g;
					img_data[offset_rgba + bb] = s.b;
					img_data[offset_rgba + aa] = s.a;
				} 
				else if (dist_squared == radius_squared + 1) 
				{
					img_data[offset_rgba + rr] = (img_data[offset_rgba + rr] + s.r) / 2;
					img_data[offset_rgba + gg] = (img_data[offset_rgba + gg] + s.g) / 2;
					img_data[offset_rgba + bb] = (img_data[offset_rgba + bb] + s.b) / 2;
					img_data[offset_rgba + aa] = (img_data[offset_rgba + aa] + s.a) / 2;
				}
			}
		}
	}
}





///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
	unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}



