#include "main.h"
#include "ObjGL.h"
#include "FluidMat.h"
#include <getopt.h>
#include <fstream>
#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/OpenGL.h>
#else
#include <GL/gl.h>
#include <GL/glut.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

UI g_UI;
ObjGL obj;
FluidMat fluidmat;

static std::string g_output_dir_name;
static std::string g_import_ref_obj_path;
static std::string g_import_cal_obj_path;

static Eigen::Vector3d input_eyepos;
static Eigen::Vector3d input_eulerangle;
static Eigen::Quaterniond input_Quot;
static double input_perspective;
static Eigen::Vector2i input_windowsize;

static unsigned g_save_number_width;
static int g_output_frame;

static bool parseCommandLineOptions( int* argc, char*** argv )
{
	const struct option long_options[] =
	{
		{ "output_dir", required_argument, nullptr, 'o' },
		{ nullptr, 0, nullptr, 0 }
	};
	
	while( true )
	{
		int option_index = 0;
		const int c{ getopt_long( *argc, *argv, "a:b:c:d:e:f:g:h:", long_options, &option_index ) };
		if( c == -1 )
		{
			break;
		}
		switch( c )
		{
			case 'a':
			{
				std::cerr << "load out file dir " << optarg <<  std::endl;
				g_output_dir_name = optarg;
				break;
			}
			case 'b':
			{
				std::cerr << "load ref obj path " << optarg <<  std::endl;
				g_import_ref_obj_path = optarg;
				break;
			}
			case 'c':
			{
				std::cerr << "load cal obj path " << optarg <<  std::endl;
				g_import_cal_obj_path = optarg;
				break;
			}
			case 'd':
			{
				std::cerr << "load time frame " << optarg <<  std::endl;
				g_output_frame = std::atoi(optarg);
				break;
			}
			case 'e':
			{
				std::cerr << "load eyepos (by blender) " << optarg <<  std::endl;
				std::vector<std::string> v;
				std::string s = optarg;
				std::stringstream ss{s};
				std::string buf;
				while (std::getline(ss, buf, ',')) {
					v.push_back(buf);
				}
				input_eyepos.x() = std::stod(v[0]);
				input_eyepos.y() = std::stod(v[1]);
				input_eyepos.z() = std::stod(v[2]);
				break;
			}
			case 'f':
			{
				std::cerr << "load Quoternion (by blender) " << optarg <<  std::endl;
				std::vector<std::string> v;
				std::string s = optarg;
				std::stringstream ss{s};
				std::string buf;
				while (std::getline(ss, buf, ',')) {
					v.push_back(buf);
				}
				input_Quot.w() = std::stod(v[0]);
				input_Quot.x() = std::stod(v[1]);
				input_Quot.y() = std::stod(v[2]);
				input_Quot.z() = std::stod(v[3]);
				break;
			}
			case 'g':
			{
				std::cerr << "load perspective " << optarg <<  std::endl;
				
				input_perspective = std::stod(optarg);
				break;
			}
			case 'h':
			{
				std::cerr << "load window size  " << optarg <<  std::endl;
				std::vector<std::string> v;
				std::string s = optarg;
				std::stringstream ss{s};
				std::string buf;
				while (std::getline(ss, buf, ',')) {
					v.push_back(buf);
				}
				input_windowsize.x() = stoi(v[0]);
				input_windowsize.y() = stoi(v[1]);
				break;
			}
			default:
			{
				std::cerr << "This is a bug in the command line parser. Please file a report." << std::endl;
				return false;
			}
		}
	}
	
	return true;
}
static std::string generatePNGFileName( const std::string& prefix, const std::string& extension )
{
	std::stringstream ss;
	if( !g_output_dir_name.empty() )
	{
		ss << g_output_dir_name << "/";
	}
	ss << prefix << "_" << std::setfill('0') << std::right << std::setw(2) << g_output_frame << "." << extension;
	return ss.str();
}

static std::string generateOutputDiffFileName( const std::string& prefix, const std::string& extension )
{
	std::stringstream ss;
	if( !g_output_dir_name.empty() )
	{
		ss << g_output_dir_name << "/";
	}
	ss << prefix << "." << extension;
	return ss.str();
}

void saveImage_cal( const unsigned int imageWidth, const unsigned int imageHeight )
{
	fluidmat.fluid_mat_cal.release();
	glReadBuffer( GL_FRONT );
	
	int width = g_UI.window_width;
	int height = g_UI.window_height;
	
	GLubyte *data = new GLubyte[3*width*height];

	cv::Mat tempImage;
	cv::Mat tempImage2;
	
	glReadPixels(0,0, width, height, GL_RGB, GL_UNSIGNED_BYTE, data);
	fluidmat.fluid_mat_cal = cv::Mat(height, width, CV_8UC3, data);

	cv::cvtColor(fluidmat.fluid_mat_cal, tempImage, cv::COLOR_BGRA2GRAY);
	
	cv::threshold(tempImage, tempImage2, 160, 255, cv::THRESH_BINARY); //閾値160で2値画像に変換
	cv::flip(tempImage2, fluidmat.fluid_mat_cal, 0);
	tempImage.release();

	delete[] data;
	
}

void InputImage()
{
	cv::Mat	image = cv::imread(g_import_ref_obj_path);;

	if (image.empty() == true) {
		exit(1);
	}
	
	cv::Mat tempImage;
	cv::Mat tempImage2;
	
	cv::resize(image, tempImage2, cv::Size(input_windowsize.x(), input_windowsize.y()));
	cv::cvtColor(tempImage2, tempImage, cv::COLOR_BGRA2GRAY);
	
	cv::threshold(tempImage, fluidmat.fluid_mat_ref, 130, 255, cv::THRESH_BINARY); //閾値160で2値画像に変換
}

void RenderObj()
{
	obj.Draw();
	saveImage_cal(g_UI.window_width, g_UI.window_height);
}


void display(void)
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);

	
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3f( 1.0, 0.0, 0.0 );

	glViewport(0, 0, g_UI.window_width, g_UI.window_height);

	glMatrixMode(GL_PROJECTION);
	
	Eigen::Quaterniond Quot = input_Quot;
	Eigen::Vector3d eyelookATpre = Quot*Eigen::Vector3d{0.0, 0.0, -1.0};
	Eigen::Vector3d eyelookUPpre = Quot*Eigen::Vector3d{0.0, 1.0, 0.0};
	Eigen::Vector3d LA_Pos = eyelookATpre + input_eyepos;
	Eigen::Vector3d LU_Pos = eyelookUPpre;

	std::cout << "rotated x+: " << Quot * Eigen::Vector3d{ 1.0, 0.0, 0.0 } << std::endl;
    std::cout << "rotated y+: " << Quot * Eigen::Vector3d{ 0.0, 1.0, 0.0 } << std::endl;
    std::cout << "rotated z+: " << Quot * Eigen::Vector3d{ 0.0, 0.0, 1.0 } << std::endl;
	
	glLoadIdentity();
	gluPerspective(input_perspective, (double)g_UI.window_width / (double)g_UI.window_height, 0.1, 1000.0);
	gluLookAt(input_eyepos.x(), input_eyepos.y(), input_eyepos.z(), LA_Pos.x(), LA_Pos.y(), LA_Pos.z(), LU_Pos.x(), LU_Pos.y(), LU_Pos.z());

	std::cout << "syslookAtpre" << std::endl;
	std::cout << eyelookATpre << std::endl;
	std::cout << "LOOKAT" << std::endl;
	std::cout << LA_Pos << std::endl;
	std::cout << "LOOKUP" << std::endl;
	std::cout << LU_Pos << std::endl;
	std::cout << "EyePos" << std::endl;
	std::cout << input_eyepos << std::endl;
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	RenderObj();
	InputImage();

	cv::Mat temp1;
	cv::Mat temp2;
	cv::Mat temp3;
	cv::bitwise_not(fluidmat.fluid_mat_cal,temp1);
	cv::bitwise_not(fluidmat.fluid_mat_ref,temp2);
	cv::bitwise_xor(temp1, temp2, temp3);
	
	cv::Mat rfddiv = temp3;
	int width = rfddiv.cols;
	int height = rfddiv.rows;
	int channels = rfddiv.channels();
	int sum = 0;
	int count = 0;
	for(int j=0; j<height; j++)
	{
		int step = j*width;
		for(int i=0; i<width; i++)
		{
			int elm = i*rfddiv.elemSize();
			for(int c=0; c<channels; c++)
			{
				count ++;
				sum += rfddiv.data[step + elm + c];
			}
		}
	}
	const double sumdiv = (double)sum/255;
	std::cout << "sum point = " << sumdiv << std::endl;

	std::string output_diff_file_name{ generateOutputDiffFileName("diffvalues", "dat") };
	std::string output_snap_cal_name_obj{ generatePNGFileName( "snapcal", "png" ) };
	std::string output_snap_ref_name_obj{ generatePNGFileName( "snapref", "png" ) };
	std::string output_snap_diff_name_obj{ generatePNGFileName( "snapdiff", "png" ) };
	
	cv::imwrite(output_snap_ref_name_obj, fluidmat.fluid_mat_ref);
	cv::imwrite(output_snap_cal_name_obj, fluidmat.fluid_mat_cal);
	cv::imwrite(output_snap_diff_name_obj, rfddiv);

	try
	{
		std::cout << "writting diff value" << std::endl;
		double diffvalue = (double)sumdiv/(width*height);
		
		std::cout << "sumdiv = " << sumdiv << std::endl;
		std::cout << "diff width by height = " << width*height << std::endl;
		std::cout << "diff value is " << diffvalue << std::endl;
		
		std::ofstream ofs(output_diff_file_name, std::ios_base::out | std::ios_base::app);
		std::cout << output_diff_file_name << std::endl;
		if (!ofs)
		{
			std::cout << "CANT OPEN THE FILE" << std::endl;
			EXIT_FAILURE;
		}
		ofs << diffvalue << std::endl;
	}
	catch( const std::string& error )
	{
		std::cerr << error << std::endl;
		EXIT_FAILURE;
	}
		
	glFlush();
	
	exit(0);
}

void resize(int w, int h)
{

}

void init(void)
{
	glClearColor(0.0, 0.0, 1.0, 1.0);
}

void sequence(int argc, char* argv[])
{
	g_UI.window_width = input_windowsize.x();
	g_UI.window_height = input_windowsize.y();

	std::cout << "input_window_size: " << input_windowsize.x() << ", " << input_windowsize.y() << std::endl;
	
	glutInit(&argc, argv);
	glutInitWindowPosition(0, 0);
  	glutInitWindowSize(g_UI.window_width,g_UI.window_height);
	glutInitDisplayMode(GLUT_RGBA);

	glutCreateWindow("PantaRhei: (Toward) A Universal Rheology Simulator");
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	
	glutMainLoop();
}

int main(int argc, char *argv[])
{
	if( !parseCommandLineOptions( &argc, &argv ) )
	{
		std::cerr << "Failed to Parse Command " << argv[1] << std::endl;
		exit(-1);
	}
	
	obj.Load(g_import_cal_obj_path.c_str());
	sequence(argc, argv);
	return 0;
}