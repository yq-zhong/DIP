#ifndef EXPCODE_H
#define EXPCODE_H

#define FRAME_HEIGHT 800
#define FRAME_WIDTH 800

int midptvalue(int** pixelmat, int mheight, int mwidth);
int maxvalue(int** pixelmat, int mheight, int mwidth);
int minvalue(int** pixelmat, int mheight, int mwidth);
float avgvalue(int** pixelmat, int mheight, int mwidth);
float varvalue(int** pixelmat, int mheight, int mwidth);
int* histogram(int** pixelmat, int mheight, int mwidth);
int** binaryimg(int** pixelmat, int mheight, int mwidth);
int** histogramequ(int** pixelmat, int mheight, int mwidth);
int** graystretch(int** pixelmat, int mheight, int mwidth);
int** centralize(int** framemat, int** pixelmat, int mheight, int mwidth);
int** rotation(int** framemat, int** pixelmat, int mheight, int mwidth);
int** moveimage(int** framemat, int** pixelmat, int mheight, int mwidth);
int** scaling(int** framemat, int** pixelmat, int mheight, int mwidth);
int** DFT(int** pixelmat, int mheight, int mwidth);
int** DCT(int** pixelmat, int mheight, int mwidth);
int** walsh(int** pixelmat, int mheight, int mwidth);
int** haar(int** pixelmat, int mheight, int mwidth);
int** medianfit(int** pixelmat, int mheight, int mwidth);
int** averagefit(int** pixelmat, int mheight, int mwidth);
int** lowpassfit(int** pixelmat, int mheight, int mwidth);
int** sobel(int** pixelmat, int mheight, int mwidth);
int** laplace(int** pixelmat, int mheight, int mwidth);
int** highpassfit(int** pixelmat, int mheight, int mwidth);
int** randomnoise(int** pixelmat, int mheight, int mwidth);
int** impulsenoise(int** pixelmat, int mheight, int mwidth);
int** inversefit(int** pixelmat, int mheight, int mwidth);
int** wienerfit(int** pixelmat, int mheight, int mwidth);
int** jpeg(int** pixelmat, int mheight, int mwidth);


//****************************** 添加的部分 ******************************
struct Coordinates
{
	double x;
	double y;
};
struct Complex
{
	double re;
	double im;
	double r;
	double theta;
};

using namespace std;
const double PI = 3.1415926;

template<typename T>
void mynew(T** &mat, int height, int width);
template<typename T>
void mydelete(T** &mat, int height);
template<typename T>
void MatricMulti(T** &Matri_L, T** &Matri_R, T** &result, const int A, const int B, const int C);
template<>
void MatricMulti(Complex** &Matri_L, Complex** &Matri_R, Complex** &result, const int A, const int B, const int C);
template<typename T>
void SquarMatrixTrans(T** &Matrix, int N);
double haar_function(const int k, const double z, const int N);
double gaussrand();
void haar_matrix(double **&haarMatrix, int N);
Complex** MyDFT(Complex** Mat, int mheight, int mwidth);
void pait_fram_edge(int** &framemat);
void walsh_matrix(double **&walshMatrx, int N, int n);
//****************************** 添加的部分 ******************************


#endif // EXPCODE_H
