#include "expcode.h"
#include <iostream>
#include "algorithm"


//示例: 求图像中心点像素的灰度值
int midptvalue(int** pixelmat, int mheight, int mwidth)
{
	//pixelmat为指向图像像素值二维数组指针, 指向图像左上角的像素, 像素值类型为int;
	//mheight为图像的高度(行), mwidth为图像的宽度(列);
	int middlerow = mheight / 2 - 1;
	int middlecol = mwidth / 2 - 1;
	return pixelmat[middlerow][middlecol];
}

//求最大值
int maxvalue(int** pixelmat, int mheight, int mwidth)
{
	int MAX = 0;//最大值
	int i_row = 0;
	int i_col = 0;
	for (i_row = 0; i_row < mheight; i_row ++)
	{
		for (i_col = 0; i_col < mwidth; i_col++)
		{
			if (pixelmat[i_row][i_col] > MAX)
			{
				MAX = pixelmat[i_row][i_col];
			}
		}
	}
    return MAX;
}

//求最小值
int minvalue(int** pixelmat, int mheight, int mwidth)
{
	int MIN = 255;
	int i_row = 0;
	int i_col = 0;
	for (i_row = 0; i_row < mheight; i_row++)
	{
		for (i_col = 0; i_col < mwidth; i_col++)
		{
			if (pixelmat[i_row][i_col] < MIN)
			{
				MIN = pixelmat[i_row][i_col];
			}
		}
	}
    return MIN;
}

//求平均值
float avgvalue(int** pixelmat, int mheight, int mwidth)
{
	int SUM = 0;
	float AVG = 0;
	int i_row = 0;
	int i_col = 0;
	for (i_row = 0; i_row < mheight; i_row++)
	{
		for (i_col = 0; i_col < mwidth; i_col++)
		{
			SUM += pixelmat[i_row][i_col];
		}
	}
	if (mheight && mwidth)
	{
		AVG = (float)SUM / (mheight * mwidth);
	}
    return AVG;
}

//求方差
float varvalue(int** pixelmat, int mheight, int mwidth)
{
	float VAR = 0;
	float AVG = 0;
	float VAR_SUM = 0;
	int i_row = 0;
	int i_col = 0;
	float diff = 0;
	AVG = avgvalue(pixelmat, mheight, mwidth);
	for (i_row = 0; i_row < mheight; i_row++)
	{
		for (i_col = 0; i_col < mwidth; i_col++)
		{
			diff = pixelmat[i_row][i_col] - AVG;
			VAR_SUM += diff * diff;
		}
	}
	if (mheight && mwidth)
	{
		VAR = VAR_SUM / (mheight * mwidth);
	}
    return VAR;
}

//统计直方图, 返回长度为256的1维数组
int* histogram(int** pixelmat, int mheight, int mwidth)
{
	//注意:函数内分配数组必须使用动态分配;
	int* HISTO = new int[256];
	int i_row = 0;
	int i_col = 0;
	int gray = 0;
	memset(HISTO, 0, sizeof(int) * 256);
	for (i_row = 0; i_row < mheight; i_row++)
	{
		for (i_col = 0; i_col < mwidth; i_col++)
		{
			gray = pixelmat[i_row][i_col];
			HISTO[gray] ++;
		}
	}
	return HISTO;
}







//示例,将灰度图转化为二值图像,返回处理后的图像
int** binaryimg(int** pixelmat, int mheight, int mwidth)
{
	for(int i = 0; i < mheight; i++)
	{
		for(int j = 0; j < mwidth; j++)
		{
			//从左上角开始遍历整幅图像, 实现二值化;
			pixelmat[i][j] = pixelmat[i][j] > 128 ? 255 : 0;
		}
	}
	//根据实验要求返回对应的值;
	return pixelmat;
}

//直方图均衡, 返回处理后的图像
int** histogramequ(int** pixelmat, int mheight, int mwidth)
{
	int* hist;//直方图
	hist = histogram(pixelmat, mheight, mwidth);
	int map[256];//灰度映射关系
	int size = mheight * mwidth;
	int sum = 0;//比例累和
	int i = 0, j = 0;
	int gray = 0;
	for (i = 0; i < 256; i++)
	{
		sum += hist[i];
		map[i] = sum * 255 / size;
	}
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			gray = pixelmat[i][j];
			pixelmat[i][j] = map[gray];
		}
	}
	delete[] hist;
    return pixelmat;
}

//灰度拉伸, 返回处理后的图像
int** graystretch(int** pixelmat, int mheight, int mwidth)
{
	int MAX = maxvalue(pixelmat, mheight, mwidth);
	if (MAX == 0)
	{
		return pixelmat;
	}
	int i = 0, j = 0;
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			pixelmat[i][j] = pixelmat[i][j] * 100 / MAX;
		}
	}
    return pixelmat;
}

//中值滤波, 返回处理后的图像
int** medianfit(int** pixelmat, int mheight, int mwidth)
{
	int** mat_temp;
	mynew(mat_temp, mheight, mwidth);
	int mask[9];
	int i = 0, j = 0;
	for (i = mheight - 2; i > 0; i--)
	{
		for (j = mwidth - 2; j > 0; j--)
		{
			mask[0] = pixelmat[i - 1][j - 1];
			mask[1] = pixelmat[i - 1][j];
			mask[2] = pixelmat[i - 1][j + 1];
			mask[3] = pixelmat[i][j - 1];
			mask[4] = pixelmat[i][j];
			mask[5] = pixelmat[i][j + 1];
			mask[6] = pixelmat[i + 1][j - 1];
			mask[7] = pixelmat[i + 1][j];
			mask[8] = pixelmat[i + 1][j + 1];
			sort(mask, mask + 9);
			mat_temp[i][j] = mask[4];
		}
	}
	for (i = mheight - 2; i > 0; i--)
	{
		for (j = mwidth - 2; j > 0; j--)
		{
			pixelmat[i][j] = mat_temp[i][j];
		}
	}
	mydelete(mat_temp, mheight);
    return pixelmat;
}

//均值滤波, 返回处理后的图像
int** averagefit(int** pixelmat, int mheight, int mwidth)
{
	int** mat_temp;
	mynew(mat_temp, mheight, mwidth);
	int i = 0, j = 0;
	int sum;
	for (i = mheight - 2; i > 0; i--)
	{
		for (j = mwidth - 2; j > 0; j--)
		{
			sum = 0;
			sum += pixelmat[i - 1][j - 1];
			sum += pixelmat[i - 1][j];
			sum += pixelmat[i - 1][j + 1];
			sum += pixelmat[i][j - 1];
			sum += pixelmat[i][j];
			sum += pixelmat[i][j + 1];
			sum += pixelmat[i + 1][j - 1];
			sum += pixelmat[i + 1][j];
			sum += pixelmat[i + 1][j + 1];
			
			mat_temp[i][j] = sum / 9;
		}
	}
	for (i = mheight - 2; i > 0; i--)
	{
		for (j = mwidth - 2; j > 0; j--)
		{
			pixelmat[i][j] = mat_temp[i][j];
		}
	}
	mydelete(mat_temp, mheight);
	return pixelmat;
}

//理想低通滤波, 返回处理后的图像
int** lowpassfit(int** pixelmat, int mheight, int mwidth)
{
	int i, j, u, v;
	const int r = 50;//滤波半径

	//原图转化成复数存储
	Complex** pixel_cmlex;
	mynew(pixel_cmlex, mheight, mwidth);
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			pixel_cmlex[i][j] = { (double)(pixelmat[i][j]), 0, 0, 0 };
		}
	}

	//*********** 频谱图 ************
	Complex** MatF;//F(u,v)
	MatF = MyDFT(pixel_cmlex, mheight, mwidth);
	
	//低通滤波,并取共轭
	for (u = 0; u < mheight; u++)
	{
		for (v = 0; v < mwidth; v++)
		{
			if ((u > r && u < (mheight - r)) || (v > r && v < (mwidth - r)))
			{
				MatF[u][v] = { 0, 0, 0, 0 };
			}
			MatF[u][v].im *= -1;
		}
	}
	
	

	//*************** 频谱图的共轭做正变换，再取共轭 <==> 反变换 ****************
	pixel_cmlex = MyDFT(MatF, mheight, mwidth);
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			pixelmat[i][j] = pixel_cmlex[i][j].r;
			if (pixelmat[i][j] > 255)
			{
				pixelmat[i][j] = 255;
			}
		}
	}


	mydelete(MatF, mheight);
	mydelete(pixel_cmlex, mheight);
    return pixelmat;
}

//sobel算子, 返回处理后的图像
int** sobel(int** pixelmat, int mheight, int mwidth)
{
	int** mat_temp_x;
	int** mat_temp_y;
	mynew(mat_temp_x, mheight, mwidth);
	mynew(mat_temp_y, mheight, mwidth);
	int i = 0, j = 0;
	int sum;
	int sum_y;
	for (i = mheight - 2; i > 0; i--)
	{
		for (j = mwidth - 2; j > 0; j--)
		{
			sum = 0;
			sum += pixelmat[i - 1][j - 1];
			sum += 2 * pixelmat[i - 1][j];
			sum += pixelmat[i - 1][j + 1];
			sum -= pixelmat[i + 1][j - 1];
			sum -= 2 * pixelmat[i + 1][j];
			sum -= pixelmat[i + 1][j + 1];
			mat_temp_x[i][j] = sum;

			sum_y = 0;
			sum_y += pixelmat[i - 1][j - 1];
			sum_y += 2 * pixelmat[i][j - 1];
			sum_y += pixelmat[i + 1][j - 1];
			sum_y -= pixelmat[i - 1][j + 1];
			sum_y -= 2 * pixelmat[i][j + 1];
			sum_y -= pixelmat[i + 1][j + 1];
			mat_temp_y[i][j] = sum_y;
		}
	}
	for (i = mheight - 2; i > 0; i--)
	{
		for (j = mwidth - 2; j > 0; j--)
		{
			pixelmat[i][j] = (abs(mat_temp_x[i][j]) + abs(mat_temp_y[i][j]))/2;
		}
	}
	mydelete(mat_temp_x, mheight);
	mydelete(mat_temp_y, mheight);
	return pixelmat;
}

//laplace算子, 返回处理后的图像
int** laplace(int** pixelmat, int mheight, int mwidth)
{
	int** mat_temp;
	mynew(mat_temp, mheight, mwidth);
	int i = 0, j = 0;
	int sum;
	for (i = mheight - 2; i > 0; i--)
	{
		for (j = mwidth - 2; j > 0; j--)
		{
			sum = 0;
			sum += 8 * pixelmat[i][j];
			sum -= pixelmat[i - 1][j - 1];
			sum -= pixelmat[i - 1][j];
			sum -= pixelmat[i - 1][j + 1];
			sum -= pixelmat[i][j - 1];
			
			sum -= pixelmat[i][j + 1];
			sum -= pixelmat[i + 1][j - 1];
			sum -= pixelmat[i + 1][j];
			sum -= pixelmat[i + 1][j + 1];

			mat_temp[i][j] = (abs(sum) / 8) > 255 ? 255 : (abs(sum) / 8);
		}
	}
	for (i = mheight - 2; i > 0; i--)
	{
		for (j = mwidth - 2; j > 0; j--)
		{
			pixelmat[i][j] = abs(mat_temp[i][j]);
		}
	}
	mydelete(mat_temp, mheight);
	return pixelmat;
}

//理想高通滤波, 返回处理后的图像
int** highpassfit(int** pixelmat, int mheight, int mwidth)
{
	int i, j, u, v;
	const int r = 50;//滤波半径

					 //原图转化成复数存储
	Complex** pixel_cmlex;
	mynew(pixel_cmlex, mheight, mwidth);
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			pixel_cmlex[i][j] = { (double)(pixelmat[i][j]), 0, 0, 0 };
		}
	}

	//*********** 频谱图 ************
	Complex** MatF;//F(u,v)
	MatF = MyDFT(pixel_cmlex, mheight, mwidth);
	//低频中心移到图像中心
	Complex swap_temp;
	for (u = 0; u < mheight / 2; u++)
	{
		for (v = 0; v < mwidth / 2; v++)
		{
			swap_temp = MatF[u][v];
			MatF[u][v] = MatF[u + mheight / 2][v + mwidth / 2];
			MatF[u + mheight / 2][v + mwidth / 2] = swap_temp;

			swap_temp = MatF[u][v + mwidth / 2];
			MatF[u][v + mwidth / 2] = MatF[u + mheight / 2][v];
			MatF[u + mheight / 2][v] = swap_temp;
		}
	}
	//高通滤波,并取共轭
	for (u = 0; u < mheight; u++)
	{
		for (v = 0; v < mwidth; v++)
		{
			if (abs(u - mheight / 2) < r && abs(v - mwidth / 2) < r)
			{
				MatF[u][v] = { 0, 0, 0, 0 };
			}
			MatF[u][v].im *= -1;
		}
	}
	//
	for (u = 0; u < mheight / 2; u++)
	{
		for (v = 0; v < mwidth / 2; v++)
		{
			swap_temp = MatF[u][v];
			MatF[u][v] = MatF[u + mheight / 2][v + mwidth / 2];
			MatF[u + mheight / 2][v + mwidth / 2] = swap_temp;

			swap_temp = MatF[u][v + mwidth / 2];
			MatF[u][v + mwidth / 2] = MatF[u + mheight / 2][v];
			MatF[u + mheight / 2][v] = swap_temp;
		}
	}



	//*************** 频谱图的共轭做正变换，再取共轭 <==> 反变换 ****************
	pixel_cmlex = MyDFT(MatF, mheight, mwidth);
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			pixelmat[i][j] = pixel_cmlex[i][j].r;
			if (pixelmat[i][j] > 255)
			{
				pixelmat[i][j] = 255;
			}
		}
	}


	mydelete(MatF, mheight);
	mydelete(pixel_cmlex, mheight);
	return pixelmat;
}







//示例, 将图像平移到显示区域的中心
int** centralize(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	//framemat为指向显示区域(画板)的二维数组指针, 大小为FRAME_HEIGHT x FRAMEWIDTH = 800 x 800
	int xpt = (FRAME_HEIGHT - mheight) / 2;
	int ypt = (FRAME_WIDTH - mwidth) / 2;
	for(int i = 0; i < mheight; i++)
	{
		for(int j = 0; j < mwidth; j++)
		{
			framemat[i + xpt][j + ypt] = pixelmat[i][j];
		}
	}
	pait_fram_edge(framemat);

	return framemat;
}

//旋转图像, 返回显示区域(画板)指针
//画布中心视为原点,y轴向下为正
int** rotation(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	const double theta = PI / 9;//逆时针旋转theta
	const double COS_THETA = cos(theta);
	const double SIN_THETA = sin(theta);

	Coordinates Cord_LT, Cord_LB, Cord_RT, Cord_RB;//原图四个角旋转后在画布中的坐标
	Cord_LT.x = -1 * (mheight * SIN_THETA + mwidth * COS_THETA) / 2 + FRAME_WIDTH / 2;//-0.5mwidth*cos + -0.5mheight*sin
	Cord_LT.y = (mwidth * SIN_THETA - mheight * COS_THETA) / 2 + FRAME_HEIGHT / 2;//- -0.5mwidth*sin + -0.5mheight*cos
	Cord_LB.x = (mheight * SIN_THETA - mwidth * COS_THETA) / 2 + FRAME_WIDTH / 2;
	Cord_LB.y = (mwidth * SIN_THETA + mheight * COS_THETA) / 2 + FRAME_HEIGHT / 2;
	Cord_RT.x = (mwidth * COS_THETA - mheight * SIN_THETA) / 2 + FRAME_WIDTH / 2;
	Cord_RT.y = -1 *(mwidth * SIN_THETA + mheight * COS_THETA) / 2 + FRAME_HEIGHT / 2;
	Cord_RB.x = (mheight * SIN_THETA + mwidth * COS_THETA) / 2 + FRAME_WIDTH / 2;
	Cord_RB.y = (mheight * COS_THETA - mwidth * SIN_THETA) / 2 + FRAME_HEIGHT / 2;

	int X_TEMP[4] = { Cord_LT.x, Cord_LB.x, Cord_RT.x, Cord_RB.x };
	int Y_TEMP[4] = { Cord_LT.y, Cord_LB.y, Cord_RT.y, Cord_RB.y };
	sort(X_TEMP, X_TEMP + 4);
	sort(Y_TEMP, Y_TEMP + 4);
	const int MIN_X = X_TEMP[0], MAX_X = X_TEMP[3];//旋转后图像在画布中X坐标的最小值和最大值
	const int MIN_Y = Y_TEMP[0], MAX_Y = Y_TEMP[3];//旋转后图像在画布中Y坐标的最小值和最大值

	//直线方程 ax + by + c = 0
	double a_row, a_col;
	double b_row, b_col;
	double c_row, c_col;
	if (Cord_LT.y == Cord_RT.y)
	{
		a_row = 0; b_row = 1; c_row = Cord_LT.y - FRAME_HEIGHT / 2;
		a_col = 1; b_col = 0; c_col = Cord_LT.x - FRAME_WIDTH / 2;
	}
	else if (Cord_LT.x == Cord_RT.x)
	{
		a_row = 1; b_row = 0; c_row = Cord_LT.x - FRAME_WIDTH / 2;
		a_col = 0; b_col = 1; c_col = Cord_LT.y - FRAME_HEIGHT / 2;
	}
	else
	{
		a_row = ((double) (Cord_LT.y - Cord_RT.y)) / (Cord_RT.x - Cord_LT.x);
		b_row = 1;
		c_row = a_row * (FRAME_WIDTH / 2 - Cord_LT.x) + b_row * (FRAME_HEIGHT / 2 - Cord_LT.y);

		a_col = ((double) (Cord_LT.y - Cord_LB.y)) / (Cord_LB.x - Cord_LT.x);
		b_col = 1;
		c_col = a_col * (FRAME_WIDTH / 2 - Cord_LT.x) + b_col * (FRAME_HEIGHT / 2 - Cord_LT.y);
	}

	int i, j;
	int x, y;
	double x_ori, y_ori;
	int x_l, x_r;
	int y_u, y_d;
	double inter_l, inter_r;
	int pix_lu;
	int pix_ld;
	int pix_ru;
	int pix_rd;
	for (i = MIN_Y; i <= MAX_Y; i++)
	{
		y = i - FRAME_HEIGHT / 2;
		for (j = MIN_X; j <= MAX_X; j++)
		{
			x = j - FRAME_WIDTH / 2;
			if ((a_row * x + b_row * y + c_row) * (a_row * x + b_row * y - c_row) <= 0 &&
				(a_col * x + b_col * y + c_col) * (a_col * x + b_col * y - c_col) <= 0)
			{
				x_ori = x * COS_THETA - y * SIN_THETA + mwidth / 2;
				y_ori = x * SIN_THETA + y * COS_THETA + mheight / 2;
				x_l = x_ori; x_r = x_l + 1;//插值近邻的x坐标
				y_u = y_ori; y_d = y_u + 1;//插值近邻的y坐标
				if (x_r >= mwidth)
				{
					int temp = x_r - mwidth + 1;
					x_l -= temp;
					x_r -= temp;
					x_ori -= temp;
				}
				if (y_d >= mheight)
				{
					int temp = y_d - mheight + 1;
					y_u -= temp;
					y_d -= temp;
					y_ori -= temp;
				}
				//四个近邻像素值
				pix_lu = pixelmat[y_u][x_l];
				pix_ld = pixelmat[y_d][x_l];
				pix_ru = pixelmat[y_u][x_r];
				pix_rd = pixelmat[y_d][x_r];
				inter_l = pix_lu + (pix_ld - pix_lu) * (y_ori - y_u);
				inter_r = pix_ru + (pix_rd - pix_ru) * (y_ori - y_u);
				framemat[i][j] = inter_l + (inter_r - inter_l) * (x_ori - x_l);
			}
		}
	}

	

	//左上角旋转后坐标
	//const int X_LR = -1 * (mheight * SIN_THETA + mwidth * COS_THETA) / 2 + FRAME_WIDTH / 2;//-0.5mwidth*cos + -0.5mheight*sin
	//const int Y_LR = (mwidth * SIN_THETA - mheight * COS_THETA) / 2 + FRAME_HEIGHT / 2;//- -0.5mwidth*sin + -0.5mheight*cos
	//int X, Y;//待求解点在画布中的坐标
	//int X_ori, Y_ori;//旋转前的点在原图中的坐标
	//int i, j;
	//for (i = 0; i < mheight; i++)
	//{
	//	for (j = 0; j < mwidth; j++)
	//	{
	//		//待求解点在画布中的坐标，画布中心为原点
	//		X = X_LR + i * SIN_THETA + j * COS_THETA;
	//		Y = Y_LR + i * COS_THETA - j * SIN_THETA;
	//		//旋转前的点在原图中的坐标，左上角为原点
	//		//X_ori = (X - FRAME_WIDTH / 2) * COS_THETA + (Y - FRAME_HEIGHT / 2) * SIN_THETA + mwidth / 2;
	//		//Y_ori = (X - FRAME_WIDTH / 2) * SIN_THETA + (Y - FRAME_HEIGHT / 2) * COS_THETA + mheight / 2;
	//		
	//		framemat[Y][X] = pixelmat[i][j];
	//	}
	//}
	pait_fram_edge(framemat);

	return framemat;
}

//平移图像, 返回显示区域(画板)指针
int** moveimage(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	int x = 100;
	int y = 130;
	int i = 0, j = 0;
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			framemat[i + x][j + y] = pixelmat[i][j];
		}
	}
	pait_fram_edge(framemat);

    return framemat;
}

//缩放图像, 返回显示区域(画板)指针
//双线性插值
int** scaling(int** framemat, int** pixelmat, int mheight, int mwidth)
{
	const double scale = 1.7;//缩放倍数
	const double scale_rev = 1 / scale;
	const int M = mheight * scale;
	const int N = mwidth * scale;
	int i = 0, j = 0;
	double x_ori, y_ori;//后相映射到原图的坐标
	int x_l, x_r;
	int y_u, y_d;
	double inter_l, inter_r;
	int pix_lu;
	int pix_ld;
	int pix_ru;
	int pix_rd;
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			x_ori = j * scale_rev;
			y_ori = i * scale_rev;
			x_l = x_ori; x_r = x_l + 1;//插值近邻的x坐标
			y_u = y_ori; y_d = y_u + 1;//插值近邻的y坐标
			if (x_r >= mwidth)
			{
				x_l--;
				x_r--;
				x_ori--;
			}
			if (y_d >= mheight)
			{
				y_u--;
				y_d--;
				y_ori--;
			}
			//四个近邻像素值
			pix_lu = pixelmat[y_u][x_l];
			pix_ld = pixelmat[y_d][x_l];
			pix_ru = pixelmat[y_u][x_r];
			pix_rd = pixelmat[y_d][x_r];
			inter_l = pix_lu + (pix_ld - pix_lu) * (y_ori - y_u);
			inter_r = pix_ru + (pix_rd - pix_ru) * (y_ori - y_u);
			framemat[i][j] = inter_l + (inter_r - inter_l) * (x_ori - x_l);
		}
	}
	pait_fram_edge(framemat);
	return framemat;
}

//DFT变换, 返回处理后的图像, 注意缩放到0~255的整型
int** DFT(int** pixelmat, int mheight, int mwidth)
{
	int u = 0, v = 0;

	int** i_mat;// = new int*[M];
	mynew(i_mat, mheight, mwidth);

	Complex **MatF;
	mynew(MatF, mheight, mwidth);
	Complex **Matf;
	mynew(Matf, mheight, mwidth);
	for (int i = 0; i < mheight; i++)
	{
		for (int j = 0; j < mwidth; j++)
		{
			Matf[i][j] = { (double)(pixelmat[i][j]), 0, 0, 0 };
		}
	}
	MatF = MyDFT(Matf, mheight, mwidth);
	for (u = 0; u < mheight; u++)
	{
		for (v = 0; v < mwidth; v++)
		{
			i_mat[u][v] = MatF[u][v].r;
			if (i_mat[u][v] > 255)
			{
				i_mat[u][v] = 255;
			}
		}
	}

	int swap_temp;
	for (u = 0; u < mheight / 2; u++)
	{
		for (v = 0; v < mwidth / 2; v++)
		{
			swap_temp = i_mat[u][v];
			i_mat[u][v] = i_mat[u + mheight / 2][v + mwidth / 2];
			i_mat[u + mheight / 2][v + mwidth / 2] = swap_temp;

			swap_temp = i_mat[u][v + mwidth / 2];
			i_mat[u][v + mwidth / 2] = i_mat[u + mheight / 2][v];
			i_mat[u + mheight / 2][v] = swap_temp;
		}
	}

	mydelete(Matf, mheight);
	mydelete(MatF, mheight);

    return i_mat;
}

Complex** MyDFT(Complex** Mat, int mheight, int mwidth)
{
	int u = 0, v = 0;
	int M = mheight, N = mwidth;
	double sumRe = 0, sumIm = 0;
	Complex sum = { 0, 0, 0, 0 };
	double theta = 0;

	Complex **H_M, **H_N;//左乘和右乘方阵
	mynew(H_M, M, M);
	mynew(H_N, N, N);

	//计算左乘方阵
	for (u = 0; u < M; u++)
	{
		for (v = u; v < M; v++)
		{
			theta = -2 * PI * u * v / M;
			H_M[u][v] = { cos(theta), sin(theta), 1, theta };
			H_M[v][u] = H_M[u][v];
		}
	}
	//计算右乘方阵
	for (u = 0; u < N; u++)
	{
		for (v = u; v < N; v++)
		{
			theta = -2 * PI * u * v / N;
			H_N[u][v] = { cos(theta), sin(theta), 1, theta };
			H_N[v][u] = H_N[u][v];
		}
	}

	Complex **tempMat;
	mynew(tempMat, mheight, mwidth);
	Complex **rtnMat;
	mynew(rtnMat, mheight, mwidth);

	int k = 0;
	MatricMulti(H_M, Mat, tempMat, mheight, mheight, mwidth);//左乘
	MatricMulti(tempMat, H_N, rtnMat, mheight, mwidth, mwidth);//右乘
	const double coeffi = sqrt(mheight * mwidth);
	for (u = 0; u < M; u++)
	{
		for (v = 0; v < N; v++)
		{
			rtnMat[u][v].re /= coeffi;
			rtnMat[u][v].im /= coeffi;
			rtnMat[u][v].r /= coeffi;
		}
	}

	mydelete(H_M, M);
	mydelete(H_N, N);
	mydelete(tempMat, M);

	return rtnMat;
}

//DCT变换, 返回处理后的图像
int** DCT(int** pixelmat, int mheight, int mwidth)
{
	int i, j;

	//**************** 左乘和右乘矩阵 ******************
	double **Matric_L, **Matric_R;
	mynew(Matric_L, mheight, mheight);
	mynew(Matric_R, mwidth, mwidth);
	//计算左乘矩阵
	double coeffi = sqrt(2.0 / mheight);
	for (j = 0; j < mheight; j++)
	{
		Matric_L[0][j] = coeffi / 1.414;
		for (i = 1; i < mheight; i++)
		{
			Matric_L[i][j] = coeffi * cos(PI * i * (2 * j + 1) / (2 * mheight));
		}
	}
	//计算右乘矩阵
	coeffi = sqrt(2.0 / mwidth);
	for (i = 0; i < mwidth; i++)
	{
		Matric_R[i][0] = coeffi / 1.414;
		for (j = 1; j < mwidth; j++)
		{
			Matric_R[i][j] = coeffi * cos(PI * j * (2 * i + 1) / (2 * mwidth));
		}
	}
	//**************** 左乘和右乘矩阵 ******************

	double **Mat_temp, **Mat_result, **Mat_pixel;
	mynew(Mat_temp, mheight, mwidth);
	mynew(Mat_result, mheight, mwidth);
	mynew(Mat_pixel, mheight, mwidth);
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			Mat_pixel[i][j] = (double)(pixelmat[i][j]);
		}
	}
	MatricMulti(Matric_L, Mat_pixel, Mat_temp, mheight, mheight, mwidth);
	MatricMulti(Mat_temp, Matric_R, Mat_result, mheight, mwidth, mwidth);
	int temp = 0;
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			temp = abs(Mat_result[i][j]);
			pixelmat[i][j] = temp > 255 ? 255 : temp;
		}
	}

	mydelete(Mat_temp, mheight);
	mydelete(Mat_result, mheight);
	mydelete(Mat_pixel, mheight);
	mydelete(Matric_L, mheight);
	mydelete(Matric_R, mwidth);
    return pixelmat;
}

//walsh变换, 返回处理后的图像
int** walsh(int** pixelmat, int mheight, int mwidth)
{
	int i, j, u, v, x, y;
	int expo[2] = { 1, -1 };

	//******************* 计算mheight与mwidth的指数 *********************
	int n_l = log2(mheight);
	int n_r = log2(mwidth);

	//********************* 左乘和右乘矩阵 ***********************
	double **Matric_L, **Matric_R;
	mynew(Matric_L, mheight, mheight);
	mynew(Matric_R, mwidth, mwidth);
	
	walsh_matrix(Matric_L, mheight, n_l);
	walsh_matrix(Matric_R, mwidth, n_r);
	//********************* 左乘和右乘矩阵 ***********************

	double **Mat_temp, **Mat_result, **Mat_pixel;
	mynew(Mat_temp, mheight, mwidth);
	mynew(Mat_result, mheight, mwidth);
	mynew(Mat_pixel, mheight, mwidth);
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			Mat_pixel[i][j] = (double)(pixelmat[i][j]);
		}
	}
	//walsh变换
	MatricMulti(Matric_L, Mat_pixel, Mat_temp, mheight, mheight, mwidth);
	MatricMulti(Mat_temp, Matric_R, Mat_result, mheight, mwidth, mwidth);
	double coeffi = sqrt(mheight * mwidth);
	int temp;
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			temp = abs(Mat_result[i][j]) / coeffi;
			pixelmat[i][j] = temp > 255 ? 255 : temp;
		}
	}


	mydelete(Mat_temp, mheight);
	mydelete(Mat_result, mheight);
	mydelete(Mat_pixel, mheight);
	mydelete(Matric_L, mheight);
	mydelete(Matric_R, mwidth);
    return pixelmat;
}

//构造walsh矩阵, N*N, n = log2(N)
void walsh_matrix(double **&walshMatrx, int N, int n)
{
	int i, j, k;
	int i_quot, i_rem;
	int j_quot, j_rem;
	int sum;
	int temp;
	int expo[2] = { 1, -1 };

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			i_quot = i;
			i_rem = 0;
			j_quot = 0;
			j_rem = j;
			temp = N / 2;
			sum = 0;
			for (k = 0; k < n; k++)
			{
				i_rem = i_quot % 2;
				i_quot = i_quot / 2;

				j_quot = j_rem / temp;
				j_rem = j_rem % temp;
				temp /= 2;
				sum += i_rem * j_quot;
			}
			sum = sum % 2;
			walshMatrx[i][j] = expo[sum];
		}
	}
}

//haar变换, 返回处理后的图像
int** haar(int** pixelmat, int mheight, int mwidth)
{
	int **mat_temp;
	mynew(mat_temp, mheight, mwidth);
	int i, j;
	int halfHeight = mheight / 2;
	int halfWidth = mwidth / 2;
	int quarterHeight = mheight / 4;
	int quarterWidth = mwidth / 4;
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < halfWidth; j++)
		{
			mat_temp[i][j] = (pixelmat[i][2 * j] + pixelmat[i][2 * j + 1]) / 2;
			mat_temp[i][j + halfWidth] = (pixelmat[i][2 * j] - pixelmat[i][2 * j + 1]) / 2;
		}
	}
	for (j = 0; j < mwidth; j++)
	{
		for (i = 0; i < halfHeight; i++)
		{
			pixelmat[i][j] = (mat_temp[2 * i][j] + mat_temp[2 * i + 1][j]) / 2;
			pixelmat[i + halfHeight][j] = (mat_temp[2 * i][j] - mat_temp[2 * i + 1][j]) / 2;
		}
	}
	mydelete(mat_temp, mheight);
	
    return pixelmat;
}

double haar_function(const int k, const double z, const int N)
{
	//计算hk(z),haar矩阵阶数为N
	//k = 2 ^ p + q - 1
	int p, q;
	double coeffi = 1 / sqrt(N);
	double result;

	if (k == 0)
	{
		result = coeffi;
		return result;
	}
	p = log2(k);
	double exp_part = pow(2, p);
	q = k + 1 - exp_part;
	if (z >= (q - 1.0) / exp_part && z < (q - 0.5) / exp_part)
	{
		result = coeffi * sqrt(exp_part);
	}
	else if ((q - 0.5) / exp_part <= z && z < (double)q / exp_part)
	{
		result = -1 * coeffi * sqrt(exp_part);
	}
	else
	{
		result = 0;
	}
	
	return result;
}

//构造haar矩阵,N*N
void haar_matrix(double **&haarMatrix, int N)
{
	double coefficient = 1.0 / N;
	int i, j;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			haarMatrix[i][j] = haar_function(i, j * coefficient, N);
		}
	}
}

//生成随机噪声, 返回处理后的图像;
int** randomnoise(int** pixelmat, int mheight, int mwidth)
{
	int i, j;
	int noise;
	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			noise = gaussrand();
			pixelmat[i][j] += noise * 30;
			if (pixelmat[i][j] < 0)
			{
				pixelmat[i][j] = 0;
			}
			if (pixelmat[i][j] > 255)
			{
				pixelmat[i][j] = 255;
			}
		}
	}
    return pixelmat;
}

double gaussrand()
{
	static double U, V;
	static int phase = 0;
	double Z;

	if (phase == 0)
	{
		U = rand() / (RAND_MAX + 1.0);
		V = rand() / (RAND_MAX + 1.0);
		Z = sqrt(-2.0 * log(U))* sin(2.0 * PI * V);
	}
	else
	{
		Z = sqrt(-2.0 * log(U)) * cos(2.0 * PI * V);
	}

	phase = 1 - phase;
	return Z;
}

//生成椒盐噪声, 返回处理后的图像
int** impulsenoise(int** pixelmat, int mheight, int mwidth)
{
	int i, j;
	int low = RAND_MAX >> 7;
	int high = RAND_MAX - low;
	int tmp;

	for (i = 0; i < mheight; i++)
	{
		for (j = 0; j < mwidth; j++)
		{
			tmp = rand();
			if (tmp < low)
			{
				pixelmat[i][j] = 0;
			}
			if (tmp > high)
			{
				pixelmat[i][j] = 255;
			}
		}
	}
    return pixelmat;
}







//逆滤波复原
int** inversefit(int** pixelmat, int mheight, int mwidth)
{
    return NULL;
}

//维纳滤波
int** wienerfit(int** pixelmat, int mheight, int mwidth)
{
    return NULL;
}

//示例: JPEG压缩及解压缩
int** jpeg(int** pixelmat, int mheight, int mwidth)
{
    return NULL;
}









//************** 动态分配删除空间 ****************
template<typename T>
void mynew(T** &mat, int height, int width)
{
	mat = new T*[height];
	int i = 0;
	for (i = 0; i < height; i++)
	{
		mat[i] = new T[width];
	}
}

template<typename T>
void mydelete(T** &mat, int height)
{
	int i = 0;
	for (i = 0; i < height; i++)
	{
		delete[] mat[i];
	}
	delete[] mat;
}
//************** 动态分配删除空间 ****************

void pait_fram_edge(int** &framemat)
{
	int i;
	for (i = 0; i < FRAME_HEIGHT; i++)
	{
		framemat[i][0] = 0;
		framemat[i][FRAME_WIDTH - 1] = 0;
	}
	for (i = 0; i < FRAME_WIDTH; i++)
	{
		framemat[0][i] = 0;
		framemat[FRAME_HEIGHT - 1][i] = 0;
	}
}

//******************** 矩阵乘法 **********************
template<typename T>
void MatricMulti(T** &Matri_L, T** &Matri_R, T** &result, const int A, const int B, const int C)
{
	int i, j, k;
	T sum = 0;
	for (i = 0; i < A; i++)
	{
		for (int j = 0; j < C; j++)
		{
			sum = 0;
			for (k = 0; k < B; k++)
			{
				sum += Matri_L[i][k] * Matri_R[k][j];
			}
			result[i][j] = sum;
		}
	}
}

template<>
void MatricMulti(Complex** &Matri_L, Complex** &Matri_R, Complex** &result, const int A, const int B, const int C)
{
	int i, j, k;
	Complex sum = { 0, 0, 0, 0 };
	for (i = 0; i < A; i++)
	{
		for (int j = 0; j < C; j++)
		{
			sum = { 0, 0, 0, 0 };
			for (k = 0; k < B; k++)
			{
				sum.re += (Matri_L[i][k].re * Matri_R[k][j].re - Matri_L[i][k].im * Matri_R[k][j].im);
				sum.im += (Matri_L[i][k].re * Matri_R[k][j].im + Matri_L[i][k].im * Matri_R[k][j].re);
			}
			sum.r = sqrt(sum.re * sum.re + sum.im * sum.im);
			sum.theta = atan2(sum.im, sum.re);
			result[i][j] = sum;
		}
	}
}
//******************** 矩阵乘法 **********************

//方阵转置
template<typename T>
void SquarMatrixTrans(T** &Matrix, int N)
{
	int i, j;
	T temp;
	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			temp = Matrix[i][j];
			Matrix[i][j] = Matrix[j][i];
			Matrix[j][i] = temp;
		}
	}
}


