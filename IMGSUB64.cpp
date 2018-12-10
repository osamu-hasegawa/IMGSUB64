// IMGSUB64.cpp : DLL アプリケーション用にエクスポートされる関数を定義します。
//

#include "stdafx.h"


//#include "windows.h"
#include <math.h>
#include <float.h>
#include "opencv/cv.h"
#include "opencv/highgui.h"
#ifdef _DEBUG
#pragma comment(lib, ".\\/opencv/x64/vc10/lib/opencv_core2411d.lib")
#pragma comment(lib, ".\\/opencv/x64/vc10/lib/opencv_imgproc2411d.lib")
#pragma comment(lib, ".\\/opencv/x64/vc10/lib/opencv_highgui2411d.lib")
#else
#pragma comment(lib, ".\\/opencv/x64/vc10/lib/opencv_core2411.lib")
#pragma comment(lib, ".\\/opencv/x64/vc10/lib/opencv_imgproc2411.lib")
#pragma comment(lib, ".\\/opencv/x64/vc10/lib/opencv_highgui2411.lib")
#endif

//extern
//VOID THRESH_HSV(LPBYTE pIMGH, LPBYTE pIMGS, LPBYTE pIMGV, LPBYTE pIMGB, int width, int height, int stride, int hmin, int hmax, int smin, int smax, int vmin, int vmax);

IplImage	 *m_img_a = NULL;//rgb
IplImage	 *m_img_h = NULL;//hsv
IplImage	 *m_img_rgb[] = { NULL, NULL, NULL };//rgb
IplImage	 *m_img_hsv[] = { NULL, NULL, NULL };//hsv
IplImage	 *m_img_g = NULL;
IplImage	 *m_img_m = NULL;
IplImage	*m_img_b = NULL;
IplImage	*m_img_f = NULL;
IplImage	*m_img_d = NULL;
IplImage	*m_img_t = NULL;
IplImage	*m_img_u = NULL;
IplImage	*m_img_v = NULL;
IplImage	*m_img_w = NULL;
CvFont		m_fnt;


int		WID;
int		HEI;
enum {
	IMG_A = 0,
	IMG_G = 1,
	IMG_B = 2,
	IMG_M = 3,
	IMG_H = 4,
	IMG_RGB_R = 5,
	IMG_RGB_G = 6,
	IMG_RGB_B = 7,
	IMG_HSV_H = 8,
	IMG_HSV_S = 9,
	IMG_HSV_V =10,
	IMG_F = 11,
	IMG_D = 12,
	IMG_T = 13,
	IMG_U = 14,
	IMG_V = 15,
	IMG_W = 16,
};
static
void THRESH_HSV(LPBYTE pIMGH, LPBYTE pIMGS, LPBYTE pIMGV, LPBYTE pIMGB, LONG width, LONG height, LONG stride, LONG hmin, LONG hmax, LONG smin, LONG smax, LONG vmin, LONG vmax)
{
	LPBYTE	ph = pIMGH,
			ps = pIMGS,
			pv = pIMGV,
			pb = pIMGB;
	int		x, y, ofs = stride-width;
	// 170-180-30
	if (hmin > hmax) {
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				if ((*ph >= hmin || *ph <= hmax)
				  &&(*ps >= smin && *ps <= smax)
				  &&(*pv >= vmin && *pv <= vmax)) {
					*pb = 255;
				}
				else {
					*pb = 0;
				}
				ph++;
				ps++;
				pv++;
				pb++;
			}
			ph += ofs;
			ps += ofs;
			pv += ofs;
			pb += ofs;
		}
	}
	else {
		for (y = 0; y < height; y++) {
			for (x = 0; x < width; x++) {
				if ((*ph >= hmin && *ph <= hmax)
				  &&(*ps >= smin && *ps <= smax)
				  &&(*pv >= vmin && *pv <= vmax)) {
					*pb = 255;
				}
				else {
					*pb = 0;
				}
				ph++;
				ps++;
				pv++;
				pb++;
			}
			ph += ofs;
			ps += ofs;
			pv += ofs;
			pb += ofs;
		}
	}
}
static
IplImage* TO_IMG(LONG I)
{
	IplImage*	img;

	switch (I) {
	case IMG_A    :img = m_img_a; break;
	case IMG_G    :img = m_img_g; break;
	case IMG_B    :img = m_img_b; break;
	case IMG_M    :img = m_img_m; break;
	case IMG_H    :img = m_img_h; break;
	case IMG_F    :img = m_img_f; break;
	case IMG_D    :img = m_img_d; break;
	case IMG_T    :img = m_img_t; break;
	case IMG_U    :img = m_img_u; break;
	case IMG_V    :img = m_img_v; break;
	case IMG_W    :img = m_img_w; break;
	case IMG_RGB_R:img = m_img_rgb[0]; break;
	case IMG_RGB_G:img = m_img_rgb[1]; break;
	case IMG_RGB_B:img = m_img_rgb[2]; break;
	case IMG_HSV_H:img = m_img_hsv[0]; break;
	case IMG_HSV_S:img = m_img_hsv[1]; break;
	case IMG_HSV_V:img = m_img_hsv[2]; break;
	default       :img = NULL        ; break;
	}
	return(img);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_TERM(void)
{
	int	i;
	if (m_img_a != NULL) {
		cvReleaseImage(&m_img_a);
		m_img_a = NULL;
	}
	if (m_img_h != NULL) {
		cvReleaseImage(&m_img_h);
		m_img_h = NULL;
	}
	for (i = 0; i < 3; i++) {
		if (m_img_rgb[i] != NULL) {
			cvReleaseImage(&m_img_rgb[i]);
			m_img_rgb[i] = NULL;
		}
		if (m_img_hsv[i] != NULL) {
			cvReleaseImage(&m_img_hsv[i]);
			m_img_hsv[i] = NULL;
		}
	}
	if (m_img_g != NULL) {
		cvReleaseImage(&m_img_g);
		m_img_g = NULL;
	}
	if (m_img_b != NULL) {
		cvReleaseImage(&m_img_b);
		m_img_b = NULL;
	}
	if (m_img_m != NULL) {
		cvReleaseImage(&m_img_m);
		m_img_m = NULL;
	}
	if (m_img_f != NULL) {
		cvReleaseImage(&m_img_f);
		m_img_f = NULL;
	}
	if (m_img_d != NULL) {
		cvReleaseImage(&m_img_d);
		m_img_d = NULL;
	}
	if (m_img_t != NULL) {
		cvReleaseImage(&m_img_t);
		m_img_t = NULL;
	}
	if (m_img_u != NULL) {
		cvReleaseImage(&m_img_u);
		m_img_u = NULL;
	}
	if (m_img_v != NULL) {
		cvReleaseImage(&m_img_v);
		m_img_v = NULL;
	}
	if (m_img_w != NULL) {
		cvReleaseImage(&m_img_w);
		m_img_w = NULL;
	}
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_RESET_MASK(LONG x, LONG y, LONG w, LONG h)
{
	CvPoint p1, p2;
	p1.x = x;
	p1.y = y;
	p2.x = x+w;
	p2.y = y+h;
	cvZero(m_img_m);
	cvRectangle(m_img_m, p1, p2, CV_RGB(255,255,255), 1, 8, 0);
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_RESET(LONG wid, LONG hei)//, LONG mx, LONG my, LONG mw, LONG mh)
{
	int	i;
	WID = wid;
	HEI = hei;
	//---
	//img_src = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	//img_gry = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	//img_dst = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);


	OCV_TERM();
	m_img_a = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 3);
	m_img_h = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 3);
	for (i = 0; i < 3; i++) {
		m_img_rgb[i] = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
		m_img_hsv[i] = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	}
	m_img_g = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	m_img_b = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	m_img_m = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	m_img_f = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_32F, 1);
	m_img_d = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	m_img_t = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	m_img_u = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	m_img_v = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	m_img_w = cvCreateImage(cvSize(WID, HEI), IPL_DEPTH_8U, 1);
	//---
	if (TRUE) {
		double hscale = 3.2 * WID / 2592.0;
		double vscale = 3.2 * HEI / 1944.0;
		int thickness = (int)(4 * WID / 2592.0 + 0.5);
		if (thickness < 1) {
			thickness = 1;
		}
		cvInitFont(&m_fnt, CV_FONT_HERSHEY_COMPLEX, hscale, vscale, 1.0, thickness, 8);
		//Cv.InitFont(out fnt, FontFace.HersheySimplex, 4.0, 4.0, 1.0, 4);
		//Cv.InitFont(out fnt, FontFace.HersheyDuplex, 4.0, 4.0, 1.0, 4);
		//Cv.InitFont(out fnt, FontFace.HersheyPlain, 3.5, 3.5, 1.0, 4);
		//Cv.InitFont(out fnt, FontFace.HersheyTriplex, 3.5, 3.5, 1.0, 4);
		//Cv.InitFont(out fnt, FontFace.HersheyScriptComplex, 3.5, 3.5, 1.0, 4);
		/*
		CV_FONT_HERSHEY_SIMPLEX 普通サイズの sans-serif フォント
		CV_FONT_HERSHEY_PLAIN 小さいサイズの sans-serif フォント
		CV_FONT_HERSHEY_DUPLEX 普通サイズの sans-serif フォント（ CV_FONT_HERSHEY_SIMPLEX よりも複雑）
		CV_FONT_HERSHEY_COMPLEX 普通サイズの serif フォント
		CV_FONT_HERSHEY_TRIPLEX 普通サイズの serif フォント（ CV_FONT_HERSHEY_COMPLEX よりも複雑）
		CV_FONT_HERSHEY_COMPLEX_SMALL CV_FONT_HERSHEY_COMPLEX の小さいサイズ版
		CV_FONT_HERSHEY_SCRIPT_SIMPLEX 手書きスタイルのフォント
		CV_FONT_HERSHEY_SCRIPT_COMPLEX CV_FONT_HERSHEY_SCRIPT_SIMPLEX 
		*/
	}
//	reset_mask_img();
	//OCV_RESET_MASK(mx, my, mw, mh);
}

// bpp:32 ... 4byte/pixel
extern "C" __declspec(dllexport)
LONG APIENTRY OCV_SET_IMG(void*ptr, LONG wid, LONG hei, LONG str, LONG bpp)
{
	int	BPP = bpp/8;
	int	size = hei*str;//*BPP;

	if (wid != WID || hei != HEI) {
		return(FALSE);
	}
	else if (size == m_img_a->imageSize) {
		memcpy(m_img_a->imageData, ptr, size);
	}
	else if (bpp == 32) {
		LPBYTE pd = (LPBYTE)m_img_a->imageData;
		LPBYTE ps = (LPBYTE)ptr;
		int	os = (str-wid*BPP);
		int od = (m_img_a->widthStep-m_img_a->width*3);

		for (int y = 0; y < hei; y++) {
			for (int x = 0; x < wid; x++) {
				*pd++ = *ps++;
				*pd++ = *ps++;
				*pd++ = *ps++;
				ps++;
			}
			ps += os;
			pd += od;
		}
	}
	else {
		return(FALSE);
	}
	return(TRUE);
}
extern "C" __declspec(dllexport)
LONG APIENTRY OCV_GET_IMG(void*ptr, LONG wid, LONG hei, LONG str, LONG bpp)
{
	int	BPP = bpp/8;
	int	size = hei*str;//*BPP;

	if (wid != WID || hei != HEI) {
		return(FALSE);
	}
	else if (size == m_img_a->imageSize) {
		memcpy(ptr, m_img_a->imageData, size);
	}
	else {
		return(FALSE);
	}
	return(TRUE);
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_TO_GRAY(LONG I, LONG H)
{
	IplImage* I_IMG = TO_IMG(I);
	IplImage* H_IMG = TO_IMG(H);

	cvCvtColor(I_IMG, H_IMG, CV_RGBA2GRAY);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_TO_HSV(LONG I, LONG H)
{
	IplImage* I_IMG = TO_IMG(I);
	IplImage* H_IMG = TO_IMG(H);

	cvCvtColor(I_IMG, H_IMG, CV_BGR2HSV);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_MERGE(LONG H1, LONG H2, LONG H3, LONG I)
{
	IplImage* I_IMG = TO_IMG(I);
	IplImage* H1_IMG = TO_IMG(H1);
	IplImage* H2_IMG = TO_IMG(H2);
	IplImage* H3_IMG = TO_IMG(H3);

	cvMerge(H1_IMG, H2_IMG, H3_IMG, NULL, I_IMG);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_SPLIT(LONG I, LONG H1, LONG H2, LONG H3)
{
	IplImage* I_IMG = TO_IMG(I);
	IplImage* H1_IMG = TO_IMG(H1);
	IplImage* H2_IMG = TO_IMG(H2);
	IplImage* H3_IMG = TO_IMG(H3);

	cvSplit(I_IMG, H1_IMG, H2_IMG, H3_IMG, NULL);// m_img_a:BGRの順
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_SMOOTH(LONG I, LONG cof)
{
	IplImage* I_IMG = TO_IMG(I);

	cvSmooth(I_IMG, I_IMG, CV_GAUSSIAN, cof, cof, 0, 0);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_THRESH_BIN(LONG I, LONG H, LONG thval, LONG inv)
{
	IplImage* I_IMG = TO_IMG(I);
	IplImage* H_IMG = TO_IMG(H);

	if (inv != 0) {
		//白背景に黒丸の時は反転しておく
		cvThreshold(I_IMG, H_IMG, thval, 255, CV_THRESH_BINARY_INV);
	}
	else {
		//白背景に黒丸の時は反転しておく
		cvThreshold(I_IMG, H_IMG, thval, 255, CV_THRESH_BINARY);
	}
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_THRESH_HSV(LONG I1, LONG I2, LONG I3, LONG H, LONG minh, LONG maxh, LONG mins, LONG maxs, LONG minv, LONG maxv)
{
	IplImage* I1_IMG = TO_IMG(I1);
	IplImage* I2_IMG = TO_IMG(I2);
	IplImage* I3_IMG = TO_IMG(I3);
	IplImage* H1_IMG = TO_IMG(H);

	THRESH_HSV(
		(LPBYTE)I1_IMG->imageData,
		(LPBYTE)I2_IMG->imageData,
		(LPBYTE)I3_IMG->imageData,
		(LPBYTE)H1_IMG->imageData,
		H1_IMG->width,
		H1_IMG->height,
		H1_IMG->widthStep,
		minh, maxh,//160, 20
		mins, maxs,//40, 140,
		minv, maxv //90, 180
		);
}
static
double C_NAN(void)
{
	double	k = -1;
	double	f = log(k);
	_int64	*p = (_int64*)&f;
	//非数: *p = 0xfff8 0000 0000 0000	__int64

	return(f);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_CAL_HIST(LONG I, LONG bMASK, double* pval, double* pmin, double* pmax, double* pavg)
{
	CvHistogram*hist;
	int			hist_size[] = { 256 };
	float		range_0[] = { 0, 256 };
	float*		ranges[] = { range_0 };
	IplImage*	img = TO_IMG(I);
	double		fmin, fmax;//, favg;
	int cnt = 0;
	double f, sum = 0;

	hist = cvCreateHist(1, hist_size, CV_HIST_ARRAY, ranges, 1);
	if (bMASK) {
		cvCalcHist(&img, hist, 0, m_img_m);
	}
	else {
		cvCalcHist(&img, hist, 0, NULL);
	}

	fmin = fmax = C_NAN();

	for (int i = 0; i < 256; i++) {
		f = cvQueryHistValue_1D(hist, i);
		pval[i] = f;
		cnt += (int)f;
		sum += (i * f);
		if (f > 0) {
			//contrast計算用,画素値の最小と最大
			if (_isnan(fmin)) {
				fmin = i;
			}
			if (TRUE) {
				fmax = i;
			}
		}
	}
	//f = m_width * m_height;
	//if (cnt != f) {
	//	f = 0;//検算:FOR BREAK POINT
	//}
	if (pmin != NULL) {
		*pmin = fmin;
	}
	if (pmax != NULL) {
		*pmax = fmax;
	}
	if (pavg != NULL) {
		*pavg = sum / cnt;
	}
	cvReleaseHist(&hist);
}
static
CvScalar CV_COLOR(DWORD c)
{
	return(CV_RGB(GetRValue(c), GetGValue(c), GetBValue(c)));
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_PUTTEXT(LONG I, LPCSTR buf, LONG x, LONG y, DWORD c)
{
	IplImage*	I_IMG = TO_IMG(I);
	CvPoint		pt;
	pt.x = x;
	pt.y = y;
	cvPutText(I_IMG, buf, pt, &m_fnt, CV_COLOR(c));
}

CvMemStorage *storage = NULL;
CvSeq *contours = NULL;

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_FIND_FIRST(LONG I, LONG mode)
{
	IplImage*	I_IMG = TO_IMG(I);
	int find_contour_num;
	CvPoint ofs;
	ofs.x = 0;
	ofs.y = 0;
	storage = cvCreateMemStorage(0);
	contours = NULL;
	//抽出モード:
	//		0:CV_RETR_EXTERNAL:最も外側の輪郭のみ抽出
	//		1:CV_RETR_LIST	:全ての輪郭を抽出し，リストに追加
	//		2:CV_RETR_CCOMP	:全ての輪郭を抽出し，二つのレベルを持つ階層構造を構成する．
	//						:1番目のレベルは連結成分の外側の境界線，
	//						:2番目のレベルは穴（連結成分の内側に存在する）の境界線．
	//		3:CV_RETR_TREE	:全ての輪郭を抽出し，枝分かれした輪郭を完全に表現する階層構造を構成する．
	//----------------------------------
	//		CV_CHAIN_APPROX_SIMPLE	:輪郭の折れ線の端点を取得
	//		CV_CHAIN_APPROX_NONE	:輪郭の全ての点を取得,Teh-Chinチェーンの近似アルゴリズム中の一つを適用する 
	//		CV_CHAIN_APPROX_TC89_L1
	//		CV_CHAIN_APPROX_TC89_KCOS


	find_contour_num = cvFindContours(
								I_IMG,                     // 入力画像
								storage,                      // 抽出された輪郭を保存する領域
								&contours,                  // 一番外側の輪郭へのポインタへのポインタ
								sizeof (CvContour),      // シーケンスヘッダのサイズ
#if 1
								mode,	// CV_RETR_EXTERNAL,
//								CV_CHAIN_APPROX_SIMPLE
								CV_CHAIN_APPROX_NONE,
#else
								CV_RETR_TREE,
								CV_CHAIN_APPROX_NONE,
#endif
								ofs
		);

}

#define PI	(3.141592)
//double	fmax = -1e99;
//CvSeq*	pmax = NULL;

extern "C" __declspec(dllexport)
LPVOID APIENTRY OCV_FIND_NEXT(void*ppos, LONG smax, LONG smin, LONG lmax, LONG lmin, double cmax, double cmin, double *ps, double *pl, double *pc)
{
	CvSeq*	pos = (ppos == NULL) ? contours: (CvSeq*)ppos;

	if (ppos == NULL) {
		pos = contours;
	}
	else {
		pos = ((CvSeq*)ppos)->h_next;
	}

	for (; pos != NULL; pos = pos->h_next) {
		double so = cvContourArea(pos);
		double l = cvArcLength(pos);
		double	s = fabs(so);
		//円形度＝4π×（面積）÷（周囲長）^2(1:真円,正方形:0.785,正三角形:0.604)
		double c = 4 * PI * s / pow(l, 2);
		double p= C_NAN();
		//CvRect rc;

		if (s > smax || s < smin) {
			continue;
		}
		if (l > lmax || l < lmin) {
			continue;
		}
		if (c > cmax || c < cmin) {
			continue;
		}
		*ps = s; 
		*pl = l;
		*pc = c;
		return((void*)pos);
	}
	return(0);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_FIND_TERM(void)
{
	//メモリストレージの解放
	cvReleaseMemStorage (&storage);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_DRAW_CONTOURS(LONG I, void*pos, DWORD c1, DWORD c2)
{
	IplImage*	img = TO_IMG(I);
	cvDrawContours(img, (CvSeq*)pos, CV_COLOR(c1), CV_COLOR(c2), 0, 2);//Cv.FILLED);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_DRAW_CONTOURS2(LONG I, void*pos, DWORD c1, DWORD c2, DWORD thickness)
{
	IplImage*	img = TO_IMG(I);
	cvDrawContours(img, (CvSeq*)pos, CV_COLOR(c1), CV_COLOR(c2), 0, thickness);
}
extern "C" __declspec(dllexport)
LONG APIENTRY OCV_CONTOURS_CNT(void*pos)
{
	return((LONG)((CvSeq*)pos)->total);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_CONTOURS_PTS(void*pos, LONG idx, LPPOINT p)
{
	CvPoint *pc;
	pc = CV_GET_SEQ_ELEM(CvPoint, (CvSeq*)pos, idx);
//	;
	p->x = pc->x;
	p->y = pc->y;
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_FIT_LINE(void*pos, float* pf)
{
	//それぞれの距離関数における数値パラメータ（C）．0を指定した場合，最適な値が選択される． 
	double	param = 0;
	//reps, aeps
	//半径（座標原点と線の距離）と角度に対する精度．それぞれ0.01が初期値として適している． 
	double reps = 0.01, aeps = 0.01;

	//CV_DIST_L2 (L2):
	// ρ(r)=r2/2 （最も単純で高速な最小二乗法）


	cvFitLine((CvSeq*)pos, CV_DIST_L2, param, reps, aeps, pf);

	//line
	// 出力される線のパラメータ．2次元フィッティングの場合，
	// 四つの浮動小数点型数(vx, vy, x0, y0)の配列で，
	// (vx, vy)は正規化された方向ベクトル，
	// (x0, y0)は線上の点を意味する．
}



#define MAX_POINT	(1024)
POINT pts[MAX_POINT];

extern "C" __declspec(dllexport)
LONG APIENTRY OCV_APPROX_PTS(void*pos, LONG bSIGNE, LONG PREC)
{
	CvSeq* tmp;
	int n, i;

	tmp = cvApproxPoly((CvSeq*)pos, sizeof(CvContour), NULL, CV_POLY_APPROX_DP, PREC);

	n = tmp->total;
	//n = n;
	if (n >= 4) {
		CvPoint	cp;
		if (n > MAX_POINT) {
			n = MAX_POINT;
		}
		for (i = 0; i < n; i++) {
			if (bSIGNE) {//左回り時は順番の入れ替え
				cp = *CV_GET_SEQ_ELEM( CvPoint, tmp, n - 1 - i);
			}
			else {
				cp = *CV_GET_SEQ_ELEM( CvPoint, tmp, i );
			}
			pts[i].x = cp.x;
			pts[i].y = cp.y;
		}
	}
	return(n);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_GET_PTS(LONG idx, LPPOINT p)
{
	p->x = pts[idx].x;
	p->y = pts[idx].y;
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_DRAW_LINE(LONG I, LPPOINT p1, LPPOINT p2, DWORD c, DWORD thick)
{
	IplImage*	img = TO_IMG(I);
	CvPoint	cp1, cp2;
	cp1.x = p1->x;
	cp1.y = p1->y;
	cp2.x = p2->x;
	cp2.y = p2->y;

	cvLine(img, cp1, cp2, CV_COLOR(c), thick);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_DRAW_RECT(LONG I, LPRECT pr, DWORD c, DWORD thickness)
{
	IplImage*	img = TO_IMG(I);
	CvPoint	cp1, cp2;
	cp1.x = pr->left  ;
	cp1.y = pr->top   ;
	cp2.x = pr->right ;
	cp2.y = pr->bottom;

	cvRectangle(img, cp1, cp2, CV_COLOR(c), thickness);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_BOUNDING_RECT(void*pos, LPRECT pr)
{
	CvRect	cr = cvBoundingRect(pos);

	pr->left   = cr.x;
	pr->top    = cr.y;
	pr->right  = cr.x+cr.width;
	pr->bottom = cr.y+cr.height;
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_DRAW_TEXT(LONG I, LONG x, LONG y, LPCSTR buf, DWORD c)
{
	IplImage*	img = TO_IMG(I);
	int baseline;
	CvSize tsize;
	cvGetTextSize(buf, &m_fnt, &tsize, &baseline);

	CvPoint pnt;
	pnt.x = x - tsize.width / 2;
	pnt.y = y - tsize.height / 2;
	if (pnt.x < 0) {
		pnt.x = 0;
	}
	if (pnt.y < 0) {
		pnt.y = 0;
	}
	if ((pnt.x + tsize.width) >= WID) {
		pnt.x = WID - tsize.width;
	}
	if ((pnt.y + tsize.height) >= HEI) {
		pnt.y = HEI - tsize.height;
	}
	cvPutText(img, buf, pnt, &m_fnt, CV_COLOR(c));
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_MIN_AREA_RECT2(void*ppos, LPPOINT p1, LPPOINT p2, LPPOINT p3, LPPOINT p4)
{
	CvSeq*	pos = (CvSeq*)ppos;
	
	CvBox2D box = cvMinAreaRect2(pos);
	CvPoint2D32f ptf[4];// = new CvPoint2D32f[4];
	//CvPoint[] pts = new CvPoint[4];

	cvBoxPoints(box, ptf);
	p1->x = (int)(ptf[0].x+0.5f);
	p1->y = (int)(ptf[0].y+0.5f);
	p2->x = (int)(ptf[1].x+0.5f);
	p2->y = (int)(ptf[1].y+0.5f);
	p3->x = (int)(ptf[2].x+0.5f);
	p3->y = (int)(ptf[2].y+0.5f);
	p4->x = (int)(ptf[3].x+0.5f);
	p4->y = (int)(ptf[3].y+0.5f);
}
static
CvPoint	cpts[MAX_POINT];

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_FILL_POLY(LONG I, LPPOINT p, LONG n, DWORD c)
{
	IplImage*	img = TO_IMG(I);
	if (n > MAX_POINT) {
		n = MAX_POINT;
	}
	for (int i = 0; i < n; i++, p++) {
		cpts[i].x = p->x;
		cpts[i].y = p->y;
	}
	cvFillConvexPoly(img, cpts, n, CV_COLOR(c));
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_ZERO(LONG I)
{
	IplImage*	img = TO_IMG(I);
	cvZero(img);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_SOBEL(LONG I, LONG H, DWORD xorder, DWORD yorder, DWORD apert_size)
{
	IplImage*	src_img = TO_IMG(I);
	IplImage*	dst_img = TO_IMG(H);
	IplImage*	tmp_img = TO_IMG(IMG_F);
	cvSobel(src_img, tmp_img, xorder, yorder, apert_size);
	cvConvertScaleAbs(tmp_img, dst_img, 1, 0);
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_LAPLACE(LONG I, LONG H, DWORD apert_size)
{
	IplImage*	src_img = TO_IMG(I);
	IplImage*	dst_img = TO_IMG(H);
	IplImage*	tmp_img = TO_IMG(IMG_F);
	cvLaplace(src_img, tmp_img, apert_size);
	cvConvertScaleAbs(tmp_img, dst_img, 1, 0);
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_CANNY(LONG I, LONG H, double th1, double th2, DWORD apert_size)
{
	IplImage*	src_img = TO_IMG(I);
	IplImage*	dst_img = TO_IMG(H);
	cvCanny(src_img, dst_img, th1, th2, apert_size);
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_MINMAX(LONG I, double* pmin, double* pmax)
{
	IplImage*	src_img = TO_IMG(I);
	cvMinMaxLoc(src_img, pmin, pmax);
}

// Y = scale * X + shift
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_SCALE(LONG I, LONG H, double scale, double shift)
{
	IplImage*	src_img = TO_IMG(I);
	IplImage*	dst_img = TO_IMG(H);
	cvConvertScaleAbs(src_img, dst_img, scale, shift);
}

extern "C" __declspec(dllexport)
VOID APIENTRY OCV_SMOOTH2(LONG I, LONG cof, double sig1, double sig2)
{
	IplImage* I_IMG = TO_IMG(I);

	cvSmooth(I_IMG, I_IMG, CV_GAUSSIAN, cof, cof, sig1, sig2);
}
//
// J=I-H
//
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_DIFF(LONG I, LONG H, LONG J)
{
	IplImage*	src1_img = TO_IMG(I);
	IplImage*	src2_img = TO_IMG(H);
	IplImage*	dst_img = TO_IMG(J);

	cvSub(src1_img, src2_img, dst_img, NULL);
}


static
void thinningIte(IplImage* img, IplImage* tmp, LONG pattern)
{
    int v9,v2,v3;
    int v8,v1,v4;
    int v7,v6,v5;
	int	x, y;
	int v9i, v2i, v3i;
	int v8i, v1i, v4i;
	int v7i, v6i, v5i;

	cvSet(tmp, cvScalar(1));

	v1i = ( 0) * img->widthStep + ( 0) ;
    v2i = (-1) * img->widthStep + ( 0) ;
    v3i = (-1) * img->widthStep + (+1) ;
    v4i = ( 0) * img->widthStep + (+1) ;
    v5i = (+1) * img->widthStep + (+1) ;
    v6i = (+1) * img->widthStep + ( 0) ;
    v7i = (+1) * img->widthStep + (-1) ;
    v8i = ( 0) * img->widthStep + (-1) ;
    v9i = (-1) * img->widthStep + (-1) ;

	for (y = 1; y < img->height-1; ++y) {
		v1i = y * img->widthStep;
        for (x = 1; x < img->width-1; ++x, v1i++) {
			//---
			v1=img->imageData[ v1i ];
            v2=img->imageData[ v1i+v2i];
            v3=img->imageData[ v1i+v3i];
            v4=img->imageData[ v1i+v4i];
            v5=img->imageData[ v1i+v5i];
            v6=img->imageData[ v1i+v6i];
            v7=img->imageData[ v1i+v7i];
            v8=img->imageData[ v1i+v8i];
            v9=img->imageData[ v1i+v9i];
			//---
            int S  =
				(v2 == 0 && v3 != 0) + (v3 == 0 && v4 != 0) +
				(v4 == 0 && v5 != 0) + (v5 == 0 && v6 != 0) +
				(v6 == 0 && v7 != 0) + (v7 == 0 && v8 != 0) +
				(v8 == 0 && v9 != 0) + (v9 == 0 && v2 != 0);
            
            int N  = v2 + v3 + v4 + v5 + v6 + v7 + v8 + v9;
            
            int m1=0, m2=0;
            
            if(pattern==0) m1 = (v2 * v4 * v6);
            if(pattern==1) m1 = (v2 * v4 * v8);
            
            if(pattern==0) m2 = (v4 * v6 * v8);
            if(pattern==1) m2 = (v2 * v6 * v8);
            
            if (S == 1 && (N >= 2 && N <= 6) && m1 == 0 && m2 == 0) {
				tmp->imageData[v1i]=0;//ピクセル削除
			}
        }
    }
    cvAnd(img, tmp, img, NULL);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_TO_01(LONG I, LONG ZERO_VAL, LONG NONZERO_VAL)
{
	IplImage*	img = TO_IMG(I);
	LPBYTE pb;
	BYTE	zero = (BYTE)ZERO_VAL;
	BYTE	none = (BYTE)NONZERO_VAL;
	int x, y;
	for (y = 0; y < img->height; y++) {
		pb = (LPBYTE)(img->imageData + y*img->widthStep);
        for (x = 0; x < img->width; x++, pb++) {
			if (!*pb) {
				*pb = zero;
			}
			else {
				*pb = none;
			}
		}
	}
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_THINNING(LONG I, LONG H, LONG cnt)
{
	IplImage*	src = TO_IMG(I);
	IplImage*	dst = TO_IMG(H);
	IplImage*	pre = TO_IMG(IMG_U);
	IplImage*	dif = TO_IMG(IMG_V);
	IplImage*	tmp = TO_IMG(IMG_W);

	cvCopy(src, dst);
	OCV_TO_01(H, 0, 1); // 0は0 , 1以上は1に変換される
    
	cvZero(pre);
    
    for ( ;cnt != 0; cnt--) {
        thinningIte(dst, tmp, 0);
        thinningIte(dst, tmp, 1);
        cvAbsDiff(dst, pre, dif);
        cvCopy(dst, pre);
		if (cvCountNonZero(dif) <= 0) {
			break;
		}
    }
    
	OCV_TO_01(H, 0, 255); 
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_COPY(LONG I, LONG H)
{
	IplImage*	src = TO_IMG(I);
	IplImage*	dst = TO_IMG(H);

	cvCopy(src, dst);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_NOT(LONG I, LONG H)
{
	IplImage*	src = TO_IMG(I);
	IplImage*	dst = TO_IMG(H);

	cvNot(src, dst);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_ERODE(LONG I, LONG H, LONG kernel_size, LONG cnt)
{
	IplImage*	src = TO_IMG(I);
	IplImage*	dst = TO_IMG(H);
	IplConvKernel*
				elem = NULL;
	int			size = kernel_size;
	int			offs = size/2;

	elem = cvCreateStructuringElementEx(size, size, offs, offs, CV_SHAPE_RECT, NULL);

	cvErode(src, dst, elem, cnt);

	cvReleaseStructuringElement(&elem);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_DILATE(LONG I, LONG H, LONG kernel_size, LONG cnt)
{
	IplImage*	src = TO_IMG(I);
	IplImage*	dst = TO_IMG(H);
	IplConvKernel*
				elem = NULL;
	int			size = kernel_size;
	int			offs = size/2;

	elem = cvCreateStructuringElementEx(size, size, offs, offs, CV_SHAPE_RECT, NULL);

	cvDilate(src, dst, elem, cnt);

	cvReleaseStructuringElement(&elem);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_MINMAX_ROI(LONG I, LONG x, LONG y, LONG w, LONG h, LONG* pmin, LONG* pmax)
{
	IplImage*	img = TO_IMG(I);
	LPBYTE		ps = (LPBYTE)img->imageData;
	LPBYTE		p;
	int			stride = img->widthStep;
	int			j, k;
	int			min = 255, max = 0;

	for (j = 0; j < h; j++, y++) {

		p = ps + x + y*stride;

		for (k = 0; k < w; k++, p++) {
			if (min > *p) {
				min = *p;
			}
			if (max < *p) {
				max = *p;
			}
		}
	}
	*pmin = min;
	*pmax = max;
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_AND(LONG I, LONG H, LONG J)
{
	IplImage*	src1 = TO_IMG(I);
	IplImage*	src2 = TO_IMG(H);
	IplImage*	dst  = TO_IMG(J);

	cvAnd(src1, src2, dst, NULL);
}
extern "C" __declspec(dllexport)
VOID APIENTRY OCV_SAVE(LONG I, LPCSTR file)
{
	IplImage*	I_IMG = TO_IMG(I);
	::cvSaveImage(file, I_IMG, NULL);
}
