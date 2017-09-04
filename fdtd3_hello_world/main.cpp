#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

const double c = 2.9979246e8; //! 光速
const double eps0 = 8.8541878e-12; //! 真空の誘電率
const double mu0 = 1.256637e-6; //! 真空の透磁率
const double z0 = 376.73031; //! 真空のインピーダンス

const std::string filename = "../gnuplot_src/fdtd_data.dat";

const int nx0 = 120; //! 解析領域x軸の分割数
const int ny0 = 120; //! 解析領域y軸の分割数
const double dx = 0.005; //! x軸セルサイズ
const double dy = 0.005; //! y軸セルサイズ
const int nstep = 100; //! 計算のステップ数
const int lpml = 8; //! PMLの総数
const int order = 4; //! PMLの次数
const double rmax = -120; //! PMLの要求精度dB
const double copml = -1.5280063e-4; // ln(10)/(40Z0)
const int nx = nx0 + 2 * lpml; //! x軸全分割数
const int ny = ny0 + 2 * lpml; //! y軸全分割数

std::vector<std::vector<double>> ex(nx + 1, std::vector<double>(ny + 1, 0.0)); //! x軸電界を記憶する配列
std::vector<std::vector<double>> ey(nx + 1, std::vector<double>(ny + 1, 0.0)); //! y軸電界を記憶する配列
std::vector<std::vector<double>> ez(nx + 1, std::vector<double>(ny + 1, 0.0)); //! z軸電界を記憶する配列
std::vector<std::vector<double>> hx(nx + 1, std::vector<double>(ny + 1, 0.0)); //! x軸磁界を記憶する配列
std::vector<std::vector<double>> hy(nx + 1, std::vector<double>(ny + 1, 0.0)); //! y軸磁界を記憶する配列
std::vector<std::vector<double>> hz(nx + 1, std::vector<double>(ny + 1, 0.0)); //! z軸磁界を記憶する配列

std::vector<std::vector<double>> aex(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> aey(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> aez(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列

std::vector<std::vector<double>> bexy(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> beyx(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> bezx(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> bezy(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列

std::vector<std::vector<double>> amx(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> amy(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> amz(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列

std::vector<std::vector<double>> bmxy(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> bmyx(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> bmzx(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列
std::vector<std::vector<double>> bmzy(nx + 1, std::vector<double>(ny + 1, 0.0)); //! 係数配列

std::vector<std::vector<double>> epsd(nx + 2, std::vector<double>(ny + 2, 0.0)); //! 比誘電率
std::vector<std::vector<double>> sgmed(nx + 2, std::vector<double>(ny + 2, 0.0)); //! 導電率
std::vector<std::vector<double>> mud(nx + 2, std::vector<double>(ny + 2, 0.0)); //! 比誘電率
std::vector<std::vector<double>> sgmmd(nx + 2, std::vector<double>(ny + 2, 0.0)); //! 導電率

const double epsbk = 1.0; //!背景誘電率
const double mubk = 1.0; //!背景透磁率
const double sigebk = 0.0; //!背景伝導率
const double sigmbk = 0.0; //!背景伝磁率

const int ic = nx / 2; //! 四角柱の中心x
const int jc = ny / 2; //! 四角柱の中心y
const int lx2 = 20; //! 四角柱1/2長さx
const int ly2 = 20; //! 四角柱1/2長さx
const double epsr = 3.0; //! 誘電体の比誘電率

const double duration = 0.1e-9; //! 励起電流源 パルス幅
const double t0 = 4.0*duration; //! 励起電流源 ピーク時刻
const int ifed = ic - lx2 - 20; //! 給電位置x
const int jfed = jc; //! 給電位置y
double befed; //! 係数

const double v = c / sqrt(epsbk*mubk);
const double dt = 0.99999 / (v*sqrt(1.0 / (dx*dx) + 1.0 / (dy*dy))); //! 時間ステップ [s]
double t; //! 時刻

class pml {
public:
	int i0, i1, j0, j1;
	std::vector<std::vector<double>> expml, eypml, ezx, ezy;
	std::vector<std::vector<double>> hxpml, hypml, hzx, hzy;
	std::vector<std::vector<double>> aexpml, aeypml;
	std::vector<std::vector<double>> beypml, bexpml;
	std::vector<std::vector<double>> amxpml, amypml;
	std::vector<std::vector<double>> bmypml, bmxpml;
};

pml pml_l, pml_r, pml_d, pml_u;

void setup(void);

void e_cal(void);
void h_cal(void);

void feed(void);

void initpml(void);
void init_pml(pml p,int x0, int x1, int y0, int y1);
void epml(void);
void e_pml(pml p);
void hpml(void);
void h_pml(pml p);

int main() {

	std::ofstream outputfile(filename);

	setup();
	initpml();

	t = dt;

	for (int n = 1; n <= nstep; n++) {
		e_cal();
		feed();
		epml();
		t += 0.5*dt;
		h_cal();
		hpml();
		t += 0.5*dt;

		
		if (n % 10 == 0) {
			for (int j = 0; j <= ny; j++) {
				for (int i = 0; i <= nx; i++) {
					outputfile << i*dx << " " << j*dy << " " << ez[i][j] << std::endl;
				}
			}
			outputfile << std::endl << std::endl;

			std::cout << n << "/" << nstep << std::endl;
		}
		
	}

	outputfile.close();

	return 0;
}

void feed(void) {
	double iz;

	iz = exp(-1 * pow((t-0.5*dt-t0)/duration,2));
	ez[ifed][jfed] = ez[ifed][jfed] - befed*iz / (dx*dy);
}

void setup(void) {

	double epsx, epsy, epsz;
	double sgex, sgey, sgez;
	double mux, muy, muz;
	double sgmx, sgmy, sgmz;
	double a;

	//! 背景媒質設定
	for (int j = 0; j <= ny; j++) {
		for (int i = 0; i <= nx; i++) {
			epsd[i][j] = epsbk;
			mud[i][j] = mubk;
			sgmed[i][j] = sigebk;
			sgmmd[i][j] = sigmbk;
		}
	}

	//! 誘電体設定
	for (int j = jc - ly2; j <= jc + ly2 - 1; j++) {
		for (int i = ic - lx2; i <= ic + lx2 - 1; i++) {
			epsd[i][j] = epsr;
			mud[i][j] = 1.0;
			sgmed[i][j] = 0.0;
			sgmmd[i][j] = 0.0;
		}
	}

	//! 電流源の係数決定
	befed = dt / (eps0*0.25*(epsd[ifed][jfed]+ epsd[ifed-1][jfed]+ epsd[ifed][jfed-1]+ epsd[ifed-1][jfed-1]));

	//! 係数配列の決定
	for (int j = 0; j <= ny; j++) {
		for (int i = 0; i <= nx; i++) {
			epsx = 0.5*(epsd[i+1][j+1] + epsd[i+1][j])*eps0;
			sgex = 0.5*(sgmed[i+1][j+1] + sgmed[i+1][j]);
			a = 0.5*sgex*dt / epsx;
			aex[i][j] = (1.0-a) / (1.0+a);
			bexy[i][j] = dt / epsx / (1.0 + a) / dy;

			epsy = 0.5*(epsd[i + 1][j + 1] + epsd[i][j+1])*eps0;
			sgey = 0.5*(sgmed[i + 1][j + 1] + sgmed[i][j+1]);
			a = 0.5*sgey*dt / epsy;
			aey[i][j] = (1.0 - a) / (1.0 + a);
			beyx[i][j] = dt / epsy / (1.0 + a) / dx;

			epsz = 0.25*(epsd[i + 1][j + 1] + epsd[i][j+1] + epsd[i+1][j] + epsd[i][j])*eps0;
			sgez = 0.25*(sgmed[i + 1][j + 1] + sgmed[i][j + 1] + sgmed[i+1][j] + sgmed[i][j]);
			a = 0.5*sgez*dt / epsz;
			aez[i][j] = (1.0 - a) / (1.0 + a);
			bezy[i][j] = dt / epsz / (1.0 + a) / dy;
			bezx[i][j] = dt / epsz / (1.0 + a) / dx;

			mux = 0.5*(mud[i + 1][j + 1] + mud[i + 1][j])*mu0;
			sgmx = 0.5*(sgmmd[i + 1][j + 1] + sgmmd[i + 1][j]);
			a = 0.5*sgmx*dt / mux;
			amx[i][j] = (1.0 - a) / (1.0 + a);
			bmxy[i][j] = dt / mux / (1.0 + a) / dy;

			muy = 0.5*(mud[i + 1][j + 1] + mud[i][j + 1])*mu0;
			sgmy = 0.5*(sgmmd[i + 1][j + 1] + sgmmd[i][j + 1]);
			a = 0.5*sgmy*dt / muy;
			amy[i][j] = (1.0 - a) / (1.0 + a);
			bmyx[i][j] = dt / muy / (1.0 + a) / dx;

			muz = 0.25*epsd[i + 1][j + 1]*mu0;
			sgmz = 0.25*sgmmd[i + 1][j + 1];
			a = 0.5*sgmz*dt / muz;
			amz[i][j] = (1.0 - a) / (1.0 + a);
			bmzx[i][j] = dt / muz / (1.0 + a) / dx;
			bmzy[i][j] = dt / muz / (1.0 + a) / dy;
		}
	}
}

void e_cal(void) {

	//! Ex
	for (int j = 1; j <= ny-1; j++)
		for (int i = 0; i <= nx-1; i++)
			ex[i][j] = aex[i][j] * ex[i][j] + bexy[i][j] * (hz[i][j] - hz[i][j-1]);

	//! Ey
	for (int j = 0; j <= ny - 1; j++)
		for (int i = 1; i <= nx - 1; i++)
			ey[i][j] = aey[i][j] * ey[i][j] - beyx[i][j] * (hz[i][j] - hz[i-1][j]);
	
	//! Ez
	for (int j = 1; j <= ny - 1; j++)
		for (int i = 1; i <= nx - 1; i++)
			ez[i][j] = aez[i][j] * ez[i][j] + bezx[i][j] * (hy[i][j] - hy[i - 1][j]) - bezy[i][j] * (hx[i][j] - hx[i][j-1]);
}

void h_cal(void) {

	//! Hx
	for (int j = 0; j <= ny - 1; j++)
		for (int i = 1; i <= nx - 1; i++)
			hx[i][j] = amx[i][j] * hx[i][j] - bmxy[i][j] * (ez[i][j+1] - ez[i][j]);

	//! Hy
	for (int j = 1; j <= ny - 1; j++)
		for (int i = 0; i <= nx - 1; i++)
			hy[i][j] = amy[i][j] * hy[i][j] + bmyx[i][j] * (ez[i+1][j] - ez[i][j]);

	//! Hz
	for (int j = 0; j <= ny - 1; j++)
		for (int i = 0; i <= nx - 1; i++)
			hz[i][j] = amz[i][j] * hz[i][j] - bmzx[i][j] * (ey[i+1][j] - ey[i][j]) + bmzy[i][j] * (ex[i][j+1] - ex[i][j]);
}

void initpml(void) {
	init_pml(pml_l, 0, lpml, 0, ny);
	init_pml(pml_r, nx-lpml, nx, 0, ny);
	init_pml(pml_d, 0, nx, 0, lpml);
	init_pml(pml_u, 0, nx, ny-lpml, ny);
}

void init_pml(pml p, int x0, int x1, int y0, int y1){

	double smax0x, smax0y; //! x,y方向の導電率の最大値
	double epspml, mupml; //! PMLの比誘電率, 比透磁率

	double sigmxm, sigmxe, sigmym, sigmye;
	double a;

	p.i0 = x0;
	p.i1 = x1;
	p.j0 = y0;
	p.j1 = y1;

	p.expml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.eypml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.ezx = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.ezy = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.hxpml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.hypml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.hzx = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.hzy = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));

	p.aeypml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.aexpml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.amypml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.amxpml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.beypml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.bexpml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.bmypml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));
	p.bmxpml = std::vector<std::vector<double>>(x1 - x0 + 1, std::vector<double>(y1 - y0 + 1, 0.0));

	smax0x = copml * rmax * (order + 1) / (lpml * dx);
	smax0y = copml * rmax * (order + 1) / (lpml * dy);

	mupml = mubk * mu0;
	epspml = epsbk * eps0;

	for (int i = x0; i <= x1; i++) {
		for (int j = y0; j <= y1; j++) {
			if (i < lpml) {
				sigmxm = pow((((double)(lpml - i) - 0.5) / (double)(lpml)), order)*smax0x;
				sigmxe = pow(((double)(lpml - i) / (double)(lpml)), order)*smax0x;
			}
			else if (i >= nx - lpml) {
				sigmxm = pow((((double)(i-nx+lpml) - 0.5) / (double)(lpml)), order)*smax0x;
				sigmxe = pow(((double)(i-nx+lpml) / (double)(lpml)), order)*smax0x;
			}
			else {
				sigmxm = 0.0;
				sigmxe = 0.0;
			}

			if (j < lpml) {
				sigmym = pow((((double)(lpml - j) - 0.5) / (double)(lpml)), order)*smax0y;
				sigmye = pow(((double)(lpml - j) / (double)(lpml)), order)*smax0y;
			}
			else if (j >= nx - lpml) {
				sigmym = pow((((double)(j - ny + lpml) - 0.5) / (double)(lpml)), order)*smax0y;
				sigmye = pow(((double)(j - ny + lpml) / (double)(lpml)), order)*smax0y;
			}
			else {
				sigmym = 0.0;
				sigmye = 0.0;
			}

			sigmxe = sigmxe*epsbk;
			a = 0.5*sigmxe*dt / epspml;
			p.aexpml[i - x0][j - y0] = (1.0 - a) / (1.0 + a);
			p.bexpml[i - x0][j - y0] = dt/epspml/(1.0+a)/dx;

			sigmye = sigmye*epsbk;
			a = 0.5*sigmye*dt / epspml;
			p.aeypml[i - x0][j - y0] = (1.0 - a) / (1.0 + a);
			p.beypml[i - x0][j - y0] = dt / epspml / (1.0 + a) / dy;

			sigmxm = sigmxm*epsbk;
			a = 0.5*sigmxm*dt / epspml;
			p.amxpml[i - x0][j - y0] = (1.0 - a) / (1.0 + a);
			p.bmxpml[i - x0][j - y0] = dt / mupml / (1.0 + a) / dx;

			sigmym = sigmym*epsbk;
			a = 0.5*sigmym*dt / epspml;
			p.amypml[i - x0][j - y0] = (1.0 - a) / (1.0 + a);
			p.bmypml[i - x0][j - y0] = dt / mupml / (1.0 + a) / dy;
		}
	}

}

void epml(void) {
	e_pml(pml_l);
	e_pml(pml_r);
	e_pml(pml_d);
	e_pml(pml_u);
}

void e_pml(pml p) {
	
	//! Ex
	for (int j = p.j0 + 1; j <= p.j1 - 1; j++) {
		for (int i = p.i0 + 1; i <= p.i1 - 1; i++) {
			p.expml[i - p.i0][j - p.j0] = p.aeypml[i - p.i0][j - p.j0] * p.expml[i - p.i0][j - p.j0] + p.beypml[i - p.i0][j - p.j0] * (hz[i][j] - hz[i][j-1]);
			ex[i][j] = p.expml[i - p.i0][j - p.j0];
		}
	}

	//! Ey
	for (int j = p.j0; j <= p.j1 - 1; j++) {
		for (int i = p.i0 + 1; i <= p.i1 - 1; i++) {
			p.eypml[i - p.i0][j - p.j0] = p.aexpml[i - p.i0][j - p.j0] * p.eypml[i - p.i0][j - p.j0] - p.bexpml[i - p.i0][j - p.j0] * (hz[i][j] - hz[i - 1][j]);
			ey[i][j] = p.eypml[i - p.i0][j - p.j0];
		}
	}

	//! Ez
	for (int j = p.j0 + 1; j <= p.j1 - 1; j++) {
		for (int i = p.i0 + 1; i <= p.i1 - 1; i++) {
			p.ezx[i - p.i0][j - p.j0] = p.aexpml[i - p.i0][j - p.j0] * p.ezx[i - p.i0][j - p.j0] + p.bexpml[i - p.i0][j - p.j0] * (hy[i][j] - hy[i - 1][j]);
			p.ezy[i - p.i0][j - p.j0] = p.aeypml[i - p.i0][j - p.j0] * p.ezy[i - p.i0][j - p.j0] - p.beypml[i - p.i0][j - p.j0] * (hx[i][j] - hx[i][j - 1]);
			ez[i][j] = p.ezx[i - p.i0][j - p.j0] + p.ezy[i - p.i0][j - p.j0];
		}
	}

}

void hpml(void) {
	h_pml(pml_l);
	h_pml(pml_r);
	h_pml(pml_d);
	h_pml(pml_u);
}

void h_pml(pml p) {

	//! Hx
	for (int j = p.j0; j <= p.j1 - 1; j++) {
		for (int i = p.i0 + 1; i <= p.i1 - 1; i++) {
			p.hxpml[i - p.i0][j - p.j0] = p.amypml[i - p.i0][j - p.j0] * p.hxpml[i - p.i0][j - p.j0] - p.bmypml[i - p.i0][j - p.j0] * (ez[i][j+1] - ez[i][j]);
			hx[i][j] = p.hxpml[i - p.i0][j - p.j0];
		}
	}

	//! Hy
	for (int j = p.j0+1; j <= p.j1 - 1; j++) {
		for (int i = p.i0; i <= p.i1 - 1; i++) {
			p.hypml[i - p.i0][j - p.j0] = p.amxpml[i - p.i0][j - p.j0] * p.hypml[i - p.i0][j - p.j0] + p.bmxpml[i - p.i0][j - p.j0] * (ez[i+1][j] - ez[i][j]);
			hy[i][j] = p.hypml[i - p.i0][j - p.j0];
		}
	}

	//! Hz
	for (int j = p.j0; j <= p.j1 - 1; j++) {
		for (int i = p.i0; i <= p.i1 - 1; i++) {
			p.hzx[i - p.i0][j - p.j0] = p.amxpml[i - p.i0][j - p.j0] * p.hzx[i - p.i0][j - p.j0] - p.bmxpml[i - p.i0][j - p.j0] * (ey[i+1][j] - ey[i][j]);
			p.hzy[i - p.i0][j - p.j0] = p.amypml[i - p.i0][j - p.j0] * p.hzy[i - p.i0][j - p.j0] + p.bmypml[i - p.i0][j - p.j0] * (ex[i][j+1] - ex[i][j]);
			hz[i][j] = p.hzx[i - p.i0][j - p.j0] + p.hzy[i - p.i0][j - p.j0];
		}
	}

}
