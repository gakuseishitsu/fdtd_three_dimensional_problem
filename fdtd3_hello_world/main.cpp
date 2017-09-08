#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cmath>

const double c = 2.9979246e8; //! speed of light
const double eps0 = 8.8541878e-12; //! permittivity of vacuum
const double mu0 = 1.256637e-6; //! magnetic permeability of vacuum
const double z0 = 376.73031; //! characteristic impedance of vacuum

const std::string filename = "../gnuplot_src/fdtd_data.dat";

const int nx0 = 120; //! x axis division of calc scope
const int ny0 = 120; //! y axis division of calc scope
const int nz0 = 120; //! z axis division of calc scope
const double dx = 0.005; //! x axis cell size 
const double dy = 0.005; //! y axis cell size
const double dz = 0.005; //! y axis cell size
const int nstep = 700; //! the number of steps of the calculation
const int pml_layer_length = 8;
const int pml_order = 4;
const double pml_accuracy = -120; //! dB
const double copml = -1.5280063e-4; // ln(10)/(40Z0)

const int nx = nx0 + 2 * pml_layer_length; //! area size of x axis
const int ny = ny0 + 2 * pml_layer_length; //! area size of y axis
const int nz = nz0 + 2 * pml_layer_length; //! area size of z axis

double ex[nx + 1][ny + 1][nz + 1]; //! electric field of x axis
double ey[nx + 1][ny + 1][nz + 1]; //! electric field of y axis
double ez[nx + 1][ny + 1][nz + 1]; //! electric field of z axis

double hx[nx + 1][ny + 1][nz + 1]; //! magnetic field of x axis
double hy[nx + 1][ny + 1][nz + 1]; //! magnetic field of x axis
double hz[nx + 1][ny + 1][nz + 1]; //! magnetic field of x axis

double aex[nx + 1][ny + 1][nz + 1];
double aey[nx + 1][ny + 1][nz + 1];
double aez[nx + 1][ny + 1][nz + 1];

double bex[nx + 1][ny + 1][nz + 1];
double bey[nx + 1][ny + 1][nz + 1];
double bez[nx + 1][ny + 1][nz + 1];

double amx[nx + 1][ny + 1][nz + 1];
double amy[nx + 1][ny + 1][nz + 1];
double amz[nx + 1][ny + 1][nz + 1];

double bmx[nx + 1][ny + 1][nz + 1];
double bmy[nx + 1][ny + 1][nz + 1];
double bmz[nx + 1][ny + 1][nz + 1];

double eps_field[nx + 1][ny + 1][nz + 1];
double sgm_e_field[nx + 1][ny + 1][nz + 1];
double mu_field[nx + 1][ny + 1][nz + 1];
double sgm_m_field[nx + 1][ny + 1][nz + 1];

const double eps_background = 1.0; //! relative permittivity
const double mu_background = 1.0; //! relative magnetic permeability 
const double sig_e_background = 0.0; //! electric conductivety
const double sig_m_background = 0.0; //! magnetic conductivety

//! configration of scatterers
const int ic = nx / 2; //! x position of center
const int jc = ny / 2; //! y position of center 
const int kc = nz / 2; //! z position of center
const int lx2 = 20; //! (x length of one side) /2
const int ly2 = 20; //! (y length of one side) /2
const int lz2 = 20; //! (z length of one side) /2
const double eps_scatterer = 3.0; //! relative permittivity of scatterers

//! configration of excitation current source
const double duration = 0.1e-9; //! pulse width
const double t0 = 4.0 * duration; //! peak time
const int ifed = ic - lx2 - 20; //! x feeding position
const int jfed = jc; //! y feeding position
const int kfed = kc; //! z feeding position
double befed; //! coefficient of soft feeding
double dl = 0.001; //! gap of hard feeding

const double v = c / sqrt(eps_background * mu_background);
const double dt = 0.99999 / (v*sqrt(1.0 / (dx*dx) + 1.0 / (dy*dy) + 1.0 / (dz*dz))); //! time step [sec]
double t = 0.0; //! current time

/*
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
*/

void setup(void);
void electric_field_step(void);
void magnetic_field_step(void);
void feed(void);

/*
void initpml(void);
void init_pml(pml p,int x0, int x1, int y0, int y1);
void epml(void);
void e_pml(pml p);
void hpml(void);
void h_pml(pml p);
*/

int main() {

	std::ofstream outputfile(filename);

	setup();
	//initpml();

	t = dt;

	for (int n = 1; n <= nstep; n++) {
		electric_field_step();
		//e_pml();
		feed();
		t += 0.5 * dt;
		magnetic_field_step();
		//h_pml();
		t += 0.5 * dt;

		if (n % 10 == 0) {
			for (int j = 0; j <= ny; j++) {
				for (int i = 0; i <= nx; i++) {
					outputfile << i*dx << " " << j*dy << " " << ez[i][j][nz/2] << std::endl;
				}
			}
			outputfile << std::endl << std::endl;
			std::cout << n << "/" << nstep << std::endl;
			std::cout << t/t0 << "," << ez[ifed-10][jfed][kfed] << std::endl;
		}
		
	}

	outputfile.close();

	return 0;
}

void feed(void) {
	double iz;

	//! soft feeding
	//iz = exp(-1 * pow((t-0.5*dt-t0)/duration,2));
	//ez[ifed][jfed][kfed] = ez[ifed][jfed][kfed] - befed*iz / (dx*dy);

	//! hard feeding
	ez[ifed][jfed][kfed] = exp(-1 * pow((t - t0) / duration, 2)) / dl;
}

void setup(void) {

	double eps, sgm, mu, sgmm, a;

	//! set the mudium of background
	for (int k = 0; k <= nz; k++) {
		for (int j = 0; j <= ny; j++) {
			for (int i = 0; i <= nx; i++) {
				eps_field[i][j][k] = eps_background;
				mu_field[i][j][k] = mu_background;
				sgm_e_field[i][j][k] = sig_e_background;
				sgm_m_field[i][j][k] = sig_m_background;
			}
		}
	}

	//! set scatterer
	for (int k = kc - lz2; k <= kc + lz2 - 1; k++) {
		for (int j = jc - ly2; j <= jc + ly2 - 1; j++) {
			for (int i = ic - lx2; i <= ic + lx2 - 1; i++) {
				eps_field[i][j][k] = eps_scatterer;
				mu_field[i][j][k] = 1.0;
				sgm_e_field[i][j][k] = 0.0;
				sgm_m_field[i][j][k] = 0.0;
			}
		}
	}

	//! set coefficient of excitation current source
	befed = dt / (eps0*0.25*(eps_field[ifed][jfed][kfed] + eps_field[ifed-1][jfed][kfed] + eps_field[ifed][jfed-1][kfed] + eps_field[ifed-1][jfed-1][kfed]));

	//! ŒW””z—ñ‚ÌŒˆ’è
	for (int k = 0; k <= nz; k++) {
		for (int j = 0; j <= ny; j++) {
			for (int i = 0; i <= nx; i++) {
				eps = 0.25 * (eps_field[i + 1][j + 1][k + 1] + eps_field[i + 1][j][k + 1] + eps_field[i + 1][j + 1][k] + eps_field[i + 1][j][k])*eps0;
				sgm = 0.25 * (sgm_e_field[i + 1][j + 1][k + 1] + sgm_e_field[i + 1][j][k + 1] + sgm_e_field[i + 1][j + 1][k] + sgm_e_field[i + 1][j][k]);
				a = 0.5 * sgm * dt / eps;
				aex[i][j][k] = (1.0 - a) / (1.0 + a);
				bex[i][j][k] = dt / eps / (1.0 + a);

				eps = 0.25 * (eps_field[i + 1][j + 1][k + 1] + eps_field[i][j + 1][k + 1] + eps_field[i + 1][j + 1][k] + eps_field[i][j + 1][k])*eps0;
				sgm = 0.25 * (sgm_e_field[i + 1][j + 1][k + 1] + sgm_e_field[i][j + 1][k + 1] + sgm_e_field[i + 1][j + 1][k] + sgm_e_field[i][j + 1][k]);
				a = 0.5 * sgm * dt / eps;
				aey[i][j][k] = (1.0 - a) / (1.0 + a);
				bey[i][j][k] = dt / eps / (1.0 + a);

				eps = 0.25 * (eps_field[i + 1][j + 1][k + 1] + eps_field[i][j + 1][k + 1] + eps_field[i + 1][j][k + 1] + eps_field[i][j][k + 1])*eps0;
				sgm = 0.25 * (sgm_e_field[i + 1][j + 1][k + 1] + sgm_e_field[i][j + 1][k + 1] + sgm_e_field[i + 1][j][k + 1] + sgm_e_field[i][j][k + 1]);
				a = 0.5 * sgm * dt / eps;
				aez[i][j][k] = (1.0 - a) / (1.0 + a);
				bez[i][j][k] = dt / eps / (1.0 + a);

				mu = 0.5 * (mu_field[i + 1][j + 1][k + 1] + mu_field[i][j + 1][k + 1])*mu0;
				sgmm = 0.5 * (sgm_m_field[i + 1][j + 1][k + 1] + sgm_m_field[i][j + 1][k + 1]);
				a = 0.5 * sgmm * dt / mu;
				amx[i][j][k] = (1.0 - a) / (1.0 + a);
				bmx[i][j][k] = dt / mu / (1.0 + a);

				mu = 0.5 * (mu_field[i + 1][j + 1][k + 1] + mu_field[i + 1][j][k + 1])*mu0;
				sgmm = 0.5 * (sgm_m_field[i + 1][j + 1][k + 1] + sgm_m_field[i + 1][j][k + 1]);
				a = 0.5 * sgmm * dt / mu;
				amy[i][j][k] = (1.0 - a) / (1.0 + a);
				bmy[i][j][k] = dt / mu / (1.0 + a);

				mu = 0.5 * (mu_field[i + 1][j + 1][k + 1] + mu_field[i + 1][j + 1][k])*mu0;
				sgmm = 0.5 * (sgm_m_field[i + 1][j + 1][k + 1] + sgm_m_field[i + 1][j + 1][k]);
				a = 0.5 * sgmm * dt / mu;
				amz[i][j][k] = (1.0 - a) / (1.0 + a);
				bmz[i][j][k] = dt / mu / (1.0 + a);
			}
		}
	}
}

void electric_field_step(void) {

	//! Ex
	for (int k = 1; k <= nz - 1; k++)
		for (int j = 1; j <= ny - 1; j++)
			for (int i = 0; i <= nx - 1; i++)
				ex[i][j][k] = aex[i][j][k] * ex[i][j][k] + bex[i][j][k] * ((hz[i][j][k] - hz[i][j - 1][k])/dy - (hy[i][j][k] - hy[i][j][k - 1])/dz);

	//! Ey
	for (int k = 1; k <= nz - 1; k++)
		for (int j = 0; j <= ny - 1; j++)
			for (int i = 1; i <= nx - 1; i++)
				ey[i][j][k] = aey[i][j][k] * ey[i][j][k] + bey[i][j][k] * ((hx[i][j][k] - hx[i][j][k - 1])/dz - (hz[i][j][k] - hz[i - 1][j][k])/dx);
	
	//! Ez
	for (int k = 0; k <= nz - 1; k++)
		for (int j = 1; j <= ny - 1; j++)
			for (int i = 1; i <= nx - 1; i++)
				ez[i][j][k] = aez[i][j][k] * ez[i][j][k] + bez[i][j][k] * ((hy[i][j][k] - hy[i - 1][j][k])/dx - (hx[i][j][k] - hx[i][j - 1][k])/dy);
}

void magnetic_field_step(void) {

	//! Hx
	for (int k = 0; k <= nz - 1; k++)
		for (int j = 0; j <= ny - 1; j++)
			for (int i = 1; i <= nx - 1; i++)
				hx[i][j][k] = amx[i][j][k] * hx[i][j][k] - bmx[i][j][k] * ((ez[i][j + 1][k] - ez[i][j][k])/dy - (ey[i][j][k + 1] - ey[i][j][k])/dz);

	//! Hy
	for (int k = 0; k <= nz - 1; k++)
		for (int j = 1; j <= ny - 1; j++)
			for (int i = 0; i <= nx - 1; i++)
				hy[i][j][k] = amy[i][j][k] * hy[i][j][k] - bmy[i][j][k] * ((ex[i][j][k + 1] - ex[i][j][k])/dz - (ez[i + 1][j][k] - ez[i][j][k])/dx);

	//! Hz
	for (int k = 1; k <= nz - 1; k++)
		for (int j = 0; j <= ny - 1; j++)
			for (int i = 0; i <= nx - 1; i++)
				hz[i][j][k] = amz[i][j][k] * hz[i][j][k] - bmz[i][j][k] * ((ey[i + 1][j][k] - ey[i][j][k])/dx - (ex[i][j + 1][k] - ex[i][j][k])/dy);
}

/*
void initpml(void) {
	init_pml(pml_l, 0, lpml, 0, ny);
	init_pml(pml_r, nx-lpml, nx, 0, ny);
	init_pml(pml_d, 0, nx, 0, lpml);
	init_pml(pml_u, 0, nx, ny-lpml, ny);
}

void init_pml(pml p, int x0, int x1, int y0, int y1){

	double smax0x, smax0y; //! x,y•ûŒü‚Ì“±“d—¦‚ÌÅ‘å’l
	double epspml, mupml; //! PML‚Ì”ä—U“d—¦, ”ä“§Ž¥—¦

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

*/
