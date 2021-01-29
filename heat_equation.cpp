#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

const double h = 0.005, tau = 0.0025,
			 time_0 = 10, x_0 = 1, alpha_coef = 1;
const int M = (int) (x_0 / h), N = (int) (time_0 / tau),
		  steps = 100, timefile = N < steps ? N : N / steps;

double y_max = 0, y_min = 10000;

double u(double t, double x) {
	//return exp(5.0 - M_PI * M_PI * t) * sin(x * M_PI);
	//return M_PI * x * tan(t) * t + 2 * cos(10 * t) * sin(5 * x);
	return (1 - cos(t)) * sin(2 * x);
}

double u_0(double x) {
	//return exp(5.0) * sin(x * M_PI);
	return 0;
	//return 2 * sin(5 * x);
}

double phi1(double t) {
	//return 0;
	//return 0; // 1
	return sin(t);
}

double phi2(double t) {
	//return 0;
	return 0;
	//return M_PI * tan(t) * t + 2 * cos(10 * t) * sin(5);
}

double psi1(double t) {
	//return t;
	return 2 * (1 - cos(t));
	//return M_PI * tan(t) * t + 10 * cos(10 * t);
}

double psi2(double t) {
	return 0;
	//return 	-exp(5.0 - M_PI * M_PI * t) * M_PI;
	//return M_PI * tan(t) * t + 10 * cos(10 * t) * cos(5);
}


double f(double t, double x) {
	//return 0;
	return 0;
	//return M_PI * x * tan(t) + M_PI * x * t /cos(t) / cos(t) + sin(5 * x) * (50 * cos(10 * t) - 20 * sin(10 * t));
}

void explicit_method(int mode, int left_border, int right_border, ofstream &out_file) {
	const double coef = tau * alpha_coef * alpha_coef / h / h;
	double *layer_1 = new double [(int) M + 1],
		   *layer_2 = new double [(int) M + 1],
		   *tmp;

	for (int i = 0; i <= M; i++) {
		layer_1[i] = u_0(i * h);
	}

	for (int i = 1; i <= N; i++) {
		for (int j = 1; j <= M - 1; j++) {
			layer_2[j] = layer_1[j] + tau * f(i * tau, j * h)
					   + coef * (layer_1[j + 1] - 2 * layer_1[j] + layer_1[j - 1]);
		}
		if (left_border == 1) {
			layer_2[0] = phi1(tau * i);
		} else if (left_border == 2) {
			layer_2[0] = (4 * layer_2[1] - layer_2[2] - psi1(tau * i) * 2 * h) / 3;
		} else {
			cout << "Left border error" << endl;
			return;
		}

		if (right_border == 1) {
			layer_2[M] = phi2(tau * i);
		} else if (right_border == 2) {
			layer_2[M] = (2 * h * psi2(tau * i)
						  + 4 * layer_2[M - 1]
						  - layer_2[M - 2])
					   / 3;
		} else {
			cout << "Right border error" << endl;
			return;
		}

		tmp = layer_2;
		layer_2 = layer_1;
		layer_1 = tmp;

		if ((i % timefile == 0) && mode == 2) {
			for (int j = 0; j <= M; j++) {
				out_file << j * h << ' ' << layer_1[j] << endl;
			}

			out_file << endl << endl;
		}

		for (int j = 0; j <= M; j++) {
			if (abs(layer_1[j]) > y_max) {
				y_max = abs(layer_1[j]);
			}

			if (layer_1[j] < y_min) {
				y_min = layer_1[j];
			}
		}
	}

	if (mode == 1) {
		double abs_err_c_h = 0, abs_err_l_2h = 0, u_abs_err_c_h = 0,
			   relative_err_c_h, relative_err_l_2h, u_abs_err_l_2h = 0;

		for (int i = 0; i <= M; i++) {
			if (abs(u(time_0, i * h) - layer_1[i]) > abs_err_c_h) {
				abs_err_c_h = abs(u(time_0, i * h) - layer_1[i]);
			}
		}

		for (int i = 0; i <= M; i++) {
			if (abs(u(time_0, i * h)) > u_abs_err_c_h) {
				u_abs_err_c_h = abs(u(time_0, i * h));
			}
		}

		for (int i = 0; i <= M; i++) {
			abs_err_l_2h += (u(time_0, i * h) - layer_1[i])
						 * (u(time_0, i * h) - layer_1[i]);
			u_abs_err_l_2h = u(time_0, i * h) * u(time_0, i * h);
		}

		abs_err_l_2h *= h;
		abs_err_l_2h = sqrt(abs_err_l_2h);

		u_abs_err_l_2h *= h;
		u_abs_err_l_2h = sqrt(u_abs_err_l_2h);

		relative_err_c_h = abs_err_c_h / u_abs_err_c_h;
		relative_err_l_2h = abs_err_l_2h / u_abs_err_l_2h;

		cout << "Absolute C_h error: " << abs_err_c_h << endl;
		cout << "Relative C_h error: " << relative_err_c_h << endl;
		cout << "Absolute l_2,h error: " << abs_err_l_2h << endl;
		cout << "Relative l_2,h error: " << relative_err_l_2h << endl;
	}

	delete[] layer_1;
	delete[] layer_2;

	return;
}

void implicit_method(int mode, int left_border, int right_border, ofstream &out_file) {
	const double B = alpha_coef * alpha_coef * tau / h / h,
				 C = (1 + 2 * B), A = B;
	double chi_1, mju_1, chi_2, mju_2,
		   alpha[M + 1], beta[M + 1], layer[M + 1];

	for (int i = 0; i <= M; i++) {
		layer[i] = u_0(i * h);
	}

	if (mode == 2) {
		for (int j = 0; j <= M; j++) {
			out_file << j * h << ' ' << layer[j] << endl;
		}
		out_file << endl << endl;
	}

	for (int i = 1; i <= N; i++) {
		if (left_border == 1) {
			chi_1 = 0;
			mju_1 = phi1(tau * i);
		} else if (left_border == 2) {
			chi_1 = 1;
			mju_1 = -h * psi1(tau * i);
		} else {
			cout << "Left border error" << endl;
			return;
		}

		if (right_border == 1) {
			chi_2 = 0;
			mju_2 = phi2(tau * i);
		} else if (right_border == 2) {
			chi_2 = 1;
			mju_2 = h * psi2(tau * i);
		} else {
			cout << "Right border error" << endl;
			return;
		}

		alpha[1] = chi_1;
		beta[1] = mju_1;

		for (int j = 2; j <= M; j++) {
			alpha[j] = B / (C - alpha[j - 1] * A);
			beta[j] = (
						layer[j - 1] + A * beta[j - 1]
						+ tau * f(i * tau, (j - 1) * h)
					  )
					/ (C - alpha[j - 1] * A);
		}

		if (right_border == 1) {
			layer[M] = phi2(tau * i);
		} else {
			layer[M] = (mju_2 + chi_2 * beta[M])
					 / (1 - chi_2 * alpha[M]);
		}

		for (int j = M - 1; j != -1; j--) {
			layer[j] = alpha[j + 1] * layer[j + 1] + beta[j + 1];
		}


		if ((i % timefile == 0) && mode == 2) {
			for (int j = 0; j <= M; j++) {
				out_file << j * h << ' ' << layer[j] << endl;
			}

			out_file << endl << endl;
		}

		for (int j = 0; j <= M; j++) {
			if (abs(layer[j]) > y_max) {
				y_max = abs(layer[j]);
			}

			if (layer[j] < y_min) {
				y_min = layer[j];
			}
		}
	}

	if (mode == 1) {
		double abs_err_c_h = 0, abs_err_l_2h = 0, u_abs_err_c_h = 0,
			   relative_err_c_h, relative_err_l_2h, u_abs_err_l_2h = 0;

		for (int i = 0; i <= M; i++) {
			if (abs(u(time_0, i * h) - layer[i]) > abs_err_c_h) {
				abs_err_c_h = abs(u(time_0, i * h) - layer[i]);
			}
		}

		for (int i = 0; i <= M; i++) {
			abs_err_l_2h += (u(time_0, i * h) - layer[i])
				  * (u(time_0, i * h) - layer[i]);
			u_abs_err_l_2h += u(time_0, i * h) * u(time_0, i * h);
		}

		abs_err_l_2h *= h;
		abs_err_l_2h = sqrt(abs_err_l_2h);

		u_abs_err_l_2h *= h;
		u_abs_err_l_2h = sqrt(u_abs_err_l_2h);

		for (int i = 0; i <= M; i++) {
			if (abs(u(time_0, i * h)) > u_abs_err_c_h) {
				u_abs_err_c_h = abs(u(time_0, i * h));
			}
		}

		relative_err_c_h = abs_err_c_h / u_abs_err_c_h;
		relative_err_l_2h = abs_err_l_2h / u_abs_err_l_2h;

		cout << "Absolute C_h error: " << abs_err_c_h << endl;
		cout << "Relative C_h error: " << relative_err_c_h << endl;
		cout << "Absolute l_2,h error: " << abs_err_l_2h << endl;
		cout << "Relative l_2,h error: " << relative_err_l_2h << endl;
	}
}

int main() {
	int mode, method, left_border, right_border;

	cout << "Enter an operating mode:" << endl;
	cin >> mode;
	cout << "Enter a method:" << endl;
	cin >> method;
	cout << "Enter a left border type: " << endl;
	cin >> left_border;
	cout << "Enter a right border type: " << endl;
	cin >> right_border;

	ofstream out_file("out.txt");

	if (method == 1) {
		explicit_method(mode, left_border, right_border, out_file);
	} else if (method == 2) {
		implicit_method(mode, left_border, right_border, out_file);
	} else {
		cout << "Method error" << endl;
		return 1;
	}

	ofstream omain("main.gn"),
			 oplotter("plotter.gn");
	omain << "set xrange [0:1]\nset yrange [" << (int) y_min - 1
		  << ": " <<  (int) y_max + 1 << "]\niter = 0\nload\"plotter.gn\"";
	oplotter << "iter = iter + 1\nplot \"out.txt\" i iter u 1:2 w l lt 6 notitle\npause 0.1\nif (iter < "
			 << steps << ") reread\n";
	out_file.close();
	return 0;
}
