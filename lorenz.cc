#include <boost/numeric/ublas/vector.hpp>
#include "matplotlib.hpp"

namespace ub = boost::numeric::ublas;

template <class T, class F>
void rk(F f, ub::vector<T>& init, T start, T end) {

	ub::vector<T> x, dx1, dx2, dx3, dx4, x1, x2, x3;

	T t = start;
	T dt = end - start;

	x = init;
	dx1 = f(x, t) * dt;
	x1 = x + dx1 / 2.;
	dx2 = f(x1, t + dt/2.) * dt;
	x2 = x + dx2 / 2.;
	dx3 = f(x2, t + dt/2.) * dt;
	x3 = x + dx3;
	dx4 = f(x3, t + dt) * dt;

	init = x + dx1 / 6. + dx2 / 3. + dx3 / 3. + dx4 / 6.;
}

struct Lorenz {
	template <class T> ub::vector<T> operator() (const ub::vector<T>& x, T t){
		ub::vector<T> y(3);

		y(0) = 10. * ( x(1) - x(0) );
		y(1) = 28. * x(0) - x(1) - x(0) * x(2);
		y(2) = (-8./3.) * x(2) + x(0) * x(1);

		return y;
	}
};

int main()
{
	int i, j;
	ub::vector<double> x, x_new;
	double t, t_new, dt;
	const char *colors[3] = {"blue", "red", "green"};

	matplotlib g;

	g.open();

	g.screen(0., -30., 10., 50.);

	x.resize(3);
	x(0) = 15.; x(1) = 15.; x(2) = 36.;
	t = 0.;
	dt = pow(2., -5);

	for (i=0; i < (1 << 5) * 10; i++) {
		x_new = x;
		t_new = t + dt;
		rk(Lorenz(), x_new, t, t_new);
		for (j=0; j<3; j++) {
			g.line(t, x[j], t_new, x_new[j], colors[j]);
		}
		x = x_new;
		t = t_new;
	}

	getchar();
	g.close();
}
