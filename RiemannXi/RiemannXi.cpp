
#include <iostream>

class Curve
{
	private:

		long double _pi;

		long double oneTerm(long double t, long double n);
		long double jensenG(long double t);
		long double simposen(long double y, long double t1, long double t2);
		long double func(long double t, long double y);

	public:
		Curve() { _pi = 4L * atanl(1L); }
		~Curve() {};
		long double xiFunc(long double y);
};


int main()
{
    std::cout << "Hello Riemann!\n";
	Curve c;

	const long num = 200;
	long double y[num + 1], u[num + 1];
	long double Ymin = 65L;
	long double Ymax = 85L;
	long double dy = (Ymax - Ymin) / num;

	char out[3000];
	std::string output = "";
	for (long i = 0; i <= num; i++)
	{
		y[i] = Ymin + i * dy;
		u[i] = c.xiFunc(y[i]);

		sprintf_s(out, "%10.5f\t%20.10e\n", y[i], u[i]);
		std::cout << "i = " << i << " y, u = " << out;
		output += out;
	}

	return 1;
}

long double
Curve::xiFunc(long double y)
{
	long double N = 10000;

	long double xi = 0L;
	for (long double i = 0; i < N; i++)
	{
		long double t1 = (i - 0.5L) * _pi / y;
		long double t2 = (i + 0.5L) * _pi / y;
		if (t1 < 0) t1 = 0;

		long double u = simposen(y, t1, t2);
		xi += u;
		if ((t2) >= 2L && fabs(u) < 1E-30) break;
	}

	return xi;
}


long double
Curve::simposen(long double y, long double t1, long double t2)
{
	long double TOR = 1E-30;

	long double u1 = func(t1, y);
	long double u2 = func(t2, y);

	long double t = 0.5L * (t1 + t2);
	long double sigma = func(t, y), tau = 0L;

	long double h = (t2 - t1) / 2L;
	long double value = (u1 + 4L * sigma + 2L * tau + u2) * h / 3L;

	long double n = 1L;
	while (n < 1000000L)
	{
		tau += sigma;
		sigma = 0L;

		h /= 2L;
		n *= 2L;

		for (long double i = 0L; i < n; i++)
		{
			t = t1 + (2L * i + 1L) * h;
			long double u = func(t, y);
			sigma += u;
		}

		long double val = value;
		value = (u1 + 4L * sigma + 2L * tau + u2) * h / 3L;

		long double err = fabs(value - val);
		if ((n >= 10000L) && err < TOR) break;
	}

	return value;
}

long double
Curve::func(long double t, long double y)
{
	long double f = jensenG(t) * cos(y * t);
	return f;
}

long double
Curve::jensenG(long double t)
{
	long double TOR = 1E-30;
	long double n1 = 1L, n2 = 10000L;

	long double value = 0L;
	for (long double n = n1; n <= n2; n++)
	{
		long double v = oneTerm(t, n);
		value += v;

		if (fabs(v) < TOR) break;
	}

	return value;
}

long double
Curve::oneTerm(long double t, long double n)
{
	long double alpha = _pi * n * n;
	long double exp4t = expl(4L * t);

	long double expfn = expl(9L * t - alpha * exp4t);
	long double other = 16L * alpha * alpha - 24L * alpha / exp4t;

	long double oneTerm = expfn * other;
	return oneTerm;
}


