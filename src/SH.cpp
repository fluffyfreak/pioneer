#include <math.h>
#include "SH.h"

/**********************************************************************\
* AUTHOR : HILLAIRE Sébastien
*
* MAIL   : hillaire_sebastien@yahoo.fr
* SITE   : sebastien.hillaire.free.fr
*
*	You are free to totally or partially use this file/code.
* If you do, please credit me in your software or demo and leave this
* note.
*	Share your work and your ideas as much as possible!
\*********************************************************************/
//some code from Spherical Harmonic Lighting: The Gritty Details by Robin Green

namespace SH
{

	namespace
	{
		// Precomputed data used by the factorial
		const int nbPrecomputedFactorial = 16;							  // 0    1    2    3    4     5      6      7       8        9         10         11          12           13            14             15
		static const double precomputedFactorial[nbPrecomputedFactorial] = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0, 3628800.0, 39916800.0, 479001600.0, 6227020800.0, 87178291200.0, 1307674368000.0 };

		const double sqrt2 = 1.41421356237;
	}

	// Compute the factorial of x.
	// NB: if using float, errors occur at factorial(14) and greater values.
	// Therefore, use double precision if you plan to use SH with more than 3 bands.
	inline double Factorial(const int x)
	{
		if (x < nbPrecomputedFactorial)
			return precomputedFactorial[x];	//return precomputed value

		//return non precomputed value starting using the last precomputed value in the array
		double result = precomputedFactorial[nbPrecomputedFactorial - 1];
		for (int i = nbPrecomputedFactorial; i <= x; ++i)
			result = result *= i;
		return result;
	}

	// Evaluate an Associated Legendre polynomial at position x
	//
	// @param l : band
	// @param m : parameter
	// @param x : abscissa coordinate
	// @return the value at abscissa x
	double Polynomial(const int l, const int m, const double x)
	{
		// Generate the value of P(m, m) at x
		double pmm = 1.0;
		if (m > 0)
		{
			double sqrtOneMinusX2 = sqrt(1.0 - x * x);
			double fact = 1.0;
			for (int i = 1; i <= m; ++i)
			{
				pmm *= (-fact) * sqrtOneMinusX2;
				fact += 2.0;
			}
		}

		// If l==m, P(l, m)==P(m, m)
		if (l == m)
			return pmm;

		// Use rule 3 to calculate P(m+1, m) from P(m, m)
		double pmp1m = x * (2.0 * m + 1.0) * pmm;

		// If l==m+1, P(l, m)==P(m+1, m)
		if (l == m + 1)
			return pmp1m;

		// Otherwise, l>m+1.
		// Iterate rule 1 to get the result
		double plm = 0.0;
		for (int i = m + 2; i <= l; ++i)
		{
			plm = ((2.0 * i - 1.0) * x * pmp1m - (i + m - 1.0) * pmm) / (i - m);
			pmm = pmp1m;
			pmp1m = plm;
		}
		return plm;
	}

	//	Compute the scaling factor for Legendre polyniomal corresponding to band l and parameter m.
	//
	// @param l : band
	// @param m : parameter
	// @return the scale factor
	double PolynomialScalingFactor(const int l, const int m)
	{
		double temp = ((2.0*l + 1.0)*Factorial(l - m)) / ((4.0*M_PI)*Factorial(l + m));
		return sqrt(temp);
	}

	//Sample the SH func at position (theta, phi) with band l and parameter m.
	//
	// @param l : band
	// @param m : parameter
	// @param theta : spherical coordinate
	// @param phi : spherical coordinate
	// @return the function result
	double SampleSHFunc(const int l, const int m, const double theta, const double phi)
	{
		if (m == 0)
			return PolynomialScalingFactor(l, 0) * Polynomial(l, m, cos(theta));
		else if (m > 0)
			return sqrt2 * PolynomialScalingFactor(l, m) * cos(m * phi) * Polynomial(l, m, cos(theta));
		//m<0
		return sqrt2 * PolynomialScalingFactor(l, -m) * sin(-m * phi) * Polynomial(l, -m, cos(theta));
	}
}
