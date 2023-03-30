// TODO: remove
#include <stdio.h>

#include "itpsolv.h"
#include "geom.h"
#include <cmath>

bool itp::ensureSafe(double (*func)(double), SOLVEBOUND bound, double eps)
{
	// If already verified, return
	if(bound.verified) { return bound.safe; }
	double low = func(bound.lower);
	double up = func(bound.upper);
	// Check for zeroes on edges, do not verify in this edge case
	if(geom::eq(0.0,low,eps) || geom::eq(0.0,up,eps)) { return true; }
	// Otherwise check sign bits, verify
	bound.safe = std::signbit(low) != std::signbit(up);
	bound.verified = true;
	return bound.safe;
}
// Is safe, will do checks and throw errors if needed
double itp::ITPsolvSafe(double (*func)(double),SOLVEBOUND bound, double eps)
{
	if(!ensureSafe(func, bound, eps))
	{
		throw new itp::BoundError(
			{
				bound.verified,
				bound.safe,
				func(bound.lower),
				func(bound.upper)
			});
	}
	return ITPsolv(func, bound, eps);
}
// Use the Interpolate, Truncate, Project algorithm to solve the given function
// on the given interval. Is not safe, will not do checks
double itp::ITPsolv(double (*function)(double),SOLVEBOUND inbound, double eps)
{
	double window = inbound.upper - inbound.lower;
	double nbisection = std::log2(window/(2.0*eps));
	double nmax = nbisection + 1.0; // n0 is one, based on prev work
	int iteration = 0;
	SOLVEBOUND bound = inbound;
	while(bound.upper - bound.lower > 2*eps)
	{
		// Interpolation: bisection and regula falsi points
		double fa = function(bound.lower);
		double fb = function(bound.upper);
		double bisection = (bound.lower + bound.upper) / 2;
		double regulafalsi = (bound.upper*fa - bound.lower*fb)/(fa-fb);
		// Truncation, perturb towards center
		bool alph = std::signbit(bisection - regulafalsi);
		double d1 = k1*std::pow((bound.upper - bound.lower),k2);
		double d2 = std::abs(bisection - regulafalsi);
		double del = d1 < d2 ? d1 : d2;
		double xt = alph ? regulafalsi + del : regulafalsi - del;
		// Projection, project the estimator to new minmax interval
		double win = (bound.upper - bound.lower);
		double pr = eps*std::pow(2.0,nmax - iteration)-(win/2.0);
		double pf = std::abs(xt - bisection);
		double proj = pr < pf ? pr : pf;
		// Find the new bound!
		double xITP = alph ? bisection - proj : bisection + proj;
		double yITP = function(xITP);
		if(yITP > 0) { bound.upper = xITP; }
		else if(yITP < 0) { bound.lower = xITP; }
		else { bound.lower = xITP; bound.upper = xITP; }
		iteration++;
	}
	return (bound.upper + bound.lower) / 2.0;
}
