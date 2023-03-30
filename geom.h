#ifndef GEOM_H
#define GEOM_H
namespace geom
{
	double eq(double a ,double b, double eps)
	{
		return a > b ? (a-b) < eps : (b-a) < eps;
	}
}
#endif
