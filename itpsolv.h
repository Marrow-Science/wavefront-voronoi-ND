#include <exception>
#include <sstream>

#ifndef ITPSOLV_H
#define ITPSOLV_H
namespace itp
{
	typedef struct
	{
		bool verified, safe;
		double lower, upper;
	} SOLVEBOUND;
	
	class BoundError : public std::exception
	{
		private:
		SOLVEBOUND bind;
		public:
		BoundError(SOLVEBOUND bnd) { bind = bnd; }
		const char* what()
		{
			std::ostringstream err;
			err << "Invalid boundary, not opposite sign";
			err << bind.lower << " " << bind.upper;
			return err.str().c_str();
		}
	};

	const double k1 = 0.1;
	const double k2 = 2.56567330897; // (0.98*(phi+1)), adjusted phi+1
	bool ensureSafe(double (*f)(double), SOLVEBOUND bound, double eps);
	double ITPsolv(double (*f)(double), SOLVEBOUND bound, double eps);
	double ITPsolvSafe(double (*f)(double),SOLVEBOUND bound, double eps);
};
#endif
