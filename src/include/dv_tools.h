#ifndef _DV_TOOLS_
#define _DV_TOOLS_
#include <vector>
#include <cmath>

namespace dvtools {
	inline bool isInt(const std::string & string){
		const char firstchar = string[0];
		if(string.empty() || ( (!isdigit(firstchar)) && (firstchar != '-') && (firstchar != '+') ) ) return false;

		char * p;
		strtol(string.c_str(), &p, 10);

		return (*p == 0);
	}

	inline bool isPosInt(const std::string & string){
		const char firstchar = string[0];
		if(string.empty() || ( (!isdigit(firstchar)) && (firstchar != '+') ) ) return false;

		char * p;
		strtol(string.c_str(), &p, 10);

		return (*p == 0);
	}

	inline double cov(double n, double sum_xy, double sum_x, double sum_y ) {
		return( (sum_xy - sum_x*sum_y/n) / (n-1)  ); //Cov = \frac{\sum xy - \frac{1}{n}\sum x \sum y}{n-1}
	}
	
	inline double sd(int n, double sum, double sumsquared){
		return( (n>1)?std::sqrt(std::abs(n*sumsquared - sum*sum) / (n*(n-1))):0.0 ); //SD = \sqrt ( \frac {| n \sum x - \sum^2 x |}{n(n-1)} )
	}


	inline int Rmod(int a, int b) {
		int r = a%b; 
		return r >= 0 ? r : r + b;
	}

	inline double fracpart(double x){
		return x - (int) x;
	}

	int brokenStickVals(double *values, int noChoices, double sum, double random); 
        int brokenStickVals(std::vector<double> &values, double random);

	class quickPosVals {
		public:
			double last;

			//Constructor
			quickPosVals(int _max, double (*_f)(int));
			quickPosVals();

			//Copy Constructor
			quickPosVals(const quickPosVals &obj);

			//Destructor
			~quickPosVals();

			//indexing
			double &operator[](int i){
				if(i <= max) last = vals[i];
				else last = f(i);
				return last;
			}

			//give a new function
			void setFunc(double (*_f)(int) );

			//calculate for the given range
			void setMax(int newmax);

			//get max
			int getMax();

		private:
			double *vals;
			int max;
			double (*f)(int);
			
	};

}

#endif
