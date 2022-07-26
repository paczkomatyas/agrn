#include "dv_tools.h"

#include <iostream>

namespace dvtools {
	int brokenStickVals(std::vector<double> &values, double random) {
               /* torott palca egyszeru bemeneti ertekekkel
                * visszateresi ertekek:
                *      0-(noChoices-1): melyiket valassztottuk
                *      -1: egyik sem (hiba)
                * values: a pointer a kumulalt ertekekhez
                * noChoices: hany tagja van a values-nak
                * random: a random szam
                * 
                * choice: szamlalo, vegigmegy a values-on
                * cumulate: ebbe kumulaljuk az ertkekeket
                */
               
               if(values.size() < 1) return(-3);

               double cumulate=0.0, sum=0.0;
               
               for(auto val = values.begin(); val != values.end(); val++) {
                       sum += *val;
               }
               
               if(sum <= 0) return -1;
               
               random *= sum;
               for (unsigned int choice = 0; choice < values.size(); choice++) {                
                       if (random <= (cumulate += values[choice]) ) {
//                             std::cout << cumulate << "/" << sum << ">" << random << ", so choice is: " << choice << std::endl;
                               return(choice);
                       }
//                     else std::cout << cumulate << "/" << sum << "<=" << random << ", so choice is NOT" << choice << std::endl;
               }
//             std::cerr << "brokenStickVals: some error, maybe sum (" << sum << ") is not valid" << std::endl;
//             for(choice=0; choice < noChoices; choice++) std::cerr << choice << ". choice: " << values[choice] << std::endl;

               return(-2);
       }
 
	
	int brokenStickVals(double *values, int noChoices, double sum, double random) {
		/* torott palca egyszeru bemeneti ertekekkel
		 * visszateresi ertekek:
		 * 	0-(noChoices-1): melyiket valassztottuk
		 * 	-1: egyik sem (hiba)
		 * values: a pointer a kumulalt ertekekhez
		 * noChoices: hany tagja van a values-nak
		 * random: a random szam
		 * 
		 * choice: szamlalo, vegigmegy a values-on
		 * cumulate: ebbe kumulaljuk az ertkekeket
		 */
		
		int choice=0;
		double cumulate=0.0;
		
		if(sum <= 0) {
			for(sum = choice = 0; choice < noChoices; choice++) {
				sum += *(values + choice);
//				std::cout << "calculating sum: + " << *(values + choice) << " = " << sum << std::endl;			
			}
			if(sum < 0) return -1;
		}
		

		random *= sum;
		
		for (choice = 0; choice < noChoices; choice++) {		
			if (random < (cumulate += *(values + choice))  ) {
//				std::cout << cumulate << "/" << sum << ">" << random << ", so choice is: " << choice << std::endl;
				return(choice);
			}
//			else std::cout << cumulate << "/" << sum << "<=" << random << ", so choice is NOT" << choice << std::endl;
		}
//		std::cerr << "brokenStickVals: some error, maybe sum (" << sum << ") is not valid" << std::endl;
//		for(choice=0; choice < noChoices; choice++) std::cerr << choice << ". choice: " << values[choice] << std::endl;

		return(-2);
	}
 

	//Constructor for quickPosVals
	quickPosVals::quickPosVals(){
		max=-1;
	}

	quickPosVals::quickPosVals(int _max, double (*_f)(int)) : max{_max}, f{_f}{
		if(max > -1) {
			vals = new double [max + 1];
			for(int i = 0; i <= max; i++){
				vals[i] = f(i);
			}
			last = vals[0];
		} 

	}

	//Copy Constructor
	quickPosVals::quickPosVals(const quickPosVals &obj){
		max = obj.max;
		f = obj.f;

		if(max > -1) {
			last = obj.last;
			vals = new double [max + 1];
			for(int i = 0; i <= max; i++){
				vals[i] = obj.vals[i];
			}
		}
	}

	//Destructor
	quickPosVals::~quickPosVals(){
		if(max > -1) delete [] (vals);
	}

	//calculate for the given range
	void quickPosVals::setMax(int newmax){
		if(max != newmax){
			if(max > -1) delete [] (vals);
			max = newmax;
			if (max > -1) {
				vals = new double (max + 1);
				for(int i = 0; i <= max; i++){
					vals[i] = f(i);
				}
			}
		}
	}

	void quickPosVals::setFunc( double (*_f)(int) ){
		f = _f;
	}

	//get max
	int quickPosVals::getMax() {return max;}


}
