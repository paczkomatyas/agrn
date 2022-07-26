#ifndef _EXPR_MOD_
#define _EXPR_MOD_

//add RCPP library to wrapper provided MACRO _EXPR_MOD_RCPP_
#ifdef _EXPR_MOD_RCPP_

#include <Rcpp.h>

#else

//#include "randomgen.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
extern "C" {
#include <time.h>
#include <stdio.h>
}

#endif 

#include <vector>
#include <string>
#include <bitset>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstring>
#include <sstream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>
#include "dv_tools.h"
#include <map>
#include <list>


namespace dv_expr {

	struct RegTo{
		RegTo(const unsigned int _to, const unsigned int _type, const unsigned int _trigger = -1, const bool is_fork = false): to(_to), trigger(_trigger), type(_type) {}

		const unsigned int to;
		const int trigger; //-1 denotes "I have no trigger" 
		const unsigned int type : 2; //0(00): hetero/autoasszoc; 1(01): conditional; 2(10): fork-default; 3(11): fork-triggered
	};

	struct RegFrom{
		RegFrom(const unsigned int _from, const unsigned int _to, const unsigned int _type = 0, const unsigned int _trigger = -1);
		const unsigned int from;
		unsigned int no_default;
		unsigned int no_conditional;
		bool autoassoc;
		std::vector< dv_expr::RegTo > tos; 

		void add(const unsigned int to, const unsigned int type = 0, const int trigger = -1);
	
	};

	struct RegUnit{
		RegUnit(const unsigned int _id_M, const unsigned int _from, const unsigned int _to, const unsigned int _type = 0, const unsigned int _trigger = -1): id_M(_id_M){
			add(_from, _to, _type, _trigger);
		}

		const unsigned int id_M;
		std::vector< dv_expr::RegFrom > froms;

		void add(const unsigned int from, const unsigned int to, const unsigned int type = 0, const int trigger = -1);
		std::vector<RegFrom>::iterator search(const unsigned int from);
	};

	class MMatrices{
		private:
			//stored values
			std::vector< double* > matrices; //Ms
			unsigned int last_size; //store last size to auto invalidate

			//refs / pointers to modells parameters
			class ExpressionModell *parent;
			unsigned int * size;
			double ** p;
			double ** M_actual; //to write over

		public:
			//constructor
			MMatrices(ExpressionModell *_parent); 

			//destructor
			~MMatrices(){
				for(auto m = matrices.begin(); m != matrices.end(); m++) delete [] (*m);
			}

			//copy constructor
			MMatrices( const MMatrices &obj); 

			//Functions
			//set values
			void setParent(ExpressionModell *_parent);

			//add new matrices
			int addM();

			//build M-s
			void force_build();
			void build(); //only if neccessary

			//set M to value
			void set(const int mid);
			void set_back();
			void set_front();

			//check sizes
			inline unsigned int no_M() {return matrices.size();}
			//inline unsigned int no_triggers() {return triggers.size();}
	};

	//class Action; //forward declaration

	class ExpressionModell;

	enum action_types {enone, eActionChangeM, eActionSwitch, eActionSet, eActionRun, eActionSetDecay};
	action_types action2enum(std::string input);

	class Action{
		public:
			Action() {};
			virtual ~Action() {};

			virtual void apply(ExpressionModell* obj) = 0;
			virtual Action* Clone() = 0;
			virtual void print(){ std::cout << "No info" << std::endl;}
	};

	class ActionChangeM : public Action{
		public:
			ActionChangeM(): no_target_M(0){}
			ActionChangeM(unsigned int _targetM): no_target_M(_targetM){}
			~ActionChangeM(){};

			void apply(ExpressionModell* obj);
			virtual Action* Clone(){return new ActionChangeM(*this); }
		private:
			unsigned int no_target_M; //the number of the appropriate M
	};

	class ExpressionModell{
	
		public:
			//exposed variables
			double t_output;

			//Constructor
			ExpressionModell();

			//Copy constructor
			ExpressionModell( const ExpressionModell &obj);
				
			//Destructor
			~ExpressionModell();

			class Trigger{
				public:
					Trigger(const boost::dynamic_bitset<> *_inner_trigger, double _timelag, double _treshold): 
						inner_trigger(_inner_trigger), 
						treshold(_treshold), 
						timelag(_timelag),
						time_elapsed(0.0){}

					Trigger(const boost::dynamic_bitset<> *_inner_trigger, double _treshold = 0.95): 
						inner_trigger(_inner_trigger), 
						treshold(_treshold), 
						timelag(0.0),
						time_elapsed(0.0){}

					Trigger( const Trigger &obj): 
						inner_trigger(obj.inner_trigger),
						treshold(obj.treshold),
						timelag(obj.timelag),
						time_elapsed(obj.time_elapsed) {
//							std::cout << "Trigger copy constructor called on " << *inner_trigger << "\n";
							// copy actions
							//actions.clear();
							for(auto act = obj.actions.begin(); act != obj.actions.end(); act++){
								actions.push_back( (*act)->Clone() );
							}
					}

					~Trigger(){
//						std::cout << "Trigger Destructor called on" << *inner_trigger <<  ", destroying actions...\n";
						for(auto act : actions) delete act;
					}


					bool operator== (const Trigger& comp){ //ignores actions and time_elapsed
//						std::cout << "Trigger: comparison of triggers\n";
						if(comp.inner_trigger == inner_trigger && comp.treshold == treshold && comp.timelag == timelag) return true;
						return false;
					}

					void addAction(Action* new_action){
//						std::cout << "Action added to trigger " << *inner_trigger << " ";
//						new_action->print();
						actions.push_back( new_action );
					}

					void apply(ExpressionModell* obj){
						for(auto act : actions) act->apply(obj);
					}
				//private:
					std::vector< Action* > actions;
					const boost::dynamic_bitset<> *inner_trigger; //the address of triggering state
					const double treshold; //switchng will be applied if pearson r value of p-inner_trigger exceeds this limit
					const double timelag; //given minimal time period the trigger value has to exceed the treshold, defaults to 0
					double time_elapsed; //how much time has elapsed dsince continously exceeding treshold
			};

			//Functions
			///initialisation
			int addPattern(const std::vector<int> &state);
			int addPattern(const std::string &state);

			int addPatternFromFile(std::ifstream &file);
			int addPatternFromFile(const std::string &filename);
			int addPatternFromFile(char *filename);
			
			void setDeath(double d); //inic it to a single value
			void setDeath(const std::vector<double> &d); //inic it to a vector
			void setDeath(int id, double val, char block); //set degradation of given gene to val
			void setDeath(const std::string &id, double val); //set degradation of given gene to val

			template <typename T>
			int addNewM(T id, const double timelag, const double treshold){
				auto itr = addInnerTrigger(id, timelag, treshold); 
				if(itr == inner_triggers.end()) return -1;
				else{
					int m = Ms.addM();
					Action *p = (Action*) new ActionChangeM(m);
					itr->addAction( p );
					return(m);
				}
			}

			template <typename T>
			int addNewM(T id, const double treshold = 0.95){
				return( addNewM(id, 0.0, treshold) );
			}
			
			int addNewM(){
				return Ms.addM();
			}
			
			//handle register

			template <typename Tfrom, typename Tto>
			void addToReg(Tfrom _from, Tto _to, std::string trigger, int id_M = 0){
				unsigned int trig_type = 0;
				//test from and to (in range)
				int from = asStateInt(_from), to = asStateInt(_to);
				if(from < 0 || to < 0) {
					std::cerr 
						<< "ERROR: ExpressionModell::addToReg: from(" 
						<< from 
						<< ") or to(" 
						<< to 
						<< ") are not in range of existing patterns (0 - " 
						<< states.size() 
						<< "). Please add pattern first!" << std::endl;
					return ;
				}

				//test trigger 
				//first character
				if(!trigger.empty()) switch(trigger[0]){
					default:
						std::cerr << "WARNING: addToReg: unrecognised first character of trigger value: " << trigger << ". Trigger ignored" << std::endl;
						break;
					case 'F':
						if(trigger.length() < 2) std::cerr << "WARNING: addToReg: too short trigger of type fork: " << trigger  << ". Trigger ignored" << std::endl;
						else { switch(trigger[1]){
								case '-':
									//if(trigger.length() < 3) {
									//	std::cerr << "WARNING: addToReg: too short trigger of type fork: " 
									//		<< trigger  << ". Trigger ignored" << std::endl;
									//} else {
										trig_type = 2; //10 - default branch of fork
										trigger.erase(0,2);
									//}
									break;

								case '+':
									if(trigger.length() < 3) {
										std::cerr << "WARNING: addToReg: too short trigger of type fork: " 
											<< trigger  << ". Trigger ignored" << std::endl;
									} else {
										trig_type = 3; //11 - triggered branch of fork
										trigger.erase(0,2);
									}
									break;
								default:
									trig_type = 2; //10 - default branch of fork
									trigger.erase(0,1);
							}
						}
						break;
					case 'C':
						if(trigger.length() < 2) std::cerr 
							<< "WARNING: addToReg: too short trigger of type conditinal: "
							<< trigger  << ". Trigger ignored" << std::endl;
						else{
							trig_type = 1; //01 - conditional
							trigger.erase(0,1);
						}
						break;
					case '0':
						break;
				}
				const int int_trigger = trig_type?asTriggerInt(trigger):-1;
				if(int_trigger < -1) return;
				

				// test id_M (<0)
				if( id_M < 0) {
					std::cerr 
						<< "ERROR: ExpressionModell::addToReg: id_M(" 
						<< id_M 
						<< ") is not valid (<0)!" << std::endl;
					return ;
				}

				//search for preexisting reg entry
				std::vector< RegUnit >::iterator m_it = searchRegM(id_M);
				if(m_it == reg.end()){ //it does not exist
					reg.emplace_back(id_M, from, to, trig_type, int_trigger);
				} else { //it exists
					m_it->add(from, to, trig_type, int_trigger);
				}
				
				//make sure endpoint has autoassotiative tag
				if(from != to) addToReg(to, to, std::string("0"), id_M);
			}


			/*template <typename Tfrom, typename Tto, typename Ttr = int>
			void addToReg(Tfrom _from, Tto _to, Ttr _trigger = -1, int id_M = 0){
				//test from and to (in range)
				int from = asStateInt(_from), to = asStateInt(_to);
				if(from < 0 || to < 0) {
					std::cerr 
						<< "ERROR: ExpressionModell::addToReg: from(" 
						<< from 
						<< ") or to(" 
						<< to 
						<< ") are not in range of existing patterns (0 - " 
						<< states.size() 
						<< "). Please add pattern first!" << std::endl;
					return ;
				}

				//test trigger 
				const int trigger = asTriggerInt(_trigger);
				if(trigger < -1) return;
				
				// test id_M (<0)
				if( id_M < 0) {
					std::cerr 
						<< "ERROR: ExpressionModell::addToReg: id_M(" 
						<< id_M 
						<< ") is not valid (<0)!" << std::endl;
					return ;
				}

				//search for preexisting reg entry
				std::vector< RegUnit >::iterator m_it = searchRegM(id_M);
				if(m_it == reg.end()){ //it does not exist
					reg.emplace_back(id_M, from, to, trigger);
				} else { //it exists
					m_it->add(from, to, trigger);
				}
			}*/


			bool buildFromReg();

			void printReg();


			//overwrite names
			void addPatternName(int id, char* name);
			void addPatternName(int id, std::string &name);
			//void add_gene_name(int id, char* name);
			//void add_gene_name(int id, std::string &name);
			//void add_factor_name(int id, char* name);
			//void add_factor_name(int id, std::string &name);
			//void add_trigger_name(int id, char* name);
			//void add_trigger_name(int id, std::string &name);
			
			///M matrix building blocks

			template <typename T1, typename T2>
			void heteroassoc(T1 id_from, T2 id_to){
				const boost::dynamic_bitset<> *pfrom = asState(id_from);
				const boost::dynamic_bitset<> *pto = asState(id_to);
				if(pfrom != NULL && pto != NULL){
					Ms.build();
					addToM( *pto, *pfrom );
					addToM( *pto, *pto);
				}
			}
			
			template <typename T>
			void homoassoc(T id){
				auto ptr = asState(id);
				if(ptr != NULL) {
					Ms.build();
					addToM(*ptr, *ptr);
				}
			}

			///M matrix user functions
			
			template <typename T>
			void transition(T id_from, T id_to){ //homoassoc + heteroasszoc
				heteroassoc(id_from, id_to);
			}

			void fork(int id_from, int id_to1, int id_to2, int trigger);
			bool _fork(int id_from, int id_to1, int id_to2, int trigger);
			void _heteroassoc(const boost::dynamic_bitset<> *from, const boost::dynamic_bitset<> *to);
			
			int fork(int id_from, int id_to1, int id_to2); //add new trigger automatically, return id of trigger
			
			void transition(int id_from, int id_to, int trigger);	//conditional transition, given trigger
			bool _conditional(int id_from, int id_to, int trigger);	//conditional transition, given trigger
			
			void condTransition(int id_from, int id_to, int trigger); //same as above
			int condTransition(int id);	//conditional transition, get automatically trigger, return its id

			///auto-assemble variables
			void addFactors();
			bool addFactor(const int id);
			bool hasFactor(const int id) const;
			int addTrigger(std::string label = "");
			void buildM();
			void perturbM(double prob, double sd);
			void switchM(int id);
			
			///simulation
			void run(double number_of_timesteps);
			void reset();

			bool setState(const std::vector<double> &state);
			void setState(const double value); //inic state to the same value overall
			bool setState(double mu, double sd); //inic state to random values from a normal distribution
			void setState(const int id); //inic current state to given pattern
			bool setState(const std::string &id, double value = -1); //inic current state to given pattern
			void setMatrixValue(const int row, const int col, const double value);
			void boostPattern(int state, double value);
			
			//triggers
			template <typename T>
			void switchIt(T _trigger){
				int trigger = asTriggerInt(_trigger);
				if(trigger < 0 || ((unsigned int) trigger) >= no_triggers){
					std::cerr << "WARNING! ExpressionModell::switchOn: invalid trigger: " << trigger << std::endl;
				} else {
					if(p[no_target_factors + no_genes + trigger]){ //trigger is not 0
						p[no_target_factors + no_genes + trigger] = 0;
					} else { //trigger is 0 
						p[no_target_factors + no_genes + trigger] = E;
					};
				}
			}

			template <typename T>
			void switchValue(T _trigger, const double value){
				int trigger = asTriggerInt(_trigger);
				if(trigger < 0 || ((unsigned int) trigger) >= no_triggers){
					std::cerr << "WARNING! ExpressionModell::switchOn: invalid trigger: " << trigger << std::endl;
				} else {
					p[no_target_factors + no_genes + trigger] = value;
				}
			}

			void switchOn(const int trigger);
			void switchOff(const int trigger);
			//void switchOnAtState(const int trigger, const int state, const double above_r = 0.95);
			template <typename T>
			void switchOnAtState(const int trigger, const T _state, const double timelag=0, const double above_r = 0.95){
				int state = asStateInt(_state);
				if(state < 0) {
					std::cerr 
						<< "ERROR: ExpressionModell::switchOnAtState: state(" 
						<< state 
						<< ") is not in range of existing patterns (0 - " 
						<< states.size() 
						<< "). Please add pattern first!" << std::endl;
					return ;
				}
				switchValueAtState(trigger, state, E, timelag, above_r);
			}

			void switchValueAtState(const int trigger, const int state, double value, const double above_r = 0.95);
			void switchValueAtState(const int trigger, const int state, double value, const double timelag, const double above_r );
			void switchOffAtState(const int trigger, const int state, const double above_r = 0.95);
			//void addMask(const std::string &mask_condition, const std::string &mask_target, double diff); 

			///outputting
			bool changeOutput(const char type);

			bool load(const std::string &filename);
			bool load(char *filename);

			bool loadActions(const std::string &filename);
			bool loadActions(char *filename);
			bool loadActions(std::ifstream &file);
			
			bool loadRules(char *filename);
			bool loadRules(const std::string &filename);
			bool loadRules(std::ifstream &file);

			bool saveRules(char *filename);
			void printPatterns() const;
			void printState() ;
			void printStatePearson() ;
			double getPearson(const int pattern_id) ;
			void printStateHeader() ;
			void printStatePearsonHeader() ;
			void printM() const;
			void printDeath() const;
			void printTriggers() const;
			void setSequence(const std::vector<int> &_seq); 
			void setSequenceTimes(const std::vector<double> &_times); 
			template <typename T>
			void setSequence2(const std::vector<T> &_seq){
				seq.clear();
				seq.reserve(_seq.size());
				for(auto s = _seq.begin(); s != _seq.end(); s++) seq.push_back( asStateInt(*s) );
			}
			double getRating() const;
			void outRanking(); //function for ranking
			void outRanking2(); //function for ranking
			void outRanking3(); //function for ranking
			void inicRanking(); //function for ranking
			std::vector<double> getStartVec(unsigned int no, double multiply);


#ifdef _EXPR_MOD_RCPP_
			Rcpp::DataFrame getRdf();
			Rcpp::DataFrame getPattern_df();
			Rcpp::NumericVector getState_Rvec();
#endif

			template <typename T>
			std::vector<Trigger>::iterator addInnerTrigger(T watch_id, const double timelag, const double treshold){   
				//check r value
				if(treshold * treshold > 1){
					std::cerr << "ERROR! ExpressionModell::addInnerTrigger: Pearson correlation r treshold value must be between 0 and 1!" << std::endl;
					return inner_triggers.end();
				}

				const boost::dynamic_bitset<> *watch = asState(watch_id);
				if(watch == NULL){
					std::cerr << "ERROR: ExpressionModell::addInnerTrigger: invalid watch id! " << watch_id << std::endl;
					return inner_triggers.end();
				}

				//create it
				Trigger new_tr(watch, timelag, treshold);

				//find out if it already has been registered
				auto tr_it = inner_triggers.begin();
				while(tr_it != inner_triggers.end() && !(new_tr == *tr_it) ) tr_it++;
				if(tr_it == inner_triggers.end()) {
					inner_triggers.push_back(new_tr); // it has not been issued so far, so push it back
					tr_it = --(inner_triggers.end());
				}

				return tr_it;
			}

			template <typename T>
			std::vector<Trigger>::iterator addInnerTrigger(T watch_id, const double treshold = 0.95){   
				const boost::dynamic_bitset<> *watch = asState(watch_id);
				if(watch == NULL){
					std::cerr << "ERROR: ExpressionModell::addInnerTrigger: invalid watch id! " << watch_id << std::endl;
					return inner_triggers.end();
				}

				//create it
				Trigger new_tr(watch, treshold);

//				std::cout << "addInnerTrigger temp trigger created with watch " << watch_id << "\n";

				//find out if it already has been registered
				auto tr_it = inner_triggers.begin();
				while(tr_it != inner_triggers.end() && !(new_tr == *tr_it) ) tr_it++;
//				std::cout << "addInnerTrigger comparison finished\n";
				if(tr_it == inner_triggers.end()) {
//					std::cout << "there is no thigger like this";	
					//inner_triggers.emplace_back(watch, treshold); // it has not been issued so far, so push it back
					inner_triggers.push_back(new_tr); // it has not been issued so far, so push it back
//					std::cout << "pushed back\n";
					tr_it = --(inner_triggers.end());
				}

//				std::cout << "addInnerTrigger finishes\n";
				return tr_it;
			}

			template <typename T>
			bool addInnerTriggerAction(Action *action, T watch_id, const double timelag, const double treshold){   
				auto itr = addInnerTrigger(watch_id, timelag, treshold);
				if (itr == inner_triggers.end()) return false;
				itr->addAction(action);
				return(true);
			}

			template <typename T>
			bool addInnerTriggerAction(Action *action, T watch_id, const double treshold = 0.95){   
				auto itr = addInnerTrigger(watch_id, treshold);
				if (itr == inner_triggers.end()) return false;
				itr->addAction(action);
				return(true);
			}

			template <typename T>
			bool addMTrigger(T watch_id, const int mid, const double timelag, const double treshold){   
				return addInnerTriggerAction( new ActionChangeM(mid), watch_id, timelag, treshold);
			}

			/*
			template <typename T>
			bool addMTrigger(T watch_id, const int mid = 0, const double treshold = 0.95){   
				const boost::dynamic_bitset<> *watch = asState(watch_id);
				if(watch == NULL){
					std::cerr << "ERROR: ExpressionModell::addTrigger: invalid watch id! " << watch_id << std::endl;
					return false;
				}

				Ms.addTrigger(watch, mid, treshold);

				return true;
			}*/

			//getters
			int get_no_genes() const;
			double get_time() const;
			int get_no_target_factors() const;
			int get_no_triggers() const;
			int get_size() const;
			int get_no_states() const;
			double get_tau() const;
			double get_death() const;
			std::vector<double> get_deathvec() const;
			double get_death(const unsigned int no) const;
			double get_shift() const;
			double get_E() const;
			double get_slope_tanh() const;
			double get_delta_t() const;

			// setters
			bool set_time(double _time);
#ifndef _EXPR_MOD_RCPP_
			void set_seed(unsigned long int seed);
#endif

			std::vector< dv_expr::RegUnit > reg; //to register M changes for better fool stability

		//private:
			friend class MMatrices;

			//modell parameters
			double tau;
			double tau_per2;
			double death;
			double shift;
			double E;
			double E2; //square of E
			double slope_tanh;
			double delta_t;
			std::vector<double> deathvec;
			bool use_death_vector;
#ifndef _EXPR_MOD_RCPP_
			gsl_rng * r;
#endif
			//other variables 
			double time;
			unsigned int size; //no_genes + no_target_factors + no_triggers = length ot state
			unsigned int no_genes; //number of "non-unique" genes
			unsigned int no_target_factors; //number of "unique" targeting factors -> they aid in finding right state
			unsigned int no_triggers; //number of triggers (at forks or conditional transitions)
			boost::dynamic_bitset<> mask_genes; //mask used to access/modify genes
			boost::dynamic_bitset<> mask_target_factors; //mask used to access/modify targeting factors
			boost::dynamic_bitset<> mask_triggers; //mask used to access/modify triggers
			
			//ranking
			std::vector<int> seq;
			std::vector<double> seq_times;
			std::vector<double> ranks;
			double rank;
			double seq_time_elapsed;
			int last_seq;

			//names
			std::vector<std::string> state_names; //string labels to states (e.g. juvenile - adultus, etc.)
			std::vector<std::string> gene_names; //string labels to genes (gene names)
			std::vector<std::string> trigger_names; //string labels to triggers (e.g. heat, adrenaline, etc.)
			std::vector<std::string> factor_names; //string labels to factors (e.g. miRNA1, etc.) 

			//model variables
			std::vector< boost::dynamic_bitset<> > states; //this stores the gene expression patterns of the landscape
			std::vector< Trigger > inner_triggers; //vector storing triggers
			MMatrices Ms; //masks for hierarchial masking
			double *p; //current state of modell
			double *M; //matrix describing transitions (its dimensions are M[size, size] )

			//function pointer used by run
			void (ExpressionModell::*inic_out)();
			void (ExpressionModell::*do_out)();
			void (ExpressionModell::*finish_out)();

			//FUNCTIONS
			
			///basic operations
			bool addToM(const boost::dynamic_bitset<> &to, boost::dynamic_bitset<> from);
			bool substractFromM(const boost::dynamic_bitset<> &to, boost::dynamic_bitset<> from);

			//internal memory allocation
			void invalidate();
			void update(); //it does not update except for p, checks or anything...
			void update_masks();
			double getPearson(const boost::dynamic_bitset<> &st) const;
			void checkInnerTriggers();
			void resetInnerTriggers();
			
			//get pointers to values
			const int asTriggerInt(const int id);
			const int asTriggerInt(const std::string &id) ;
			const int asStateInt(const int id) const;
			const int asStateInt(const std::string &id) ;
			const boost::dynamic_bitset<>* asState(const int id) const;
			boost::dynamic_bitset<>* asState(const std::string &id) ;

			//register handling
			std::vector< RegUnit >::iterator searchRegM(const unsigned int search);

			//defining varaibles/methods for Rcpp wrapper compatible output
#ifdef _EXPR_MOD_RCPP_
			std::vector< std::vector<double> > df_out;
			Rcpp::DataFrame df;
			void inicRPearsonOutput(); 
			void outRPearson();
			void inicRComposite(); 
			void outRComposite();
			void outRComposite2();
#endif
	};



	class ActionSwitch : public Action{
		public: 
			ActionSwitch(): target("0"), value(-1) {};
			ActionSwitch(const double _target, const double _value): target(std::to_string(_target) ), value(_value) {};
			ActionSwitch(const std::string & _target, const double _value): target(_target), value(_value) {};
			ActionSwitch(const double _target): target( std::to_string(_target) ), value(-1) {};
			ActionSwitch(const std::string & _target): target(_target), value(-1) {};
			~ActionSwitch(){
//				std::cout << "ActionSwitch destroyed" << std::endl;
			};

			virtual Action* Clone(){return new ActionSwitch(*this); }
			void apply(ExpressionModell* obj){
//				std::cout << "Action: switch target " << target << " to value " << value << std::endl;
				if(value < 0){ //value not set -> if it is on, then off or the other way around
					obj->switchIt(target);	
				} else { //value set
					obj->switchValue(target, value);
				}
			}
			void print(){
				std::cout << "ActionSwitch with target" << target << " and value " << value << std::endl;
			}
		private:
			const std::string target;  //which value will be altered in case of inner_triggering
			const double value;   //the assigned value
	};

	class ActionRun : public Action{
		public: 
			ActionRun(): length(-1) {};
			ActionRun(const double _length): length( _length ) {};
			~ActionRun(){};

			virtual Action* Clone(){return new ActionRun(*this); }
			void apply(ExpressionModell* obj){
				if(length < 0){ //value not set -> minimal
					obj->run(obj->delta_t);	
				} else { //value set
					obj->run(length);
				}
			}
			void print(){
				std::cout << "ActionRun: run for additional " << length << " time." << std::endl;
			}
		private:
			double length;   //the assigned value
	};

	class ActionSet : public Action{
		public: 
			ActionSet(): target("0"), value(-1) {};
			ActionSet(const double _target, const double _value): target(std::to_string(_target) ), value(_value) {};
			ActionSet(const std::string & _target, const double _value): target(_target), value(_value) {};
			ActionSet(const double _target): target( std::to_string(_target) ), value(-1) {};
			ActionSet(const std::string & _target): target(_target), value(-1) {};
			~ActionSet(){
//				std::cout << "ActionSet destroyed" << std::endl;
			};

			virtual Action* Clone(){return new ActionSet(*this); }
			void apply(ExpressionModell* obj){
//				std::cout << "Action: switch target " << target << " to value " << value << std::endl;
				obj->setState(target, value);	
			}
			void print(){
				std::cout << "ActionSet with target" << target << " and value " << value << std::endl;
			}
		private:
			const std::string target;  //which value will be altered in case of inner_triggering
			const double value;   //the assigned value
	};

	class ActionSetDecay : public Action{
		public: 
			ActionSetDecay(){};
			ActionSetDecay(const double _value){value.push_back(_value);};
			ActionSetDecay(const std::vector<double> &_value): value(_value) {};
			~ActionSetDecay(){};

			virtual Action* Clone(){return new ActionSetDecay(*this); }
			void apply(ExpressionModell* obj){
				if(value.empty()){
					obj->setDeath(0.2);
				} else{
					for(auto val = value.begin(); val != value.end(); val++) if(*val < 0) *val = 0.2; //set minuses to default value
					obj->setDeath(value);
				}
			}
			void print(){
				std::cout << "ActionSetDecay with value(s): ";
				for(auto val = value.begin(); val != value.end(); val++) std::cout << *val << " ";
				std::cout << std::endl;
			}
		private:
			std::vector<double> value;   //the assigned value
	};

#ifndef _EXPR_MOD_RCPP_
	ExpressionModell* MutateModell(const ExpressionModell &basis, const unsigned int startstate, const double length_of_sim = 300, const bool change_unique = true, const unsigned int no_steps = 100000, const unsigned int to_keep = 100, const double sd = 0.4);
#endif

}

//#include "r_wrapper.h"

#endif

