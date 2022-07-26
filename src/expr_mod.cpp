#include "expr_mod.h"
#include "include/expr_mod.h"

namespace dv_expr{

	//////////////////////
	// ActionChangeM /////
	//////////////////////
	
	action_types action2enum(std::string input){
		boost::algorithm::to_lower(input); 
		if(input == "actionchangem") return(eActionChangeM);
		if(input == "actionswitch") return(eActionSwitch);
		if(input == "actionset") return(eActionSet);
		if(input == "actionrun") return(eActionRun);
		if(input == "actionsetdecay") return(eActionSetDecay);
		return(enone);
	}
	
	void ActionChangeM::apply(ExpressionModell* obj){
		obj->switchM(no_target_M);
	}


	///////////////////////
	// RegUnit, RegFrom ///
	///////////////////////
	
	//constructor RegFrom
	RegFrom::RegFrom(const unsigned int _from, const unsigned int _to, const unsigned int _type, const unsigned int _trigger ):
		from(_from),
		no_default(0),
		no_conditional(0),
		autoassoc(false){
			add(_to, _type, _trigger);
		}

	/// searching
	std::vector<dv_expr::RegFrom>::iterator RegUnit::search(const unsigned int search){
		auto found = froms.begin();
		while(found != froms.end() && found->from != search) found++;
		return(found);
	}

	/// adding values
	void RegUnit::add(const unsigned int from, const unsigned int to, const unsigned int type, const int trigger ){
		std::vector< RegFrom >::iterator from_it = search(from);
		if(from_it == froms.end()){ //it does not exist
			froms.emplace_back(from, to, type, trigger);
		} else { //it exists
			from_it->add(to, type, trigger);
		}
	}

	/*void RegFrom::add(const unsigned int to, const int trigger ){
		if(from == to){ //homoassoc
			//warnings
			if(trigger >= 0) std::cerr 
				<< "WARNING: RegUnit: trigger (" 
				<< trigger
				<< ") added to homoassotiation of pattern "
				<< from
				<< ". It will be ignored!" << std::endl;
			if(no_homo) std::cerr
				<< "WARNING: RegUnit: there is more than one ("
				<< (no_homo + 1)
				<< ") homoassotiation given to pattern "
				<< from
				<< " , the model`s behaviour will be unpredictable!" << std::endl;
			
			//add data
			//tos.emplace_back(to, -2);
			no_homo++;
		} else { //heteroasszoc
			tos.emplace_back(to, trigger);
			if(trigger >= 0){ //... and there is a trigger
				no_conditional++;
				//issue warnings
				if(no_conditional > 1 && no_default) std::cerr
					<< "WARNING: RegUnit: additional conditinal transition added to node "
					<< from 
					<< " containing a fork already. To prevent unpredictable behaviour default branch will be ragarded as simple heteroassotiation and other branches as simple conditional transitions!" << std::endl;
			} else { //... and there is no trigger
				//issue warnings
				if(no_default) std::cerr
					<< "WARNING: RegUnit: more than one (" 
					<< (no_default + 1) 
					<< ") heteroassotiations given to pattern " 
					<< from 
					<< ". The model`s behaviour may become unpredictable!" << std::endl; 

				//do stuff
				no_default++; //and it is a default transition
			}
		} //heteroasszoc
	}*/

	void RegFrom::add(const unsigned int to, const unsigned int type, const int trigger ){
		if(from == to){ //homoassoc
			//warnings
			if(trigger >= 0) std::cerr 
				<< "WARNING: RegUnit: trigger (" 
				<< trigger
				<< ") added to homoassotiation of pattern "
				<< from
				<< ". It will be ignored!" << std::endl;
			
			//add data
			//tos.emplace_back(to, -2);
			autoassoc=true;
		} else { //heteroasszoc
			tos.emplace_back(to, type, trigger);
			if(type & 1){ //... and there is a trigger
				no_conditional++;
				//issue warnings
				if(no_conditional > 1 && no_default) std::cerr
					<< "WARNING: RegUnit: additional conditinal transition added to node "
					<< from 
					<< " containing a fork already. To prevent unpredictable behaviour default branch will be ragarded as simple heteroassotiation and other branches as simple conditional transitions!" << std::endl;
			} else if(!type){ //... and there is no trigger
				//issue warnings
				if(no_default) std::cerr
					<< "WARNING: RegUnit: more than one (" 
					<< (no_default + 1) 
					<< ") heteroassotiations given to pattern " 
					<< from 
					<< ". The model`s behaviour may become unpredictable!" << std::endl; 

				//do stuff
				no_default++; //and it is a default transition
			}
		} //heteroasszoc
	}

	//////////////////
	// MMatrices /////
	//////////////////
	
	MMatrices::MMatrices(ExpressionModell *_parent): 
		last_size(_parent->size),
		parent(_parent), 
		size( &(parent->size) ),
		p( &(parent->p) ),
		M_actual( &(parent->M) )
	{
		addM();
		*M_actual = NULL;
	}

	MMatrices::MMatrices( const MMatrices &obj): 
		last_size(obj.last_size), 
		parent(nullptr), 
		size(nullptr), 
		p(nullptr), 
		M_actual(nullptr)  {
		for(auto original_matrix = obj.matrices.begin(); original_matrix != obj.matrices.end(); original_matrix++) {
			const unsigned int msize = last_size * last_size;
			matrices.push_back(new double [msize]);
			double *new_matrix = matrices.back();

			for(unsigned int cell=0; cell < msize; cell++){
				new_matrix[cell] = (*original_matrix)[cell];
			}
		}
	}

	void MMatrices::setParent(ExpressionModell *_parent) {
		last_size = _parent->size;
		parent = _parent;
		size = &(parent->size) ;
		p = &(parent->p) ;
		M_actual = &(parent->M) ;
		if(matrices.empty()) {
			addM();
			*M_actual = NULL;
		}
		else *M_actual = matrices[0];
	}

	int MMatrices::addM(){
		build(); //check if size has been changed

		if(*size > 0) {
			int msize = *size;
			msize *= msize;

			//allocate
			double *new_M = new double [msize];
			
			//store it
			matrices.push_back( new_M );

			//null it
			while( msize-- ) {
				*new_M = 0.0;
				new_M++;
			}

		} else {
			//std::cerr << "WARNING: MMatrices::newM: size is null!" << std::endl;
			matrices.push_back(NULL);
		}

		return(matrices.size() -1 );
	}

	void MMatrices::force_build(){
		last_size = *size;
		const unsigned int msize = last_size * last_size;

		auto m_old = matrices.begin(); //have to set back M_actuel
		for(auto m = matrices.begin(); m != matrices.end(); m++){
			if(*m != NULL) {
				if(*M_actual == *m) m_old = m;
				delete [] (*m);
			}
			*m = new double [msize]; 
			for(unsigned int cell = 0; cell < msize; cell++) (*m)[cell] = 0.0;
		}
		*M_actual = *m_old;
	}

	void MMatrices::build(){
		if(*size != last_size) {
			force_build();
		}
	}

	void MMatrices::set(const int mid){
		build();
		if(mid < 0 || ((unsigned int ) mid) >= matrices.size()) {
			std::cerr << "WARNING: MMatrices::set: invalid id (" << mid << "). It should be between 0 and " << matrices.size() << std::endl;
			return ;
		}

		*M_actual = matrices[mid];
	}

	void MMatrices::set_back() {build(); *M_actual = matrices.back();}

	void MMatrices::set_front() {build(); *M_actual = matrices.front();}
	
	///////////////////////
	// ExpressionModell ///
	///////////////////////
	
	// Register methods
	/// building
	bool ExpressionModell::buildFromReg(){
		//create first M
		Ms.force_build();

		//check p
		unsigned int last_size= size;
		
		//build it
		for(auto reg_m = reg.begin(); reg_m != reg.end(); reg_m++){//loopoing tru M-s
			while(reg_m->id_M >= Ms.no_M()) Ms.addM(); //if not enough M -> make them
			Ms.set(reg_m->id_M); //set target

			for(auto cf = reg_m->froms.begin(); cf != reg_m->froms.end(); cf++) { //looping tru froms (cf: current from struct)
				std::vector<unsigned int> types{0, 0, 0, 0};
				for(auto ct = cf->tos.begin(); ct != cf->tos.end(); ct++) types[ct->type]++;

				//forks
				if(types[2] != types[3]){//not every fork is closed
					std::cerr << "ERROR: buildFromReg: not every fork is defined!" << std::endl;
					return (false); 
				} else if(types[2] == 1 ){ // it is a good fork
					int to_def=-1, to_trig=-1, fork_trigger = -1; 
//					std::cerr << "in node " << cf->from << " number of output edges:" << cf->tos.size() << std::endl;
					for(auto ct = cf->tos.begin(); ct != cf->tos.end() && (to_def < 0 || to_trig < 0); ct++) {
//						std::cerr << "of type " << ct->type << " to node " << ct->to << std::endl;
						switch(ct->type){
							case(2):
								to_def = ct->to;
								break;
							case(3):
								to_trig = ct->to; 
								fork_trigger = ct->trigger;
							default:
								break;
						}
					}
					if(to_def == to_trig || to_def < 0 || to_trig < 0 || fork_trigger < 0){
						std::cerr << "ERROR: buildFromReg: fork not defined properly (targets are equal or not found: " 
							<< to_def << " / " << to_trig << ", trigger: " << fork_trigger << ") from node " << cf->from << std::endl;
						return false;
					}
					_fork(cf->from, to_def, to_trig, fork_trigger);
				} else if(types[2]){ //too much forks
					std::cerr << "ERROR: buildFromReg: multiple forks (" << types[2] << ") from one node (" << cf->from << ") is prohibited!" << std::endl;
					return(false);
				}

				//hetero and conditional transitions
				if(types[0] + types[1]){
					if(types[0] > 1) {
						std::cerr << "ERROR: buildFromReg: from node " << cf->from << ": multiple heteroassociations is prohibited!" << std::endl;
						return false;
					}
					for(auto ct = cf->tos.begin(); ct != cf->tos.end() ; ct++) {
						switch(ct->type){
							case 0:
								_heteroassoc(asState(cf->from), asState(ct->to));
								break;
							case 1:
								_conditional(cf->from, ct->to, ct->trigger);
								break;
							default:
								break;
						}
					}
				}

				//autoassotiations
				if(cf->autoassoc || types[1] || types[3]) homoassoc(cf->from);
			}
		}

		Ms.set(0);

		if(last_size != size || p == NULL) setState(0.0);
		return true;
	}

	std::vector<RegUnit>::iterator ExpressionModell::searchRegM(const unsigned int search){
		auto found = reg.begin();
		for(; found != reg.end(); found++){
			if(found->id_M == search) {
				break;
			}
		}
		return found;
	}

	/// print
	void ExpressionModell::printTriggers() const {
		for(auto tr = inner_triggers.begin(); tr != inner_triggers.end(); tr++){
			std::cout << *(tr->inner_trigger) << " " << tr->actions.size() << std::endl;
			for(auto act = tr->actions.begin(); act != tr->actions.end(); act++) {
				(*act)->print();
			}
		}
	}
	void ExpressionModell::printReg(){
		for(auto obj = reg.begin(); obj != reg.end(); obj++) {
			for(auto from = obj->froms.begin(); from != obj->froms.end(); from++ ){
				for(auto to = from->tos.begin(); to != from->tos.end(); to++ ){
					std::cout << "transition from " << from->from << " to " << to->to << " in M " << obj->id_M;
					if(to->trigger >= 0) std::cout << " in case of trigger " << to->trigger;
					std::cout << std::endl;
				}	
			}
		}
	}
	//constuctor
	ExpressionModell::ExpressionModell(): 	t_output (0.1),
				tau(1.0),
				death(0.2),
				shift(0.05),
				slope_tanh(50),
				delta_t(0.1),
				use_death_vector(false),
				time(0.0), 
				size(0),
				no_genes(0),
				no_target_factors(0),
				no_triggers(0),
				Ms(this)
				//valid_M(false),
		{
			E = tau/death;
			E2 = E * E;
			tau_per2 = tau * 0.5;
			//Ms.emplace_back(1, (boost::dynamic_bitset<> *) NULL, 0);
			M = NULL;
			p = NULL;

			inic_out = &ExpressionModell::printStatePearsonHeader;
			do_out = &ExpressionModell::printStatePearson;

#ifndef _EXPR_MOD_RCPP_
			//initialise rng
			time_t timer;
			r = (gsl_rng *) gsl_rng_alloc (gsl_rng_mt19937);
			gsl_rng_set(r, std::time(&timer));
			//randomszam_inic(std::time(&timer) , r);
#endif
		}
	//copy constructor
	ExpressionModell::ExpressionModell( const ExpressionModell &obj): 	
				t_output (obj.t_output),
				tau(obj.tau),
				tau_per2(obj.tau_per2),
				death(obj.death),
				shift(obj.shift),
				E(obj.E),
				E2(obj.E2),
				slope_tanh(obj.slope_tanh),
				delta_t(obj.delta_t),
				deathvec(obj.deathvec),
				use_death_vector(obj.use_death_vector),
				time(0.0), 
				size(obj.size),
				no_genes(obj.no_genes),
				no_target_factors(obj.no_target_factors),
				no_triggers(obj.no_triggers),
				mask_genes(obj.mask_genes),
				mask_target_factors(obj.mask_target_factors),
				mask_triggers(obj.mask_triggers),
				seq(obj.seq),
				seq_times(obj.seq_times),
				ranks(obj.ranks),
				rank(0),
				seq_time_elapsed(0),
				last_seq(-1),
				state_names(obj.state_names),
				gene_names(obj.gene_names),
				trigger_names(obj.trigger_names),
				factor_names(obj.factor_names),
				inner_triggers(obj.inner_triggers),
				Ms(obj.Ms),
				p(NULL)
		{
			Ms.setParent(this);

			// copy p
			if(obj.p != NULL){
				std::vector<double> p_original;
				p_original.reserve(size);
				for(unsigned int pval = 0; pval < size; pval++) p_original.push_back(obj.p[pval]);
				setState(p_original);
			}

			// copy states
			for(auto ost = obj.states.begin(); ost != obj.states.end(); ost++) states.push_back(*ost);

			inic_out = &ExpressionModell::printStatePearsonHeader;
			do_out = &ExpressionModell::printStatePearson;

#ifndef _EXPR_MOD_RCPP_
			//initialise rng
			time_t timer;
			r = (gsl_rng *) gsl_rng_alloc (gsl_rng_mt19937);
			gsl_rng_set(r, std::time(&timer));
			//randomszam_inic(std::time(&timer) , r);
#endif
		}

	//destructor
	ExpressionModell::~ExpressionModell(){
		invalidate();

#ifndef _EXPR_MOD_RCPP_
		//close rng
		gsl_rng_free(r);
#endif
	}
	
	//checking and convert
	const int ExpressionModell::asStateInt(const int id) const{
		if(id < 0 || id >= (long int) states.size()) {
			std::cerr << "ERROR: ExpressionModell::asState: invalid id (" << id << ") returned -1" << std::endl;
			return -1;
		}
		return( id );
	}

	const int ExpressionModell::asStateInt(const std::string &id) {
		std::vector< std::string >::iterator found = std::find(state_names.begin(), state_names.end(), id);
		if(found != state_names.end()){
			return (const int) (found - state_names.begin()) ;
		} else {
			std::cerr << "ERROR: ExpressionModell::asState: invalid id (" << id << ") returned -1" << std::endl;
			return -1;
		}
	}
	
	const boost::dynamic_bitset<>* ExpressionModell::asState(const int id) const{
		if(id < 0 || id >= (long int) states.size()) {
			std::cerr << "ERROR: ExpressionModell::asState: invalid id (" << id << ") returned NULL" << std::endl;
			return NULL;
		}
		return( &(states[id]) );
	}

	boost::dynamic_bitset<>* ExpressionModell::asState(const std::string &id) {
		std::vector< std::string >::iterator found = std::find(state_names.begin(), state_names.end(), id);
		if(found != state_names.end()){
			return &(states[ (int) (found - state_names.begin()) ]);
		} else {
			std::cerr << "ERROR: ExpressionModell::asState: invalid id (" << id << ") returned NULL" << std::endl;
			return NULL;
		}
	}
	
	const int ExpressionModell::asTriggerInt(const int id){
		if(id < 1 ) {
			if(id != -1) std::cerr << "ERROR: ExpressionModell::asTrigger: invalid id (" << id << ") returned -1" << std::endl;
			return -1;
		}
		while( ((int) no_triggers) < id){ //have to create new triggers
			addTrigger();
		}
		return( id-1 );
	}

	const int ExpressionModell::asTriggerInt(const std::string &id) {
		if(dvtools::isPosInt(id)) {
			const int int_id = std::stoi(id);
			return( asTriggerInt(int_id) );

		} else {
			std::vector< std::string >::iterator found = std::find(trigger_names.begin(), trigger_names.end(), id);
			if(found != trigger_names.end()){
				return (const int) (found - trigger_names.begin()) ;
			} else {
				return addTrigger(id);
			}
		}
	}
	
	//death vector handling
	void ExpressionModell::setDeath(double d){
		death = d;
		use_death_vector = false;
	}

	void ExpressionModell::setDeath(const std::vector<double> &d){
		deathvec.clear();
		deathvec.assign(d.begin(), d.end());
		use_death_vector = true;
	}

	void ExpressionModell::setDeath(const std::string &id, double val){ //set degradation of given gene to val
		if(deathvec.size() > size) deathvec.clear();
		while(deathvec.size() < size){
			deathvec.push_back(death); //in case there has been a deathvec but shorter
		}

		//try to find id
		if( std::find(gene_names.begin(), gene_names.end(), id) != gene_names.end()) { //it is a gene name
			deathvec[ no_target_factors + (int) (std::find(gene_names.begin(), gene_names.end(), id) - gene_names.begin()) ] = val; 
		} else if( std::find(trigger_names.begin(), trigger_names.end(), id) != trigger_names.end()) { //it is a trigger name
			deathvec[ no_target_factors + no_genes + (int) (std::find(trigger_names.begin(), trigger_names.end(), id) - trigger_names.begin()) ] = val; 
		} else if( std::find(factor_names.begin(), factor_names.end(), id) != factor_names.end()) { //it is a factor name
			deathvec[ (int) (std::find(factor_names.begin(), factor_names.end(), id) - factor_names.begin()) ] = val; 
		} else {
			std::cerr << "WARNING: ExpressionModell::setDeath: invalid id " << id << std::endl;
			return ;
		}

		use_death_vector = true;
	}

	void ExpressionModell::setDeath(int id, double val, char block = 'g'){
		if(deathvec.size() > size) deathvec.clear();
		while(deathvec.size() < size){
			deathvec.push_back(death); //in case there has been a deathvec but shorter
		}
		switch(block){
			case 'a':
				if(id < 0 || ((unsigned int) id) >= size) {
					std::cerr << "WARNING: ExpressionModell::setDeath: invalid id (" << id << ") given for block " << block << std::endl;
					return ;
				}
				deathvec[id] = val;
				break;
			case 'g':
				if(id < 0 || ((unsigned int) id) >= no_genes) {
					std::cerr << "WARNING: ExpressionModell::setDeath: invalid id (" << id << ") given for block " << block << std::endl;
					return ;
				}
				deathvec[no_target_factors + id] = val;
				break;
			case 'G':
				if(id < 0 || ((unsigned int) id) >= (no_genes + no_target_factors) ) {
					std::cerr << "WARNING: ExpressionModell::setDeath: invalid id (" << id << ") given for block " << block << std::endl;
					return ;
				}
				deathvec[id] = val;
				break;
			case 't':
				if(id < 0 || ((unsigned int) id) >= no_triggers) {
					std::cerr << "WARNING: ExpressionModell::setDeath: invalid id (" << id << ") given for block " << block << std::endl;
					return ;
				}
				deathvec[no_target_factors + no_genes + id] = val;
				break;
			case 'f':
				if(id < 0 || ((unsigned int) id) >= no_target_factors) {
					std::cerr << "WARNING: ExpressionModell::setDeath: invalid id (" << id << ") given for block " << block << std::endl;
					return ;
				}
				deathvec[id] = val;
				break;
			default:
				std::cerr << "WARNING: ExpressionModell::setDeath: invalid block: " << block << ". Expected values: a(absolute position), g(genes), t(triggers), f(target factors)" << std::endl;
				return ;
		}
		use_death_vector = true;
	}
	
	//adding values to modell
	void ExpressionModell::perturbM(double prob, double sd){
		int size2 = size*size;
		for(int m_counter = 0; m_counter < size2; m_counter++){
#ifdef _EXPR_MOD_RCPP_
			if(R::runif(0,1) < prob) {
				//M[m_counter] += R::rnorm(0,sd) ;
				M[m_counter] *= std::abs(R::rnorm(1, sd)) ;
				//M[m_counter] *= R::rnorm(1, sd) ;
			}
#else
			if(gsl_rng_uniform(r) < prob) {
				// pertrb with a normal distribution (negative values considered nulling)
				double rand = gsl_ran_gaussian(r, sd)+1 ;
				M[m_counter] *= (rand<0)?0:rand ;

				// perturb with a normal distribution (absolute values)
				//M[m_counter] *= std::abs( gsl_ran_gaussian(r, sd)+1 ) ;
				
				// perturb with a normal distribution
				//M[m_counter] *= gsl_ran_gaussian(r, sd)+1 ;

				// perturb with nulling elements
				//M[m_counter] = 0.0 ;
			}
#endif
		}
	}

	void ExpressionModell::boostPattern(int state, double value){
		auto st= states[state];
		for(unsigned int gene_no = 0; gene_no < size; gene_no++){
			 if(st[gene_no]) {
				 //std::cout << "gene " << gene_no << " changed" << std::endl;
				 p[gene_no] += value;
			 }
		}
	}

	int ExpressionModell::addTrigger(std::string label){
		//add it
		for(auto st = states.begin(); st != states.end(); st++){
			st->push_back(false);
		}
		no_triggers++;
		size++;
		invalidate();

		//add label
		if(label == ""){ //no label set
			trigger_names.push_back("tr" + std::to_string(size-1) );
		} else {
			trigger_names.push_back(label);
		}

		update_masks();

		return(no_triggers-1);
	}

	void ExpressionModell::addFactors(){
		int no_states = states.size(), no_new_factors = no_states - no_target_factors;

		//add factors to all the states
		for(unsigned int no_st = 0; no_st < (unsigned int) no_states; no_st++){
			//add places to states
			for(int i = no_new_factors; i--;){
				states[no_st].push_back(false);
			}

			//shift them
			states[no_st] <<= no_states;
			
			//add unique identifier
			if(no_st >= no_target_factors) states[no_st][no_st-no_target_factors] = true;
		}

		//add new names
		for(int i = 0; i < no_new_factors; i++) factor_names.push_back("f" + std::to_string(i));
		
		no_target_factors +=  no_new_factors;
		size +=  no_new_factors;
		invalidate();
		update_masks();

	}

	int ExpressionModell::addPattern(const std::vector<int> &state){
		//check if new state is different in size than previous ones
		if(states.size()){
			if(state.size() != no_genes ){
				std::cerr << "WARNING: state added to ExpressionModell has different size than previous element(s)!" << std::endl;
				return (-1);
			}
		} else {
			size = no_genes = state.size(); //in case it is the first set number of genes
			for(unsigned int newgenename=0; newgenename < no_genes; newgenename++) gene_names.push_back( "g" + std::to_string(newgenename)); //give name to genes
			invalidate(); //invalidate M, as dimensions have changed
		}

		//create new state
		states.emplace_back(size, 0 );
		auto current = states.end();

		for(unsigned int st = 0; st < no_genes; st++){
			(*current)[st] = (bool) state[st];
		}

		//shift to make place for target factors
		if(no_target_factors){
			*current <<= no_target_factors; // shifted left to make place unique target factors (current not included!!)
		}

		//give name
		state_names.push_back("state" + std::to_string(states.size()));

		update_masks();

		return( (int) states.size() -1);
		
	}

	int ExpressionModell::addPattern(const std::string &state){
		//check if new state is different in size than previous ones
		if(states.size()){
			if(state.size() != no_genes ){
				std::cerr << "WARNING: state added to ExpressionModell has different size than previous element(s)!" << std::endl;
				return (-1);
			}
		} else {
			size = no_genes = state.size(); //in case it is the first set number of genes
			for(unsigned int newgenename=0; newgenename < no_genes; newgenename++) gene_names.push_back( "g" + std::to_string(newgenename)); //give name to genes
			invalidate(); //invalidate M, as dimensions have changed
			
		}

		//create new state
		states.emplace_back(state) ;

		if(no_target_factors || no_triggers){ //resize its
			states.end()->resize(size, false);
		}

		//shift to make place for target factors
		if(no_target_factors){
			*states.end() <<= no_target_factors; // shifted left to make place unique target factors (current not included!!)
		}

		//give name
		state_names.push_back("state" + std::to_string(states.size()));

		update_masks();

		return( (int) states.size() -1);
	}
		
	int ExpressionModell::addPatternFromFile(std::ifstream &file){
		int no_new_genes = 0, no_new_states=0;
		std::string line, word;
		
		//check for first line: is it a header?
		if(file.peek() == '%'){ //it is - than jump!
			getline(file, line);
		}

		while(file.peek() == '#') getline(file, line); //step tru coments between section header and header

		//read header
		if( getline(file, line) ){
			std::istringstream linestream(line);

			//add state names
			linestream >> word; //ignore first word
			while(linestream >> word){
				no_new_states++;

				if(word.size()){
					state_names.push_back(word);
				} else {
					state_names.push_back("state" + std::to_string(no_new_states-1) );
				}
			}
		} else {
			std::cerr << "WARNING: ExpressionModell::addPatternFromFile: file is empty!" << std::endl;
			return( -2 );
		}

		//allocate states
		for(int no_st = no_new_states; no_st--;) states.emplace_back();

		//read in rest of file (containing gene names and expression patterns)
		while( file.peek() != '%' && getline(file, line)){
			if(line.size() == 0 || line[0] == '#') break; //ignore empty lines and lines starting with #

			std::istringstream linestream(line);
			bool genestate = false;

			no_new_genes++;

			//add gene name
			linestream >> word;
			if(word.size()){
				gene_names.push_back(word);
			} else {
				gene_names.push_back("g" + std::to_string(no_new_genes-1) );

			}

			//read in expression of gene
			for(int no_st=0; no_st < no_new_states; no_st++){
				linestream >> genestate;
				//if(word == ) 
				states[no_st].push_back(genestate);
			}
		}

		no_genes += no_new_genes;
		size += no_new_genes;

		update_masks();

		// closing remarks
		if(!no_target_factors) addFactors();
		if(p == NULL || no_new_genes) setState(0.0);
		return(states.size() -1);
	}

	int ExpressionModell::addPatternFromFile(const std::string &filename){
		return (addPatternFromFile( (char *) filename.c_str()));
	}

	int ExpressionModell::addPatternFromFile(char *filename){
		std::ifstream file(filename);

		if(!file.is_open()) {
			std::cerr << "ERROR: ExpressionModell::addPatternFromFile: file (" << filename << ") not found!" << std::endl;
			return( -1 );
		}

		return addPatternFromFile(file);
	}
	
	bool ExpressionModell::loadRules(std::ifstream &file){
		std::string from, to, line, trigger;
		int matrix;
		
		//check for first line: is it a header?
		if(file.peek() == '%'){ //it is - than jump!
			getline(file, line);
		}

		while(file.peek() == '#') getline(file, line); //step tru coments between section header and header

		//read header
		getline(file, line);

		//read in rest of file (containing: from to trigger matrix)
		while(file.peek() != '%' && getline(file, line)){
			if(line.size() == 0 || line[0] == '#') continue; //ignore empty lines and lines starting with #

			std::istringstream linestream(line);

			//read in data
			linestream >> from;	//string
			if(linestream.eof()) {std::cerr << "ERROR: loadRules: line should contain at least two words!" << std::endl; continue;}
			linestream >> to;	//string
			if(linestream.eof()) trigger = "";
			else linestream >> trigger;	//string
			if(linestream.eof()) matrix = 0;
			else linestream >> matrix;	//int

			//add data
			addToReg(from, to, trigger, matrix);
			if(file.peek() == '%') {
				if(reg.size()) buildFromReg();
				return true;
			}
		}

		//closing remarks
		if(reg.size()) buildFromReg();
		return(true);

	}

	bool ExpressionModell::loadRules(char *filename){
		std::ifstream file(filename);

		if(!file.is_open()) {
			std::cerr << "ERROR: ExpressionModell::loadRules: file (" << filename << ") not found!" << std::endl;
			return( false );
		}
	
		return(loadRules(file));
	}

	bool ExpressionModell::loadRules(const std::string &filename){
		return (loadRules( (char *) filename.c_str()));
	}

	bool ExpressionModell::load(const std::string &filename){
		return load( (char *) filename.c_str());
	}
	
	bool ExpressionModell::load(char *filename){
		std::ifstream file(filename);
		std::string line;
		
		std::vector<unsigned int> datas, topologies, settings; 

		if(!file.is_open()) {
			std::cerr << "ERROR: ExpressionModell::load: file (" << filename << ") not found!" << std::endl;
			return( false );
		}
		
		//read in rest of file (containing: from to trigger matrix)
		while( getline(file, line)){
			//if(line.size() == 0 || line[0] == '#') continue; //ignore empty lines and lines starting with #
			if(line[0] == '%' && line.size() != 0){ //if it seems a section header
				std::istringstream linestream(line);
				std::string sectiontype;
				linestream >> sectiontype; //jump the %
				linestream >> sectiontype; //read in second word

				if( sectiontype == "data" ){
					//addPatternFromFile(file);
					datas.push_back(file.tellg());
					//break;
				} else if ( sectiontype == "topology" ){
					//loadRules(file);
					topologies.push_back(file.tellg());
				} else if ( sectiontype == "settings" ){
					//loadActions(file);
					settings.push_back(file.tellg());
				} else {
					std::cerr << "WARNING: load: unknown section header (" << sectiontype << ")! Please choose from the followings: data, topology or settings. Section ignored." << std::endl;
				}
			}
		}

		//loading stuff in right order 
		switch(datas.size()){
			case 0: 
				std::cerr << "ERROR: load: no data were found. The use of single state representation models with biuldFromReg() is not yet supported. The program will crash." << std::endl;
				//return false;
				break;
			default:
				std::cerr << "Warning: load: more than one data sections were found. For safety reasons you should arrange your data to one unique block. Nevermind, trying to read them..." << std::endl;
			case 1:
				file.clear();
				for(auto pos = datas.begin(); pos != datas.end(); pos++) {
					file.seekg(*pos);
					addPatternFromFile(file);
				}
		}

		switch(topologies.size()){
			case 0: 
				std::cerr << "ERROR: load: no expression pathways (topology) were found. The program may crash." << std::endl;
				//return false;
				break;
			default:
				std::cerr << "Warning: load: more than one topology sections were found. For safety reasons you should arrange your data to one unique block. Nevermind, trying to read them..." << std::endl;
			case 1:
				file.clear();
				for(auto pos = topologies.begin(); pos != topologies.end(); pos++) {
					file.seekg(*pos);
					loadRules(file);
				}
		}

		switch(settings.size()){
			case 0: 
				break;
			default:
				std::cerr << "Warning: load: more than one settings sections were found. For safety reasons you should arrange your settings to one unique block. Nevermind, trying to read them..." << std::endl;
			case 1:
				file.clear();
				for(auto pos = settings.begin(); pos != settings.end(); pos++) {
					file.seekg(*pos);
					loadActions(file);
				}
		}

		return true;
	}

	bool ExpressionModell::loadActions(std::ifstream &file){
		std::string line;

		//check for first line: is it a header?
		if(file.peek() == '%'){ //it is - than jump!
			getline(file, line);
		}

		//read in rest of file (containing: from to trigger matrix)
		while( file.peek() != '%' && getline(file, line)){
			if(line.size() == 0 || line[0] == '#') {
				continue; //ignore empty lines and lines starting with #
			}

			//initialise
			std::string action, temp;
			std::vector<std::string> argoments;
			std::istringstream linestream(line);

			//read in data
			linestream >> action;
			while(linestream >> temp) argoments.push_back(temp);

			//add data
			/* actions: 
			 * 	- changeM
			 * 	- switch
			 * types of trigger:
			 * 	- by state
			 * 	- by state, delayed
			 * 	- by state, delayed, treshold
			 * 	- at time 
			 * */

			// make sense of data
			Action *todo = nullptr;
			auto timer_from = argoments.end();

			for(auto a = argoments.begin(); a != argoments.end(); a++) if(*a == "@") {timer_from = ++a; break;}
			int no_args = timer_from==argoments.end()?argoments.size():(std::distance( argoments.begin(), timer_from)-1);
			int no_timers = std::distance(timer_from, argoments.end());

			switch(action2enum(action)){
				case enone:
					std::cerr << "WARNING: loadActions: could not recognise action " << action << ". Ignoring this line." << std::endl;
					break;
				case eActionChangeM:
					if( no_args ){
						todo = (Action*) new ActionChangeM(std::stoi(argoments[0]));
					} else { //it is just a reset to 0 immediately order
						todo = (Action*) new ActionChangeM();
					}
					break;
				case eActionSwitch:
					switch( no_args ){
						case 0:
							todo = (Action*) new ActionSwitch();
							break;
						case 1:
							todo = (Action*) new ActionSwitch( argoments[0]);
							break;
						case 2:
							todo = (Action*) new ActionSwitch( argoments[0], std::stof(argoments[1]) );
							break;
						default: 
							std::cerr << "WARNING: loadActions: incorrect number of argoments to " << action << ". Allowed numbers: 0, 1 or 2. Ignoring line." << std::endl;
					}
					break;
				case eActionSet:
					switch( no_args ){
						case 0:
							todo = (Action*) new ActionSet();
							break;
						case 1:
							todo = (Action*) new ActionSet( argoments[0]);
							break;
						case 2:
							todo = (Action*) new ActionSet( argoments[0], std::stof(argoments[1]) );
							break;
						default: 
							std::cerr << "WARNING: loadActions: incorrect number of argoments to " << action << ". Allowed numbers: 0, 1 or 2. Ignoring line." << std::endl;
					}
					break;
				case eActionSetDecay:
					{
						std::vector<double> doublevec;
						auto aend = timer_from;
						if(no_timers) aend--;
						for(auto a = argoments.begin(); a != aend; a++) doublevec.push_back(std::stof(*a));
						todo = (Action*) new ActionSetDecay( doublevec );
					}
					break;
				case eActionRun:
					if(no_args) todo = (Action*) new ActionRun( std::stof(argoments[0]) );
					else todo = (Action*) new ActionRun();
					break;
			}

			if(todo != nullptr){ //if there is a valid action
				if(no_timers){
					if(no_timers == 1 && dvtools::isPosInt(*timer_from)){ //at give time
						std::cout << "ERROR: loadActions: doing actions at given time is not yet implemented!" << std::endl;
					} else { //at state
						std::string watch = *timer_from++;
						double delay = no_timers>1?std::stod(*timer_from++):0, treshold = no_timers>2?std::stod(*timer_from):0.95;

						if( !addInnerTriggerAction(todo, watch, delay, treshold) ){
							std::cerr << "WARNING: loadActions: wrong timing syntax. Try the followings instead:" << std::endl << "\t@ state" << std::endl 
								<< "\t@ state delay" << std::endl 
								<< "\t@ state delay treshold" << std::endl;
						}
					}
				} else{
					//do it immediately
					todo->apply(this);
				}
			}
		} // while in lines

		return(true);

	}
	
	bool ExpressionModell::loadActions(const std::string &filename){
		return loadActions( (char *) filename.c_str());
	}
	
	bool ExpressionModell::loadActions(char *filename){
		std::ifstream file(filename);

		if(!file.is_open()) {
			std::cerr << "ERROR: ExpressionModell::loadActions: file (" << filename << ") not found!" << std::endl;
			return( false );
		}
		
		return loadActions(file);
	}

	bool ExpressionModell::addFactor(const int id){
		if(id < 0 || ((unsigned int) id) >= states.size() ){
			std::cerr << "WARNING! ExpressionModell::add_factor: invalid id (" << id << ")! Returning false" << std::endl;
			return(false);
		}

		//add empty factors
		for(auto st = states.begin(); st != states.end(); st++){
			st->push_back(false);
			*st <<= 1; //shift if forward
		}
		no_target_factors++;
		size++;

		//assign it
		states[id][0] = true;

		invalidate(); //as new values added to states dimensions of M are invalid
		update_masks();

		return true;
	}

	void ExpressionModell::addPatternName(int id, char* name){
		std::string name_string(name);
		ExpressionModell::addPatternName(id, name_string);
	}

	void ExpressionModell::addPatternName(int id, std::string &name){
		state_names[id] = name;
	}
	

	//basic topology
	void ExpressionModell::fork(int id_from, int id_to1, int id_to2, int trigger){
		if(_fork(id_from, id_to1, id_to2, trigger)){
			homoassoc(id_to1); //stabilize end states
			homoassoc(id_to2); 
		}
	}

	void ExpressionModell::_heteroassoc(const boost::dynamic_bitset<> *from, const boost::dynamic_bitset<> *to){
		if(from != NULL && to != NULL){
			Ms.build();
			addToM( *to, *from );
			//addToM( *to, *to);
		}
	}
	
	bool ExpressionModell::_fork(int id_from, int id_to1, int id_to2, int trigger){
		//check trigger (ids will be checked in hetero/homoasszoc)
		if(trigger < 0 || ((unsigned int) trigger) >= no_triggers){
			std::cerr << "WARNING! ExpressionModell::_fork: invalid trigger (" << trigger << ")! Operation aborted" << std::endl;
			return false;
		} else {
			//check M
			 Ms.build();
			
			//make trigger
			boost::dynamic_bitset<> tvec(size);
			tvec[no_target_factors + no_genes + trigger] = true;
			
			//add to M
			addToM(states[id_to1], states[id_from]); //transition to to_1 as default (with trigger 0)
			substractFromM(states[id_to1], tvec); //in case of trigger 0 it does neot go to to_1
			addToM(states[id_to2], tvec); //in case of trigger 1 it goes to to_2

			return true;
		}
	}

	int ExpressionModell::fork(int id_from, int id_to1, int id_to2){ //add new trigger automatically, return id of trigger
		//get new trigger 
		int trigger = addTrigger(); //addTrigger invalidetes M! Please use ther function!  

		//check M
		Ms.build();
		
		//make trigger
		boost::dynamic_bitset<> tvec(size);
		tvec[no_target_factors + no_genes + trigger] = true;
		
		//add to M
		homoassoc(id_to1); //stabilize end states
		homoassoc(id_to2); 
		//heteroassoc(id_from, id_to1); //transition to to_1 as default (with trigger 0)
		addToM(states[id_to1], states[id_from]); //transition to to_1 as default (with trigger 0)
		substractFromM(states[id_to1], tvec); //in case of trigger 0 it does neot go to to_1
		addToM(states[id_to2], tvec); //in case of trigger 1 it goes to to_2

		return(trigger);

	}

	bool ExpressionModell::_conditional(int id_from, int id_to, int trigger){	//conditional transition, given trigger
		//check trigger (ids will be checked in hetero/homoasszoc)
		if(trigger < 0 || ((unsigned int) trigger) >= no_triggers){
			std::cerr << "WARNING! ExpressionModell::transition: invalid trigger (" << trigger << ")! Operation aborted" << std::endl;
			return false;
		} else if(id_from < 0 || ((unsigned int) id_from) >= states.size()){ //check id_from (id_to will be invetigated in homoasszoc)
			std::cerr << "WARNING! ExpressionModell::transition: invalid id_from (" << id_from << ")! Operation aborted" << std::endl;
			return false;
		} else {
			//check M
			Ms.build();

			//make trigger
			boost::dynamic_bitset<> tvec(size);
			tvec[no_target_factors + no_genes + trigger] = true;

			//add to M
			substractFromM(states[id_from], tvec); //it does not goes backwards in case trigger 1
			addToM(states[id_to], tvec); //in case of trigger 1 it goes to to 

			return true;
		}
	}

	void ExpressionModell::transition(int id_from, int id_to, int trigger){	//conditional transition, given trigger
		if(_conditional(id_from, id_to, trigger)){
			//add to M
			homoassoc(id_to); //stabilises in to
		}
	}

	void ExpressionModell::condTransition(int id_from, int id_to, int trigger){ //same as above
		transition(id_from, id_to, trigger);
	}

	//int ExpressionModell::condTransition(int id){return -1;}	//conditional transition, get automatically trigger, return its id

	bool ExpressionModell::hasFactor(const int id = -1) const {
		if(id < 0 || ((unsigned int) id) >= states.size() ){
			std::cerr << "WARNING! ExpressionModell::has_factor: invalid id (" << id << ")! Returning false" << std::endl;
			return(false);
		}

		//checking existing
		for(unsigned int no_tf = 0; no_tf < no_target_factors; no_tf++){
			 if( states[id][no_tf] ){
				 return( true ); //if other methods work properly, it has to be unique factor
			 }
		}

		return false;
	}


	void ExpressionModell::invalidate(){
		if(p != NULL) {
			delete [] p;
			p = NULL;
		}
	}
	
	///simulation
	void ExpressionModell::run(double number_of_timesteps){
		//check
		Ms.build();
		
		if(p == NULL){
			setState(0.0);
		}

		//check deathvec
		if(use_death_vector){
			if(deathvec.size() != size){
				std::cerr << "ERROR: ExpressionModell::run: invalid deathvector. Set to use individual deathvector instead" << std::endl;
				use_death_vector = false;
			}
		}


		//updating
		if(t_output != 0) (*this.*inic_out)();
		//below is neccessary due to double inprecision
		for(int until = number_of_timesteps/delta_t, timecount = 0, print_time = t_output/delta_t; timecount < until; timecount++, time += delta_t ){
			//Ms.set();
			if( timecount % print_time == 0) {
				(*this.*do_out)();
			}
			update();
			checkInnerTriggers();
		}
		if(t_output != 0) (*this.*do_out)();
	}

	void ExpressionModell::reset(){
		time=0;
		ranks.clear();
		ranks.reserve( seq_times.size() );
		for(int i = seq_times.size(); i--; ) ranks.push_back(0.0);
		rank = 0;
		seq_time_elapsed = 0;
		last_seq = -1;
		resetInnerTriggers();
	}

	void ExpressionModell::update(){
		/* $$\frac{\Delta \underline p}{\Delta t} = \frac{\tau}{2} \left \{ tanh \left [ slope \left ( \underline {\underline M} \cdot \underline p + shift \right ) \right ] +1 \right \} - \delta \underline p$$
		 *
		 * 1. take M*p
		 * 2. take its tangens hyperbolicus (to transform it to 0-1)
		 * 3. apply death part
		 *
		 * */

		//updateMactual();
		double *M_ptr = M;
		for(unsigned int row = 0; row < size; row++){
			double derivate = 0.0;

			//sum M[row, ] * p[] 
			for(unsigned int col = 0; col < size; col++)	{
				derivate += (*M_ptr * p[col]);
//				std::cout << "M " << *M_ptr << " p " << p[col] << std::endl;;
				M_ptr++;
			}

			//apply tanh and substract death
//			std::cout << "reriv_bef " << derivate;
			if(use_death_vector) derivate = tau_per2 * (1.0 + std::tanh(slope_tanh * ( derivate + shift))) - deathvec[row] * p[row];
			else derivate = tau_per2 * (1.0 + std::tanh(slope_tanh * ( derivate + shift))) - death * p[row];
//			std::cout << " after " << derivate << std::endl;
			//apply numeric derivate
			p[row] += (derivate * delta_t);
//			if(p[row] < 0.0) p[row] = 0.0;
		}
	}

	void ExpressionModell::setMatrixValue(const int row, const int col, const double value){
		M[row * size + col] = value;
	}
	
	void ExpressionModell::setState(const int id){ //inic current state to given pattern
		const boost::dynamic_bitset<> *pattern = asState(id);
		if(pattern != NULL){
			if(p != NULL) { // p is set, free it!
				delete [] (p);
			}

			//allocate p
			p = new double [size];

			//set new values
			for(unsigned int val = 0; val < size; val++){ 
				if( (*pattern)[val] )
					p[val] = E;
				else 
					p[val] = 0;
			}
		}
	}
	bool ExpressionModell::setState(const std::string &id, double value){ //inic current state to given pattern
		const boost::dynamic_bitset<> *pattern = asState(id);
		if(value < 0) value = E;

		if(pattern != NULL){
			if(p != NULL) { // p is set, free it!
				delete [] (p);
			}

			//allocate p
			p = new double [size];

			//set new values
			for(unsigned int val = 0; val < size; val++){ 
				if( (*pattern)[val] )
					p[val] = value;
				else 
					p[val] = 0;
			}
		} else return false;
		return true;
	}

	bool ExpressionModell::setState(const std::vector<double> &state){
		//check input data
		if(state.size() != size){
			std::cerr << "ERROR! ExpressionModell::set_state: new state has incorrect size: " << state.size() << " instead of " << size << ". Operation aborted." << std::endl;
			return false;
		}

		if(p == NULL) { // p is not set, allocate
			p = new double [size];
		}

		//write data
		for(unsigned int expr = 0; expr < size; expr++){
			p[expr] = state[expr];
		}

		return true;
			
	}

	void ExpressionModell::setState(const double value){ //inic state to the same value overall
		if(p == NULL) { // p is not set, allocate
			p = new double [size];
		}

		for(unsigned int expr = 0; expr < size; expr++){
			p[expr] = value;
		}
	}

	bool ExpressionModell::setState(double mu, double sd){return false;} //inic state to random values from a normal distribution
	
	///outputting
	bool ExpressionModell::saveRules(char *filename){return false;}

	// Ranking
	double ExpressionModell::getRating() const{
		return(rank);
	}

	void ExpressionModell::setSequenceTimes(const std::vector<double> &_times){
		if(_times.size() != seq.size()){
			std::cerr << "ExpressionModell::setSequenceTimes: set appropriate seq before seq_times!" << std::endl;
		} else {
			seq_times.clear();
			seq_times = _times;

			ranks.clear();
			ranks.reserve( seq_times.size() );
			for(int i = seq_times.size(); i--; ) ranks.push_back(0.0);

			seq_time_elapsed = 0.0;
		}
	}

	void ExpressionModell::setSequence(const std::vector<int> &_seq){
		seq.clear();
		seq = _seq;
	}

	void ExpressionModell::inicRanking() {
		last_seq = -1;
		rank = 0;
	}

	void ExpressionModell::outRanking(){
		//check if we are at the end
		if(rank != 1){

			// look for maximum pearson (if it is above treshold)
			double max = -2.0;
			int max_it= -1;
			for(unsigned int st = 0; st < states.size(); st++){
				double pr = getPearson( states[st] );
				if( pr > max){
					max = pr;
					if(pr >= 0.95){
					//if(pr >= 0.8){
						max_it = st;
					}
				}
			}

			//check stuff
			if(max_it >= 0){ //only check if we a maximum
				if(last_seq >= 0)  { //if we have already started a sequence 
					if( max_it == seq[last_seq] ) return; //if the current is the dominant

					// if we are at the next step: step forward
					if(max_it == seq[last_seq + 1]){
						if(++last_seq == (int)(seq.size()-1) ){
							rank = 1;
						}
						return;
					} else { // it is something else -> delete progress
						last_seq = -1;
					}
				}
				//if we are here last_seq is -1!
				if(max_it == seq[0]) last_seq = 0; //check if it is first state
			}
		}

	}

	void ExpressionModell::outRanking3(){
		// look for maximum pearson (if it is above treshold)
		double max = -2.0;
		int max_it= -1;
		for(unsigned int st = 0; st < states.size(); st++){
			double pr = getPearson( states[st] );
			if( pr > max){
				max = pr;
				if(pr >= 0.95){
					max_it = st;
				}
			}
		}

		if(max_it >= 0){ //only check if we have a maximum
			//check stuff
			if(last_seq >= 0)  { //if we have already started a sequence 
				//if the current is the dominant
				if( max_it == seq[last_seq] ) {
					ranks[last_seq] += delta_t;					
					return; 
				} else {  
					//if we were at the final lap
					if( last_seq == (int)(seq.size()-1) ){
						rank = 0.0;
						for(auto rank_it = ranks.begin(); rank_it != ranks.end(); rank_it++) rank += *rank_it;
						rank /= ranks.size();
						last_seq = -1;
					} else { //we were NOT in the final lap
						// if we are at the next step: step forward
						if(max_it == seq[last_seq + 1]){
							ranks[++last_seq] = delta_t;
							return;
						} else { // it is something else -> delete progress
							last_seq = -1;
						}
					}
				}
			} 

			//if we are here last_seq is -1!
			if(max_it == seq[0]) {
				last_seq = 0; //check if it is first state
				ranks[0] = delta_t;
			}
			
		}

	}

	void ExpressionModell::outRanking2(){
		if(last_seq >= 0)  { //if we have already started a sequence 
			if( seq_time_elapsed >= seq_times[last_seq] ){ // if we have ended our observance window
				last_seq++;
				seq_time_elapsed = 0.0;

				if((unsigned int)last_seq == seq.size()) {//if we have ended our observance window
					rank=0.0;

					for(unsigned int s = 0; s < seq.size(); s++){
						rank += ranks[s] / seq_times[s];
						ranks[s] = 0.0;
					}
					rank /= (double) seq.size() / delta_t;

					last_seq = -1;
					return;
				}
			}

			ranks[last_seq] += getPearson( states[seq[last_seq]] );
			seq_time_elapsed += delta_t;
		} else { // should i start a ranking sequence?
			if( getPearson( states[seq[0]] ) > 0.95) {
				last_seq = 0;
				ranks[0] = getPearson( states[seq[0]] );
			}
		}

	}
	
#ifdef _EXPR_MOD_RCPP_

	void ExpressionModell::inicRComposite() {
		inicRPearsonOutput();
		inicRanking();
	}

	void ExpressionModell::outRComposite2(){
		outRPearson();
		outRanking2();
	}

	void ExpressionModell::outRComposite(){
		outRPearson();
		outRanking();
	}

	void ExpressionModell::inicRPearsonOutput() {
		//if its not empty, clear it
		if(df_out.size() != (states.size()+1) ) {
			df_out.clear();
//			std::cout << "Cleared" << std::endl;
		}
		//pushing vectors for columns containing time and pearson r values for (p, states)
		if(df_out.empty()) for( int n = states.size() + 1; n--;){
			df_out.emplace_back();
			df_out.back().reserve(100);
		}
	}

	void ExpressionModell::outRPearson(){
		auto df_col = df_out.begin();
		//assign time
		df_col->push_back(time);
		df_col++;

		//adding stat Pearson vals
		for(auto st = states.begin(); st != states.end(); st++, df_col++){
			df_col->push_back( getPearson(*st) );
		}
	}

	Rcpp::DataFrame ExpressionModell::getRdf(){
		if(time == 0){
			df = Rcpp::DataFrame::create();
			return df;
		}
		if(df.size() > 0) { //df has content
			Rcpp::NumericVector firstcol = df[0];
			if(df_out[0].back() == (double) firstcol[0] ) return df; //if it hasnot changed, return it
			df.erase(df.begin(), df.end()); //clear it
		} 
		auto df_col = df_out.begin();
		auto st_name = state_names.begin();

		//add first column
		df = Rcpp::DataFrame::create(Rcpp::Named("time") = *df_col);
		
		//adding other rows
		for(df_col++ ;df_col != df_out.end(); df_col++, st_name++){
			df.push_back( *df_col, *st_name);
		}

		//Rcpp::Rcout << df.length() << std::endl;

		return(df);
	}

	Rcpp::DataFrame ExpressionModell::getPattern_df(){
		Rcpp::DataFrame df_patterns;
		auto st_name = state_names.begin();

		for(auto st = states.begin(); st != states.end(); st++, st_name++){
			Rcpp::LogicalVector logvec;
			int expr = 0;
			for(auto name = factor_names.begin(); name != factor_names.end(); name++, expr++){
				logvec.push_back( (*st)[expr], *name );
			}
			for(auto name = gene_names.begin(); name != gene_names.end(); name++, expr++){
				logvec.push_back( (*st)[expr], *name );
			}
			for(auto name = trigger_names.begin(); name != trigger_names.end(); name++, expr++){
				logvec.push_back( (*st)[expr], *name );
			}
			df_patterns.push_back( Rcpp::clone(logvec), *st_name);
		}

		return(df_patterns);

	}

	Rcpp::NumericVector ExpressionModell::getState_Rvec(){
		Rcpp::NumericVector vec;
		double *pval = p;
		for(auto name = factor_names.begin(); name != factor_names.end(); name++, pval++){
			vec.push_back( *pval, *name );
		}
		for(auto name = gene_names.begin(); name != gene_names.end(); name++, pval++){
			vec.push_back( *pval, *name );
		}
		for(auto name = trigger_names.begin(); name != trigger_names.end(); name++, pval++){
			vec.push_back( *pval, *name );
		}
		return(vec);
	}
#else
	void ExpressionModell::set_seed(unsigned long int seed){
		gsl_rng_set(r, seed);
	}
#endif

	void ExpressionModell::printStatePearsonHeader() {
		std::cout << "time";

		//print state names
		for( auto f = state_names.begin(); f != state_names.end(); f++ ){
			std::cout << "\t" << *f;
		}
		std::cout << std::endl;
	}

	void ExpressionModell::printStateHeader() {
		std::cout << "time";

		//print factor names
		for( auto f = factor_names.begin(); f != factor_names.end(); f++ ){
			std::cout << "\t" << *f;
		}
		std::cout << "\t|";

		//print gene names
		for( auto gene = gene_names.begin(); gene != gene_names.end(); gene++ ){
			std::cout << "\t" << *gene;
		}
		std::cout << "\t|";

		//print trigger names
		for( auto tr = trigger_names.begin(); tr != trigger_names.end(); tr++ ){
			std::cout << "\t" << *tr;
		}
		std::cout << std::endl;
	}

	void ExpressionModell::printState() {
		std::cout << time;
		for(unsigned int expr = 0; expr < size ; expr++){
			if(expr == no_target_factors || expr == (no_target_factors + no_genes) ) {
				std::cout << "\t|";
			}
			std::cout << "\t" << p[expr];
		}
		std::cout << std::endl;
	}

	std::vector<double> ExpressionModell::getStartVec(unsigned int no, double multiply){
		std::vector<double> vec;
		vec.reserve(size);

		if(no >= states.size() ) {
			std::cerr << "ERROR: getStartVec: unknown state requested: " << no << ". Number of patterns is " << states.size() << std::endl;
		} else {
			auto st = states[no];
			for(unsigned int val = 0; val < st.size(); val++){
				vec.push_back( ((double) st[val]) * multiply );
			}
		}
		return(vec);
	}

	double ExpressionModell::getPearson(int pattern_id) {
		if(pattern_id < 0 || ((unsigned int) pattern_id) >= states.size()){
			std::cerr << "ERROR: ExpressionModell::getPearson: not valid id! " << pattern_id << std::endl;
			return(100);
		} else {
			return getPearson(states[pattern_id]);
		}
	}

	double ExpressionModell::getPearson(const boost::dynamic_bitset<> &st) const{
		//ATTENTION!!:
		// it does not calculate paerson r value for the whole bitset, only for it first N character!!!
		int N = no_genes + no_target_factors;

		//Calculating sum for pearson correlation
		// sum_xy = sum(  p[i]*st[i]  )
		// sum_x  = sum(  p[i]   )
		// sum_y  = sum(  st[i]  )
		double sum_xy = 0.0, sum_x = 0.0, sum_y = 0.0, sum_x2 = 0.0, sum_y2 = 0.0;
		for(int i = 0; i < N; i++){
			if(st[i]) {
				sum_y2 += E2;
				sum_y += E;
				sum_xy += (E * p[i]);
			}
			sum_x += p[i];
			sum_x2 += (p[i] *p [i]);
		}

		//Calculate covariance ans SD -> r = Cov_xy / SD_x / SD_y
		double rvalue = dvtools::cov((double) N, sum_xy, sum_x, sum_y);  //temporalily it holds Cov_x,y
		double sdx = dvtools::sd(N, sum_x, sum_x2);
		//rvalue /= dvtools::sd(N, sum_x, sum_x2);
		double sdy = dvtools::sd(N, sum_y, sum_y2);
		//rvalue /= dvtools::sd(N, sum_y, sum_y2);
		
		if(sdx != 0) rvalue /= sdx;
		else rvalue = (rvalue>0)?1:-1;
			
		rvalue /=  sdy;

		return(rvalue);
	
	}

	void ExpressionModell::printStatePearson(){
		std::cout << time;

		//going tru states
		for(auto st = states.begin(); st != states.end(); st++){
			std::cout << '\t' << getPearson(*st);
		}
		std::cout << std::endl;

	}

	void ExpressionModell::printPatterns() const{
//		std::cout << gene_names.size() << std::endl;
	
		//print header
		std::cout << "state";

		//print factor names
		for( auto f = factor_names.begin(); f != factor_names.end(); f++ ){
			std::cout << "\t" << *f;
		}
		std::cout << "\t|";

		//print gene names
		for( auto gene = gene_names.begin(); gene != gene_names.end(); gene++ ){
			std::cout << "\t" << *gene;
		}
		std::cout << "\t|";

		//print trigger names
		for( auto tr = trigger_names.begin(); tr != trigger_names.end(); tr++ ){
			std::cout << "\t" << *tr;
		}
		std::cout << std::endl;

		//print patterns
		for(auto st = states.begin(); st != states.end(); st++){
			std::cout << state_names[st - states.begin()];
			for(unsigned int expr = 0; expr < size ; expr++){
				if(expr == no_target_factors || expr == (no_target_factors + no_genes) ) {
					std::cout << "\t|";
				}
				std::cout << "\t" << (*st)[expr];
			}
			std::cout << std::endl;
		}
	}

	void ExpressionModell::printDeath() const{
		std::cout << "name\tdeath" << std::endl;

		auto n = factor_names.begin();
		for(auto d = deathvec.begin(); d != deathvec.end(); d++, n++){
			if(n == factor_names.end()) n = gene_names.begin();
			if(n == gene_names.end()) n = trigger_names.begin();
			std::cout << *n << '\t' << *d << std::endl;
		}
	}

	void ExpressionModell::printM() const{
		if(M != NULL){
			double *M_ptr = M;
			for(int row = size; row--; ){
				for(int col = size; col--;) {
					std::cout << '\t' << (double) *(M_ptr);
					M_ptr++;
				}
				std::cout << std::endl;
			}
		} else {
			std::cerr << "ERROR: ExpressionModell::printM: There is no valid M matrix!" << std::endl;
		}
	}

	bool ExpressionModell::addToM(const boost::dynamic_bitset<> &to, boost::dynamic_bitset<> from){
		Ms.build();
		if(M != NULL){
			double *M_ptr = M;

			//nullify non-unique genes from from vector (as it is not consistent part of targeting, should it have no effect in targeting)
			if(no_target_factors) from &= ~mask_genes;
			/* from:	101010
			 * mask_genes:	001100
			 * ~mask_genes:	110011
			 * 
			 * 	 10 10 10
			 * 	&11 00 11
			 * 	---------
			 *	 10 00 10
			 * */

			for(unsigned int col = 0; col < size; col++){
				for(unsigned int row = 0; row < size; row++){
					if(from[row]) {
						*M_ptr += (to[col])?1:-1;
						M_ptr++;
					} else {
					M_ptr ++; //advance pointer size times
					}
				} //in rows
			} //in cols
		} else return false;

		return true;
	}

	bool ExpressionModell::substractFromM(const boost::dynamic_bitset<> &to, boost::dynamic_bitset<> from){
		Ms.build();
		if(M != NULL){
			double *M_ptr = M;

			//nullify non-unique genes from from vector (as it is not consistent part of targeting, should it have no effect in targeting)
			if(no_target_factors) from &= ~mask_genes;

			for(unsigned int col = 0; col < size; col++){
				for(unsigned int row = 0; row < size; row++){
					if(from[row]) {
						*M_ptr -= (to[col])?1:-1;
						M_ptr++;
					} else {
					M_ptr ++; //advance pointer size times
					}
				} //in rows
			} //in cols
		} else return false;

		return true;
	}

	void ExpressionModell::update_masks(){
		//realloc gene mask
		mask_genes.clear();
		mask_genes.resize(size, false);
		//realloc facotr mask
		mask_target_factors.clear();
		mask_target_factors.resize(size, false);
		//reallocate trigger mask
		mask_triggers.clear();
		mask_triggers.resize(size, false);

		//rewrite them
		boost::dynamic_bitset<> *bitset_to_write = &mask_target_factors;
		//std::cout << "mask_genes(" << no_genes << ") = " << mask_genes << ", mask_target_factors(" << no_target_factors << ") = " << mask_target_factors << ", mask_triggers(" << no_triggers << ") = " << mask_triggers << std::endl;		
		for(unsigned int b = 0; b < size; b++){
			if(b == no_target_factors) bitset_to_write = &mask_genes;
			if(b == (no_genes + no_target_factors) ) bitset_to_write = &mask_triggers;
			(*bitset_to_write)[b] = true; //there are more elegant solutions for this, but i am busy
		}
		//std::cout << "mask_genes = " << mask_genes << ", mask_target_factors = " << mask_target_factors << ", mask_triggers = " << mask_triggers << std::endl;		
	}


	void ExpressionModell::switchOff(const int trigger){
		if(trigger < 0 || ((unsigned int) trigger) >= no_triggers){
			std::cerr << "WARNING! ExpressionModell::switchOn: invalid trigger: " << trigger << std::endl;
		} else {
			p[no_target_factors + no_genes + trigger] = 0;
		}
	}

	void ExpressionModell::switchOn(const int trigger){
		if(trigger < 0 || ((unsigned int) trigger) >= no_triggers){
			std::cerr << "WARNING! ExpressionModell::switchOn: invalid trigger: " << trigger << std::endl;
		} else {
			p[no_target_factors + no_genes + trigger] = E;
		}
	}

	/*void ExpressionModell::switchOnAtState(const int trigger, const int state, const double above_r){
		switchValueAtState(trigger, state, E, above_r);
	}*/

	void ExpressionModell::switchValueAtState(const int trigger, const int state, double value, const double above_r){
		switchValueAtState(trigger, state, value, 0, above_r);
	}

	void ExpressionModell::switchValueAtState(const int trigger, const int state, double value, const double timelag, const double above_r){
		/* trigger: it is the trigger to be assigned value
		 * state: the innerTrigger. If it exceeds a given value, trigger will be switched on/off
		 * value: the value to be assigned to the trigger
		 * above_r: the treshold value for state
		 * */
//		std::cout << "switchValueAtState called\n";
		//check trigger
		if(trigger < 0 || ((unsigned int) trigger) >= no_triggers) {
			std::cerr << "ERROR! ExpressionModell::switchValueAtState: there is no such trigger as " << trigger << std::endl;
			return;
		}

		//check state
		/*if(state < 0 || ((unsigned int) state) >= states.size()){
			std::cerr << "ERROR! ExpressionModell::switchValueAtState: there is no such state as " << state << std::endl;
			return;
		}*/

		//check r value
		/*if(above_r * above_r > 1){
			std::cerr << "ERROR! ExpressionModell::switchValueAtState: Pearson correlation r treshold value must be between 0 and 1!" << std::endl;
			return;
		}*/

		//get a trigger
		auto tr = addInnerTrigger(state, timelag, above_r );

		//add action to trigger
//		std::cout << "switchValueAtState: got trigger\n";
		tr->addAction( (Action*) new ActionSwitch(trigger, value) );
//		std::cout << "switchValueAtState: Action added to trigger " << *(tr->inner_trigger) << std::endl;
	}

	void ExpressionModell::switchOffAtState(const int trigger, const int state, const double above_r){
		switchValueAtState(trigger, state, 0.0, above_r);
	}

       void ExpressionModell::resetInnerTriggers(){
               for(auto tr = inner_triggers.begin(); tr != inner_triggers.end(); tr++){
                                       tr->time_elapsed = 0.0;
               }
       }


	void ExpressionModell::checkInnerTriggers(){
//		std::cout << "Triggers at time " << time << ": ";
		for(auto tr = inner_triggers.begin(); tr != inner_triggers.end(); tr++){
			//if(getPearson( *(tr->inner_trigger) ) > tr->treshold) *(tr->target) = tr->value;
			if(getPearson( *(tr->inner_trigger) ) > tr->treshold) {
//				std::cout << "T ";
//				std::cout << "Trigger " << *(tr->inner_trigger) << "called at " << time << std::endl;
				if(tr->time_elapsed >= tr->timelag) {
					tr->apply(this);
					tr->time_elapsed = 0.0;
				}
				else {
					tr->time_elapsed += delta_t;
				}
			}
//			else std::cout << "N ";
		}
//		std::cout << std::endl;
	}
	bool ExpressionModell::changeOutput(const char type){
		switch(type){
			case 'S':
				inic_out = &ExpressionModell::printStateHeader;
				do_out = &ExpressionModell::printState;
				return true;
			case 'P':
				inic_out = &ExpressionModell::printStatePearsonHeader;
				do_out = &ExpressionModell::printStatePearson;
				return true;
			case 'm':
				inic_out = &ExpressionModell::inicRanking;
				do_out = &ExpressionModell::outRanking2;
				return true;
			case '3':
				inic_out = &ExpressionModell::inicRanking;
				do_out = &ExpressionModell::outRanking3;
				return true;
			case 'M':
				inic_out = &ExpressionModell::inicRanking;
				do_out = &ExpressionModell::outRanking;
				return true;
#ifdef _EXPR_MOD_RCPP_
			case 'R':
				inic_out = &ExpressionModell::inicRPearsonOutput;
				do_out = &ExpressionModell::outRPearson;
				return true;
			case 'c':
				inic_out = &ExpressionModell::inicRComposite;
				do_out = &ExpressionModell::outRComposite2;
				return true;
			case 'C':
				inic_out = &ExpressionModell::inicRComposite;
				do_out = &ExpressionModell::outRComposite;
				return true;
#endif
			default:
				std::cerr << "WARNING: ExpressionModell::changeOutput: type not correct: " << type << std::endl;
				return false;
		}
	}

	void dv_expr::ExpressionModell::buildM() {Ms.force_build();}

	void dv_expr::ExpressionModell::switchM(int id) {Ms.set(id);}

	int ExpressionModell::get_no_genes() const {return no_genes;}
	int ExpressionModell::get_no_target_factors() const {return no_target_factors;}
	int ExpressionModell::get_no_triggers() const {return no_triggers;}
	int ExpressionModell::get_size() const {return size;}
	double ExpressionModell::get_time() const {return time;}
	int ExpressionModell::get_no_states() const {return states.size();}
	double ExpressionModell::get_tau() const {return tau;}
	double ExpressionModell::get_death() const {return death;}
	std::vector<double> ExpressionModell::get_deathvec() const {
		std::vector<double> output;
		output.reserve(size);
		if(use_death_vector) output = deathvec;
		else while(output.size() < size){
			output.push_back(death); 
		}
		return output;
	}
	double ExpressionModell::get_death(const unsigned int no) const {
		if(no >= size) return(-1.0);
		else if(use_death_vector) return deathvec[no];
		else return(death);
	}
	double ExpressionModell::get_shift() const {return shift;}
	double ExpressionModell::get_E() const {return E;}
	double ExpressionModell::get_slope_tanh() const {return slope_tanh;}
	double ExpressionModell::get_delta_t() const {return delta_t;}

	// setters
	bool ExpressionModell::set_time(double _time){
		if(_time < 0) return false;
		time = _time;
		return true;
	}
}


#ifdef _EXPR_MOD_RCPP_
//Rcpp staff


bool check_doublevector_valid(SEXP* args, int nargs){
     if( nargs != 1 ) return false ;
     if( TYPEOF(args[0]) != REALSXP ) return false ;
     //return ( LENGTH(args[0]) == 1 ) ;
     return ( true ) ;
}

/*template<typename... Targs>
bool check_args(SEXP* args, int nargs){
	//types
	// REALSXP
	// STRSXP
	// INTSXP
	
	if ( nargs != sizeof...(Targs) ) return false ;
	for (int a = 0; a < nargs, a++) {
		( TYPEOF(args[a] != ) return false ;
	}
}*/

bool check_intvector_valid(SEXP* args, int nargs){
     if( nargs != 1 ) return false ;
     if( TYPEOF(args[0]) != INTSXP ) return false ;
     //return ( LENGTH(args[0]) == 1 ) ;
     return ( true ) ;
}

bool check_string_valid(SEXP* args, int nargs){
     //if( nargs != 1 ) return false ;
     if( TYPEOF(args[0]) != STRSXP ) return false ;
     //return ( LENGTH(args[0]) == 1 ) ;
     return ( true ) ;
}

bool check_1double_valid(SEXP* args, int nargs){
     if( nargs != 1 ) return false ;
     if( TYPEOF(args[0]) != REALSXP ) return false ;
     //return ( LENGTH(args[0]) == 1 ) ;
     return ( true ) ;
}

RCPP_MODULE(expr_mod) {
	using namespace Rcpp;
	class_<dv_expr::ExpressionModell>("ExpressionModell")

	.constructor()
	.property("time", &dv_expr::ExpressionModell::get_time, "Generation time of the modell. Initialised to 0.")
	.field("t_output", &dv_expr::ExpressionModell::t_output, "Frequency of outputting in run().")
	.property("size", &dv_expr::ExpressionModell::get_size, "Full length of state vector.")
	.property("no_genes", &dv_expr::ExpressionModell::get_no_genes, "Number of unique core genes.")
	.property("no_triggers", &dv_expr::ExpressionModell::get_no_triggers, "Number of triggers.")
	.property("no_target_factors", &dv_expr::ExpressionModell::get_no_target_factors, "Number of target factors.")
	.property("delta_t", &dv_expr::ExpressionModell::get_delta_t)
	.property("tau", &dv_expr::ExpressionModell::get_tau, "time scaling property")
	.property("death_rate", &dv_expr::ExpressionModell::get_death, "Degradation rate")
	.property("shift", &dv_expr::ExpressionModell::get_shift, "property to scale tanh transformation")
	.property("slope_tanh", &dv_expr::ExpressionModell::get_slope_tanh, "property to scale tanh transformation")
	.property("no_patterns", &dv_expr::ExpressionModell::get_no_states, "Number of patterns set")
	.property("E", &dv_expr::ExpressionModell::get_E, "default value of gene expression")
	//.property("p", )
	.method("run", &dv_expr::ExpressionModell::run)
	.method("loadRules",
		(bool (dv_expr::ExpressionModell::*) (const std::string &) )
		(&dv_expr::ExpressionModell::loadRules)
	)
	.method("addPatternFromFile",
		(int (dv_expr::ExpressionModell::*) (const std::string &) )
		(&dv_expr::ExpressionModell::addPatternFromFile)
	)
	.method("addPattern",								//name in R
		( int (dv_expr::ExpressionModell::*)(const std::vector<int> &) )	//c++ method signature
		( &dv_expr::ExpressionModell::addPattern) ,				//c++ method pointer
		"text",									//description in R show() - it has to be added in case of a validator
		&check_intvector_valid							//validator. if false, then goes to next method
	)
	.method("addPattern",
		( int (dv_expr::ExpressionModell::*)(const std::string &) )
		(&dv_expr::ExpressionModell::addPattern)
	)
	.method("heteroassoc", 
		( void (dv_expr::ExpressionModell::*)(int, int) )
		( &dv_expr::ExpressionModell::heteroassoc )
	)
	.method("homoasszoc",
		( void (dv_expr::ExpressionModell::*)(const std::string &) )
		( &dv_expr::ExpressionModell::homoassoc ),
		"homoassotiation: it makes a pattern to pull points to itself",
		&check_string_valid
	)
	.method("homoasszoc",
		( void (dv_expr::ExpressionModell::*)(int) )
		( &dv_expr::ExpressionModell::homoassoc)
	)
	.method("transition",
		( void (dv_expr::ExpressionModell::*)(const std::string &, const std::string &))
		(&dv_expr::ExpressionModell::transition),
		"transition from id_from to id_to",
		&check_string_valid 
	)
	.method("transition",
		( void (dv_expr::ExpressionModell::*)(int, int))
		(&dv_expr::ExpressionModell::transition),
		"transition from id_from to id_to"
	)
	.method("transition",
		( void (dv_expr::ExpressionModell::*)(int, int, int))
		(&dv_expr::ExpressionModell::transition),
		"transition from id_from to id_to"
	)
	.method("fork",
		( int (dv_expr::ExpressionModell::*)(int, int, int))
		( &dv_expr::ExpressionModell::fork)
	)
	.method("fork",
		( void (dv_expr::ExpressionModell::*)(int, int, int, int))
		( &dv_expr::ExpressionModell::fork)
	)
	.method("condTransition",
		(void (dv_expr::ExpressionModell::*)(int, int, int))
		(&dv_expr::ExpressionModell::condTransition)
	)
	.method("buildM", &dv_expr::ExpressionModell::buildM)
	/*.method("buildM",
		(void (dv_expr::ExpressionModell::*)())
		(&dv_expr::ExpressionModell::buildM)
	)
	.method("buildM",
		(void (dv_expr::ExpressionModell::*)(double))
		(&dv_expr::ExpressionModell::buildM)
	)*/
	.method("addFactors", &dv_expr::ExpressionModell::addFactors)
	.method("addFactor", &dv_expr::ExpressionModell::addFactor)
	.method("hasFactor", &dv_expr::ExpressionModell::hasFactor)
	.method("addTrigger", &dv_expr::ExpressionModell::addTrigger)
	//.property("p", &dv_expr::ExpressionModell::getState_Rvec, (bool (dv_expr:ExpressionModell::*)(const std::vector<double> &) )(&dv_expr::ExpressionModell::setState) )
	.property("p", &dv_expr::ExpressionModell::getState_Rvec )
	.method("setState",
		( bool (dv_expr::ExpressionModell::*)(const std::vector<double> &))
		(&dv_expr::ExpressionModell::setState),
		"set state of modell",
		&check_doublevector_valid 
	)
	.method("setState",
		( void (dv_expr::ExpressionModell::*)(const double))
		(&dv_expr::ExpressionModell::setState),
		"set state of modell",
		&check_1double_valid 
	)
	.method("setState",
		( bool (dv_expr::ExpressionModell::*)(double, double))
		(&dv_expr::ExpressionModell::setState)
	)
	//.method("switchValue", &dv_expr::ExpressionModell::switchValue)
	//.method("switchOn", &dv_expr::ExpressionModell::switchOn)
	.method("switchOff", &dv_expr::ExpressionModell::switchOff)
	.method("switchOnAtState",
			(void (dv_expr::ExpressionModell::*)(const int, const std::string &, const double, const double) )
			(&dv_expr::ExpressionModell::switchOnAtState)
	)
	.method("switchOffAtState", &dv_expr::ExpressionModell::switchOffAtState)
	.method("switchValueAtState", 
			(void (dv_expr::ExpressionModell::*)(const int, const int, double, const double) )
			(&dv_expr::ExpressionModell::switchValueAtState)
	)
	.method("switchValueAtState2", 
			(void (dv_expr::ExpressionModell::*)(const int, const int, double, const double, const double) )
			(&dv_expr::ExpressionModell::switchValueAtState)
	)
	.method("printPatterns", &dv_expr::ExpressionModell::printPatterns)
	.method("printState", &dv_expr::ExpressionModell::printState)
	.method("printStatePearson", &dv_expr::ExpressionModell::printStatePearson)
	.method("printM", &dv_expr::ExpressionModell::printM)
	.method("changeOutput", &dv_expr::ExpressionModell::changeOutput)
	.method("df", &dv_expr::ExpressionModell::getRdf)
	.method("patterns", &dv_expr::ExpressionModell::getPattern_df)
	.method("setDeathvec", 
		( void (dv_expr::ExpressionModell::*) (const std::vector<double> &))
		(&dv_expr::ExpressionModell::setDeath)
	)
	.method("setDeath", 
		( void (dv_expr::ExpressionModell::*)(const std::string &, double))
		(&dv_expr::ExpressionModell::setDeath),
		"set death to gene name"
	)
	.method("setDeathToVal", 
		( void (dv_expr::ExpressionModell::*)(const double))
		(&dv_expr::ExpressionModell::setDeath),
		"set death to on given val"
	)
	.method("addNewM", 
		(int (dv_expr::ExpressionModell::*)(int, const double, const double))
		(&dv_expr::ExpressionModell::addNewM)
	)
	.method("switchM", &dv_expr::ExpressionModell::switchM)
	.method("printDeath", &dv_expr::ExpressionModell::printDeath)
	.method("addMTrigger", 
		(bool (dv_expr::ExpressionModell::*)(const int, const int, const double, const double))
		( &dv_expr::ExpressionModell::addMTrigger )
	)
	.method("addToReg",
		(void (dv_expr::ExpressionModell::*)( int, int, std::string, int ))
		( &dv_expr::ExpressionModell::addToReg ),
		"argoments in order: from, to, trigger, id of M"
	)
	.method("printReg", &dv_expr::ExpressionModell::printReg)
	.method("buildFromReg", &dv_expr::ExpressionModell::buildFromReg)
	.method("perturbM", &dv_expr::ExpressionModell::perturbM)
	.method("boostPattern", &dv_expr::ExpressionModell::boostPattern)
	.method("printTriggers", &dv_expr::ExpressionModell::printTriggers)
	.method("getRating", &dv_expr::ExpressionModell::getRating)
	.method("setSequence", &dv_expr::ExpressionModell::setSequence)
	.method("setSequence2", 
			(void (dv_expr::ExpressionModell::*)(const std::vector<std::string> &))
			(&dv_expr::ExpressionModell::setSequence2)
	)
	.method("setSequenceTimes", &dv_expr::ExpressionModell::setSequenceTimes)
	.method("setMatrixValue", &dv_expr::ExpressionModell::setMatrixValue)
	.method("load", 
			(bool (dv_expr::ExpressionModell::*) (const std::string &))
			(&dv_expr::ExpressionModell::load)
	)
	;
}
#endif

