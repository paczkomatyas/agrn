// [[Rcpp::depends(BH)]]

#include<Rcpp.h>
using namespace Rcpp;

RCPP_MODULE(unif_module) {
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
  //.method("addPatternFromFile", &dv_expr::ExpressionModell::addPatternFromFile)
	//.method("addPattern", &dv_expr::ExpressionModell::addPattern)
	;
}
