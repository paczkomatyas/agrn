#include "expr_mod.h"

using namespace std;

int main(int argc, char *argv[]){

	dv_expr::ExpressionModell hema;
	hema.load( std::string(argv[1]) );

	return 0;
}
