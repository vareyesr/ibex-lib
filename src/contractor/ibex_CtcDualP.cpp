/*
 * ibex_CtcDualP.cpp
 *
 *  Created on: Apr 13, 2021
 *      Author: victor
 */

#include "ibex_CtcDualP.h"
#include "ibex_LinearizerFixed.h"
#include "ibex_Interval.h"

using namespace std;

namespace ibex {


CtcDualP::CtcDualP(LinearizerXTaylor& lr, IntervalVector initial_box) :
		Ctc(lr.nb_var()), lr(lr),
		mylinearsolver(nb_var),
		last_box(initial_box),dual_sols(0,0),b_input(0),A_input(0,0) {

	update_dual_sols(initial_box); //sets up for the initial box
}

CtcDualP::~CtcDualP() {
}

void CtcDualP::update_dual_sols(IntervalVector box){
	last_box = box;
	int cont = lr.linearize(box, mylinearsolver);
	A_input = lr.A_input;
	b_input = lr.b_input;
	dual_sols.resize(2*A_input.nb_cols(),A_input.nb_rows());
	compute_dual(A_input,b_input,last_box);
}

/**
 * A_dual= input
 * b_dual= coefs
 * c_dual= obj func
 */
void CtcDualP::compute_dual(Matrix A, Vector b, IntervalVector box){

	LPSolver dual_solver(A.nb_cols()*2+A.nb_rows());
	Vector lambdas(b.size());
	IntervalVector dual_bounds(A.nb_cols()*2+A.nb_rows());
	//for min-max, i=0 min ; i=1 max
	for (int i = 0 ; i < 2 ; i++){
		for (int j = 0 ; j < A.nb_cols() ; j++){
			dual_solver.clear_bounds();
			dual_solver.clear_constraints();
			//setting bounds
			for(int k = 0 ; k < 2*A.nb_cols() ; k++){
				if (i == 0) dual_bounds[k] = Interval(NEG_INFINITY,0);
				else dual_bounds[k] = Interval(0,POS_INFINITY);
			}
			for (int k = 2*A.nb_cols() ; k <  2*A.nb_cols()+A.nb_cols() ; k++)
				dual_bounds[k] = Interval(NEG_INFINITY,POS_INFINITY);
			Vector c_dual(A.nb_cols()*2+A.nb_rows());
			Vector b_dual(A.nb_cols());
			Matrix A_dual(A.nb_cols(),2*A.nb_cols()+A.nb_rows());
			//c_dual
			for (int k = 0 ; k < A.nb_cols() ; k++){
				if (k != j){
					if (i == 0){//min
						c_dual[k] = box[k].lb();
					    c_dual[k+A.nb_cols()] = -box[k].ub();
					}
					else if (i == 1){//max
			            c_dual[k] = -box[k].lb();
			            c_dual[k+A.nb_cols()] = box[k].ub();
					}
				}
			}
			//b_dual
			b_dual[A_dual.nb_rows()-1] = 1;
			//A_dual
			int l = 0;
			for (int k = 0 ; k < A_dual.nb_rows()-1 ; k++){
				if (l == j) l++;
				if (k < A.nb_cols()){
					A_dual[k][l] = 1;
					A_dual[k][l+A.nb_cols()] = -1;
				}
				for (int m = 0 ; i < A.nb_rows() ;m++)
					A_dual[k][2*A.nb_cols()+m] = A[m][l];
				l++;
			}
			for (int k = 0 ; k < A.nb_rows() ; k++)
				A_dual[A_dual.nb_rows()-1][k+2*A.nb_cols()] = A[k][j];

			//set the dual solver
			dual_solver.set_bounds(dual_bounds);
			dual_solver.set_cost(c_dual);
			dual_solver.add_constraints(A_dual,EQ,b_dual);


			//solve
			LPSolver::Status stat = dual_solver.minimize();
			Vector v(A.nb_cols()*2+A.nb_rows());
			v = dual_solver.not_proved_primal_sol();
			for (int k = 0 ; k < dual_sols.nb_cols() ; k++){
				dual_sols[(i+1)*j][k] = v[2*A.nb_cols()+k];
			}
		}

	}
}

void CtcDualP::contract(IntervalVector& box){
	//call optimized fwd_bwd
}

}

