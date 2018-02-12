/*
 * ibex_CtcGauss.cpp
 *
 *  Created on: 06-02-2018
 *      Author: victor
 */

#include "ibex_CtcGauss.h"

using namespace std;

namespace ibex {

	GaussContractor::GaussContractor (const System& sys, IntervalVector & initial_box) : sys(sys), Ctc(sys.ctrs.size()), A(1,1), b(1){
		/*need try-catch?*/
		init_system(initial_box,sys);
	}

	void GaussContractor::contract(IntervalVector & box){
		vector <vector <pair <int,int> > > proj_vars;
		/*cleaning of permutation, PA,Pb lists*/
		perm_list.clear();
		An.clear();
		bn.clear();
		IntervalVector xn = box-box.mid();
		IntervalVector box_aux = xn;
		/*Perform gauss Jordan on the matrix A in order to create the permutation list*/
		best_gauss_jordan (A, box, perm_list, proj_vars,1e-8);
		while(box_aux != xn){
			for (int i = 0 ; i < proj_vars.size() ; i++){
				for (int j = 0 ; j < proj_vars[i].size() ; j++){
					int var = proj_vars[i][j].first;
					int eq = proj_vars[i][j].second;
					Interval value(b[eq]);
					for (int k = 0; k < A.nb_cols(); k++){
						if (k != var){
							value = value +(-A[eq][k]*xn[k]);
						}
					}
					value = value/A[eq][var];
					xn[var] = value&=xn[var];
					}
				}
			if (xn.is_empty()){
				box.set_empty(); /*???*/
				return;
			}
		}
		for (int i = 0 ; i < box.size() ; i++)
			box[i] = box[i].mid()+xn[i];
	}

	void GaussContractor::init_system(IntervalVector initial_box, const System& sys){
		A = sys.ctrs_jacobian(initial_box);
		b = -sys.ctrs_eval(initial_box.mid());
	}


} /* namespace ibex */

