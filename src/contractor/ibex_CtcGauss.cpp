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
		init_system(box,sys);
		vector <vector <pair <int,int> > > proj_vars;

		/*cleaning of permutation, PA,Pb lists*/
		perm_list.clear();
		bool box_size_change = false;
		IntervalVector xn = box-box.mid();
		IntervalVector box_aux = xn;
		/*Perform gauss Jordan on the matrix A in order to create the permutation list*/
		best_gauss_jordan (A, xn, perm_list, proj_vars,1e-8);
		bool do_contraction = true;
		IntervalVector box_aux_aux =box;
		while(do_contraction){
			box_aux = xn;
//			cout << perm_list.size() << endl;
			for (int i = 0 ; i < perm_list.size() ; i++){
				IntervalMatrix An = perm_list[i]*A;
				IntervalVector bn = perm_list[i]*b;
				for (int j = 0 ; j < proj_vars[i].size() ; j++){
					int var = proj_vars[i][j].first;
					int eq = proj_vars[i][j].second;
					for (int k = 0 ; k < An.nb_rows(); k++){
						if (k != eq) An[k][var] = 0;
						else An[k][var] = 1;
					}
				}
				for (int j = 0 ; j < proj_vars[i].size() ; j++){
					int var = proj_vars[i][j].first;
					int eq = proj_vars[i][j].second;
					Interval value(bn[eq]);
					for (int k = 0; k < An.nb_cols(); k++){
						if (k != var){
							value = value +(-An[eq][k]*xn[k]);
						}
					}
					xn[var] = xn[var]&=value;
				}
			}
			if (xn.is_empty()){
				box.set_empty();
				return;
			}
			if (box_aux == xn){
				do_contraction = false;
			}
			else {
				perm_list.clear(); proj_vars.clear();
				best_gauss_jordan (A, xn, perm_list, proj_vars,1e-8);
				box_size_change = true;
			}
		}

		if (box_size_change){
			box = box.mid()+xn;
			return;
		}
	}

	void GaussContractor::init_system(IntervalVector initial_box, const System& sys){
		A = sys.ctrs_jacobian(initial_box);
		b = -sys.ctrs_eval(initial_box.mid());
	}


} /* namespace ibex */

