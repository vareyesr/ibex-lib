/*
 * ibex_CtcGauss.cpp
 *
 *  Created on: 06-02-2018
 *      Author: victor
 */

#include "ibex_CtcGauss.h"

using namespace std;

namespace ibex {

void GaussContractor::write_ext_box(const IntervalVector& box, IntervalVector& ext_box) {
	int i2=0;
	for (int i=0; i<box.size(); i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		ext_box[i2]=box[i];
	}
}

void GaussContractor::read_ext_box(const IntervalVector& ext_box, IntervalVector& box) {
	int i2=0;
	for (int i=0; i<box.size(); i++,i2++) {
		if (i2==goal_var) i2++; // skip goal variable
		box[i]=ext_box[i2];
	}
}


	GaussContractor::GaussContractor (const System& sys, int goal_var) : sys(sys), Ctc(sys.ctrs.size()), A(1,1), b(1), goal_var(goal_var) {
		/*only for contrained problems*/
//		init_system(initial_box,sys);

	}

	void GaussContractor::contract(IntervalVector & ext_box){
		IntervalVector box(sys.nb_var);
		read_ext_box(ext_box,box);
		init_system(box,sys);
//		IntervalVector xn = box-box.mid();
//		IntervalVector last_box= xn;
//		bwd_mul(b,A,xn,0.01);
		vector <vector <pair <int,int> > > proj_vars;
		/*cleaning of permutation, PA,Pb lists*/
		perm_list.clear();
		bool box_size_change = false;
		IntervalVector xn = box-box.mid();
		IntervalVector box_aux = xn;
		/*Perform gauss Jordan on the matrix A in order to create the permutation list*/
		best_gauss_jordan (A, xn, perm_list, proj_vars,1e-8);
//		all_gauss_jordan (A, perm_list,proj_vars, 1e-8);
//			bool do_contraction = true;
//			IntervalVector last_box= xn;
//		while(do_contraction){
		IntervalVector last_box= xn;
		for (int i = 0 ; i < perm_list.size() ; i++){
			IntervalMatrix P = perm_list[i];
			IntervalMatrix An=P*A;
			IntervalVector bn= P*b;

				for (int j = 0 ; j < proj_vars[i].size() ; j++){
					int var = proj_vars[i][j].first;
					int eq = proj_vars[i][j].second;
					for (int k = 0 ; k < An.nb_rows(); k++){
						if (k != eq) An[k][var] = 0;
						else An[k][var] = 1;
					}
				}
				bwd_mul(bn,An,xn,0.01);
//				for (int k = 0 ; k  < proj_vars[i].size() ; k++){
//					int var = proj_vars[i][k].first;
//					int eq = proj_vars[i][k].second;
//					Interval value(bn[eq]);
//					for (int l = 0; l < An.nb_cols(); l++){
//						if (l != var){
//							value = value +(-An[eq][l]*xn[l]);
//						}
//					}
//					value = value/An[eq][var];
//					xn[var] = value&=xn[var];
//					if(xn.is_empty()){
//						box.set_empty();
//						return;
//					}
//				}
			}
//			if (xn == last_box) do_contraction = false;
//			else{
//
//				perm_list.clear(); proj_vars.clear();
//				best_gauss_jordan (A, xn, perm_list, proj_vars,1e-8);
//				box_size_change = true;
//			}
//		}
		if (xn.is_empty()){
			ext_box.set_empty();
			return;
		}
		if (xn != last_box){
			box = (box.mid()+xn);
			write_ext_box(box,ext_box);
		}
	}

	void GaussContractor::init_system(IntervalVector initial_box, const System& sys){
//		A.resize(sys.ctrs.size(),sys.box.size());
//		sys.f_ctrs.hansen_matrix(initial_box,A);
		b = -sys.ctrs_eval(initial_box.mid());
		A = sys.ctrs_jacobian(initial_box);
	}


} /* namespace ibex */

