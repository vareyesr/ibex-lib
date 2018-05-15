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
		counter = 0;

	}

	void GaussContractor::contract(IntervalVector & ext_box){
		if (counter > ext_box.size()) return;
		IntervalVector box(sys.nb_var);
		IntervalVector xn(sys.nb_var);
		read_ext_box(ext_box,box);
		init_system(box,sys);

		/*just the linearization*/
//		xn = box-box.mid();
//		IntervalVector last_box= xn;
//		bwd_mul(b,A,xn,0.01);
//		if (xn.is_empty()){
//			ext_box.set_empty();
//			return;
//		}
//		if (xn != last_box){
//			box = (box.mid()+xn);
//			write_ext_box(box,ext_box);
//		}

		/*Gauss+linearization*/
		vector <vector <pair <int,int> > > proj_vars;
//		/*cleaning of permutation, PA,Pb lists*/
		perm_list.clear();
		bool box_size_change = false;
		xn = box-box.mid();
		IntervalVector box_aux = xn;
		/*Perform gauss Jordan on the matrix A in order to create the permutation list*/
		best_gauss_jordan (A, xn, perm_list, proj_vars,1e-8);
		IntervalVector last_box= xn;
		if (perm_list.size() > 0){
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
				if (xn.is_empty()){
					ext_box.set_empty();
					return;
				}
			}
			if (xn != last_box){
				box = (box.mid()+xn);
				if (box.is_subset(ext_box)){
					write_ext_box(box,ext_box);
					counter = 0;
				}
			}
			else counter++;
		}
		else {
			bwd_mul(b,A,xn,0.01);
			if (xn.is_empty()){
				ext_box.set_empty();
				return;
			}
			if (xn != last_box){
				box = (box.mid()+xn);
				if (box.is_subset(ext_box)){
					write_ext_box(box,ext_box);
					counter = 0;
				}
			}
			else counter++;
		}
	}

	void GaussContractor::init_system(IntervalVector initial_box, const System& sys){
		b.resize(sys.ctrs.size());
		A.resize(sys.ctrs.size(),sys.box.size());
		sys.f_ctrs.hansen_matrix(initial_box,A);
		for (int i = 0; i < b.size() ; i++){
			if (sys.ops[i] == EQ) b[i] = -sys.ctrs_eval(initial_box.mid())[i];
			else if (sys.ops[i] == LEQ || sys.ops[i] == LT) b[i] = Interval(NEG_INFINITY,-sys.ctrs_eval(initial_box.mid()).ub()[i]);
			else if (sys.ops[i] == GEQ || sys.ops[i] == GT) b[i] = Interval(-sys.ctrs_eval(initial_box.mid()).lb()[i],POS_INFINITY);
		}
//		A = sys.ctrs_jacobian(initial_box);
	}


} /* namespace ibex */

