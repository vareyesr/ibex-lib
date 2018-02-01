/*
 * ibex_Conditioners.cpp
 *
 *  Created on: 01-02-2018
 *      Author: victor
 */

#include "ibex_Conditioners.h"

using namespace std;

namespace ibex {

	bool compare(const std::pair<double, std::pair<int,int> >&i, const std::pair<double, std::pair<int,int> >&j){
		return i.first < j.first;
	}


	pair<int,int> find_next_pivot(IntervalMatrix A, IntervalVector x,set<int> ban_rows, set<int> ban_cols){

		vector<pair <double, pair<int,int> > > order_cols; /*impact value,variable,equation*/
		double impact_values[A.nb_rows()][A.nb_cols()];
		for (int j = 0 ; j < A.nb_cols() ; j++){
			if (ban_cols.count(j) != -1){
				order_cols.push_back(make_pair(-1,make_pair(j,-1)));
				for (int i = 0 ; i < A.nb_rows() ; i++){
					if (ban_rows.count(i) != -1)
						impact_values[i][j] = A[i][j].mag()*x[j].diam();
				}
			}
		}
		for (int j  = 0; j < order_cols.size() ; j++){
				for (int i = 0 ; i < A.nb_rows() ; i++){
					if (ban_rows.count(i) != -1){
						if (impact_values[i][j] > order_cols[j].first){
							order_cols[j].first = impact_values[i][j];
							order_cols[j].second.second = i;
						}
					}
				}
		}
		std::sort(order_cols.begin(), order_cols.end(), compare);
		if (std::abs(order_cols[0].first) < 1e-8) return make_pair(-1,-1);
		else return order_cols[0].second;
	}

	void best_gauss_jordan (IntervalMatrix& A, IntervalVector x, vector<Matrix> & perm_list,
				vector <vector <pair <int,int> > > proj_vars, double prec=1e-8){

		Matrix B(1,1);
		B.resize(A.nb_rows(),A.nb_cols());
		Matrix perm(1,1);
		set<int> ban_rows;
		set<int> ban_cols;
		perm.resize(B.nb_rows(),B.nb_rows());
		pair<int,int> max_values;

		while (ban_cols.size() != A.nb_cols()){
			ban_rows.clear();
			/*Initialize B*/
			B = A.mid();
			/*Initialize the permutation matrix*/
			for (int i = 0; i<A.nb_rows() ; i++)
				for (int j = 0; j<A.nb_rows() ; j++){
					if (i == j) perm[i][j] = 1;
					else perm[i][j] = 0;
				}

			while (ban_rows.size() != A.nb_rows()){
				pair<int,int> var_eq = find_next_pivot(A, x, ban_rows, ban_cols);
				if (var_eq.first !=-1){
					double coef = B[var_eq.second][var_eq.first];
					Matrix aux_perm(1,1);
					aux_perm.resize(A.nb_rows(),A.nb_rows());
					for (int k = 0; k<A.nb_rows() ; k++)
						for (int l = 0; l<A.nb_rows() ; l++){
							if (k == l) aux_perm[k][l] = 1;
							else aux_perm[k][l] = 0;
						}
					for (int k = 0 ; k < B.nb_rows() ; k++){
						if ((k != var_eq.second) &&( (B[k][var_eq.first] < -prec) || (B[k][var_eq.first] > prec))) {
							double factor = B[k][var_eq.first];
							aux_perm[k][var_eq.second] = -factor/coef;
							for (int l = 0 ; l < B.nb_cols() ; l++)
								B[k][l]	= B[k][l]-(B[var_eq.second][l]*factor/coef);
						}
					}
					/*make the pivot position 1*/
					for (int i = 0 ; i < B.nb_cols() ; i++)
						B[var_eq.second][i] = B[var_eq.second][i]/coef;
					aux_perm[var_eq.second][var_eq.second] = 1/coef;
					perm = aux_perm*perm;
				}
			}
			perm_list.push_back(perm);
		}
	}

	void all_gauss_jordan (IntervalMatrix& A, IntervalVector x, vector<Matrix> & perm_list,
					vector <vector <pair <int,int> > > proj_vars){



	}

} /* namespace ibex */
