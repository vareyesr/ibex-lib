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
		return i.first > j.first;
	}

	void combinatorial(IntervalMatrix A, int cols,int rows,std::vector< std::vector <int> > & comb_piv){
		vector<int> pivots;
		/*Initialize the possible pivot combination*/
		for (int i = 0; i < rows ; i++)
			pivots.push_back(i);
		int last = pivots.size()-1;
		bool end = true;
		while (end){
			comb_piv.push_back(pivots);
			if (pivots[last] < cols-1)
				pivots[last]++;
			else{
				int k=0;
				for (int i = last-1 ; i>=0 ; i--){
					if (pivots[i] != cols-2-k){
						end = true;
						pivots[i]++;
						int aux = pivots[i];
						k=1;
						for (int j = i+1 ; j < last+1 ; j++){
							pivots[j] = aux + k;
							k++;
						}
						break;
					}
					else
						end =false;
					k++;
				}
			}
		}
	}

	pair<int,int> find_next_pivot(IntervalMatrix & PA, IntervalVector x,set<int> & ban_rows, set<int> & ban_cols){

		vector<pair <double, pair<int,int> > > order_cols; /*impact value,variable,equation*/
		Matrix impact_values(PA.nb_rows(),PA.nb_cols());
		for (int j = 0 ; j < PA.nb_cols() ; j++){
			if (ban_cols.count(j) != 1){
				order_cols.push_back(make_pair(-1,make_pair(j,-1)));
				for (int i = 0 ; i < PA.nb_rows() ; i++){
					if (ban_rows.count(i) != 1)
						impact_values[i][j] = PA[i][j].mag()*x[j].diam();
				}
			}
		}

		for (int j  = 0; j < order_cols.size() ; j++){
				for (int i = 0 ; i < PA.nb_rows() ; i++){
					if (ban_rows.count(i) != 1){
						int col = order_cols[j].second.first;
						if (impact_values[i][col] > order_cols[j].first){
							order_cols[j].first = impact_values[i][col];
							order_cols[j].second.second = i;
						}
					}
				}
		}
		std::sort(order_cols.begin(), order_cols.end(), compare);
		if (std::abs(order_cols[0].first) < 1e-8) return make_pair(-1,-1);
		else{
			ban_rows.insert(order_cols[0].second.second);
			ban_cols.insert(order_cols[0].second.first);
			return order_cols[0].second;
		}
	}

	void best_gauss_jordan (IntervalMatrix A, IntervalVector x, vector<Matrix> & perm_list,
				vector <vector <pair <int,int> > > & proj_vars, double prec=1e-8){

		vector <pair<int,int> > aux_list;
		Matrix B(1,1);
		B.resize(A.nb_rows(),A.nb_cols());
		Matrix perm(1,1);
		set<int> ban_rows;
		set<int> ban_cols;
		perm.resize(B.nb_rows(),B.nb_rows());
		pair<int,int> max_values;
		bool available_cols = true;
		while (available_cols){
			aux_list.clear();
			ban_rows.clear();
			/*Initialize B*/
			B = A.mid();
			IntervalMatrix PA = A;
			/*Initialize the permutation matrix*/
			for (int i = 0; i<A.nb_rows() ; i++)
				for (int j = 0; j<A.nb_rows() ; j++){
					if (i == j) perm[i][j] = 1;
					else perm[i][j] = 0;
				}
			while (ban_rows.size() != A.nb_rows()){
				if (ban_cols.size() == A.nb_cols()){
					ban_cols.clear();
					available_cols = false;
				}
				pair<int,int> var_eq = find_next_pivot(PA, x, ban_rows, ban_cols);
				if (var_eq.first !=-1){
					aux_list.push_back(make_pair(var_eq.first,var_eq.second));
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
					PA = aux_perm*PA;
					perm = aux_perm*perm;
					proj_vars.push_back(aux_list);
				}
			}
			perm_list.push_back(perm);
			if(A.nb_cols()==ban_cols.size()) available_cols = false;
		}
	}

	void gauss_jordan_all (IntervalMatrix& A, vector<Matrix>& permutations,vector < vector < pair <int, int> > > & pair_contr_all , double prec){
			int temp_piv;
			set <int> rows_checked;
			std::vector< std::vector <int> > comb_piv;
			/*get all possible combinations of pivots*/
			combinatorial(A,A.nb_cols(),A.nb_rows(),comb_piv);
			vector <pair<int,int> > aux_list;
			Matrix B(1,1);
			B.resize(A.nb_rows(),A.nb_cols());
			Matrix perm(1,1);
			perm.resize(B.nb_rows(),B.nb_rows());

			/*perform the gauss elimination for each comb_piv element*/
			while(comb_piv.size() > 0){
				aux_list.clear();
				/*Initialize the permutation matrix*/
				for (int i = 0; i<A.nb_rows() ; i++)
					for (int j = 0; j<A.nb_rows() ; j++){
						if (i == j) perm[i][j] = 1;
						else perm[i][j] = 0;
				}
				/*Initialize B*/
				B = A.mid();
				rows_checked.clear();
				for (int k = 0 ; k < comb_piv[comb_piv.size()-1].size() ; k++){
					int actual_col = comb_piv[comb_piv.size()-1][k];
					temp_piv = -1;
					for (int i = 0; i < B.nb_rows() ; i++)
						if (( (B[i][actual_col] < -prec) || (B[i][actual_col] > prec)) && (rows_checked.count(i) != 1)){
							rows_checked.insert(i);
							temp_piv = i;
							aux_list.push_back(make_pair(temp_piv,actual_col));
							break;
						}
					if (temp_piv==-1)
						break;
					else{
						double coef = B[temp_piv][actual_col];
						Matrix aux_perm(1,1);
						aux_perm.resize(A.nb_rows(),A.nb_rows());
						for (int m = 0; m<A.nb_rows() ; m++)
							for (int l = 0; l<A.nb_rows() ; l++){
								if (m == l) aux_perm[m][l] = 1;
								else aux_perm[m][l] = 0;
							}
						for (int m = 0 ; m < B.nb_rows() ; m++){
							if (m != temp_piv){
								double factor = B[m][actual_col];
								aux_perm[m][temp_piv] = -factor/coef;
								for (int l = 0 ; l < B.nb_cols() ; l++)
									B[m][l]	= B[m][l]-(B[temp_piv][l]*factor/coef);
							}
						}
						/*new:dejar con 1 test*/
						for (int i = 0 ; i < B.nb_cols() ; i++)
							B[temp_piv][i] = B[temp_piv][i]/coef;
						aux_perm[temp_piv][temp_piv] = 1/coef;
						/*end: dejar con 1 test*/
						perm = aux_perm*perm;
					}
				}

				/*If gauss is perform complete, add the matrix perm to the list permutation*/
				if (temp_piv!=-1){
					pair_contr_all.push_back(aux_list);
					permutations.push_back(perm);
				}
				comb_piv.pop_back();
			}
		}

} /* namespace ibex */
