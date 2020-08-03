/*
 * ibex_Conditioners.cpp
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#include "ibex_Conditioners.h"

using namespace std;

namespace ibex {

	double impact_heuristic_new(int row, int column,IntervalMatrix PA, IntervalVector x,int heuristic){

			if (heuristic == 1){
				return std::abs(PA[row][column].mag());
			}

			else if (heuristic == 2) {
				return std::abs(PA[row][column].mag()*x[column].diam());
			}

			else if (heuristic == 3) {
				return std::abs(PA[row][column].mag());
			}
			else if (heuristic == 4) {
				return std::abs(PA[row][column].mag()*x[column].diam());
			}


			else return 0;
		}

	pair<int,int> find_pivot_new(IntervalMatrix & PA, IntervalVector x,set<int> & ban_rows, set<int> & ban_cols,int heuristic){
		pair<int,int> position;
		position.first = -1;
		position.second = -1;
		Matrix impact_values(PA.nb_rows(),PA.nb_cols());
		for (int i = 0 ; i < impact_values.nb_rows(); i++)
			for (int j = 0 ; j < impact_values.nb_cols() ; j++)
				impact_values[i][j] = 0;
		for (int j = 0 ; j < PA.nb_cols() ; j++){
			if (ban_cols.count(j) != 1){
				for (int i = 0 ; i < PA.nb_rows() ; i++){
					if ((ban_rows.count(i) != 1) &&  !(PA[i][j].contains(0))){
						impact_values[i][j] = impact_heuristic_new(i,j,PA,x,heuristic);
					}
				}
			}
		}
		/*search for the maximum value*/
		double max_value = 1e-8;
		if ((heuristic == 1) || (heuristic == 2))
			for (int i = 0 ; i < impact_values.nb_rows(); i++)
				for (int j = 0 ; j < impact_values.nb_cols() ; j++){
					if (impact_values[i][j] > max_value){
						max_value = impact_values[i][j];
						position.first = j;
						position.second = i;

					}
				}
		if ((heuristic == 3) || (heuristic == 4)){
			int best_col = -1;
			for (int i = 0 ; i < impact_values.nb_rows();i++){
				double sum = 0;
				for (int j =0 ; j < impact_values.nb_cols();j++)
					sum = sum + impact_values[i][j];
				if (sum != 0){
					for (int j =0 ; j < impact_values.nb_cols();j++)
						impact_values[i][j] = (double)(impact_values[i][j]/sum);
				}
			}
			Matrix impact_relative(1,PA.nb_cols());
			for (int j = 0 ; j < PA.nb_cols() ; j++)
				for (int i = 0 ; i < PA.nb_rows() ; i++)
					impact_relative[0][j] = impact_relative[0][j]+impact_values[i][j];
			for (int j = 0 ; j < PA.nb_cols() ; j++){
				if (impact_relative[0][j] > max_value){
					max_value = impact_relative[0][j];
					best_col = j;
				}
			}

			if (best_col == -1) return position;
			else{
				max_value = 1e-8 ;
				for (int i = 0 ; i < impact_values.nb_rows() ; i++){
					if (impact_values[i][best_col] > max_value){
						max_value = impact_values[i][best_col];
						position.first = best_col;
						position.second = i;

					}
				}
			}
		}




		if(position.first != -1){
			ban_rows.insert(position.second);
			ban_cols.insert(position.first);
		}
		return position;
	}


	void best_gauss_jordan_new (IntervalMatrix A, IntervalVector x, vector<IntervalMatrix> & perm_list,
						vector <vector <pair <int,int> > > & proj_vars, double prec,int heuristic){

		vector <pair<int,int> > aux_list;
		IntervalMatrix B(1,1);
		B.resize(A.nb_rows(),A.nb_cols());
		IntervalMatrix perm(1,1);
		set<int> ban_rows;
		set<int> ban_cols;
		perm.resize(B.nb_rows(),B.nb_rows());
		pair<int,int> max_values;
		bool available_cols = true;
		B = A;
		while (available_cols){

			aux_list.clear();
			ban_rows.clear();
			/*Initialize B*/
			B = A;
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
//				pair<int,int> var_eq = find_next_pivot(B, x, ban_rows, ban_cols,heuristic);
				pair<int,int> var_eq = find_pivot_new(B, x, ban_rows, ban_cols,heuristic);
				if (var_eq.first !=-1){
					/*dddd*/
					aux_list.push_back(make_pair(var_eq.first,var_eq.second));
					Interval coef = B[var_eq.second][var_eq.first];
					IntervalMatrix aux_perm(1,1);
					aux_perm.resize(A.nb_rows(),A.nb_rows());
					for (int k = 0; k<A.nb_rows() ; k++)
						for (int l = 0; l<A.nb_rows() ; l++){
							if (k == l) aux_perm[k][l] = 1;
							else aux_perm[k][l] = 0;
						}
					for (int k = 0 ; k < B.nb_rows() ; k++){
						if ((k != var_eq.second) && (  (B[k][var_eq.first].mag() > prec))) {
							Interval factor = B[k][var_eq.first];
							aux_perm[k][var_eq.second] = -factor/coef;
							for (int l = 0 ; l < B.nb_cols() ; l++){
								if (l == var_eq.first) B[k][l] = 0;
								else B[k][l] = B[k][l]-(B[var_eq.second][l]*factor/coef);

							}
						}
					}
					/*make the pivot position 1*/
					for (int i = 0 ; i < B.nb_cols() ; i++){
						if (i == var_eq.first) B[var_eq.second][i] = 1;
						else B[var_eq.second][i] = B[var_eq.second][i]/coef;
					}
					aux_perm[var_eq.second][var_eq.second] = 1/coef;
					perm = aux_perm*perm;

				}
				else{
					available_cols = false;
					break;
				}
			}
			if (aux_list.size()>0){
				proj_vars.push_back(aux_list);
				perm_list.push_back(perm);
				 available_cols = false;
//					cout << perm.nb_cols() << " " << perm.nb_rows() << endl;
			}
			if(A.nb_cols()==ban_cols.size()) available_cols = false;
		}
	}

	Matrix best_P (IntervalMatrix& A, IntervalVector x, int size_b, double prec, int extended){

		if (extended) extended = size_b;

		Matrix perm(A.nb_cols()-extended,A.nb_rows());

		double max =0;
		bool found;
		IntervalVector box2(2*A.nb_cols()+A.nb_rows());
		vector <pair<double,Vector > > P_rows;
		//Vector row(box2.size());
		double suma[A.nb_cols()];
		int nb_eq=0;
		int lol = 0;
		for (int i = 0 ; i < A.nb_cols()-extended ; i++){
//			LPSolver lp_solver(box2.size(), 1000); // p, s y w
			LPSolver lp_solver(box2.size());
			lp_solver.reset(box2.size());
			lp_solver.set_max_iter(10000000000000000);
			//lp_solver.set_bounds(IntervalVector(row.size()));

			//initializing domains and objective function
			for (int j = 0 ; j < box2.size() ; j++){
				if (j < A.nb_rows()){
					//initialize the variables p_i
					box2[j] = Interval(-1e7,1e7);
					lp_solver.set_cost(j,0);
				}
				else if (j >= A.nb_rows() && j < A.nb_rows() + A.nb_cols()){
					//initialize variables s_i
					box2[j]=Interval(-1e7, 1e7);
					lp_solver.set_cost(j,0);
				}
				else{
					//initialize auxiliary variables w_i
					box2[j]=Interval(-1e7, 1e7);
//					lp_solver.set_obj_var(j,0);
					lp_solver.set_cost(j,(x[j- A.nb_rows() -A.nb_cols()].diam()));
//					lp_solver.set_obj_var(j,(suma[j- A.nb_rows() -A.nb_cols()])*(x[j- A.nb_rows() -A.nb_cols()].diam()));
				}
			}
			lp_solver.set_bounds(box2);
			Vector row(box2.size());
			//s_i = 1
			row[A.nb_rows()+i] = 1;
			lp_solver.add_constraint(row,GEQ,1);
			lp_solver.add_constraint(row,LEQ,1);

			// w_i - s_i > ; w_i + s_i > 0
			for (int j = 0 ; j < A.nb_cols() ; j++){
				Vector row(box2.size());
				row[A.nb_rows() + A.nb_cols()+j] = 1;
				row[A.nb_rows()+j]=-1;
				lp_solver.add_constraint(row,GEQ,0);

				Vector row2(box2.size());
				row2[A.nb_rows() + A.nb_cols()+j] = 1;
				row2[A.nb_rows()+j]=1;
				lp_solver.add_constraint(row2,GEQ,0);
			}
			//s_j - sum_{j=1} A[i][j] p_i = 0
			for (int j = 0 ; j < A.nb_cols() ; j++){
				Vector row(box2.size());
				for (int ii = 0 ; ii < A.nb_rows() ; ii++)
					row[ii] = -A[ii][j].mid();
				row[A.nb_rows()+j]=1;
				lp_solver.add_constraint(row,GEQ,0);
				lp_solver.add_constraint(row,LEQ,0);
			}

			LPSolver::Status stat = lp_solver.minimize();

			Vector v(box2.size());
			v = lp_solver.not_proved_primal_sol();
			v.resize(A.nb_rows());

			found = false;
			for (int j = 0 ; j < nb_eq ; j++)
				if (v==perm[j]) {found = true; break;}

			if (found == false){
				perm.put(nb_eq,0,v,true);
				nb_eq++;
			}

		}

		for (int i = 0 ; i < perm.nb_rows() ; i++)
				for (int j = 0 ; j < perm.nb_cols() ; j++)
					if (std::abs(perm[i][j])<1e-7)  perm[i][j]=0;
		perm.resize(nb_eq,A.nb_rows());
		for (int i = 0 ; i < A.nb_rows() ; i++)
			for (int j = 0 ; j < A.nb_cols() ; j++)
				if (std::abs(A[i][j].mid())<1e-7)  A[i][j]=0;
		return perm;
	}

	Matrix best_P_bound (IntervalMatrix& A, IntervalVector x, int size_b, double prec, int extended,int wbound){

		if (extended) extended = size_b;

		Matrix perm(A.nb_cols()-extended,A.nb_rows());

		double max =0;
		bool found;
		IntervalVector box2(2*A.nb_cols()+A.nb_rows());
		vector <pair<double,Vector > > P_rows;
		//Vector row(box2.size());
		double suma[A.nb_cols()];

		int nb_eq=0;

		for (int i = 0 ; i < A.nb_cols()-extended ; i++){
			LPSolver lp_solver(box2.size()); // p, s y w

			lp_solver.reset(box2.size());
			lp_solver.set_max_iter(10000000000000000);
			//lp_solver.set_bounds(IntervalVector(row.size()));

			//initializing domains and objective function
			if (wbound == 0)
				for (int j = 0 ; j < box2.size() ; j++){
					if (j < A.nb_rows()){
						//initialize the variables p_i
						box2[j] = Interval(-1e7,1e7);
						lp_solver.set_cost(j,0);
					}
					else if (j >= A.nb_rows() && j < A.nb_rows() + A.nb_cols()){
						//initialize variables s_i
						box2[j]=Interval(-1e7, 1e7);
						lp_solver.set_cost(j,0);
					}
					else{
						//initialize auxiliary variables w_i
						box2[j]=Interval(-1e7, 1e7);
						lp_solver.set_cost(j,-1);
	//					lp_solver.set_obj_var(j,-(x[j- A.nb_rows() -A.nb_cols()].diam()));
	//					lp_solver.set_obj_var(j,(suma[j- A.nb_rows() -A.nb_cols()])*(x[j- A.nb_rows() -A.nb_cols()].diam()));
					}
				}
			else if (wbound == 1)
				for (int j = 0 ; j < box2.size() ; j++){
					if (j < A.nb_rows()){
						//initialize the variables p_i
						box2[j] = Interval(-1e7,1e7);
						lp_solver.set_cost(j,0);
					}
					else if (j >= A.nb_rows() && j < A.nb_rows() + A.nb_cols()){
						//initialize variables s_i
						box2[j]=Interval(-1e7, 1e7);
						lp_solver.set_cost(j,0);
					}
					else{
						//initialize auxiliary variables w_i
						box2[j]=Interval(-1e7, 1e7);
						lp_solver.set_cost(j,1);
	//					lp_solver.set_obj_var(j,-(x[j- A.nb_rows() -A.nb_cols()].diam()));
	//					lp_solver.set_obj_var(j,(suma[j- A.nb_rows() -A.nb_cols()])*(x[j- A.nb_rows() -A.nb_cols()].diam()));
					}
				}
			lp_solver.set_bounds(box2);
			Vector row(box2.size());
			//s_i = 1
			row[A.nb_rows()+i] = 1;
			lp_solver.add_constraint(row,GEQ,1);
			lp_solver.add_constraint(row,LEQ,1);

			// w_j-s_j*lb(x_j) <= 0 ; w_j-s_j*ub(x_j) <= 0
			if (wbound == 0){
				for (int j = 0 ; j < A.nb_cols() ; j++){
					Vector row(box2.size());
					row[A.nb_rows() + A.nb_cols()+j] = 1;
					row[A.nb_rows()+j]=-1*x[j].lb();
					/*nuevo*/
//						for (int ii = 0 ; ii < A.nb_rows() ; ii++)
//							row[ii] = -x[j].lb()*A[ii][j].mid();
					lp_solver.add_constraint(row,LEQ,0);
					Vector row2(box2.size());
					row2[A.nb_rows() + A.nb_cols()+j] = 1;
					row2[A.nb_rows()+j]=-1*x[j].ub();
					/*nuevo*/
//						for (int ii = 0 ; ii < A.nb_rows() ; ii++)
//							row2[ii] = -x[j].ub()*A[ii][j].mid();
					lp_solver.add_constraint(row2,LEQ,0);
				}
			}
			// w_j-s_j*lb(x_j) >= 0 ; w_j-s_j*ub(x_j) >= 0
			else if (wbound == 1){
				for (int j = 0 ; j < A.nb_cols() ; j++){
					Vector row(box2.size());
					row[A.nb_rows() + A.nb_cols()+j] = 1;
					row[A.nb_rows()+j]=-1*x[j].lb();
					/*nuevo*/
//						for (int ii = 0 ; ii < A.nb_rows() ; ii++)
//							row[ii] = -x[j].lb()*A[ii][j].mid();
					lp_solver.add_constraint(row,GEQ,0);

					Vector row2(box2.size());
					row2[A.nb_rows() + A.nb_cols()+j] = 1;
					row2[A.nb_rows()+j]=-1*x[j].ub();
					/*nuevo*/
//						for (int ii = 0 ; ii < A.nb_rows() ; ii++)
//							row2[ii] = -x[j].ub()*A[ii][j].mid();
					lp_solver.add_constraint(row2,GEQ,0);
				}
			}
//				s_j - sum_{j=1} A[i][j] p_i = 0
			for (int j = 0 ; j < A.nb_cols() ; j++){
				Vector row(box2.size());
				for (int ii = 0 ; ii < A.nb_rows() ; ii++)
					row[ii] = -A[ii][j].mid();
				row[A.nb_rows()+j]=1;
				lp_solver.add_constraint(row,GEQ,0);
				lp_solver.add_constraint(row,LEQ,0);
			}
			/*nueva parte*/
			Vector rowk(box2.size());
			for (int ii = 0 ; ii < A.nb_cols() ; ii++)
				rowk[ii+A.nb_rows()+A.nb_cols()] = 1;
			if (wbound == 0)
				lp_solver.add_constraint(rowk,LT,0);
			else
				lp_solver.add_constraint(rowk,GT,0);


			LPSolver::Status stat = lp_solver.minimize();
			Vector v(box2.size());
			v = lp_solver.not_proved_primal_sol();
			found = false;

			v.resize(A.nb_rows());
//				cout << lp_solver.get_obj_value() << endl;

			for (int j = 0 ; j < nb_eq ; j++)
				if (v==perm[j]) {found = true; break;}

			if (found == false){
				perm.put(nb_eq,0,v,true);
				nb_eq++;
			}

		}

		for (int i = 0 ; i < perm.nb_rows() ; i++)
				for (int j = 0 ; j < perm.nb_cols() ; j++)
					if (std::abs(perm[i][j])<1e-7)  perm[i][j]=0;
		perm.resize(nb_eq,A.nb_rows());
		for (int i = 0 ; i < A.nb_rows() ; i++)
			for (int j = 0 ; j < A.nb_cols() ; j++)
				if (std::abs(A[i][j].mid())<1e-7)  A[i][j]=0;
		return perm;
	}
} /* namespace ibex */

