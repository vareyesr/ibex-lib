#ifndef SRC_NUMERIC_IBEX_CONDITIONERS_H_
#define SRC_NUMERIC_IBEX_CONDITIONERS_H_

#include <list>
#include <set>
#include "ibex_IntervalMatrix.h"
#include "ibex_Matrix.h"
#include "ibex_Linear.h"
#include "ibex_Expr.h"
#include <algorithm>
#include "ibex_LPSolver.h"
#include "ibex_System.h"
#include <random>
using namespace std;



namespace ibex {

	double impact_heuristic_new(int row, int column,IntervalMatrix PA, IntervalVector x,int heuristic);
	pair<int,int> find_pivot_new(IntervalMatrix & PA, IntervalVector x,set<int> & ban_rows, set<int> & ban_cols,int heuristic);
	void best_gauss_jordan_new (IntervalMatrix A, IntervalVector x, vector<IntervalMatrix> & perm_list,
							vector <vector <pair <int,int> > > & proj_vars, double prec,int heuristic);

	Matrix best_P (IntervalMatrix& A, IntervalVector x, int size_b, double prec, int exteded);
	Matrix best_P_bound (IntervalMatrix& A, IntervalVector x, int size_b, double prec, int extended,int wbound);


}



#endif /* SRC_NUMERIC_IBEX_CONDITIONERS_H_ */

