/*
 * ibex_Conditioners.h
 *
 *  Created on: 01-02-2018
 *      Author: victor
 */

#ifndef SRC_NUMERIC_IBEX_CONDITIONERS_H_
#define SRC_NUMERIC_IBEX_CONDITIONERS_H_

#include <list>
#include <set>
#include <algorithm>
#include "ibex_IntervalMatrix.h"
#include "ibex_Matrix.h"
#include "ibex_Linear.h"
#include "ibex_IntervalVector.h"
#include "ibex_Interval.h"

using namespace std;

namespace ibex {


	/*TODO: add comments combinatorial function/*

	/*
    * \brief .
    */

	void combinatorial(IntervalMatrix A, int cols,int rows,std::vector< std::vector <int> > & comb_piv);

	/*TODO: add comments compare function/*

	/*
    * \brief .
    */

	bool compare(const std::pair<double, std::pair<int,int> >&i, const std::pair<double, std::pair<int,int> >&j);

	/*TODO: explicar bien la regla de seleccion/*

	/*
    * \brief This function finds the next variable from A , not included in ban_cols, to be the next variable to be pivot.
    * We use the set ban_rows to avoid the selection of the same row (as can be selected once).
    * This selection is perform by using the rule..
    */

	pair<int,int> find_next_pivot(IntervalMatrix & A, IntervalVector x,set<int> & ban_rows, set<int> & ban_cols);

    /*
     * \brief This function performs (n/m) Gauss-Jordan eliminations to the matrix A. All the row operations are stored
     * on the list perm_list. The pivot variables and the corresponding equation are stored in proj_vars in order to be
     * contracted by the projection operator.
     */
	void best_gauss_jordan (IntervalMatrix A, IntervalVector x, vector<Matrix> & perm_list, vector <vector <pair <int,int> > > proj_vars, double prec);

	 /*
	  * \brief This function performs all the possible Gauss-Jordan eliminations to the matrix A. All the row operations are stored
	  * on the list perm_list. The pivot variables and the corresponding equation are stored in proj_vars in order to be
	  * contracted by the projection operator.
	 */
	void all_gauss_jordan (IntervalMatrix A, IntervalVector x, vector<Matrix> & perm_list,vector <vector <pair <int,int> > > proj_vars , double prec);
}



#endif /* SRC_NUMERIC_IBEX_CONDITIONERS_H_ */
