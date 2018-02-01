/*
 * ibex_Conditioners.h
 *
 *  Created on: 13-10-2017
 *      Author: victor
 */

#ifndef SRC_NUMERIC_IBEX_CONDITIONERS_H_
#define SRC_NUMERIC_IBEX_CONDITIONERS_H_

#include <list>
#include <set>
#include "ibex_IntervalMatrix.h"
#include "ibex_Matrix.h"
#include "ibex_Linear.h"

using namespace std;

namespace ibex {

    /*
     * \brief This function performs (n/m) Gauss-Jordan eliminations to the matrix A. All the row operations are stored
     * on the list perm_list. The pivot variables and the corresponding equation are stored in proj_vars in order to be
     * contracted by the projection operator.
     */
	void gauss_jordan (IntervalMatrix& A, IntervalVector x, vector<Matrix> & perm_list,
			vector <vector <pair <int,int> > > proj_vars , double prec=1e-7);



}



#endif /* SRC_NUMERIC_IBEX_CONDITIONERS_H_ */
