/*
 * ibex_CtcGauss.h
 *
 *  Created on: 06-02-2018
 *      Author: victor
 */

#ifndef SRC_CONTRACTOR_IBEX_CTCGAUSS_H_
#define SRC_CONTRACTOR_IBEX_CTCGAUSS_H_

using namespace std;
#include <list>
#include <vector>
#include <set>
#include "ibex_IntervalMatrix.h"
#include "ibex_System.h"
#include "ibex_Ctc.h"
#include "ibex_Conditioners.h"

namespace ibex {

	class GaussContractor : public Ctc {
	public:

		GaussContractor (const System& sys, int goal_var);

		/*
	    * \brief This function
	    */
		void contract(IntervalVector & ext_box);


		/*
	    * \brief This function
	    */
		void init_system(IntervalVector box, const System& sys);

		void write_ext_box(const IntervalVector& box, IntervalVector& ext_box);

		void read_ext_box(const IntervalVector& ext_box, IntervalVector& box);

		const System& sys;
		vector<IntervalMatrix> perm_list;
		IntervalMatrix A;
		IntervalVector b;
		/**
		 * \brief Index of the goal variable.
		 *
		 * See #ExtendedSystem.goal_var().
		 */
		const int goal_var;
	};


}

#endif /* SRC_CONTRACTOR_IBEX_CTCGAUSS_H_ */
