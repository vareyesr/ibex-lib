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

		GaussContractor (const System& sys, IntervalVector & initial_box);

		/*
	    * \brief This function
	    */
		void contract(IntervalVector & box);


		/*
	    * \brief This function
	    */
		void init_system(IntervalVector box, const System& sys);



		const System& sys;
		vector<IntervalMatrix> perm_list;
		IntervalMatrix A;
		IntervalVector b;
	};


}

#endif /* SRC_CONTRACTOR_IBEX_CTCGAUSS_H_ */
