/*
 * ibex_CtcDualP.h
 *
 *  Created on: Apr 13, 2021
 *      Author: victor
 */

#ifndef __IBEX_CTC_DUAL_P_H__
#define __IBEX_CTC_DUAL_P_H__

#include "ibex_Linearizer.h"
#include "ibex_Ctc.h"
#include "ibex_LPSolver.h"

namespace ibex {
class CtcDualP : public Ctc {
public:

	CtcDualP(Linearizer& lr, int max_iter=LPSolver::default_max_iter,
			int time_out=LPSolver::default_timeout, double eps=LPSolver::default_tolerance);

	/**
	 * \brief Delete this.
	 */
	virtual ~CtcDualP();

	/**
	 * \brief Contract the box.
	 *
	 * Linearize the system and performs 2n calls to Simplex in order to reduce
	 * the 2 bounds of each variable
	 */
	virtual void contract(IntervalVector& box);

	/**
	 * \brief Add linearizer properties to the map
	 */
	virtual void add_property(const IntervalVector& init_box, BoxProperties& map);
	/*
	 * \brief update the current preconditioner matrix
	 */
	void update_dual_sols();

protected:

	/**
	 * \brief The linearization technique
	 */
	Linearizer& lr;
	/*
	 * \brief box used in the last updated
	 */
	IntervalVector last_box;
	/**
	 * \brief  The linear solver that will be used
	 */
	LPSolver mylinearsolver;
	/*
	 * The dual solution, which is updated throw the search.
	 */
	Matrix dual_sols;
}
}

#endif // __IBEX_CTC_DUAL_P_H__
