/*
 * ibex_CtcDualP.h
 *
 *  Created on: Apr 13, 2021
 *      Author: victor
 */

#ifndef __IBEX_CTC_DUAL_P_H__
#define __IBEX_CTC_DUAL_P_H__

#include "ibex_LinearizerXTaylor.h"
#include "ibex_Ctc.h"
#include "ibex_LPSolver.h"

namespace ibex {
class CtcDualP : public Ctc {
public:

	CtcDualP(LinearizerXTaylor& lr, IntervalVector initial_box);

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
	virtual void add_property(const IntervalVector& init_box, BoxProperties& map){};
	/*
	 * \brief update the current preconditioner matrix
	 */
	void update_dual_sols(IntervalVector box);
	void compute_dual(Matrix A, Vector b, IntervalVector box);

protected:

	/**
	 * \brief The linearization technique
	 */
	LinearizerXTaylor& lr;
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
	Matrix A_input;
	Vector b_input;
};
}

#endif // __IBEX_CTC_DUAL_P_H__
