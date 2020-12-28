/* ============================================================================
 * I B E X - AbsTaylor linealization
 * ============================================================================
 *
 * Author(s)   : Ignacio Araya, Victor Reyes
 * Created     : April 2018
 * Updated     : December 2020
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_LOUP_FINDER_ABS_TAYLOR_H__
#define __IBEX_LOUP_FINDER_ABS_TAYLOR_H__

#include "ibex_LinearizerAbsTaylor.h"
#include "ibex_LoupFinder.h"

namespace ibex {


class LoupFinderAbsTaylor : public LoupFinder {
public:

	/**
	 * \ingroup optim
	 *
	 * \brief Upper-bounding algorithm based on AbsTaylor restriction.
	 *
	 * The algorithm builds an inner polytope inside the
	 * current box by using a Taylor form with absolute values.
	 * Then, it minimizes a linear approximation of the goal function
	 * on this polytope via a LP solver.
	 */
	LoupFinderAbsTaylor(const System& sys,bool trace = false);

	/**
	 * \brief Find a new loup in a given box.
	 *
	 * \see comments in LoupFinder.
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& exp_point, double current_loup);

	/**
	 * \brief The NLP problem.
	 */
	const System& sys;

protected:

	/** Linearization technique. */
	Linearizer* lr;
	/** A trace of the linealization and the loup (if exists) founded by the technique**/
	bool trace;
	/** linear solver */
	LPSolver lp_solver;

};

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_ABS_TAYLOR_H__ */