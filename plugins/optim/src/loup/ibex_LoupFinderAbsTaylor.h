/* ============================================================================
 * I B E X - AbsTaylor linear restriction
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
	 * \brief Create the algorithm for a given system.
	 *
	 * \param sys         - The NLP problem.
	 */
	LoupFinderAbsTaylor(const System& sys);

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

	/** linear solver */
	LPSolver lp_solver;

};

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_ABS_TAYLOR_H__ */
