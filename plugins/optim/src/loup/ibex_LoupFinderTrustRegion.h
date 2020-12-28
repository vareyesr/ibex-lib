#ifndef __IBEX_LOUP_FINDER_TRUST_REGION_H__
#define __IBEX_LOUP_FINDER_TRUST_REGION_H__

#include "ibex_LoupFinder.h"
#include "ibex_System.h"
#include "ibex_LoupFinderXTaylor.h"
#include "ibex_LoupFinderAbsTaylor.h"
#include "ibex_Vector.h"

namespace ibex {

class LoupFinderTrustRegion : public LoupFinder {

public:
	/**
	 * \ingroup optim
	 *
	 * \brief A trust region algorithm based on XTaylor and/or AbsTaylor.
	 *
	 * The algorithm builds an inner (feasible) polytope inside the
	 * current box by using either AbsTaylor or XTaylor; and then minimizes a
	 * linear approximation of the goal function on this polytope via
	 * a LP solver. If the algorithm success, then it construct a new box (inside the
	 * search space), continuing the search of better upperbounds.
	 */

	LoupFinderTrustRegion(const System& sys, const IntervalVector& initial_box,double alpha);

	/**
	 * \brief Delete this.
	 */
	virtual ~LoupFinderTrustRegion();

	/**
	 * \brief Find a new loup in a given box and the neighborhood by using
	 * XTaylor and/or AbsTaylor
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);

private:
	/**LoupFinder AbsTaylor (if needed)**/
	LoupFinderAbsTaylor finder_abs_taylor;
	/**LoupFinder XTaylor (if needed)**/
	LoupFinderXTaylor finder_x_taylor;
	/**The initial box (search space)**/
	const IntervalVector& initial_box;
	/**the goal function**/
	const Function* f_goal;
	/**User parameter for convergence purposes**/
	double alpha;

};

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_DEFAULT_H__ */
