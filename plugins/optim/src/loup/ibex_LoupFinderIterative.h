#ifndef __IBEX_LOUP_FINDER_ITERATIVE_H__
#define __IBEX_LOUP_FINDER_ITERATIVE_H__

#include "ibex_LoupFinder.h"
#include "ibex_System.h"
#include "ibex_LoupFinderXTaylor.h"
#include "ibex_LoupFinderAbsTaylor.h"
#include "ibex_Vector.h"

namespace ibex {

class LoupFinderIterative : public LoupFinder {

public:


	typedef enum  { XT, ABST, BOTH } loup_finders;

	/**
	 * \ingroup optim
	 *
	 * \brief An iterative algorithm based on XTaylor and/or AbsTaylor.
	 *
	 * The algorithm builds an inner (feasible) polytope inside the
	 * current box by using either AbsTaylor or XTaylor; and then minimizes a
	 * linear approximation of the goal function on this polytope via
	 * a LP solver. If the algorithm success, then it construct a new box (inside the
	 * search space), continuing the search of better upperbounds.
	 */
	LoupFinderIterative(const System& sys, const IntervalVector& initial_box,double alpha=0.9,loup_finders lfinders=BOTH,int max_iter=10);

	/**
	 * \brief Delete this.
	 */
	virtual ~LoupFinderIterative();

	/**
	 * \brief Find a new loup in a given box and the neighborhood by using
	 * a LoupFinder (e.g. AbsTaylor/XTaylor).
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);
	const System& sys;
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
	/*the loupfinder to be used*/
	loup_finders lfinders;
	/*maximum number of iterations*/
	int max_iter;

};

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_ITERATIVE_H__ */
