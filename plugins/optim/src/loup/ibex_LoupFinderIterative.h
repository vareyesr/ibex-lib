#ifndef __IBEX_LOUP_FINDER_ITERATIVE_H__
#define __IBEX_LOUP_FINDER_ITERATIVE_H__

#include "ibex_LoupFinder.h"
#include "ibex_System.h"
#include "ibex_LoupFinderIP.h"
#include "ibex_Vector.h"

namespace ibex {

class LoupFinderIterative : public LoupFinder {

public:

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
	LoupFinderIterative(const System& sys, const IntervalVector& initial_box,std::vector<LoupFinderIP*> loup_finders, double alpha=0.9,int max_iter=10, double prec=1e-3);

	/**
	 * \brief Delete this.
	 */
	virtual ~LoupFinderIterative();

	/**
	 * \brief Find a new loup in a given box and the neighborhood by using
	 * a LoupFinder (e.g. AbsTaylor/XTaylor).
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);
	/**
	 * \brief To use for checking the progress of the upperbounds
	 */
	void set_trace(bool trace){this->trace = trace;}
	/**
	 * \brief For printing the ub, in case the trace function is true
	 */
	void print_ub(std::pair<IntervalVector,double> p);
	/**
	 * \brief Changes the size and position of the current search box. Always inside
	 * the root box.
	 * reduce and move
	 */
	void change_box_size(IntervalVector& box_aux, Vector old_exp);
	/**The system**/
	const System& sys;
private:

	/**LoupFinders**/
	std::vector<LoupFinderIP*> loup_finders;
	/**The initial box (search space)**/
	const IntervalVector& initial_box;
	/**User parameter for convergence purposes**/
	double alpha;
	/*maximum number of iterations*/
	int max_iter;
	/*the precision */
	double prec;
	/*trace, just for testing purposes*/
	bool trace;
};

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_ITERATIVE_H__ */
