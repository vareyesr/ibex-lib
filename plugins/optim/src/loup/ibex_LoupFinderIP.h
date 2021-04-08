
#ifndef __IBEX_LOUP_FINDER_IP_H__
#define __IBEX_LOUP_FINDER_IP_H__

#include "ibex_LinearizerXTaylor.h"
#include "ibex_LinearizerAbsTaylor.h"


#include "ibex_LoupFinder.h"
#include "ibex_LPSolver.h"

namespace ibex {
/**
 * \ingroup optim
 *
 * \brief Upper-bounding algorithm based on XTaylor restriction.
 *
 * The algorithm builds an inner (feasible) polytope inside the
 * current box (see #LinearizerXTaylor) and then minimizes a
 * linear approximation of the goal function on this polytope via
 * a LP solver. The resulting point is verified a posteriori to
 * be feasible (wrt nonlinear constraint) and a new "loup".
 *
 * \note Only works with inequality constraints.
 */
class LoupFinderIP : public LoupFinder {

public:

	typedef enum  {XT, ABST} loup_finder;
	/**
	 * \brief Create the algorithm for a given system.
	 *
	 * \param sys         - The NLP problem.
	 */
	LoupFinderIP(const System& sys, loup_finder loup=XT);

	/**
	 * \brief Find a new loup in a given box.
	 *
	 * \see comments in LoupFinder.
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);

	/**
	 * \brief Find a new loup in a given box.
	 *
	 * \see comments in LoupFinder.
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup, BoxProperties& prop);

	/**
	 * \brief Add properties required by this loup finder.
	 */
	virtual void add_property(const IntervalVector& init_box, BoxProperties& prop);
	/**
	 * \brief Sets the expansion point for the AbsTaylor linearization.
	 */
	void set_expansion_point(Vector point){ exp_point=point; }
	/**
	 * \brief The NLP problem.
	 */
	const System& sys;
private:
	/*the loupfinder to be used*/
	loup_finder loup;
	/*the expansion point for abstaylor, by default the midpoint*/
	IntervalVector exp_point;
protected:

	/** Linearization technique. */
	Linearizer* lr;
	/** linear solver */
	LPSolver lp_solver;
};

/*============================================ inline implementation ============================================ */

inline std::pair<IntervalVector, double> LoupFinderIP::find(const IntervalVector& box, const IntervalVector& loup_point, double loup) {
	BoxProperties prop(box);
	return find(box, loup_point, loup, prop);
}

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_IP_H__ */
