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


	LoupFinderTrustRegion(const System& sys, const IntervalVector& initial_box,double alpha);

	/**
	 * \brief Delete this.
	 */
	virtual ~LoupFinderTrustRegion();

	/**
	 * \brief Find a new loup in a given box.
	 *
	 * \see comments in LoupFinder.
	 */
	virtual std::pair<IntervalVector, double> find(const IntervalVector& box, const IntervalVector& loup_point, double loup);

	/**
	 * Loup finder using inner polytopes.
	 */

private:

	LoupFinderAbsTaylor finder_abs_taylor;
	LoupFinderXTaylor finder_x_taylor;
	const IntervalVector& initial_box;
	const Function* f_goal;
	double alpha;

};

} /* namespace ibex */

#endif /* __IBEX_LOUP_FINDER_DEFAULT_H__ */
