/* ============================================================================
 * I B E X - AbsTaylor linear restriction
 * ============================================================================
 *
 * Author(s)   : Ignacio Araya, Victor Reyes
 * Created     : April 2018
 * Updated     : December 2020
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_LINEARIZER_ABS_TAYLOR__
#define __IBEX_LINEARIZER_ABS_TAYLOR__

#include "ibex_Linearizer.h"
#include "ibex_System.h"

namespace ibex {

/**
 * \ingroup numeric
 *
 * \brief Abs-Taylor linearization technique.
 *
 */
class LinearizerAbsTaylor : public Linearizer {

public:


	/**
	 * \brief Creates the AbsTaylor linearizer.
	 *
	 * \param sys The system
	 */
	LinearizerAbsTaylor(const System& sys, bool trace = false);

	/**
	 * \brief Deletes this.
	 */
	~LinearizerAbsTaylor();

	/**
	 * \brief Generation of the linear inequalities
	 */
	virtual int linearize(const IntervalVector& box, LPSolver& lp_solver);
	/**
	 * \brief Sets the expansion point of the method
	 *
	 * \param point a vector of real values
	 */
	void set_expansion_point(Vector point){ exp_point=point; }

private:

	/**
	 * \brief Linearization (RESTRICT mode)
	 */
	int linear_restrict(const IntervalVector& box);

	/**
	 * \brief Linearize a constraint g(x)<=0 inside a box, from the midpoint.
	 *
	 * \param dg_box:   dg([box])
	 * \param g_mid: g(mid)
	 */
	int linearize_leq_mid(const IntervalVector& box, const Vector& point, const IntervalVector& dg_box, const Interval& g_mid);

	/**
	 * \brief Add the constraint ax<=b in the LP solver.
	 */
	int check_and_add_constraint(const IntervalVector& box, const Vector& a, double b);

	/**
	 * \brief The system
	 */
	const System& sys;

	/**
	 * \brief Number of (real-valued) constraints
	 */
	int m;

	/**
	 * \brief Goal constraint (in case of extended system, -1 otherwise).
	 */
	const int goal_ctr;

	Vector exp_point;
	bool trace;

	/**
	 * Current LP solver
	 */
	LPSolver* lp_solver;

};

} // end namespace ibex

#endif /* __IBEX_LINEARIZER_ABS_TAYLOR__ */
