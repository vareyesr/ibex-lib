/* ============================================================================
 * I B E X - X-Taylor linear relaxation
 * ============================================================================
 * Copyright   : IMT Atlantique (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file COPYING.LESSER.
 *
 * Author(s)   : Gilles Chabert, Ignacio Araya, Bertrand Neveu
 * Created     : July 01th, 2012
 * Updated     : June 23th, 2017
 * ---------------------------------------------------------------------------- */

#ifndef __IBEX_LINEAR_RELAX_X_TAYLOR__
#define __IBEX_LINEAR_RELAX_X_TAYLOR__

#include "ibex_NormalizedSystem.h"
#include "ibex_LinearRelax.h"

#include <vector>

namespace ibex {

/**
 * \ingroup numeric
 *
 * \brief X-Taylor linear relaxation technique
 *
 * This class is an implementation of the X-Taylor algorithm
 */

class LinearRelaxXTaylor : public LinearRelax {

public:

	/**
	 * \brief Slope matrices computation mode.
	 *
	 * TAYLOR: one gradient evaluation for all variables simultaneously (less accurate but faster)
	 * HANSEN: recursive evaluation with successive variable instantiation (more accurate but longer)
	 */
	typedef enum  { TAYLOR, HANSEN } linear_mode;

	/**
	 * \brief Strategies for selecting corner points.
	 *
	 * INF_X - Lower-left corner
	 * SUP_X - Upper-right corner
	 * RANDOM - One random point
	 * RANDOM_INV - The opposite of the random point (to be placed after)
	 */
	typedef enum  {INF_X, SUP_X, RANDOM, RANDOM_INV } corner_point;

	/** Default max_diam_deriv value, set to 1e6  **/
	static const double default_max_diam_deriv;

	/**
	 * \brief Creates the X_Taylor contractor.
	 *
	 * \param sys             - The system (extended or not).
	 * \param cpoints         - The vector of corner selection in linearization (X_INF, X_SUP, RANDOM, RANDOM_INV)
	 * \param lmode           - TAYLOR | HANSEN : linear relaxation method.
	 * \param max_diam_deriv  - The maximum diameter of the box for for the linear solver (default value 1.e6).
	 * 	  				        Soplex may lose solutions when it is called with "big" domains.
	 */
	LinearRelaxXTaylor(const System& sys, std::vector<corner_point>& cpoints,
			linear_mode lmode=HANSEN, double max_diam_deriv=default_max_diam_deriv);

	/**
	 * \brief Deletes this instance.
	 */
	~LinearRelaxXTaylor();


	/**
	 * \brief Generation of the linear inequalities
	 */
	int linearization(const IntervalVector& box, LinearSolver& lp_solver);

private:
	/**
	 * \brief The system
	 */
	const System& sys;

	/**
	 * \brier Number of variables
	 */
	int n;

	/**
	 * \brief Number of (real-valued) constraints
	 */
	int m;

	/**
	 * \brief Goal constraint (in case of extended system, -1 otherwise).
	 */
	const int goal_ctr;

	/**
	 * \brief The vector of corner selection in linearization
	 *
	 * Can be either X_INF, X_SUP, RANDOM, or RANDOM_INV.
	 */
	std::vector<corner_point> cpoints;

	/** Maximum diameter of the derivatives for calling linear solver (default value 1.e5) */
	double max_diam_deriv;

	/** TAYLOR | HANSEN : the linear relaxation method */
	linear_mode lmode;

	/**
	 * \brief Get the corner information "inf".
	 *
	 * \param inf - (output) whether the jth variable is set to lower bound (true)
	 *              or upper bound (false)
	 */
	void get_corner(corner_point cpoint, bool* inf);

	/**
	 * \brief Get the point corresponding to "inf"
	 */
	IntervalVector get_corner_point(const IntervalVector& box, bool* inf);

	/**
	 * \brief Add the constraint ax<=b in the LP solver.
	 */
	bool check_and_add_constraint(const IntervalVector& box, const Vector& a, double b, LinearSolver& lp_solver);

};

} // end namespace ibex

#endif /* __IBEX_LINEAR_RELAX_X_TAYLOR__ */

