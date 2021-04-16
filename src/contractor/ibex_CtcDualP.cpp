/*
 * ibex_CtcDualP.cpp
 *
 *  Created on: Apr 13, 2021
 *      Author: victor
 */

#include "ibex_CtcDualP.h"

#include "ibex_LinearizerFixed.h"

using namespace std;

namespace ibex {


CtcDualP::CtcDualP(Linearizer& lr, int max_iter, int time_out, double eps) :
		Ctc(lr.nb_var()), lr(lr),
		mylinearsolver(nb_var, LPSolver::Mode::Certified, eps, time_out, max_iter),
		last_box(nb_var),dual_sols(0,0) {

}

}

