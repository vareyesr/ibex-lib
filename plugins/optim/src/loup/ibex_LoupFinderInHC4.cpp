//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderInHC4.cpp
// Author      : Gilles Chabert, Ignacio Araya, Bertrand Neveu
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : May 14, 2012
// Last Update : Jul 09, 2017
//============================================================================

#include "ibex_LoupFinderInHC4.h"
#include "ibex_LoupFinderProbing.h"

using namespace std;

namespace ibex {

LoupFinderInHC4::LoupFinderInHC4(const System& sys) : sys(sys) {
	mono_analysis_flag=true;
//	nb_inhc4=0;
//	diam_inhc4=0;
//	nb_rand=0;
//	diam_rand=0;
}

//void LoupFinderInHC4::report() {
//    if (trace) {
//      cout << " nbrand " << nb_rand << " nb_inhc4 " << nb_inhc4 << " nb simplex " << nb_simplex << endl;
//      cout << " diam_rand " << diam_rand << " diam_inhc4 " << diam_inhc4 << " diam_simplex " << diam_simplex << endl;
//    }
//}

std::pair<Vector, double> LoupFinderInHC4::find(const IntervalVector& box, const Vector& loup_point, double loup) {

	BitSet active=sys.active_ctrs(box);

	IntervalVector gx = sys.active_ctrs_eval(box);

	IntervalVector inbox=box;

	int c;
	bool inner_found=true;
	for (unsigned int i=0; i<active.size(); i++) {
		c=(i==0? active.min() : active.next(c));

		// Quick infeasibility check
		if (gx[i].lb()>0) throw NotFound();

		Interval right_cst;
		switch(sys.ops[c]) {
		case LT :
		case LEQ : right_cst=Interval::NEG_REALS; break;
		case EQ  : right_cst=Interval::ZERO;      break;
		case GEQ :
		case GT : right_cst=Interval::POS_REALS;  break;
		}

		sys.f_ctrs[c].ibwd(right_cst, inbox);

		if (box.is_empty()) {
			inner_found=false;
			break;
		}
	}

//	if (inner_found) {
//		cout <<  " inner box found:" << inbox << endl;
//		nb_inhc4++;
//		diam_inhc4 = ((nb_inhc4-1) *diam_inhc4 + inbox.max_diam()) / nb_inhc4;
//	} else {
//		nb_rand++;
//		diam_rand = ((nb_rand-1) *diam_rand + box.max_diam()) / nb_rand;
//	}

	// if we reach this point: "inbox" is an inner box
	if (mono_analysis_flag)
		monotonicity_analysis(sys, inbox, inner_found);

	return LoupFinderProbing(sys).find(inner_found? inbox : box,loup_point,loup);

}

} /* namespace ibex */
