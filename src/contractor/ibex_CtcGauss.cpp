/*
 * ibex_CtcGauss.cpp
 *
 *  Created on: 06-02-2018
 *      Author: victor
 */

#include "ibex_CtcGauss.h"

using namespace std;

namespace ibex {

	GaussContractor::GaussContractor (const System& sys) : sys(sys), Ctc(sys.ctrs.size()), A(1,1), b(1){



	}

	void GaussContractor::contract(IntervalVector& box){



	}

	void GaussContractor::linearization(IntervalVector box, const System sys){



	}


} /* namespace ibex */

