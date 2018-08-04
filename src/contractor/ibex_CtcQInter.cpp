//============================================================================
//                                  I B E X                                   
// File        : Q-intersection contractor
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Apr 30, 2012
// Last Update : Apr 30, 2012
//============================================================================

#include "ibex_CtcQInter.h"
#include "ibex_QInter.h"

namespace ibex {

CtcQInter::CtcQInter(const Array<Ctc>& list, int q) : Ctc(list), list(list), q(q), boxes(list.size(), nb_var) { }

void CtcQInter::add_property(const IntervalVector& init_box, BoxProperties& map) {
	for (int i=0; i<list.size(); i++)
		list[i].add_property(init_box, map);
}

void CtcQInter::contract(IntervalVector& box) {
	ContractContext context;
	contract(box,context);
}

void CtcQInter::contract(IntervalVector& box, ContractContext& context) {
	Array<IntervalVector> refs(list.size());

	// --------------------- context ------------------
	ContractContext c_context;
	c_context.set_impact(context.impact());
	// ------------------------------------------------

	for (int i=0; i<list.size(); i++) {
		boxes[i]=box;

		if (context.prop())
			c_context.set_properties(new BoxProperties(boxes[i], *context.prop()));

		list[i].contract(boxes[i], c_context);

		refs.set_ref(i,boxes[i]);

		if (context.prop())
			delete c_context.prop();
	}

	box = qinter(refs,q);

	if (context.prop()) {
		context.prop()->update(BoxEvent(box,BoxEvent::CONTRACT));
	}
}

} // end namespace ibex
