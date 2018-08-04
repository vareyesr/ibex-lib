//============================================================================
//                                  I B E X                                   
// File        : Composition of contractors
// Author      : Gilles Chabert
// Copyright   : Ecole des Mines de Nantes (France)
// License     : See the LICENSE file
// Created     : Apr 25, 2012
// Last Update : Apr 25, 2012
//============================================================================

#include "ibex_CtcCompo.h"

namespace ibex {

/*! Default composition ratio. */
const double CtcCompo::default_ratio = 0.1;

void CtcCompo::init_impacts() {

//	impacts = new BoolMask[list.size()];
//
//	for (int i=0; i<list.size(); i++) {
//		assert(list[i].nb_var==nb_var);
//		impacts[i].resize(nb_var);
//		impacts[i].set_all();
//	}
}

CtcCompo::CtcCompo(const Array<Ctc>& list, bool incremental, double ratio) :
		Ctc(list), list(list), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, bool incremental, double ratio) :
		Ctc(c1.nb_var),list(Array<Ctc>(c1,c2)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, bool incremental, double ratio) :
		Ctc(c1.nb_var),list(Array<Ctc>(c1,c2,c3)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}


CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, Ctc& c12, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, Ctc& c12, Ctc& c13, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, Ctc& c12, Ctc& c13, Ctc& c14, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, Ctc& c12, Ctc& c13, Ctc& c14, Ctc& c15, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, Ctc& c12, Ctc& c13, Ctc& c14, Ctc& c15, Ctc& c16, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, Ctc& c12, Ctc& c13, Ctc& c14, Ctc& c15, Ctc& c16, Ctc& c17, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, Ctc& c12, Ctc& c13, Ctc& c14, Ctc& c15, Ctc& c16, Ctc& c17, Ctc& c18, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, Ctc& c12, Ctc& c13, Ctc& c14, Ctc& c15, Ctc& c16, Ctc& c17, Ctc& c18, Ctc& c19, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::CtcCompo(Ctc& c1, Ctc& c2, Ctc& c3, Ctc& c4, Ctc& c5, Ctc& c6, Ctc& c7, Ctc& c8, Ctc& c9, Ctc& c10, Ctc& c11, Ctc& c12, Ctc& c13, Ctc& c14, Ctc& c15, Ctc& c16, Ctc& c17, Ctc& c18, Ctc& c19, Ctc& c20, bool incremental, double ratio) :
		Ctc(c1.nb_var),	list(Array<Ctc>(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20)), incremental(incremental), ratio(ratio) {
	assert(check_nb_var_ctc_list(list));

	init_impacts();
}

CtcCompo::~CtcCompo() {
//	delete[] impacts;
}

void CtcCompo::add_property(const IntervalVector& init_box, BoxProperties& map) {
	for (int i=0; i<list.size(); i++)
		list[i].add_property(init_box, map);
}

void CtcCompo::contract(IntervalVector& box) {
	ContractContext context;
	contract(box,context);
}

void CtcCompo::contract(IntervalVector& box, ContractContext& context) {

	bool inactive = true;

	BitSet* flags = context.output_flags();

	BitSet impact(BitSet::all(nb_var)); // always set to "all" for the moment (to be improved later)

	BitSet* input_impact = context.impact();
	context.set_impact(&impact);

	// TODO: a more clever impact handling could be
	// done here
	for (int i=0; i<list.size(); i++) {

		if (flags && inactive) {

			flags->clear();

			list[i].contract(box, context);

			if (!(*flags)[ContractContext::INACTIVE]) {
				inactive=false;
				// now, no need for asking sub-contractors
				// to calculate the impact:
				context.set_output_flags(NULL);
			}
		} else {
			list[i].contract(box, context);
		}

		if (box.is_empty()) {
			if (flags) {
				flags->clear();
				flags->add(ContractContext::FIXPOINT);
				context.set_output_flags(flags); // restore!
			}
			context.set_impact(input_impact); // restore!
			return;
		}
	}

	if (flags) {
		if (inactive) flags->add(ContractContext::INACTIVE);
		context.set_output_flags(flags); // restore!
	}
	context.set_impact(input_impact); // restore!
}

} // end namespace ibex
