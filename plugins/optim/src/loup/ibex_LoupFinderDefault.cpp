//============================================================================
//                                  I B E X
// File        : ibex_LoupFinderDefault.cpp
// Author      : Gilles Chabert
// Copyright   : IMT Atlantique (France)
// License     : See the LICENSE file
// Created     : Jul 09, 2017
// Last update : Aug 01, 2018
//============================================================================

#include "ibex_LoupFinderDefault.h"
#include "ibex_LoupFinderInHC4.h"
#include "ibex_LoupFinderFwdBwd.h"
#include "ibex_BxpLinearRelaxArgMin.h"
#include "ibex_LoupFinderProbing.h"
#include <random>
#include "ibex_LinearizerAbsTaylor.h"
#include "ibex_LinearizerXTaylor.h"


using namespace std;

namespace ibex {

LoupFinderDefault::LoupFinderDefault(const System& sys, bool inHC4,int opcion) :
	finder_probing(inHC4? (LoupFinder&) *new LoupFinderInHC4(sys) : (LoupFinder&) *new LoupFinderFwdBwd(sys)),
	opcion(opcion) {
	LinearizerXTaylor* lr = new LinearizerXTaylor(sys,LinearizerXTaylor::RESTRICT,LinearizerXTaylor::RANDOM);
	LinearizerAbsTaylor* lr2 = new LinearizerAbsTaylor(sys);
	finder_ip_abst = new LoupFinderIP(sys,lr2);	finder_ip_xt = new LoupFinderIP(sys,lr);
	vector<LoupFinderIP*> finders;
	finders.push_back(finder_ip_abst); finders.push_back(finder_ip_xt);
	finder_trust = new LoupFinderIterative(sys,sys.box,finders,0.95,20,1e-3);
}


void LoupFinderDefault::add_property(const IntervalVector& init_box, BoxProperties& prop) {
	finder_probing.add_property(init_box,prop);
	finder_ip_xt->add_property(init_box,prop);

	//--------------------------------------------------------------------------
	/* Using line search from LP relaxation minimizer seems not interesting. */
//	if (!prop[BxpLinearRelaxArgMin::get_id(finder_x_taylor.sys)]) {
//		prop.add(new BxpLinearRelaxArgMin(finder_x_taylor.sys));
//	}
	//--------------------------------------------------------------------------

}

std::pair<IntervalVector, double> LoupFinderDefault::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup, BoxProperties& prop) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;

	try {
		p=finder_probing.find(box,p.first,p.second,prop);
		found=true;
	} catch(NotFound&) { }

	if(opcion==2 || opcion==4){
	try {
			// TODO
			// in_x_taylor.set_inactive_ctr(entailed->norm_entailed);
			p=finder_ip_abst->find(box,p.first,p.second,prop);
			found=true;
		} catch(NotFound&) { }
	}
	if(opcion==3 ){
	try {
		// TODO
		// in_x_taylor.set_inactive_ctr(entailed->norm_entailed);
		p=finder_trust->find(box,box.mid(),p.second);
		found=true;
	} catch(NotFound&) { }
	}
	if(opcion==1 || opcion==4){
	try {
		// TODO
		// in_x_taylor.set_inactive_ctr(entailed->norm_entailed);
		p=finder_ip_xt->find(box,p.first,p.second,prop);
		found=true;
	} catch(NotFound&) { }
	}

	if (found) {
		//--------------------------------------------------------------------------
		/* Using line search from LP relaxation minimizer seems not interesting. */
		//	BxpLinearRelaxArgMin* argmin=(BxpLinearRelaxArgMin*) prop[BxpLinearRelaxArgMin::get_id(finder_x_taylor.sys)];
		BxpLinearRelaxArgMin* argmin=NULL;
		//--------------------------------------------------------------------------

		if (argmin && argmin->argmin()) {
			Vector loup_point = p.first.lb();
			double loup = p.second;
			LoupFinderProbing(finder_ip_xt->sys).dichotomic_line_search(loup_point, loup, *argmin->argmin(), false);
			//cout << "better loup found! " << loup << endl;
			p=make_pair(loup_point,loup);
		}
		return p;
	} else
		throw NotFound();

}

LoupFinderDefault::~LoupFinderDefault() {
	delete &finder_probing;
}

} /* namespace ibex */
