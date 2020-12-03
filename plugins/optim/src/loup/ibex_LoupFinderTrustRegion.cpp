

#include "ibex_LoupFinderTrustRegion.h"

using namespace std;

namespace ibex {


LoupFinderTrustRegion::LoupFinderTrustRegion(const System& sys,const IntervalVector& initial_box,double alpha) :
	finder_abs_taylor(sys),finder_x_taylor(sys),initial_box(initial_box),f_goal(sys.goal),alpha(alpha) {

}

std::pair<IntervalVector, double> LoupFinderTrustRegion::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;
	pair<IntervalVector,double> new_ub = p;
	pair<IntervalVector,double> old_ub = p;
 	bool flag = true;

 	double pp;

		try{
			pair<IntervalVector,double> new_ub=finder_abs_taylor.find(box,box.mid(),p.second);
			if(new_ub.second < p.second){
				found = true;
				p = new_ub;
			}
			else throw NotFound();
		} catch(NotFound&) { }
		try{
			pair<IntervalVector,double> new_ub=finder_x_taylor.find(box,box.mid(),p.second);
			if(new_ub.second < p.second){
				found = true;
				p = new_ub;
			}
			else throw NotFound();
		} catch(NotFound&) { }

		if (!found){
			throw NotFound();
		}
	 	IntervalVector box_aux(box.size());
	 	box_aux = box;
	 	while ((old_ub.second-p.second>1e-6) || (flag)){
	 		if (old_ub.second-p.second < 1e-6) flag = false;
	 		else flag = true;

	 		Vector old_exp = p.first.mid();
	 		for (int i = 0 ; i < box.size() ; i++){
	 			if ((std::abs(box[i].lb()-p.first[i].lb())<1e-6) || (std::abs(box[i].ub()-p.first[i].ub())<1e-6)){
	 				pp = 1;
	 			}
	 			else{
	 				pp = 1;
	 			}
				if ((old_exp[i]-(alpha/pp)*box_aux[i].diam()/2>=initial_box[i].lb()) && (old_exp[i]+(alpha/pp)*box_aux[i].diam()/2<=initial_box[i].ub()))
					box_aux[i] = Interval(old_exp[i]-(alpha/pp)*box_aux[i].diam()/2,old_exp[i]+(alpha/pp)*box_aux[i].diam()/2);
				else if ((old_exp[i]-(alpha/pp)*box_aux[i].diam()/2>=initial_box[i].lb()) && (old_exp[i]+(alpha/pp)*box_aux[i].diam()/2>initial_box[i].ub()))
					box_aux[i] = Interval(old_exp[i]-(alpha/pp)*box_aux[i].diam()/2,initial_box[i].ub());
				else if ((old_exp[i]-(alpha/pp)*box_aux[i].diam()/2<initial_box[i].lb()) && (old_exp[i]+(alpha/pp)*box_aux[i].diam()/2<=initial_box[i].ub()))
					box_aux[i] = Interval(initial_box[i].lb(),old_exp[i]+(alpha/pp)*box_aux[i].diam()/2);
				else if ((old_exp[i]-alpha*box_aux[i].diam()/2<initial_box[i].lb()) && (old_exp[i]+alpha*box_aux[i].diam()/2>initial_box[i].ub()))
					box_aux[i] = Interval(initial_box[i].lb(),initial_box[i].ub());
	 		}
		 	old_ub = p;
	 		try {
		 		new_ub=finder_abs_taylor.find(box_aux,p.first.mid(),p.second);
	 			if(new_ub.second < p.second){
	 		 		p = new_ub;
	 		 	}
	 		} catch(NotFound&) {}
	 		try {
				new_ub=finder_x_taylor.find(box_aux,p.first.mid(),p.second);
	 			if(new_ub.second < p.second){
	 		 		p = new_ub;
	 		 	}

			} catch(NotFound&) {}

	 	}
	 	if (found){
	 		return p;
	 	}
	 	else
	 		throw NotFound();
	 }

LoupFinderTrustRegion::~LoupFinderTrustRegion() {
	delete &finder_abs_taylor;
	delete &finder_x_taylor;
}

} /* namespace ibex */
