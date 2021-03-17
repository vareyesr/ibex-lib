

#include "ibex_LoupFinderIterative.h"

using namespace std;

namespace ibex {


LoupFinderIterative::LoupFinderIterative(const System& sys,const IntervalVector& initial_box,double alpha,loup_finders lfinders) :
	finder_abs_taylor(sys),finder_x_taylor(sys),initial_box(initial_box),f_goal(sys.goal),alpha(alpha),sys(sys),lfinders(lfinders) {

}

std::pair<IntervalVector, double> LoupFinderIterative::find(const IntervalVector& box, const IntervalVector& old_loup_point, double old_loup) {

	pair<IntervalVector,double> p=make_pair(old_loup_point, old_loup);

	bool found=false;
	pair<IntervalVector,double> new_ub = p;
	pair<IntervalVector,double> old_ub = p;
 	bool flag = true;



 		if ((lfinders == BOTH) || (lfinders == ABST)){
			try{
				pair<IntervalVector,double> new_ub=finder_abs_taylor.find(box,box.mid(),p.second);
				if(new_ub.second < p.second){
					found = true;
					p = new_ub;
				}
				else throw NotFound();
			} catch(NotFound&) { }
 		}
 		if ((lfinders == BOTH) || (lfinders == XT)){
			try{
				pair<IntervalVector,double> new_ub=finder_x_taylor.find(box,box.mid(),p.second);
				if(new_ub.second < p.second){
					found = true;
					p = new_ub;
				}
				else throw NotFound();
			} catch(NotFound&) { }
 		}

		if (!found){
			throw NotFound();
		}
	 	IntervalVector box_aux(box.size());
	 	box_aux = box;
	 	while ((old_ub.second-p.second > 1e-6) || (flag)){
	 		if (old_ub.second-p.second < 1e-6) flag = false;
	 		else flag = true;

	 		Vector old_exp = p.first.mid();
	 		for (int i = 0 ; i < box.size() ; i++){

				if ((old_exp[i]-(alpha)*box_aux[i].diam()/2>=initial_box[i].lb()) && (old_exp[i]+(alpha)*box_aux[i].diam()/2<=initial_box[i].ub()))
					box_aux[i] = Interval(old_exp[i]-(alpha)*box_aux[i].diam()/2,old_exp[i]+(alpha)*box_aux[i].diam()/2);
				else if ((old_exp[i]-(alpha)*box_aux[i].diam()/2>=initial_box[i].lb()) && (old_exp[i]+(alpha)*box_aux[i].diam()/2>initial_box[i].ub()))
					box_aux[i] = Interval(old_exp[i]-(alpha)*box_aux[i].diam()/2,initial_box[i].ub());
				else if ((old_exp[i]-(alpha)*box_aux[i].diam()/2<initial_box[i].lb()) && (old_exp[i]+(alpha)*box_aux[i].diam()/2<=initial_box[i].ub()))
					box_aux[i] = Interval(initial_box[i].lb(),old_exp[i]+(alpha)*box_aux[i].diam()/2);
				else if ((old_exp[i]-alpha*box_aux[i].diam()/2<initial_box[i].lb()) && (old_exp[i]+alpha*box_aux[i].diam()/2>initial_box[i].ub()))
					box_aux[i] = Interval(initial_box[i].lb(),initial_box[i].ub());
	 		}
		 	old_ub = p;

		 	if ((lfinders == BOTH) || (lfinders == ABST)){
				try {
					new_ub=finder_abs_taylor.find(box_aux,p.first.mid(),p.second);
					if(new_ub.second < p.second){
						p = new_ub;
					}
				} catch(NotFound&) {}
		 	}
		 	if ((lfinders == BOTH) || (lfinders == XT)){
				try {
					new_ub=finder_x_taylor.find(box_aux,p.first.mid(),p.second);
					if(new_ub.second < p.second){
						p = new_ub;
					}

				} catch(NotFound&) {}
		 	}
	 	}
	 	if (found){
	 		return p;
	 	}
	 	else
	 		throw NotFound();
	 }

LoupFinderIterative::~LoupFinderIterative() {

}

} /* namespace ibex */
