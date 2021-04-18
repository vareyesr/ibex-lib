#include "ibex_LoupFinderIterative.h"

using namespace std;

namespace ibex {


LoupFinderIterative::LoupFinderIterative(const System& sys,const IntervalVector& initial_box,double alpha,int max_iter,double prec) :
	initial_box(initial_box),alpha(alpha),sys(sys),max_iter(max_iter),prec(prec),trace(false) {
	LinearizerXTaylor* lr = new LinearizerXTaylor(sys,LinearizerXTaylor::RESTRICT,LinearizerXTaylor::RANDOM);
	LinearizerAbsTaylor* lr2 = new LinearizerAbsTaylor(sys);
	LoupFinderIP* finder_ip_abst = new LoupFinderIP(sys,lr2); LoupFinderIP* finder_ip_xt = new LoupFinderIP(sys,lr);
	finders.push_back(*finder_ip_abst); finders.push_back(*finder_ip_xt);
}


void LoupFinderIterative::change_box_size(IntervalVector& box_aux, Vector old_exp){
	for (int i = 0 ; i < box_aux.size() ; i++){
		if ((old_exp[i]-(alpha)*box_aux[i].diam()/2>=initial_box[i].lb()) && (old_exp[i]+(alpha)*box_aux[i].diam()/2<=initial_box[i].ub()))
			box_aux[i] = Interval(old_exp[i]-(alpha)*box_aux[i].diam()/2,old_exp[i]+(alpha)*box_aux[i].diam()/2);
		else if ((old_exp[i]-(alpha)*box_aux[i].diam()/2>=initial_box[i].lb()) && (old_exp[i]+(alpha)*box_aux[i].diam()/2>initial_box[i].ub()))
			box_aux[i] = Interval(old_exp[i]-(alpha)*box_aux[i].diam()/2,initial_box[i].ub());
		else if ((old_exp[i]-(alpha)*box_aux[i].diam()/2<initial_box[i].lb()) && (old_exp[i]+(alpha)*box_aux[i].diam()/2<=initial_box[i].ub()))
			box_aux[i] = Interval(initial_box[i].lb(),old_exp[i]+(alpha)*box_aux[i].diam()/2);
		else if ((old_exp[i]-alpha*box_aux[i].diam()/2<initial_box[i].lb()) && (old_exp[i]+alpha*box_aux[i].diam()/2>initial_box[i].ub()))
			box_aux[i] = Interval(initial_box[i].lb(),initial_box[i].ub());
	}
}

void LoupFinderIterative::print_ub(std::pair<IntervalVector,double> p){
	cout << "The point :    ";
	cout << p.first.ub() << endl;
	cout << "corresponds to an upperbound of the problem with a cost of ";
	cout << p.second << endl << endl;
}

std::pair<IntervalVector, double> LoupFinderIterative::find(const IntervalVector& box, const IntervalVector& exp_point, double old_loup) {

	pair<IntervalVector,double> p=make_pair(exp_point, old_loup);

	bool found=false;
	pair<IntervalVector,double> new_ub = p;
	pair<IntervalVector,double> old_ub = p;
 	bool flag = true;

//while XD
 	for (int i = 0 ; i < finders.size() ; i++){
		try{
			pair<IntervalVector,double> new_ub=finders[i].find(box,exp_point,p.second);
			if(new_ub.second < p.second){
				found = true;
				p = new_ub;
				if (trace) print_ub(p);
			}
			else throw NotFound();
		} catch(NotFound&) { }
 	}

	if (!found){
		throw NotFound();
	}
	IntervalVector box_aux(box.size());
	box_aux = box;
	int nb_iter = 0;
	while ((old_ub.second-p.second > prec) || (flag)){
		if (old_ub.second-p.second < prec) flag = false;
		else flag = true;

		Vector old_exp = p.first.mid();
		change_box_size(box_aux,old_exp);
		old_ub = p;

		try {
			new_ub=finders[0].find(box_aux,p.first.mid(),p.second);
			if(new_ub.second < p.second){
				p = new_ub;
				if (trace) print_ub(p);
			}
		} catch(NotFound&) {}

		nb_iter++;
		if (nb_iter >= max_iter)
			break;
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
