/* ============================================================================
 * I B E X - Inner Polytope loup finder
 * ============================================================================
 *
 * Author(s)   : Ignacio Araya, Victor Reyes
 * Created     : April 2021
 * Updated     : April 2021
 * ---------------------------------------------------------------------------- */

#include "ibex_LoupFinderIP.h"
#include <stdio.h>
#include <string>

using namespace std;

namespace ibex {

//TODO: remove this recipe for the argument of the max number of iterations of the LP solver
LoupFinderIP::LoupFinderIP(const System& sys,loup_finder loup) : sys(sys), lp_solver((loup==XT)? sys.nb_var:(2*sys.nb_var)),loup(loup),exp_point(sys.box.mid()) {
	lp_solver.set_max_iter(std::min(sys.nb_var*3, int(LPSolver::default_max_iter)));
		if (loup == XT) lr = new LinearizerXTaylor(sys,LinearizerXTaylor::RESTRICT,LinearizerXTaylor::RANDOM);
		else if (loup == ABST) lr = new LinearizerAbsTaylor(sys);
}

void LoupFinderIP::add_property(const IntervalVector& init_box, BoxProperties& prop) {
	lr->add_property(init_box,prop);
}

std::pair<IntervalVector, double> LoupFinderIP::find(const IntervalVector& box, const IntervalVector&, double current_loup, BoxProperties& prop) {


	int n=sys.nb_var;

	if (box.is_unbounded())
		throw NotFound();

	lp_solver.clear_constraints();

	if (loup == XT)
		lp_solver.set_bounds(box);
	else if (loup == ABST){
		IntervalVector box2(n*2);
		for(int i=0;i<n;i++)
			box2[i]=box[i];

		//initialize auxiliary variables u_i
		for(int i=0;i<n;i++)
			box2[n+i]=Interval(-box2[i].mag()-1, box2[i].mag()+1);
//	box2[n+i]=Interval(-box2[i].mag(), box2[i].mag());
		lp_solver.set_bounds(box2);
		LinearizerAbsTaylor* lr_abst = dynamic_cast<LinearizerAbsTaylor*>(lr);
//		lr_abst->set_expansion_point(exp_point.mid());
	}

	IntervalVector ig=sys.goal->gradient(box.mid());
	if (ig.is_empty()) // unfortunately, at the midpoint the function is not differentiable
		throw NotFound(); // not a big deal: wait for another box...

	Vector g=ig.mid();

	// set the objective coefficient
	// TODO: replace with lp_solver.set_cost(g) when implemented
	for (int j=0; j<n; j++)
		lp_solver.set_cost(j,g[j]);

	int count = lr->linearize(box,lp_solver,prop);

	if (count==-1) {
		lp_solver.clear_constraints();
		throw NotFound();
	}
	LPSolver::Status stat = lp_solver.minimize();

	if (stat == LPSolver::Status::Optimal) {
		//the linear solution is mapped to intervals and evaluated
		Vector loup_point = lp_solver.not_proved_primal_sol();
		//in case abstaylor is used
		loup_point.resize(box.size());
		// we allow finding a loup outside of the current box, but
		// not outside of the system box.
		if (!sys.box.contains(loup_point)) throw NotFound();

		/*probar si es necesario*/
		for(int i=0;i<box.size();i++){
			if(box[i].lb() > loup_point[i]) loup_point[i]=box[i].lb();
			else if(box[i].ub() < loup_point[i]) loup_point[i]=box[i].ub();
		}

		double new_loup=current_loup;

		if (check(sys,loup_point,new_loup,false)) {
			return std::make_pair(loup_point,new_loup);
		}
	}

	throw NotFound();
}

} /* namespace ibex */
