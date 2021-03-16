/* ============================================================================
 * I B E X - AbsTaylor linear restriction
 * ============================================================================
 *
 * Author(s)   : Ignacio Araya, Victor Reyes
 * Created     : April 2018
 * Updated     : December 2020
 * ---------------------------------------------------------------------------- */


#include "ibex_LoupFinderAbsTaylor.h"
#include <stdio.h>
#include <string>

using namespace std;

namespace ibex {

//TODO: remove this recipe for the argument of the max number of iterations of the LP solver
LoupFinderAbsTaylor::LoupFinderAbsTaylor(const System& sys) :
		sys(sys), lp_solver(2*sys.nb_var) {
		lr = new LinearizerAbsTaylor(sys);
}



std::pair<IntervalVector, double> LoupFinderAbsTaylor::find(const IntervalVector& box, const IntervalVector& exp_point, double current_loup) {

	int n=sys.nb_var;
	lp_solver.clear_constraints();


	IntervalVector box2(n*2);
	for(int i=0;i<n;i++)
		box2[i]=box[i];
		//initialize auxiliary variables u_i
	for(int i=0;i<n;i++)
		box2[n+i]=Interval(-box2[i].mag()-1, box2[i].mag()+1);
	lp_solver.set_bounds(box2);
	LinearizerAbsTaylor* lr_abst = dynamic_cast<LinearizerAbsTaylor*>(lr);
	lr_abst->set_expansion_point(exp_point.mid());

	IntervalVector ig=sys.goal->gradient(exp_point.mid());

	if (ig.is_empty()) // unfortunately, at the midpoint the function is not differentiable
		throw NotFound(); // not a big deal: wait for another box...

	Vector g=ig.mid();

	// set the objective coefficient
	for (int j=0; j<n; j++)
		lp_solver.set_cost(j,g[j]);

	int count = lr->linearize(box,lp_solver);

	if (count==-1) {
		lp_solver.clear_constraints();
		throw NotFound();
	}

	LPSolver::Status stat = lp_solver.minimize();

	if (stat == LPSolver::Status::Optimal) {
		//the linear solution is mapped to intervals and evaluated
		Vector loup_point = lp_solver.not_proved_primal_sol();


		loup_point.resize(box.size());

		//correction
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
