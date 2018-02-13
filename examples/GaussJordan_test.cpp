#include "ibex.h"

using namespace std;
using namespace ibex;

int main(int argc, char** argv){
/*mini example*/
//		IntervalMatrix A(1,7);
//		IntervalVector x(7);
//		IntervalVector b(1);
//		A[0][0] = Interval(4,4);
//		A[0][1] = Interval(-0.1,-0.1);
//		A[0][2] = Interval(0.2,0.2);
//		A[0][3] = Interval(4,4);
//		A[0][4] = Interval(-0.1,-0.1);
//		A[0][5] = Interval(0.2,0.2);
//		A[0][6] = Interval(0.2,0.2);
//		x[0] = Interval(-100,100);
//		x[1] = Interval(-40,30);
//		x[2] = Interval(-60,60);
//		x[3] = Interval(-90,90);
//		x[4] = Interval(-8,6);
//		x[5] = Interval(1,3);
//		x[6] = Interval(1,3);
//		b[0] = Interval(-10,30);
		vector<Matrix> perm_list;
		vector <vector <pair <int,int> > > proj_vars;
//		cout << A.mid() << endl << endl << x << endl << endl;
//		int i = 0;
//		best_gauss_jordan (A, x, perm_list,proj_vars, 1e-8);
//		cout << perm_list.size() << endl;
//		for (int i =  0; i < perm_list.size() ; i++)
//			cout <<perm_list[i] << endl << endl;

   System *sys = new System(argv[1]);
   GaussContractor* gaussctc = new GaussContractor(*sys,sys->box);
   cout << sys->box << endl;
   cout << "hola" << endl;
   gaussctc->contract(sys->box);
   cout << sys->box << endl;
   //   best_gauss_jordan (gaussctc->A, sys->box, perm_list,proj_vars, 1e-8);
//   cout << perm_list.size() << endl;
//   for (int i =  0; i < perm_list.size() ; i++)
//	   cout << perm_list[i]*gaussctc->A << endl << endl;
}
