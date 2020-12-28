/**
 *
 * This main shows how the methods AbsTaylor or TrustRegion are performed.
 *
 * As output you can expect the inner linealizations, the feasible point and the evaluation.
 */

#include "ibex.h"

using namespace std;
using namespace ibex;


int main(int argc, char** argv) {

	System *sys;
	sys = new System(argv[1]);

	LoupFinderAbsTaylor abst(*sys,true);
//	LoupFinderTrustRegion trust(*sys);

	IntervalVector point(sys->box.size());
	int opcion;

	cout << "Sistema original:" << endl;
	cout << *sys << endl << endl;
	cout << "Indique la linealizacion que desea: " << endl;
	cout << "0 para AbsTaylor, 1 para TrustRegion " << endl;
	cin >> opcion;

	double cost = 0;
	if (opcion == 0){
		try {
			abst.find(sys->box,point,cost);
			std::ifstream f("sistema.txt");
			if (f.is_open())
				std::cout << f.rdbuf();

		} catch(int e) {cout << "Upperbound not found"<<endl; }
	}
	else{
		try {
			abst.find(sys->box,point,cost);
			std::ifstream f("sistema.txt");
			if (f.is_open())
				std::cout << f.rdbuf();

		} catch(int e) {cout << "Upperbound not found"<<endl; }
	}

}
