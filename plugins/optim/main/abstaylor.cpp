/**
 *
 * This main shows how the methods AbsTaylor or TrustRegion are performed.
 *
 * As output you can expect the inner linealizations, the feasible point and the evaluation.
 */

#include "ibex.h"

using namespace std;
using namespace ibex;
#include <random>

float get_random(double low,double high)
{
	std::random_device rd;
    std::default_random_engine e(rd());
    std::uniform_real_distribution<> dis(low, high); // rage 0 - 1

    return dis(e);
}

int main(int argc, char** argv) {

	if (argc < 4){
		cout << "Use the format: ./abstaylor FILE [abst|iterative] [mid|rnd]" << endl;
		return 0;
	}
	System *sys;
	pair<IntervalVector,double> p;
	sys = new System(argv[1]);
	string loup_finder = argv[2];
	string exp_point= argv[3];
	LoupFinderAbsTaylor abst(*sys);
	LoupFinderIterative trust(*sys,sys->box);

	IntervalVector point(sys->box.size());

	if (exp_point == "mid")
		point= sys->box.mid();
	else if (exp_point == "rnd"){
		for (int i = 0 ; i < sys->box.size() ; i++){
			point[i] = get_random(sys->box[i].lb(),sys->box[i].ub());
		}
	}
	else{
		cout << "Use either mid for the midpoint or rnd for a random point inside the box" << endl;
		exit(1);
	}
	cout << endl<<"The search of upperbounds will be performed in the box : "<< sys->box << endl<<endl;
	cout << "The expansion point is: " << point.mid() << endl<<endl;
	cout << "Press a key to continue"<< endl;
	string aux; cin >> aux;

	if (loup_finder == "abst"){
		try {
			cout << "The linealization is the following:" << endl;
			p = abst.find(sys->box,point,POS_INFINITY);
			abst.lp_solver.write_to_file("system.txt");
			std::ifstream f("system.txt");
				if (f.is_open())
					std::cout << f.rdbuf();
			f.close();
			cout << "The point :    ";
			cout << p.first.ub() << endl;
			cout << "corresponds to an upperbound of the problem with a cost of ";
			cout << p.second << endl << endl;

		} catch(int e) {cout << "Upperbound not found"<<endl; }
	}
	else if (loup_finder == "iterative"){
		try {
			trust.set_trace(true);
			trust.find(sys->box,point,POS_INFINITY);
		} catch(int e) {cout << "Upperbound not found"<<endl; }
	}
	else{
		cout << "Method not found!" << endl;
		exit(1);
	}

}
