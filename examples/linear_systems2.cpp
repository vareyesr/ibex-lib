/*
 * linear_systems.cpp
 *
 *  Created on: Jul 19, 2020
 *	Author: Victor Reyes
 */

#include "ibex.h"
#include <iostream>
#include <random>
#include <fstream>
#include <cstdlib>
#include <ctime>

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */


using namespace std;
using namespace ibex;


//just for 2n-simplex (polytopehull)
class LinearSystem : public Linearizer, public Ctc {
public:
	LinearSystem(IntervalMatrix A, IntervalVector b, int nb_var) :
		Linearizer(nb_var),
		ctc(*this, LPSolver::default_max_iter,
		LPSolver::default_timeout, LPSolver::default_tolerance),
		Ctc(nb_var), A(A), b(b){

	}

	void contract(IntervalVector& box){
		ctc.contract(box);
		if(box.is_empty()) return;
	}

	int linearization(const IntervalVector& x, LPSolver& lp_solver){
		int num=0;
		for (int i=0; i<A.nb_rows(); i++) {


			Interval error = (IntervalVector(A[i].diam())* IntervalVector(x.mag())).ub();
			Vector row=A[i].mid();
			try {
				Interval eval=(row*x); // left side image
				if(eval.lb() < (b[i] - error).lb()){
				   lp_solver.add_constraint(row,GEQ, (b[i] - error).lb());
				   num++;
				}
				if(eval.ub() > (b[i] + error).ub()){
					lp_solver.add_constraint(row,LEQ, (b[i] + error).ub());
					num++;

				}
			} catch (LPException&) { }
		}
		return num;
	}

	virtual int linearize(const IntervalVector& x, LPSolver& lp_solver){
		return LinearSystem::linearization(x,lp_solver);
	}

	CtcPolytopeHull ctc;
	IntervalMatrix A;
	IntervalVector b;
};

double randnum (double a, double b)
{
  unsigned seed = 1;
  static std::default_random_engine generator (seed);
  //static std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution (a,b);
  return distribution(generator);
}

int main(int argc, char** argv) {

	if (argc<3)
		cout << "usage: ./linear_systems bench extended results"  << endl;

	/*reading from file*/
	std::ifstream file;
	file.open(argv[1]);
	int extended = atoi(argv[2]);
	int nb_rows; file >> nb_rows;
	int nb_cols; file >> nb_cols;

	int results = atoi(argv[3]);
	/*system*/
	IntervalMatrix A(nb_rows,nb_cols);
	IntervalVector b(nb_rows);
	IntervalVector box(nb_cols);
	double value;

	for (int i = 0 ; i < nb_rows ; i++){
		for (int j = 0 ;  j < nb_cols ; j++){
			file >> value;
			A[i][j] = value;
		}
		file >> value;
		b[i] = value;
	}


	/*initial box*/
	double lb_value,ub_value;
	for (int i = 0 ; i < box.size() ; i++){
		file >> lb_value; file >> ub_value;
		box[i] = Interval(lb_value,ub_value);
	}

	IntervalVector sub_box = box;
	double prop = 0.5;
	double alpha = 1-prop;
	for (int i = 0 ; i < box.size() ; i++){
		double diam_lb = -box[i].lb();
		double diam_ub = box[i].ub();
		if (diam_ub > diam_lb){
			double beta = alpha*box[i].diam();
			do{
				double new_bound = randnum (diam_ub-beta,diam_ub);
				double gamma = new_bound-box[i].diam()*prop;
				sub_box[i] = Interval(gamma,new_bound);
			}while (!sub_box[i].interior_contains(0));
		}
		else{
			do{
				double beta = alpha*box[i].diam();
				double new_bound = randnum(diam_lb-beta,diam_lb);
				double gamma = -new_bound+box[i].diam()*prop;
				sub_box[i] = Interval(-new_bound,gamma);
			}while (!sub_box[i].interior_contains(0));
		}
	}
	file.close();

	/*check if extended. if true, then transform A and b*/
	if (extended){
		IntervalMatrix AA=A;
		Matrix I = Matrix::diag(-Vector::ones(b.size()));
		A.resize(AA.nb_rows(), AA.nb_cols()+b.size());
		A.put(0,AA.nb_cols(),I);
		box.resize(AA.nb_cols()+b.size());
		int _size = b.size();
		box.put(AA.nb_cols(),b);
		b=Vector::zeros(_size);
	}




	//Contraction of the box using HC4
	IntervalVector box_hc4 = box;
	IntervalMatrix tmp(A);
	box_hc4 = sub_box;//NUEVO ! SUBCAJA
	bwd_mul(b, tmp, box_hc4, 1e-8);

	//Perimeter
	double perimeter_hc4 = 0;
	for (int i = 0 ; i < box_hc4.size() ; i++)
		perimeter_hc4+=box_hc4[i].diam();
	if (results == 1){
		cout << endl;
		cout << "HC4 Contraction"<<endl<<endl;
		cout << CYAN<<"Min diameter:  "<<RESET << box_hc4.min_diam() << endl;
		cout << YELLOW<<"Perimeter:  " << RESET<<perimeter_hc4 << endl <<endl<<endl;
	}
	//Contraction of the box using 2n-Simplex

	IntervalVector box_2nsimplex = sub_box;//NUEVO ! SUBCAJA
	LinearSystem hola(A,b,A.nb_cols());
	hola.contract(box_2nsimplex);
	double perimeter_2nsimplex = 0;
	for (int i = 0 ; i < box_2nsimplex.size() ; i++)
		perimeter_2nsimplex+=box_2nsimplex[i].diam();
	if (results == 1){
		cout << "2n Simplex Contraction"<<endl<<endl;
		cout << CYAN<<"Min diameter:  " << RESET<<box_2nsimplex.min_diam() << endl;
		cout << YELLOW<<"Perimeter:  " <<RESET<< perimeter_2nsimplex << endl <<endl<<endl;
	}

	//Contraction of the box using preconditioners
	//Gauss
	IntervalVector box_gauss = box;
	IntervalMatrix AA=A;
	vector<IntervalMatrix> perm_list;
	vector <vector <pair <int,int> > > proj_vars;
	best_gauss_jordan_new (AA,box,perm_list,proj_vars,1e-7,2);
	IntervalMatrix PA_gauss(perm_list[0]*A);
	IntervalVector Pb_gauss(perm_list[0]*b);
	box_gauss = sub_box;    //NUEVO ! SUBCAJA
	bwd_mul(Pb_gauss, PA_gauss, box_gauss, 1e-8);
	//Perimeter
	double perimeter_gauss = 0;
	for (int i = 0 ; i < box_gauss.size() ; i++)
		perimeter_gauss+=box_gauss[i].diam();
	if (results == 1){
		cout << "Gauss Contraction"<<endl<<endl;
		cout << CYAN<<"Min diameter:  " << RESET<<box_gauss.min_diam() << endl;
		cout << YELLOW<<"Perimeter:  " <<RESET<< perimeter_gauss << endl <<endl<<endl;
	}
	//Contraction using LP-1-simplex
	AA=A;
	IntervalVector LP_1_simplex = box;
	IntervalMatrix P_LP1 = best_P (AA, LP_1_simplex, b.size(), 1e-8, extended);
	IntervalMatrix PA_LP1(P_LP1*A);
	IntervalVector Pb_LP1(P_LP1*b);
	LP_1_simplex =  sub_box;    //NUEVO ! SUBCAJA
	bwd_mul(Pb_LP1, PA_LP1, LP_1_simplex, 1e-8);
	//Perimeter
	double perimeter_LP1 = 0;
	for (int i = 0 ; i < LP_1_simplex.size() ; i++)
		perimeter_LP1+=LP_1_simplex[i].diam();
	if (results == 1){
		cout << "LP 1-Simplex Contraction"<<endl<<endl;
		cout << CYAN<<"Min diameter:  " << RESET<<LP_1_simplex.min_diam() << endl;
		cout << YELLOW<<"Perimeter:  " <<RESET<< perimeter_LP1 << endl <<endl<<endl;
	}
	//Contraction using LP-2-simplex
	AA=A;
	IntervalVector LP_2_simplex = box;

	vector<IntervalMatrix> P_LP2;
	P_LP2.push_back(best_P_bound (AA, LP_2_simplex, b.size(), 1e-8, extended,0));
	P_LP2.push_back(best_P_bound (AA, LP_2_simplex, b.size(), 1e-8, extended,1));
	vector<IntervalMatrix> PA_LP2; PA_LP2.push_back(P_LP2[0]*A); PA_LP2.push_back(P_LP2[1]*A);
	vector<IntervalVector> Pb_LP2; Pb_LP2.push_back(P_LP2[0]*b); Pb_LP2.push_back(P_LP2[1]*b);
	LP_2_simplex = sub_box;    //NUEVO ! SUBCAJA

	for (int i = 0 ; i < P_LP2.size() ; i++)
		bwd_mul(Pb_LP2[i], PA_LP2[i], LP_2_simplex, 1e-8);
	//Perimeter
	double perimeter_LP2 = 0;
	for (int i = 0 ; i < LP_2_simplex.size() ; i++)
		perimeter_LP2+=LP_2_simplex[i].diam();
	if (results == 1){
		cout << "LP 2-Simplex Contraction"<<endl<<endl;
		cout << CYAN<<"Min diameter:  " << RESET<<LP_2_simplex.min_diam() << endl;
		cout << YELLOW<<"Perimeter:  " <<RESET<< perimeter_LP2 << endl <<endl<<endl;
	}


	if (results == 0){
		cout <<argv[1] <<" "<<box_hc4.min_diam() << " " << perimeter_hc4 << " " << box_2nsimplex.min_diam() << " "
				<< perimeter_2nsimplex << " " << box_gauss.min_diam() << " " << perimeter_gauss <<
				" " << LP_1_simplex.min_diam() << " " << perimeter_LP1 << " " <<
				LP_2_simplex.min_diam() << " " << perimeter_LP2 << endl;
	}
}


