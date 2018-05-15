/*
 * optimsolver.cpp
 *
 *  Created on: 21 dic. 2017
 *      Author: iaraya
 */

#include "ibex.h"
#include "args.hxx"

using namespace std;
using namespace ibex;

double get_prec(double loup, double uplo);

int main(int argc, char** argv){

	try {

	args::ArgumentParser parser("********* IbexOpt (defaultoptimizer) *********.", "Solve a Minibex file.");
	args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
	args::ValueFlag<std::string> _filtering(parser, "string", "the filtering method", {'f', "filt"});
	args::ValueFlag<std::string> _linear_relax(parser, "string", "the linear relaxation method", {"linear_relax"});
	args::ValueFlag<std::string> _bisector(parser, "string", "the bisection method", {'b', "bis"});
	args::ValueFlag<std::string> _strategy(parser, "string", "the search strategy", {'s', "search"});
	args::ValueFlag<double> _eps_x(parser, "float", "eps_x (the precision of the boxes)", {"exp_x"});
	args::ValueFlag<double> _eps(parser, "float", "eps (the precision of the objective)", {"eps"});
	args::ValueFlag<double> _timelimit(parser, "float", "timelimit", {'t',"timelimit"});
	args::ValueFlag<int> _seed(parser, "int", "seed", {"seed"});
	args::ValueFlag<std::string> _loup_mode(parser, "loup_mode", "UpperBounding mode (xt, abst or both).", {"lmode"});
	args::Flag _trace(parser, "trace", "Activate trace. Updates of loup/uplo are printed while minimizing.", {"trace"});

	args::Positional<std::string> filename(parser, "filename", "The name of the MINIBEX file.");

	try
	{
		parser.ParseCLI(argc, argv);
	}
	catch (args::Help&)
	{
		std::cout << parser;
		return 0;
	}
	catch (args::ParseError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}
	catch (args::ValidationError& e)
	{
		std::cerr << e.what() << std::endl;
		std::cerr << parser;
		return 1;
	}

	string filtering = (_filtering)? _filtering.Get() : "acidhc4";
	cout << "Filtering: " << filtering << endl;
	string linearrelaxation= (_linear_relax)? _linear_relax.Get() : "compo";
	cout << "Relax method: " << linearrelaxation << endl;
	string bisection= (_bisector)? _bisector.Get() : "lsmear";
	cout << "Bisector: " << bisection << endl;
	string strategy= (_strategy)? _strategy.Get() : "feasdiv";
	cout << "Search strategy: " << strategy << endl;
	double prec= (_eps_x)? _eps_x.Get() : 1e-7;
	cout << "prec_x: " << prec << endl;
	double goalprec= (_eps)? _eps.Get() : 1e-6 ;
	cout << "prec_goal: " << goalprec << endl;
	double timelimit = (_timelimit)? _timelimit.Get() : 3600 ;
	cout << "timelimit: " << timelimit << endl;
	int nseed= (_seed)? _seed.Get() : 1 ;
	cout << "seed: " << nseed << endl;



 	double eqeps= 1.e-8;
 	double default_relax_ratio = 0.2;
	srand(nseed);
	RNG::srand(nseed);

    System* orig_sys,*sys;
	LoupFinderDefault* loupfinder;
    std::size_t found = string(filename.Get().c_str()).find(".nl");
	if (found!=std::string::npos){
	       AmplInterface interface (argv[1]);
	       orig_sys= new System(interface);
     }else
           orig_sys = new System(argv[1]);

	// the extended system
	if(strategy=="solver") sys=orig_sys;
	else {
		sys=new ExtendedSystem(*orig_sys,eqeps);
		NormalizedSystem* norm_sys = new NormalizedSystem(*orig_sys,eqeps); //orig_sys
		loupfinder = new  LoupFinderDefault(*norm_sys,true);
	}



	// Build the bisection heuristic
	// --------------------------
	Bsc * bs;

	if (bisection=="roundrobin")
	  bs = new RoundRobin (prec);
	else if (bisection== "largestfirst")
          bs= new LargestFirst();
	else if (bisection=="smearsum")
	  bs = new SmearSum(*sys,prec);
	else if (bisection=="smearmax")
	  bs = new SmearMax(*sys,prec);
	else if (bisection=="smearsumrel")
	  bs = new SmearSumRelative(*sys,prec);
	else if (bisection=="smearmaxrel")
	  bs = new SmearMaxRelative(*sys,prec);
	else if (bisection=="lsmear")
	  bs = new LSmear(*sys,prec);

	else {cout << bisection << " is not an implemented  bisection mode "  << endl; return -1;}

	// Build the contractors

	// the first contractor called
	CtcHC4 hc4(sys->ctrs,0.01,true);
	// hc4 inside acid and 3bcid : incremental propagation beginning with the shaved variable
	CtcHC4 hc44cid(sys->ctrs,0.1,true);
	// hc4 inside xnewton loop
	CtcHC4 hc44xn (sys->ctrs,0.01,false);
	// The 3BCID contractor on all variables (component of the contractor when filtering == "3bcidhc4")
	Ctc3BCid c3bcidhc4(hc44cid);
	// hc4 followed by 3bcidhc4 : the actual contractor used when filtering == "3bcidhc4"
	CtcCompo hc43bcidhc4 (hc4, c3bcidhc4);
	// The ACID contractor (component of the contractor  when filtering == "acidhc4")
	CtcAcid acidhc4(*sys,hc44cid,1);
	// hc4 followed by acidhc4 : the actual contractor used when filtering == "acidhc4"
	CtcCompo hc4acidhc4 (hc4, acidhc4);

	string filtering2="acidhc4";
	Ctc* gauss = new GaussContractor(*orig_sys,dynamic_cast<ExtendedSystem*>(sys)->goal_var());
	Ctc* ctc;
	if (filtering == "hc4")
	  ctc= &hc4;
	else if (filtering =="acidhc4")
	  ctc= &hc4acidhc4;
	else if (filtering =="3bcidhc4")
	  ctc= &hc43bcidhc4;
	else{cout << filtering <<  " is not an implemented  contraction  mode "  << endl; return -1;}
	ctc= new CtcCompo (*gauss,*ctc);
    CtcNewton* ctcnewton=NULL;
    if(strategy=="solver"){
	    ctcnewton=new CtcNewton(sys->f_ctrs,5e8,prec,1.e-4);
		ctc =  new CtcCompo (*ctc, *ctcnewton);
	}


	Linearizer* lr;
	if (linearrelaxation=="art")
		  lr= new LinearizerCombo(*sys,LinearizerCombo::ART);
	else if  (linearrelaxation=="compo")
		  lr= new LinearizerCombo(*sys,LinearizerCombo::COMPO);
	else if (linearrelaxation=="xn")
	  lr= new LinearizerXTaylor (*sys, LinearizerXTaylor::RELAX, LinearizerXTaylor::RANDOM_OPP);

	//	else {cout << linearrelaxation  <<  " is not an implemented  linear relaxation mode "  << endl; return -1;}
	// fixpoint linear relaxation , hc4  with default fix point ratio 0.2
	CtcFixPoint* cxn;
	CtcPolytopeHull* cxn_poly;
	CtcCompo* cxn_compo;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn"|| linearrelaxation=="abs")
          {
		cxn_poly = new CtcPolytopeHull(*lr);
		cxn_compo =new CtcCompo(*cxn_poly, hc44xn);
		cxn = new CtcFixPoint (*cxn_compo, default_relax_ratio);
	  }
	//  the actual contractor  ctc + linear relaxation
	Ctc* ctcxn;
	if (linearrelaxation=="compo" || linearrelaxation=="art"|| linearrelaxation=="xn"|| linearrelaxation=="abs")
          ctcxn= new CtcCompo  (*ctc, *cxn);
	else
	  ctcxn = ctc;

	// one point probed when looking for a new feasible point (updating the loup)
	int samplesize=1;

	if(strategy=="solver"){
		// A "CellStack" means a depth-first search.
		CellStack buff;


		Solver s(*sys, *ctcxn, *bs, buff, Vector(sys->nb_var,prec), Vector(sys->nb_var,POS_INFINITY));
		s.time_limit = timelimit;;
		s.trace=_trace;  // solutions are printed as soon as found when trace=1

		// Solve the system and get the solutions
		Solver::Status state=s.solve(sys->box);


		cout << state << endl;
		s.report();


		// Display the number of boxes (called "cells")
		// generated during the search
		cout << "number of cells=" << s.get_nb_cells() << endl;
		// Display the cpu time used
		cout << "cpu time used=" << s.get_time() << "s."<< endl;

		cout << argv[1] << " " << s.get_manifold().size() << " " << s.get_time() << " " <<
		s.get_nb_cells() << " " << (s.get_time()>timelimit) << endl;

		return 1;

	}


	// the optimizer : the same precision goalprec is used as relative and absolute precision
	Optimizer* o=NULL;

	CellBufferOptim* buffer = NULL;

	if(strategy=="lbub")
		buffer = new CellDoubleHeap(*dynamic_cast<ExtendedSystem*>(sys));
	else if(strategy=="feasdiv")
		buffer = new CellBeamSearch(*new CellHeap (*dynamic_cast<ExtendedSystem*>(sys)),
								       *new CellHeap (*dynamic_cast<ExtendedSystem*>(sys)), *dynamic_cast<ExtendedSystem*>(sys));
	else {cout << strategy <<  " is not an implemented search strategy "  << endl; return -1;}

	o=new Optimizer(sys->nb_var, *ctcxn,*bs, *loupfinder, *buffer, dynamic_cast<ExtendedSystem*>(sys)->goal_var(),
	    		prec,goalprec,goalprec);



	// the trace
	o->trace=_trace;

	if (o && o->trace)	cout << " sys.box " << sys->box << endl;

	// the allowed time for search
	o->timeout=timelimit;

    vector<IntervalVector> sols;
    cout.precision(10);
	// the search itself
	o->optimize(orig_sys->box);

	// printing the results
	if (o->trace)
	  o->report();


	delete bs;
	if (linearrelaxation=="xn" ||linearrelaxation=="compo" || linearrelaxation=="art" )
	  {
		//delete lr;
	    delete ctcxn;
	  //  delete cxn;
	  }




    cout << argv[1] << " " << o->get_uplo() << "," << o->get_loup() << " " << double(o->get_time()) << " " <<
             double(o->get_nb_cells()) << " " << (o->get_time()>timelimit)<< endl;
	//loupfinderd->print_results();


	return 0;

	}


	catch(ibex::SyntaxError& e) {
	  cout << e << endl;
	}
}
