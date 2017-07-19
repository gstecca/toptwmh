/*
 * MH.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: gstecca
 */

#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <map>
#include "graphg.h"
#include "toptw.h"
#include <ilcplex/ilocplex.h>
#include "MH.h"

using namespace std;



MH::MH(IloEnv _env, IloModel _model, toptwdata _dat, params _par, map<string, IloConstraint> _mapConstr, map<string, IloNumVar> _vars, float  **_t) {
	env = _env;
	model = _model;
	dat = _dat;
	mapConstr = _mapConstr;
	mapVar = _vars;
	t = _t;
	par = _par;
	//mh_mapConstr = map<string, IloConstraint> ();
	//mh_mapVar =	map<string, IloNumVar>();

	// add variables and set free constraints
	for (int i = 1; i < dat.n+1; i++){
		string varName = "u#" + to_string(i);
		mh_mapVar[varName] = IloBoolVar(env, 1, 1, varName.c_str());
		IloExpr expr(env);
		for(int v = 0; v < par.m; v++){
            expr += mapVar["y#"+to_string(i)+"#"+to_string(v)];
		}
        expr -= mh_mapVar[varName];
		string cname = "mh_y#" + to_string(i);
        IloConstraint ic(expr==0);
        ic.setName(cname.c_str());
        mh_mapConstr[cname] = ic;
        //model.add(ic);
        expr.end();
	}

}

int MH::init_sol(){
	int n = dat.n;
	int m = par.m;

    IloCplex cplex(model);
    cplex.exportModel("model.lp");

    cplex.setParam(IloCplex::Param::TimeLimit, 5);
    cplex.solve();
    cplex.out() << "Solution status: " << cplex.getStatus() << endl;
    cplex.out() << " Solution value: " << cplex.getObjValue () << endl;



    IloNum tolerance = cplex.getParam(
       IloCplex::Param::MIP::Tolerances::Integrality);
    cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
    cplex.out() << "node, vehicle  - y[i][v]" << endl;
    for(int i = 0; i < n+2; i++) {
        for(int v = 0; v <m; v++){
        	if ( mapVar.find( "y#" + to_string(i) + "#" + to_string(v) ) == mapVar.end())
        			continue;
          if (cplex.getValue(mapVar["y#" + to_string(i) + "#" + to_string(v)]) >= 1 - tolerance) {
              cplex.out() <<  i << "," << v << endl;
          }
      }
    }

    cout << "x[i][j][k]" << endl;
    for(int i = 0; i < n+2; i++) {
        for(int j = 0; j <n+2; j++){
            for(int k = 0; k < m; k++){
              if (cplex.isExtracted(mapVar["x#" + to_string(i) + "#" + to_string(j)+ "#" + to_string(k)])
              		&& cplex.getValue(mapVar["x#" + to_string(i) + "#" + to_string(j)+ "#" + to_string(k)]) >= 1 - tolerance) {
              cout <<  i << "," << j << "," << k << "\n";
            }

          }
      }
    }

   cout << "a[i][k]" << endl;
    for(int i = 0; i < n+1; i++) {
            for(int k = 0; k < m; k++){
              if (cplex.isExtracted(mapVar["a#" + to_string(i) + "#" + to_string(k)])
              		&& cplex.getValue(mapVar["a#" + to_string(i) + "#" + to_string(k)]) >= 1 - tolerance) {
              cout <<  i << "," << k << "\n";
            }

          }
      }
    current_sol = fill_solution(cplex);
    best_sol = fill_solution(cplex);
    best_sol_value = getValue();
    cout << "BEST VALUE :::::::::::: " << best_sol_value << "\n";
	return 0;
}

Solution MH::fill_solution(IloCplex & cplex){
	int n = dat.n;
	int m = par.m;
	IloNum tolerance = cplex.getParam(
	       IloCplex::Param::MIP::Tolerances::Integrality);
	Solution * sol = new Solution();
    for(int k = 0; k < m; k++) {
    		Tour * tour = new Tour();
            for(int i = 0; i < n+1; i++){
            	float value = cplex.getValue(mapVar["a#" + to_string(i) + "#" + to_string(k)]);
            	if (value  >= 1 - tolerance) {
            	  tour->hops[value] = i;
            	  if (i!=0 && i!=dat.n+1){
            		  sol->nodes_in[i]=i;
            		  IloNumVar u = mh_mapVar["u#" + to_string(i)];
            		  u.setBounds(1,1);
            		  string cname = "mh_y#" + to_string(i);
            		  IloConstraint ic = mh_mapConstr[cname];
            		  model.add(ic);
            	  }
            }

          }
          sol->tours.push_back(*tour);
      }
    for(int i =1; i < n+1; i++){
    	if (sol->nodes_in.find(i)==sol->nodes_in.end()){
    		string cname = "mh_y#" + to_string(i);
    		IloConstraint ic = mh_mapConstr[cname];
    		model.remove(ic);
    		IloNumVar u = mh_mapVar["u#" + to_string(i)];
    		u.setBounds(0,1);
    	}
    }


    //test print
    cout << "TESTING FILL SOLUTION\n";
    for(Tour tour : sol->tours){
    	for(auto const &ent : tour.hops){
    		cout<< ent.second <<" at: " << ent.first << "\n";
    	}
    }

	return *sol;
}

int MH::update_best(){
	return best_sol_value = getValue();
	best_sol = current_sol;
}
float MH::getValue(){
	float value = 0;
	for (Tour tour : current_sol.tours){
		for(auto const & ent : tour.hops){
			value += dat.p[ent.second];
		}
	}
	return value;
}
float MH::getBestValue(){
	return best_sol_value;
}
int MH::iteration(){
	//compute probability
	vector<int> in_sol;
	for (Tour tour : current_sol.tours){
		for(auto const & ent : tour.hops){

			dyn_scores[ent.second] += dat.p[ent.second] / (max(1,(int)tour.hops.size()));
			in_sol.push_back(ent.second);
		}
	}
	// remove mh constraints

	//variable u(i)= 0 if constrained to be not in solution = 1 if in solution
	//removed (with parallel constraints if forced to be free
	for (int i = 1; i < dat.n+1; i++){
		if (mh_mapVar.find("u#"+ to_string(i)) != mh_mapVar.end()){
			cout << "ciao\n";
		}
	}


	return 0;
}
int MH::destroy(){
	return 0;
}
int MH::rebuild(){
	return 0;
}


MH::~MH() {
	// TODO Auto-generated destructor stub
}

