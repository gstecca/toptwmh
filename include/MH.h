/*
 * MH.cpp
 *
 *  Created on: Jul 3, 2017
 *      Author: gstecca
 */


#ifndef MH_H_
#define MH_H_
#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <map>
#include "graphg.h"
#include "toptw.h"
#include <ilcplex/ilocplex.h>

struct Tour{
	map<float, int> hops;

};
struct Solution {
	vector<Tour> tours;
	map<int,int> nodes_in;
	float getObj(toptwdata & dat){
		float obj = 0;
		for(Tour aT : tours){
			for (auto const & hop : aT.hops)
				obj += dat.p[hop.second];
		}
		return -1;
	}
};



class MH {
	IloEnv env;
	IloModel model;
	toptwdata dat;
	params par;
	map<string, IloConstraint> mapConstr;
	map<string, IloNumVar> mapVar;
	map<string, IloConstraint> mh_mapConstr;
	//map<string, IloConstraint *> mh_mapConstr_free;
	map<string, IloNumVar> mh_mapVar;
	Solution best_sol;
	Solution current_sol;
	vector<float> dyn_scores;
	float best_sol_value = 0;
	float  **t;
	Solution fill_solution(IloCplex & cplex);
public:
	MH(IloEnv env, IloModel _model, toptwdata dat,
			params par, map<string, IloConstraint> mapConstr,
			map<string, IloNumVar> vars, float  **t);
	int init_sol();
	int update_best();
	float getValue();
	float getBestValue();
	int iteration();
	virtual ~MH();
};



#endif
