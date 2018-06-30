#ifndef toptw_H_
#define toptw_H_

#include <string>
#include <ilcplex/ilocplex.h>


const bool debug = false;

struct coord {
  float x;
  float y;
};
struct tw {
    int l;
    int r;
};

struct params {
	double timeLimit = 2;
	int m = 1;
};

struct toptwdata {
    vector<coord> coords;
    vector<int> p;
    vector<int> id;
    int n;
    int m;
    int tau;
    float **t;
    vector<tw> atw;
    vector<int> s; //service time
    
    string to_string(){
        string ss = "";
       
        ss+= "n: " + std::to_string(n) + "\n";
        ss+= "m: " + std::to_string(m) + "\n";
        ss+= "Tmax: " + std::to_string(tau) + "\n";
        ss+= "id\tx\ty\ts\tp\tl\tr\n\n";
        for (auto i : id){
            ss += std::to_string(id[i]) + "\t" + std::to_string(coords[i].x) + "\t" 
               + std::to_string(coords[i].y) + "\t" + std::to_string(s[i]) + "\t" 
               +  std::to_string(p[i]) + "\t" + std::to_string(atw[i].l) + "\t" 
               + std::to_string(atw[i].r) + "\n";
            
        }
 
        return ss;
    }
  
    
};

int load_data(vector<int> * profits, vector<coord>* coords, int* n, int* m, float* tau, string filename );
int load_data_tw(toptwdata * data, string filename);
int init_sol(IloModel* env, IloCplex* cplex, toptwdata * data, params* par);
IloModel build_model(IloEnv& env, toptwdata & dat, params& par, map<string, IloConstraint>& mapConstr, map<string, IloNumVar>& vars, float  **t);
int solve_exactly(IloEnv& env, IloModel model, toptwdata& dat, params& par,
		map<string, IloConstraint>& mapConstr, map<string, IloNumVar>& vars, float **t );
#endif

