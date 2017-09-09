#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <map>
#include "graphg.h"
#include "toptw.h"
#include "MH.h"
#include <ilcplex/ilocplex.h>

using namespace std;
ILOSTLBEGIN

int main()
{
    bool SOLVE_EXACTLY = false;
	vector<int> p;
    vector<coord> coords;
    vector<Node> nodes;
    toptwdata dat;
    params pam;
    //string filename = "r101.txt";
    string filename = "pr01.txt";
    load_data_tw(&dat, filename);
    int n = dat.n;
    tw atw;
    atw.l = dat.atw[0].l;
    atw.r = dat.atw[0].r;
    dat.atw.push_back(atw);
    dat.p.push_back(0);
    
    string s = dat.to_string();
    cout << s << "\n";
    coords = dat.coords;
    p = dat.p;
    cout << "SIZE OF P: " << p.size() << "\n";

    
    //load_data(&p, &coords, &n, &m, &tau, "../p4.2.a.txt");

    //float t[n+2][n+2];  // in position 0 and in position n+1 there will be the depot
    float ** t;
    t = new float *[n+2];
    for(int r = 0; r < n+1; r++) {
    	t[r] = new float[n+2];
        for(int c = 0; c < n+1; c++) {
            t[r][c] = sqrt(pow(coords[r].x - coords[c].x, 2) + pow(coords[r].y - coords[c].y, 2));
        }
        //t[r][n+1] = t[r][0];
        Node aNode(r);
        // aNode.id = r;
        nodes.push_back(aNode);
        //t[n][n] = t[0][0];
    }
    t[n+1] = new float[n+2];
    for (int i=1; i<n+1; i++){
        t[i][n+1] = t[i][0];
        t[n+1][i] = t[0][i];
    }
    t[0][n+1]=0;
    t[n+1][0]=0;
    t[n+1][n+1]=0;
    t[0][0]=0;
    Node aNode(n);
    nodes.push_back(aNode);
    dat.t = t;

    //constraints container
     map<string,IloConstraint> mapConstr; // map containing all model constraints
    // variable container
     map<string, IloNumVar> mapVar;
    
    if(debug){
        cout << "PRINTING P" << "\n";
        for (int i : dat.p){
            cout << i << "\t";
        }
        cout << "\n";
        cout << "TIME WINDOWS\n";
        for (tw a : dat.atw)
            cout << a.l << " - " << a.r << "\n";
        
        cout << "\n travelling timesTIME WINDOWS\n";
        for (int i = 0; i < n+2; i++)
            for(int j = 0; j < n+2; j++)
                cout << i << "-" << j << ": " << t[i][j] << "\n";
    }
    
    try{
    // cplex model
      IloEnv env;
      IloModel model = build_model(env, dat, pam, mapConstr, mapVar, t);
      if(!SOLVE_EXACTLY){

    	  MH mh(env, model, dat, pam, mapConstr, mapVar, t);
    	  mh.init_sol();
    	  mh.update_best();
    	  for (int it=0; it <20; it++){
    		  mh.iteration();
    		  if (mh.getValue() < mh.getBestValue())
    			  mh.update_best();
    		  else
    			  mh.destroy();
    		  	  mh.rebuild();
    	  }
      } else
    	  solve_exactly(env, model, dat, pam, mapConstr, mapVar, t);

      /*
      *
      *
      *
     IloModel model(env);
      for(i = 0; i < nbClients; i++)
         model.add(IloSum(supply[i]) == 1);
      for(j = 0; j < nbLocations; j++) {
         IloExpr v(env);
         for(i = 0; i < nbClients; i++)
            v += supply[i][j];
         model.add(v <= capacity[j] * open[j]);
         v.end();
      }

      IloExpr obj = IloScalProd(fixedCost, open);
      for(i = 0; i < nbClients; i++) {
         obj += IloScalProd(cost[i], supply[i]);
      }
      model.add(IloMinimize(env, obj));
      obj.end();

      IloCplex cplex(env);
      cplex.extract(model);
      cplex.solve();

      cplex.out() << "Solution status: " << cplex.getStatus() << endl;

      IloNum tolerance = cplex.getParam(
         IloCplex::Param::MIP::Tolerances::Integrality);
      cplex.out() << "Optimal value: " << cplex.getObjValue() << endl;
      for(j = 0; j < nbLocations; j++) {
         if (cplex.getValue(open[j]) >= 1 - tolerance) {
            cplex.out() << "Facility " << j << " is open, it serves clients ";
            for(i = 0; i < nbClients; i++) {
               if (cplex.getValue(supply[i][j]) >= 1 - tolerance)
                  cplex.out() << i << " ";
            }
            cplex.out() << endl;
         }
      }

    */
      env.end();
    }
    catch (IloException& ex) {
      cerr << "Error: " << ex << endl;
   }
   catch (...) {
      cerr << "Error" << endl;
   }



   return 0;

}

IloModel build_model(IloEnv& env, toptwdata& dat, params& par,
		map<string, IloConstraint>& mapConstr, map<string, IloNumVar>& vars, float **t ){
	int n = dat.n;
	int m = par.m;
	int tau = dat.tau;
    IloModel model(env);
    int bigM = tau*2;


    // flow variable
    IloBoolVar x[n+2][n+2][m];
    for(int r = 0; r < n+2; r++) {
        for(int c = 0; c < n+2; c++) {
            for(int k=0; k<m; k++) {
                string varName = "x#" + to_string(r) + "#" + to_string(c) + "#" + to_string(k);
                x[r][c][k] = IloBoolVar(env, 0, 1, varName.c_str());
                vars[varName] = x[r][c][k];
            }

        }
    }
    IloBoolVar y[n+2][m];
    IloFloatVar a[n+2][m];
    for(int i = 0; i < n+2; i++){
        for (int v = 0; v < m; v++){
            string varName1 = "y#" + to_string(i) + "#" + to_string(v);
            y[i][v] = IloBoolVar(env, 0, 1, varName1.c_str());
            vars[varName1] = y[i][v];
        }
    }
    for(int i = 0; i < n+2; i++){
        for (int v = 0; v < m; v++){
            string varName2 = "a#" + to_string(i) + "#" + to_string(v);
            a[i][v] = IloFloatVar(env, 0, tau, varName2.c_str());
            vars[varName2] = a[i][v];
        }
    }
    
    
    
    

    //IloRange(const IloEnv env, const IloNumExprArg expr, IloNum rhs=IloInfinity, const char * name=0)
    
    // (2) link x y for nodes
    for(int i = 1; i < n+1; i++){
        for (int v = 0; v < m; v++){
            IloExpr expr(env);
            for (int j = 1; j< n+2; j++){
                if (i!=j)
                    expr += x[i][j][v];
            }
            expr -= y[i][v];
            string cname = "xtoy#" + to_string(i) + "#" + to_string(v);
            IloConstraint ic(expr==0);
            ic.setName(cname.c_str());
            mapConstr[cname] = ic;
            model.add(ic);
            expr.end();
        }
    }
    
    // (3) flow constraints
    for(int j = 1; j < n+1; j++){
        for (int v = 0; v < m; v++){
            IloExpr expr(env);
            for (int i = 0; i < n+1; i++){
                expr += x[i][j][v];
            }
            for (int i = 1; i < n+2; i++){
                expr -= x[j][i][v];
            }
            string cname = "flow#" + to_string(j) + "#" + to_string(v);
            IloConstraint ic(expr==0);
            ic.setName(cname.c_str());
            mapConstr[cname] = ic;
            model.add(ic);
            expr.end();
        }
    }
    
    // (4a) x depart m vehicles
    IloExpr expr(env);
        for (int j = 1; j < n+1; j++){
            for(int v = 0; v < m; v++){
                expr += x[0][j][v];

            }
        }
    string cname = "xdep";
    IloConstraint ic(expr==m);
    ic.setName(cname.c_str());
    mapConstr[cname] = ic;
    model.add(ic);
    // (4b) y depart m vehicles
    //expr.clear();
    expr = IloExpr(env);
    for(int v = 0; v < m; v++){
        expr += y[0][v];
    }
    cname = "ydep";
    ic = IloConstraint(expr==m);
    ic.setName(cname.c_str());
    mapConstr[cname] = ic;
    model.add(ic);

    // (4c)  x arrival m vehicles
    expr = IloExpr(env);
    for (int j = 1; j < n+1; j++){
        for(int v = 0; v < m; v++){
            expr += x[j][n+1][v];
        }
    }
    cname = "xarr";
    ic = IloConstraint(expr==m);
    ic.setName(cname.c_str());
    mapConstr[cname] = ic;
    model.add(ic);
    expr.end();
    
    // (4d) y arrival m vehicles
    expr = IloExpr(env);
    for(int v = 0; v < m; v++){
       expr += y[n+1][v];
    }
    cname = "ydep";
    ic = IloConstraint(expr==m);
    ic.setName(cname.c_str());    mapConstr[cname] = ic;
    model.add(ic);
    expr.end();
    
    //(5) max 1 node
 
    for (int i = 1; i<n+1; i++){
        expr = IloExpr(env);
        for(int v=0; v < m; v++){
            expr += y[i][v];
        }
        ic = IloConstraint(expr <= 1);
        string cname = "max1visit#" + to_string(i);
        ic.setName(cname.c_str());
        mapConstr[cname] = ic;
        model.add(ic);
        expr.end();
    }
    
    // (6) sequence
    for(int i = 0; i < n+1; i++){
        for (int j = 1; j < n+2; j++){
            if (i==j) continue;
            for (int v = 0; v < m; v++){
                IloExpr expr(env);
                cout << "t["<<i<<"]["<<j<<"]: " << t[i][j] << "\n";
                expr = a[i][v] + dat.s[i] + t[i][j] - a[j][v] + bigM*x[i][j][v];
                IloConstraint ic(expr <= bigM);
                string cname = "cseq#" + to_string(i) + "#" + to_string(j) + "#" + to_string(v);
                ic.setName(cname.c_str());
                mapConstr[cname] = ic;
                model.add(ic);
                expr.end();
            }
        }
    }
    
    //(7a) left time window ltw
    for (int i = 0; i < n+2; i++){
        for (int v = 0; v < m; v++){
            IloExpr expr(env);
            expr = dat.atw[i].l*y[i][v] - a[i][v];
            IloConstraint ic(expr <= 0);
            string cname = "cltw#" + to_string(i) + "#" + to_string(v);
            ic.setName(cname.c_str());
            mapConstr[cname] = ic;
            model.add(ic);
            expr.end();
        }
    }
    
    //(7b) right time window rtw
    for (int i = 0; i < n+2; i++){
        for (int v = 0; v < m; v++){
            IloExpr expr(env);
            expr = a[i][v] - dat.atw[i].r*y[i][v];
            IloConstraint ic(expr <= 0);
            string cname = "crtw#" + to_string(i) + "#" + to_string(v);
            ic.setName(cname.c_str());
            mapConstr[cname] = ic;
            model.add(ic);
            expr.end();
        }
    }
    
     // SOLVE AND SHOW THE SOLUTION
     IloExpr obj(env);
     for(int i = 1; i<n+1; i++){
        for(int v = 0; v < m; v++) {
            obj += dat.p[i]*y[i][v];
        }
     }

      model.add(IloMaximize(env, obj));
      obj.end();
      
	return model;
}

int solve_exactly(IloEnv& env, IloModel model, toptwdata& dat, params& par,
		map<string, IloConstraint>& mapConstr, map<string, IloNumVar>& mapVar, float **t ){
	int n = dat.n;
	int m = par.m;

    IloCplex cplex(model);
    cplex.exportModel("model.lp");

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

	return 0;
}

int load_data_tw(toptwdata * dat, string filename)
//int load_data_tw(vector<int>* profits, vector<coord>* coords, int* n, int* m, float* tau, string filename)
{
    ifstream myfile(filename);
    cout << "filename: " << filename << "\n";
    string line;
    if(myfile) {
        //* dat = new toptwdata()
        getline(myfile, line);
        string delimiter = " ";
        string token;
        int startsub = 0;
        size_t pos = line.find(delimiter);
        token = line.substr(startsub, pos);
        //

        startsub = pos + 1; //
        pos = line.find(delimiter, startsub);
        token = line.substr(startsub, pos - startsub);
            // coords[i][1] = stof(token);
        dat->m = stoi(token);
        //cout << dat->m << "\n";

        startsub = pos + 1; //
        pos = line.find(delimiter, startsub);
        token = line.substr(startsub, pos - startsub);
            // coords[i][1] = stof(token);
        dat->n  = stoi(token);

        getline(myfile, line); //skip a line

       // c->y = stof(token);
            // cout << "y coord: " << token << endl;
            // profit

        int i = 0;
        while(getline(myfile, line)) {
            /*
            * 	i x y d S f a list O C
            *
            *   Where
            *   i = vertex number
            *   x = x coordinate
            *   y = y coordinate
            *   d = service duration or visiting time
            *   S = profit of the location
            *   f = not relevant
            *   a = not relevant
            *   list = not relevant (length of the list depends on a)
            *   O = opening of time window (earliest time for start of service)
            *   C = closing of time window (latest time for start of service)
            *
            *    * REMARKS *
            *    - The first point (index 0) is the starting AND ending point.
            *    - The number of paths (P) is not included in the data file. This
            *
            *
            *
            */
            //cout << "line: " << line << "\n";

            string buf; // Have a buffer string
            stringstream ss(line); // Insert the string into a stream

            vector<string> tokens; // Create vector to hold our words


            while (ss >> buf)
                tokens.push_back(buf);
            /*cout << "LINE TOKENIZED  \n";
            for (string aaa : tokens)
                cout << aaa << "*";
            cout << "\n";
            */
            coord c = coord();
            tw ttw = tw();

            dat->id.push_back(stoi(tokens[0]));
            float x = stof(tokens[1]);
            float y = stof(tokens[2]);
            c.x = x;
            c.y = y;
            dat->coords.push_back(c);
            //visit time
            dat->s.push_back(stoi(tokens[3]));
            // profit
            dat->p.push_back(stoi(tokens[4]));
            // time windows
            ttw.l = stoi(tokens[tokens.size()-2]);
            ttw.r = stoi(tokens[tokens.size()-1]);
            dat->atw.push_back(ttw);
            if(i == 0){
                dat->tau = ttw.r;
                i = 1;
            }
        }

    } else
        cout << "error opening file\n";



    return 0;
}



int load_data(vector<int>* profits, vector<coord>* coords, int* n, int* m, float* tau, string filename)
{
    ifstream myfile(filename);

    string line;
    if(myfile) {
        getline(myfile, line);
        string ns = line.substr(2, sizeof(line) - 2);
        *n = std::stoi(ns);
        getline(myfile, line);
        string ms = line.substr(2, sizeof(line) - 2);
        *m = stoi(ms);
        getline(myfile, line);
        string ts = line.substr(5, sizeof(line) - 5);
        *tau = stof(ts);
        cout << "n: " << *n << endl;
        cout << "m: " << *m << endl;
        cout << "tau: " << *tau << endl;
        string delimiter = "\t";
        string token;
        int i = 0;
        while(getline(myfile, line)) {
            coord* c = new coord();
            int startsub = 0;
            size_t pos = line.find(delimiter);
            // x coord
            token = line.substr(startsub, pos);
            // coords[i][0] = stof(token);
            c->x = stof(token);

            // cout << "x coord: " << token << endl;
            // y coord
            startsub = pos + 1;
            pos = line.find(delimiter, startsub);
            token = line.substr(startsub, pos - startsub);
            // coords[i][1] = stof(token);
            c->y = stof(token);
            // cout << "y coord: " << token << endl;
            // profit
            startsub = pos;
            token = line.substr(pos + 1, sizeof(line) - pos + 1);
            // profits[i] = stoi(token);
            profits->push_back(stoi(token));
            coords->push_back(*c);
            // i++;
        }

    } else
        cout << "error opening file\n";

    return 0;
}
