// -------------------------------------------------------------- -*- C++ -*-
// File: ilodiet.cpp
// Version 12.6.1
// --------------------------------------------------------------------------
// Licensed Materials - Property of IBM
// 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
// Copyright IBM Corporation 2000, 2014. All Rights Reserved.
//
// US Government Users Restricted Rights - Use, duplication or
// disclosure restricted by GSA ADP Schedule Contract with
// IBM Corp.
// --------------------------------------------------------------------------
//
// A dietary model.
//
// Input data:
// foodMin[j]      minimum amount of food j to use
// foodMax[j]      maximum amount of food j to use
// foodCost[j]     cost for one unit of food j
// nutrMin[i]      minimum amount of nutrient i
// nutrMax[i]      maximum amount of nutrient i
// nutrPer[i][j]   nutrition amount of nutrient i in food j
//
// Modeling variables:
// Buy[j]          amount of food j to purchase
//
// Objective:
// minimize sum(j) Buy[j] * foodCost[j]
//
// Constraints:
// forall foods i: nutrMin[i] <= sum(j) Buy[j] * nutrPer[i][j] <= nutrMax[j]
//

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

// declaration of function which will build the model; defined later
void
buildModelByRow(IloModel           mod,
                IloNumVarArray     Buy,
                const IloNumArray  foodMin,
                const IloNumArray  foodMax,
                const IloNumArray  foodCost,
                const IloNumArray  nutrMin,
                const IloNumArray  nutrMax,
                const IloNumArray2 nutrPer,
                IloNumVar::Type type);

// main
int mymain(int argc, char **argv);
int mymain(int argc, char **argv)
{
   IloEnv env;

   try {
      const char*     filename  = "diet.dat";
      IloBool         byColumn  = IloFalse;


      //ILOFLOAT , ILOINT , or ILOBOOL for continuous, integer, or Boolean variables, respectively.
      // example: The following constructor creates an integer variable with bounds -1 and 10:
      //loNumVar myIntVar(env, -1, 10, ILOINT);

      // for the diet problem we need float
      IloNumVar::Type varType   = ILOFLOAT;  //Float variable type

      IloInt i;
	  cout << "parameters from the file diet.dat " << filename << "\n";
      ifstream file(filename);
      if ( !file ) {
         cerr << "ERROR: could not open file '" << filename
              << "' for reading" << endl;
         throw (-1);
      }

      // model data

      IloNumArray  foodCost(env), foodMin(env), foodMax(env); // parameters arrays
      IloNumArray  nutrMin(env), nutrMax(env); //parameters arrays
      IloNumArray2 nutrPer(env);  // 2 dimensional array

      file >> foodCost >> foodMin >> foodMax;
      file >> nutrMin >> nutrMax;
      file >> nutrPer;

      IloInt nFoods = foodCost.getSize();
      IloInt nNutr  = nutrMin.getSize();

      if ( foodMin.getSize() != nFoods ||
           foodMax.getSize() != nFoods ||
           nutrPer.getSize() != nNutr  ||
           nutrMax.getSize() != nNutr    ) {
         cerr << "ERROR: Data file '" << filename
              << "' contains inconsistent data" << endl;
         throw (-1);
      }

      for (i = 0; i < nNutr; i++) {
         if (nutrPer[i].getSize() != nFoods) {
            cerr << "ERROR: Data file '" << argv[0]
                 << "' contains inconsistent data" << endl;
            throw (-1);
         }
      }

      // Build model

      IloModel       mod(env);  //iloModel cplex object
      IloNumVarArray Buy(env);  // cplex object used to store variables

      buildModelByRow(mod, Buy, foodMin, foodMax, foodCost,
                         nutrMin, nutrMax, nutrPer, varType);

      // Solve model

      IloCplex cplex(mod);
      cplex.exportModel("diet.lp");  // save an exported model into lp file format

      cplex.solve();  // solve the model
      cplex.out() << "solution status = " << cplex.getStatus() << endl;  // print the status

      cplex.out() << endl;
      cplex.out() << "cost   = " << cplex.getObjValue() << endl;  // object value

      // print out the variable array Buy[]
      for (i = 0; i < foodCost.getSize(); i++) {
         cplex.out() << "  Buy" << i << " = " << cplex.getValue(Buy[i]) << endl;
      }
   }
   catch (IloException& ex) {
      cerr << "Error: " << ex << endl;
   }
   catch (...) {
      cerr << "Error" << endl;
   }

   env.end();

   return 0;
}


void
buildModelByRow(IloModel           mod,
                IloNumVarArray     Buy,
                const IloNumArray  foodMin,
                const IloNumArray  foodMax,
                const IloNumArray  foodCost,
                const IloNumArray  nutrMin,
                const IloNumArray  nutrMax,
                const IloNumArray2 nutrPer,
                IloNumVar::Type type) {

   IloEnv env = mod.getEnv();  // get the cplex environment object

   Buy.clear();  // initialize Buy[] array of variables
   IloNumVarArray tmp(env, foodMin, foodMax, type);  //declare variables as foodMin <= Buy <= foodMax; the type is ILOFLOAT
   //add the variable to Buy objext
   Buy.add(tmp);
   //free memory for tmp
   tmp.end();

   IloInt i, j;
   IloInt n = foodCost.getSize();
   IloInt m = nutrMin.getSize();

   // add the objective function  to cplex environment. the objective function is the scalar product of variable array Buy per foodCost parameters array
   mod.add(IloMinimize(env, IloScalProd(Buy,foodCost)));
   // add constraints in a for loop

   for (i = 0; i < m; i++) {
      //initialize an expression expr
      IloExpr expr(env);

      // in the loop sum up the constraint. Each loop is a multiplication of the variable buy[i] multiplied by the relative coefficient nutrPer[i][j]
      for (j = 0; j < n; j++) {
         expr += Buy[j] * nutrPer[i][j];
      }
      // add the constraint
      mod.add(nutrMin[i] <= expr <= nutrMax[i]);
      expr.end();
   }
}




/*
cost   = 14.8557
  Buy0 = 4.38525
  Buy1 = 0
  Buy2 = 0
  Buy3 = 0
  Buy4 = 0
  Buy5 = 6.14754
  Buy6 = 0
  Buy7 = 3.42213
  Buy8 = 0
*/
