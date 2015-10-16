// ------------------------------------------------------------------
// Copyright (C) 2004 INRA <eugene@ossau.toulouse.inra.fr>
//
// This program is open source; you can redistribute it and/or modify
// it under the terms of the Artistic License (see LICENSE file).
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
//
// You should have received a copy of Artistic License along with
// this program; if not, please see http://www.opensource.org
//
// $Id: OptiAlgorithm.cc,v 1.12 2008-04-23 08:36:54 tschiex Exp $
// ------------------------------------------------------------------
// File:     OptiAlgorithm.cc
// Contents: The mother class of optimization algorithms classes
// ------------------------------------------------------------------


#include <iostream>
#include <iomanip>

#include "OptiAlgorithm.h"

#include "../Param.h"

extern Parameters PAR;

//-------------------------------------------------------
// Constructor
//-------------------------------------------------------
OptiAlgorithm::OptiAlgorithm(void)
{
  int n, c;
  std::string relation;
  std::vector <int> v;
  std::string para_optimization_trace;

  para_optimization_trace = (std::string) PAR.getC("ParaOptimization.Trace");
  IsTracing = ( ( (para_optimization_trace == "1") || (para_optimization_trace == "TRUE") ) ? true : false );

  NbParaCluster = PAR.getI("ParaOptimization.NbCluster");
  for (int i = 0; i < NbParaCluster; i++) {
    ParaClusters.push_back( v );
    relation = (std::string) PAR.getC("ParaOptimization.Cluster", i);
    if (relation == "IDENTICAL")
      ParaClusterRelations.push_back(IDENTICAL); 
    else
      if (relation == "LINKED")
	ParaClusterRelations.push_back(LINKED); 
      else 	{
	std::cerr <<"ERROR: Bad cluster relation "<<relation<<"  in the parameter file."<<std::endl; 
	exit(100);
      }
  }

  n = PAR.getI("ParaOptimization.NbParameter");
  for (int i=0; i<n; i++) {
    ParaName.push_back( PAR.getC("ParaOptimization.Para.Name", i) );
    ParaMin.push_back( PAR.getD("ParaOptimization.Para.Min" , i) );
    ParaMax.push_back( PAR.getD("ParaOptimization.Para.Max" , i) );
    c = PAR.getI("ParaOptimization.Para.Cluster", i);
    ParaCluster.push_back( c );
    if (c < ParaClusters.size())
        ParaClusters[c].push_back( i );
    else
	{
	    std::cerr <<"ERROR: Bad cluster number (" << c << ") for parameter "<< ParaName.back() <<" in the parameter file."<<std::endl; 
	    exit(100);
	}
  }
  for (int i=0; i<n; i++) 
    // Name of parameters to optimize must be <name>* or <name>*[<nn>], 
    if (!( (ParaName[i][ParaName[i].size()-1] == '*') ||
	 (ParaName[i][ParaName[i].size()-1] == ']' &&  
	  (ParaName[i][ParaName[i].size()-4] == '*' || ParaName[i][ParaName[i].size()-5] == '*')) )) {
      std::cerr <<"ERROR "<<ParaName[i]<<": Parameters to optimize must have a name finishing by '*' and a number of instance < 100."<<std::endl;
      exit(100);
    }

  // Set numerical precision for displayed values
  std::cout.setf(std::ios::showpoint);
  std::cout << setiosflags (std::ios::fixed);
  std::cout.precision(4);
}

//-------------------------------------------------------
// return a string with maximum 6 signs
// if s (without [nn]) has more than 6 letters plus '*' at the end
// returns: the 2 first letters + '_' + the 3 last (without counting '*[nn]')
//-------------------------------------------------------
std::string OptiAlgorithm::ReduceName(std::string s)
{
  std::string ss, s0, s1, sn3, sn2, sn1;
  int n; // position of the last letter to consider

  if (s[s.length()-1] == ']') {
    if (s[s.length()-3] == '[') n = s.length()-4;
    else n = s.length()-5;  // assume maximum 2 digits for the number of instance
  } else n = s.length()-1; // ignore the last caracter which is always '*'
  
  if (n > 6) {
   s0 = s[0]; s1 = s[1]; sn3 = s[n-3]; sn2 = s[n-2]; sn1 = s[n-1];
   ss = s0 + s1 + "_"+ sn3 + sn2 + sn1;
 } else
   for (int k=0; k<n; k++) ss += s[k];
 
 return ss;
}

