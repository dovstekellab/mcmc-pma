//
//  mcmc_baranyi.cpp
//mcmc-pma, Bayesian parameter estimation for phenotype microarray data.
//
//Copyright 2011-2015 Matthias Gerstgrasser.
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "mcmc.h"
#include "mcmc_logistic.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include <time.h>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>


mcmc_logistic::mcmc_logistic(std::vector<double> od, std::string l) : mcmc(od, l) 


{
    dim = 3; 
    parameternames.resize(dim);
    parameternames[0] = "y0";
    parameternames[1] = "ymax";
    parameternames[2] = "r";
    //parameternames[3] = "t0";
    modelname = "Logistic";
    bool logsp[] = {0,0,1};
    logspace.assign(logsp, logsp+dim);
    double startv[] = {ody0(), odymax(), 0.05};
    double p[] = {0.25,0.25,0.00000025};
    std::vector<double> startp(startv, startv+dim);
    startparameters = lintolog(startp);
    startvariances.assign(p, p+dim);

};

bool mcmc_logistic::do_before_run(){
    
    startparameters[2] = odr(startparameters, 2, 0.001, 1);
    return TRUE;
}

double mcmc_logistic::results_lagbyparams(std::vector<double> parameters){
    return logtolin(parameters)[3]/4;
}

std::vector<double> mcmc_logistic::md(std::vector<double> parameters){
    std::vector<double> model_data(observed_data.size());
    
    std::vector<double> p = logtolin(parameters);
    
    
    int i=0;
    
    /*
    while ( i < p[3] && i < observed_data.size() ) {
        model_data[i] = p[0];
        i++;
    }
     */
    
    
    
    while (i < observed_data.size()) {
        model_data[i] = p[0] * p[1] * exp(p[2] * (i)) 
        / 
        (p[1] + (p[0]) * (exp(p[2]*(i)) - 1) ); 
        i++;
    }   

    
    return model_data;
}

bool mcmc_logistic::checkparameters(std::vector<double> parameters){
    std::vector<double> p = logtolin(parameters);
    //if (p[3] < 0) return FALSE;
    //if (p[3] > observed_data.size()) return FALSE;
    if (p[1] < p[0]) return FALSE;
    //if (p[3] < 0 || p[3] > observed_data.size()) return FALSE;
    return TRUE;
}


