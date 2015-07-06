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
#include "mcmc_gompertz.h"

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






mcmc_gomp::mcmc_gomp(std::vector<double> od, std::string l) : mcmc(od, l)


{
    dim = 4;
    parameternames.resize(dim);
    parameternames[0] = "a";
    parameternames[1] = "k";
    parameternames[2] = "b";
    parameternames[3] = "y0";
    //parameternames[3] = "t0";
    modelname = "Gompertz";
    bool logsp[] = {0,0,1,0};
    logspace.assign(logsp, logsp+dim);
    double startv[] = {2, odymax(), 0.02, ody0()};
    double p[] = {1,0.01,0.5,1};
    std::vector<double> startp(startv, startv+dim);
    startparameters = lintolog(startp);
    startvariances.assign(p, p+dim);
    
};

bool mcmc_gomp::do_before_run(){
    //startparameters[2] = odr(startparameters, 2, 0.001, 5);
    //startparameters[0] = odoptimiseparameter(startparameters, 0, 0.001, 5);
    return TRUE;
}

double mcmc_gomp::results_lagbyparams(std::vector<double> parameters){
    return logtolin(parameters)[3]/4;
}

std::vector<double> mcmc_gomp::md(std::vector<double> parameters){
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
        // model_data[i] = p[1] * exp( -exp(p[0] - i*p[2])) ;
        model_data[i] = p[3] * exp( log(p[1]/p[3]) * exp( -exp(p[0] - i*p[2]))) ;
        i++;
    }
    
    
    return model_data;
}

bool mcmc_gomp::checkparameters(std::vector<double> parameters){
    std::vector<double> p = logtolin(parameters);
    if (p[3] <= 0) return FALSE;
    //if (p[3] > observed_data.size()) return FALSE;
    //if (p[1] < p[0]) return FALSE;
    //if (p[3] < 0 || p[3] > observed_data.size()) return FALSE;
    return TRUE;
}

