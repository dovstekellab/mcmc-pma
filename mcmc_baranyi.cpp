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
#include "mcmc_baranyi.h"

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

#define T0 2500




// Class constructor. This is called in addition to the mcmc constructor.
mcmc_baranyi::mcmc_baranyi(std::vector<double> od, std::string l) : mcmc(od, l)
{
    // Set dimension of parameter space.
    dim = 6;
    // Give parameters names.
    parameternames.resize(dim);
    parameternames[0] = "lag";
    parameternames[1] = "y0";
    parameternames[2] = "ymax";
    parameternames[3] = "umax";
    parameternames[4] = "m";
    parameternames[5] = "nu";
    modelname = "Baranyi-mn";
    // Declar which parameters to work with in logarithmic scale.
    bool logsp[] = {0,0,0,0,1,1,1};
    logspace.assign(logsp, logsp+dim);
    // Guess initial parameter vector to initialise MCMC chain. Some parameters are guessed at the beginning of the MCMC function - see comment at end of this file.
    double startv[] = {odlag() ,ody0(),odymax(),0.1,1,1};

    // Guess initial variances.
    double p[] = {0.5,2,2,0.001,0.000025, 0.000001};
    std::vector<double> startp(startv, startv+dim);
    
    
    startparameters = lintolog(startp);
    startvariances.assign(p, p+dim);
    
    
};


// Function that takes a parameter vector and returns modelled data. This is "the model".
std::vector<double> mcmc_baranyi::md(std::vector<double> parameters){
    std::vector<double> model_data(observed_data.size());
    
    std::vector<double> p = logtolin(parameters);
  

    
    double q0 = (double)1/(exp(p[5]*p[0])-1);
    
    
    for (int i=0; i < observed_data.size(); i++) {


        model_data[i] = exp(
                            log(p[1])
                            +
                            p[3] * i
                            +
                            p[3]/p[5] * log( (exp(-p[5]*i) + q0)/(1+q0) )
                            -
                            1/p[4] * log( 1 + 
                                         ( exp( p[4]*p[3]*i + p[4] * p[3]/p[5] * log((exp(-p[5]*i) + q0)/(1+q0)) ) - 1 ) 
                                         / 
                                         (exp(p[4]*(log(p[2])-log(p[1])))) 
                                         )
                            
                            );
        


    }   

    
    return model_data;
}




// Function that returns true iff a parameter vector is within our (flat) prior.
bool mcmc_baranyi::checkparameters(std::vector<double> parameters){
    std::vector<double> p = logtolin(parameters);
    if (p[0] < 0) return FALSE;
    if (p[0] > observed_data.size()) return FALSE;
    if (p[2] < p[1]) return FALSE;
    if (p[1] < 0) return FALSE;
    if (p[2] < 0) return FALSE;
    if (p[3] < 0 || p[3] > 1) return FALSE;
    if (p[4] < 1 || p[4] > 100) return FALSE;
    if (p[5] < p[3] || p[5] > 100) return FALSE;
    
    
    return TRUE;
}


// Optimise initial guess of further parameters. This is done at the beginning of the MCMC chain instead of the class constructor, as this gives a slight performance benefit when using OpenMP. (MCMC chains are run in parallel, class constructors are not. When not using OpenMP this makes no difference.)
bool mcmc_baranyi::do_before_run(){
    
    startparameters[3] = odr(startparameters, 3, 0.01, 1);
    startparameters[4] = odoptimiseparameter(startparameters, 4, 1, 10);
    startparameters[5] = odoptimiseparameter(startparameters, 5, exp(startparameters[3]), 1);
    startparameters[0] = odoptimiseparameter(startparameters, 0, startparameters[0] - 5, startparameters[0] + 5);
    
    return TRUE;
    
}






