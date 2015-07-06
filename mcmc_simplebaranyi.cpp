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
#include "mcmc_simplebaranyi.h"

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


mcmc_simplebaranyi::mcmc_simplebaranyi(std::vector<double> od, std::string l) : mcmc(od, l) 


{
    dim = 5; 
    parameternames.resize(dim);
    parameternames[0] = "lag";
    parameternames[1] = "y0";
    parameternames[2] = "ymax";
    parameternames[3] = "umax";
    parameternames[4] = "nu";
    modelname = "Baranyi-nu";
    m=1;
    bool logsp[] = {0,0,0,1,1};
    logspace.assign(logsp, logsp+dim);
    double startv[] = {odlag(),ody0(),odymax(),0.1,0.1};
    double p[] = {0.5,0.1,0.1,0.0001,0.01};
    std::vector<double> startp(startv, startv+dim);
    startparameters = lintolog(startp);
    startvariances.assign(p, p+dim);

    
};

mcmc_simplebaranyi::mcmc_simplebaranyi(std::vector<double> od, std::string l, double mm) : mcmc(od, l) 


{
    dim = 5; 
    parameternames.resize(dim);
    parameternames[0] = "lag";
    parameternames[1] = "y0";
    parameternames[2] = "ymax";
    parameternames[3] = "umax";
    parameternames[4] = "nu";
    modelname = "Baranyi-nu_fixedm";
    m=mm;
    bool logsp[] = {0,0,0,1,1};
    logspace.assign(logsp, logsp+dim);
    double startv[] = {odlag(),ody0(),odymax(),0.1,0.1};
    double p[] = {0.5,0.1,0.1,0.0001,0.01};
    std::vector<double> startp(startv, startv+dim);
    startparameters = lintolog(startp);
    startvariances.assign(p, p+dim);
    
    
};

double mcmc_simplebaranyi::results_lagbyparams(std::vector<double> parameters){
    double t = 0;
    std::vector<double> p = logtolin(parameters);
    t = p[0]; 
    return t/4;
}

std::vector<double> mcmc_simplebaranyi::md(std::vector<double> parameters){
    std::vector<double> model_data(observed_data.size());
    
    std::vector<double> p = logtolin(parameters);
    
    
    double q0 = (double)1/(exp(p[4]*p[0])-1);
    
    for (int i=0; i < observed_data.size(); i++) {
        
        
        model_data[i] = exp(
                            log(p[1])
                            +
                            p[3] * i
                            +
                            (p[3] / p[4]) * log( (exp(-p[4]*i) + q0)/(1+q0) )
                            -
                            1/m * log( 1 + 
                                ( exp( m * p[3]*i + m* (p[3] / p[4]) * log((exp(-p[4]*i) + q0)/(1+q0)) ) - 1 ) 
                                / 
                                (exp(m*(log(p[2])-log(p[1])))) 
                                )
                            
                            );
    }   
    
    
    return model_data;
}

bool mcmc_simplebaranyi::checkparameters(std::vector<double> parameters){
    std::vector<double> p = logtolin(parameters);
    if (p[0] < 0) return FALSE;
    if (p[0] > observed_data.size()) return FALSE;
    if (p[2] < p[1]) return FALSE;
    if (p[1] < 0) return FALSE;
    if (p[2] < 0) return FALSE;
    if (p[3] < 0 || p[3] > 1) return FALSE;
    if (p[4] < p[3] || p[4] > 100) return FALSE; //<0.01
    
    
    return TRUE;
}


bool mcmc_simplebaranyi::do_before_run(){
    startparameters[3] = odr(startparameters, 3, 0.01, 1);
    startparameters[4] = odoptimiseparameter(startparameters, 4, exp(startparameters[3]), 1);
    return TRUE;
}

