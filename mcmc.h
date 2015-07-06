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



#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif



#ifndef mcmc_oo_mcmc_h
#define mcmc_oo_mcmc_h

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
// #include <omp.h>
#include <sys/stat.h>
#include <unistd.h>


/*
 Class: mcmc
 This abstract class contains everything required to run MCMC for some observed data, except the function md() - this is so that we can define child classes for the different models we want to look at.
 */ 



class mcmc {
public:
    mcmc(std::vector<double> od, std::string l);
    int dim;
    virtual bool run_mcmc(int iterations);
    std::vector<double> mean();
    std::vector<double> variances();
    double bic();
    double dic();
    double aic();
    virtual std::vector<double> md(std::vector<double> parameters) =0;
    gsl_matrix * cov;
    std::vector<std::string> parameternames;
    std::string modelname;
    std::vector<double> lintolog(std::vector<double> parameters);
    std::vector<double> logtolin(std::vector<double> parameters);
    std::vector<double> logvartolinvar(std::vector<double> var, std::vector<double> par);
    std::vector<double> startparameters;
    std::vector<double> maxlikelihood;
    std::string label;
    void setlabel(std::string arg);
    virtual double odmax();
    double odlag();
    double ody0();
    
    double odymax();
    double odr(std::vector<double> params, const int index, double limitlower, double limitupper);
    double odoptimiseparameter(std::vector<double> params, const int index, double limitlower, double limitupper);
    std::string textoutput();
    std::stringstream logfile;
     double results_timetoendoflag(std::vector<double> parameters);
    double results_lag(std::vector<double> parameters);
     double results_maxrate(std::vector<double> parameters);
    double results_max(std::vector<double> parameters);
    virtual double results_lag();
    virtual double results_timetoendoflag();
    virtual double results_maxrate();
    virtual double results_max();
    virtual double results_timetoendoflag_var();
    virtual double results_maxrate_var();
    virtual double results_max_var();
    virtual double results_timetoendoflag_lq();
    virtual double results_maxrate_lq();
    virtual double results_max_lq();
    virtual double results_timetoendoflag_hq();
    virtual double results_maxrate_hq();
    virtual double results_max_hq();
    int acceptance_count = 0;
    int invalid_count = 0;
    std::vector<double> observed_data;
protected:
    virtual long double logprobability(std::vector<double> parameters);
    virtual long double logprobabilityfromdata(std::vector<double> data);
    std::vector<double> proposal(std::vector<double> parameters, gsl_matrix * cov);
    std::vector<std::vector<double> > *output;
    virtual bool checkparameters(std::vector<double> parameters) {return TRUE;};
    gsl_matrix * epsilon;
    std::vector<double> startvariances;
    std::vector<bool> logspace;
    std::vector<double> results, variance;
    std::vector<double> RVresults, RVvariance, RVlowquantile, RVhighquantile;
    gsl_rng * r; 
    virtual bool do_before_run();
    virtual bool do_after_run();
    double odvariance;
};


#endif
