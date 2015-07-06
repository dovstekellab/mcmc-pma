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
#include <gsl/gsl_fit.h>
#include <math.h>
#include <time.h>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

#define T0 2500
#define WRITE_TO_DISK



extern std::string ID; // Need the global ID to know directory name if we write to disk.




// Basic setup of some cass variables etc. Usually some additional setup will need to be performed in each specific mcmc_model class.

mcmc::mcmc(std::vector<double> od, std::string l)
{
    observed_data = od;
    label = l;
    //Set up random number generator:
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    time_t seed = time(NULL);
    gsl_rng_set(r, seed);
    
    //Set up OD variance:
    double v=0;
    for (int i=4; i<od.size()-4; i++){
        v += pow(od[i]-(od[i-4]+od[i-3]+od[i-2]+od[i-1]+od[i]+od[i+1]+od[i+2]+od[i+3]+od[i+4])/9, 2);
    }
    v/=od.size()-9;
    logfile << "OD Variance: " << v << "\n";
    if (v<5) v=5;
    odvariance = v;
}








// Run the adaptive Metropolis-Hastings algorithm.

bool mcmc::run_mcmc(int iterations){
    
    
    //If this macro is defined above, we write the entire output of the chain to disk. For now, we create the directory and open a file.
#ifdef WRITE_TO_DISK
    std::string basedir = ID + "/wells/" + label;
    mkdir("output", 0750);
    mkdir(basedir.c_str(), 0750);
    
    std::stringstream timestamp;
    timestamp << (int*)(this);
    
    std::ofstream file((basedir + "/debug_" + label + "-" + modelname + "-" + timestamp.str() + ".csv").c_str());
#endif
    
    time_t starttime = time(NULL);
    
    do_before_run();    //Some models might want to perform some additional setup before running the algorithm.
    
    //Initialise variables:
    
    output = new std::vector<std::vector<double> > (1, std::vector<double>(dim)); //The states (parameters) of the MCMC chain.
    std::vector<std::vector<double> > * additionalRVs = new std::vector<std::vector<double> > (1, std::vector<double>(4));  //Some additional random variables whose posterior distribution we want to sample.
    
    //Set up the logprobability of the initial parameters.
    long double logp_current_parameters = 0;
    logp_current_parameters = logprobability(startparameters);
    long double logp_candidate = 0;
    
    //Various counters, mostly for diagnostic purposes.
    acceptance_count = 0;
    invalid_count = 0;
    int acceptance_count_2 = 0;
    int invalid_count_2 = 0;
    int iteration_count = 0;
    int nan_count = 0;
    int inf_count = 0;
    int am_index=0;
    bool accepted = 0;
    
    //Vectors for the candidate parameters, modelled data based on these.
    std::vector<double> candidate(dim);
    std::vector<double> data(observed_data.size());
    
    //Set up the first state of the chain.
    for (int d=0; d<dim; d++) {
        (*output)[0][d] = startparameters[d];
    }
    
    //Variable for the parameters with the maximum likelihood.
    maxlikelihood = startparameters;
    double logprobability_maxlikelihood = logprobability(maxlikelihood);
    
    
    //Variables for the adaptive covariance parts:
    //Set initial covariance:
    gsl_matrix * cov0 = gsl_matrix_alloc(dim, dim);
    gsl_matrix_set_zero(cov0);
    for (int j=0; j<dim; j++){
        gsl_matrix_set(cov0, j, j, startvariances[j]);

    }
    
    //Initialise iterative covariance and auxiliarry matrices:
    cov = gsl_matrix_alloc(dim, dim);
    
    gsl_matrix * tempdd = gsl_matrix_alloc(dim, dim);
    gsl_matrix * tempd1 = gsl_matrix_alloc(dim, 1);
    
    gsl_matrix * meanvec = gsl_matrix_alloc(dim, 1);
    gsl_matrix_set_zero(meanvec);
    gsl_matrix_set_zero(cov);
    for (int j=0; j<dim; j++){
        gsl_matrix_set(meanvec, j, 0, startparameters[j]);

    }
    
    epsilon = gsl_matrix_alloc(dim, dim); 
    gsl_matrix_set_zero(epsilon);
    for(int d=0; d<dim; d++){ gsl_matrix_set(epsilon, d, d, 0.0000000001); }
    
    
    gsl_matrix * covariance = gsl_matrix_alloc(dim,dim);
    gsl_matrix_set_zero(covariance);

    
    //Scaling parameter for proposal density:
    double lambda = 5.6644/dim;
    
    
    // Done with initialising various variables.
    
    
    
    
    
    //Some output to log:
    logfile << "Starting AM algorithm with " << iterations << " iterations for well " << label << " and model " << modelname << ".\n";
    
    logfile << "Initial covariance:\n";
    for (int d=0; d<dim; d++){
        for (int e=0; e<dim; e++){
            logfile << gsl_matrix_get(cov0, d, e) << "; ";
        }
        logfile << "\n";
    }
    logfile << "\n";
    
    
    
    
    
    
    //Main M-H loop.
    for (int i=1; i < iterations; i++) {

        //We keep track of whether we jumped in this iteration or not, so we know if we should update the covariance.
        accepted = 0;
        
        //Get covariance for this iteration.
        if (i<T0){
            gsl_matrix_memcpy(covariance, cov0); //Fixed for initial T0 iterations
        } else {
            gsl_matrix_memcpy(covariance, cov); //Adaptive thereafter.
        } 
        
        gsl_matrix_scale(covariance, lambda);
        
        //Generate candidate vector.
        candidate = proposal((*output)[i-1], covariance);
        

        if (checkparameters(candidate)){    //We work purely with flat priors on simple regions (usually boxes) of the parameter space. checkparameters() returns true if a parameter vector is inside this region, false otherwise.
            
            
            
            
            // Calculate acceptance probability
            data = md(candidate);
            logp_candidate = logprobabilityfromdata(data);
            long double p;
            if (isinf(logp_candidate - logp_current_parameters)){
                p = 0;
            } else {
                p = exp(logp_candidate - logp_current_parameters);
            }
            if (p>1) p=1;
            
            // Debug counters:
            if (isnan(p)) {nan_count++; p=0;}
            if (isinf(p)) {inf_count++; p=1;}
            
            
            
            
            
            // Get uniform random value
            double u = gsl_ran_flat(r, 0, 1);
            
            
            //To jump, or not to jump:
            if (u <= p) {
                //If we do jump:
                
                //Accept the current candidate parameters.
                output->push_back(candidate);
                
                //Update logprobability.
                logp_current_parameters = logp_candidate;
                
                //Debug counters
                acceptance_count++;
                acceptance_count_2++;
                
                //Calculate any additional random variables we want to sample.
                std::vector<double> RV(5);
                RV[0] = 0;
                RV[1] = results_lag(data);
                RV[2] = results_maxrate(data);
                RV[3] = results_max(data);
                RV[4] = (float)-2 * (logp_candidate);
                additionalRVs->push_back(RV);
                
                //Check if this is the highest likelihood we have seen so far.
                if(logp_candidate > logprobability_maxlikelihood){
                    maxlikelihood = candidate;
                    logprobability_maxlikelihood = logp_current_parameters;
                    
                    
                accepted = 1;
                }
            } else {
                
                //If we don't accept, stay at the current state.
                output->push_back(output->back());
                additionalRVs->push_back(additionalRVs->back());
                
            }
            
            //Update scaling factor:
            if (i>0){
                lambda = exp( log(lambda) + 1.0/(i) * (p - 0.234) );
            }
        } else {
            //If we're outside our prior, stay at the current state.
            output->push_back(output->back());
            additionalRVs->push_back(additionalRVs->back());
            
            //Debug counters.
            invalid_count++;
            invalid_count_2++;
            
            //Update scaling factor.
            if (i>0){
                lambda = exp( log(lambda) + 1.0/(i) * (0 - 0.234) );
            }
        }
        iteration_count++;
        
        
        //If this macro is set, we write the current state to disk.
#ifdef WRITE_TO_DISK
        for (int k=0; k<(output->back()).size(); k++){
            file << output->back()[k] << ",";
        }
        file << "\n";
#endif
        
        //Adapt covariance.
        //For the first T0 iterations, we only update the covariance if we've jumped. After T0 iterations, always.
        if (accepted == 1 || i > T0){
            //Generate new covariance:
            gsl_matrix_view gsl_parameters = gsl_matrix_view_array(&((*output)[i][0]), dim, 1);
            
            
            gsl_matrix_set_zero(tempdd);
            gsl_matrix_memcpy(tempd1, meanvec);
            gsl_matrix_sub(tempd1, &gsl_parameters.matrix);
            gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, tempd1, tempd1, 0, tempdd);
            
            gsl_matrix_sub(tempdd, cov);
            gsl_matrix_scale(tempdd, (double)1/(am_index+1));
            gsl_matrix_add(cov, tempdd);
            
            gsl_matrix_add(cov, epsilon);

                
            //Update mean:
            gsl_matrix_scale(meanvec, am_index);
            gsl_matrix_add(meanvec, &gsl_parameters.matrix);
            gsl_matrix_scale(meanvec, (double)1/(am_index+1));
            am_index++;
        }
        
        // We reset the scaling factor when we start using the adapted covariance.
        if (i==T0){
            lambda = 5.6644/dim;
        }
        
        
        
        //Some more output to logfile every 10k iterations:

        if (i%10000==0 || i == 1 || i == iterations -1){


            logfile << "Iterative Covariance matrix at iteration " << i << ":\n";
            for (int d=0; d<dim; d++){
                for (int e=0; e<dim; e++){
                    logfile << std::setw(6) << std::setfill(' ') << gsl_matrix_get(cov, d, e) << "; ";
                }
                logfile << "\n";
            }
            logfile << "\n";
            logfile << "Iterative mean: ";
            for (int d=0; d<dim; d++){
                logfile << parameternames[d] << ": " << gsl_matrix_get(meanvec, d, 0) << "; ";
            }

            
            logfile << "\nAcceptance ratio: " << (double)acceptance_count_2/10000 << "\n\n";
            acceptance_count_2=0;
            
            logfile << "Scaling factor: " << lambda << "\n";
            
            logfile << "inf count: " << (double)inf_count/10000 << "; nan: " << (double)nan_count/10000 << "; invalid: " << (double)invalid_count_2/10000 << ";\n";
            inf_count=0;
            nan_count=0;
            invalid_count_2=0;
            
            logfile << "\n\n";
            
        }
        
        
        
    }
    
    // Done with the main MCMC loop!
    
    
    
    
    //We discard the first half of the chain as burn-in. Also, we need the output as a C array for GSL, so we simply copy the second half of the output vector to one.
    int j = ceil(iterations/2);
    int k = floor(iterations/2);
    double * array = new double[k];
    std::vector<double> mean(dim);
    std::vector<double> var(dim);
    
    for (int d=0; d<dim; d++) {
        for (int i=0; i<k; i++) {
            array[i] = (*output)[i+j][d];
        }
        
        // Calculate mean and variance for each parameter.
        mean[d] = gsl_stats_mean(array, 1, k);
        var[d] = gsl_stats_variance_m(array, 1, k, mean[d]);
        
    }
    
    results=mean;
    variance=var;
    
    //Calculate mean, variance and quantiles for the additional random variables we sampled.
    
    RVresults.resize(5);
    RVvariance.resize(5);
    RVlowquantile.resize(5);
    RVhighquantile.resize(5);
    
    for (int d=0; d<5; d++) {
        for (int i=0; i<k; i++) {
            array[i] = (*additionalRVs)[i+j][d];
        }
        
        
        RVresults[d] = gsl_stats_mean(array, 1, k);
        RVvariance[d] = gsl_stats_variance_m(array, 1, k, RVresults[d]);
        gsl_sort(array, 1, k);
        RVlowquantile[d] = gsl_stats_quantile_from_sorted_data(array, 1, k, 0.025);
        RVhighquantile[d] = gsl_stats_quantile_from_sorted_data(array, 1, k, 0.975);
        
    }
    
    
    
    
    //A little more diagnostic output to the log.
    
    logfile << "Final covariance:\n";
    for (int d=0; d<dim; d++){
        for (int e=0; e<dim; e++){
            logfile << gsl_matrix_get(cov, d, e) << "; ";
        }
        logfile << "\n";
    }
    logfile << "\n";
    
    logfile << "Final overall acceptance ratio: " << (double)acceptance_count/iterations << "\n";
    logfile << "Total number of acceptances: " << acceptance_count << "\n";
    logfile << "Total number of invalids / outside prior support: " << invalid_count << "\n";
    
    
    
    // In case a particular model wants to do things after the chain has finished:
    do_after_run();
    
    
    
    // Free some memory.
    delete output;
    delete array;
    delete additionalRVs;
    
    gsl_matrix_free(cov0);
    gsl_matrix_free(tempdd);
    gsl_matrix_free(tempd1);
    gsl_matrix_free(meanvec);
    
    
    // And some final output to the logfile.
    time_t stoptime = time(NULL);
    
    int runtime = stoptime - starttime;
    
    logfile << "Runtime in seconds: " << runtime << " \n";
    
    
    logfile << "\n" << " \n";
    
    logfile << "DIC: " << dic() << " \n";
    
    
    logfile << "D_bar: " << RVresults[4] << " \n";
    
    
    logfile << "logl(theta_bar) " << logprobability(this->mean()) << " \n";
    logfile << "p_D_notused: " << RVresults[4] + 2*logprobability(this->mean()) << " \n";
    
    logfile << "p_D: " << 0.5 * RVvariance[4] << " \n";

 
#ifdef WRITE_TO_DISK
    file.close();
#endif
    
    return 0;
}

bool mcmc::do_after_run(){
    return TRUE;
}

bool mcmc::do_before_run(){
    return TRUE;
}


// Generate a new candidate vector based on current state and a normal distribution.
__inline std::vector<double> mcmc::proposal(std::vector<double> parameters, gsl_matrix * cov){
    
    std::vector<double> candidate(parameters.size());
    std::vector<double> norm(parameters.size());
    
    for (int i=0; i<parameters.size(); i++){ norm[i] = gsl_ran_gaussian(r, 1); }    //Sample i-dimensional normal.
    
    gsl_matrix * covariance = gsl_matrix_alloc(parameters.size(), parameters.size());
    gsl_matrix_memcpy(covariance, cov);
    int status;
    gsl_set_error_handler_off();
    
    status = gsl_linalg_cholesky_decomp(covariance);
    
    
    //In case the covariance is singular, add epsilon * I until it is not.

    if(status){
        logfile << "Covariance singular!\n";
    }
    
        while (status){
        gsl_matrix_add(covariance, epsilon);
        status = gsl_linalg_cholesky_decomp(covariance);
    }
    
    
    //Add the scaled normal sample to the candidate vector.
    
    for (int i=0; i<parameters.size(); i++){
        candidate[i] = parameters[i];
        for (int j=0; j<=i; j++){ candidate[i] += (gsl_matrix_get(covariance, i, j)) * norm[j]; }
    }
    
    gsl_matrix_free(covariance);
    return candidate;
}


// Calculate logprobability from parameter or data vector.
__inline long double mcmc::logprobability(std::vector<double> parameters){
    std::vector<double> model_data = md(parameters);    
    return logprobabilityfromdata(model_data);
}

__inline long double mcmc::logprobabilityfromdata(std::vector<double> data){
    
    long double prob = 0;
    for (int t = 0; t<observed_data.size(); t++){
        prob += (-pow((data[t]-observed_data[t]),2)/(2*odvariance)); 
    }
    return prob;
}




// Helper functions to switch between parameters in linear and logarithmic scales.

std::vector<double> mcmc::lintolog(std::vector<double> parameters){
    for (int i=0; i<dim; i++) {
        if (logspace[i]) parameters[i] = log(parameters[i]);
    }
    return parameters;
}

std::vector<double> mcmc::logtolin(std::vector<double> parameters){
    for (int i=0; i<dim; i++) {
        if (logspace[i]) parameters[i] = exp(parameters[i]);
    }
    return parameters;
}


std::vector<double> mcmc::logvartolinvar(std::vector<double> var, std::vector<double> par){
    std::vector<double> linvars(var.size());
    for (int i=0; i<linvars.size(); i++){
        if (logspace[i]){
            linvars[i] = (exp(pow(var[i],2)) - 1) * exp( 2*par[i] + pow(var[i],2) );
        } else {
            linvars[i] = var[i];
        }
    }
    return linvars;
}






//Various functions to access the results (posterior mean of lag, max, and max rate) of the MCMC run. NB, we do not check that we have actually run the MCMC chain. If we haven't these all return garbage or segfault.

std::vector<double> mcmc::mean(){
    
    return results;
}

std::vector<double> mcmc::variances(){
    
    return variance;
}




double mcmc::results_timetoendoflag(std::vector<double> parameters){
    std::vector<double> data = md(parameters);

    
    return results_lag(data);
    
    
}



double mcmc::results_timetoendoflag(){
    //return results_timetoendoflag(mean());
    return RVresults[1];
}

double mcmc::results_lag(){
    //return results_timetoendoflag(mean());
    return RVresults[1];
}

double mcmc::results_max(){
    return RVresults[3];
}

double mcmc::results_maxrate(){
    return RVresults[2];
}

double mcmc::results_timetoendoflag_var(){
    return RVvariance[1];
}

double mcmc::results_max_var(){
    return RVvariance[3];
}

double mcmc::results_maxrate_var(){
    return RVvariance[2];
}

double mcmc::results_timetoendoflag_lq(){
    return RVlowquantile[1];
}

double mcmc::results_max_lq(){
    return RVlowquantile[3];
}

double mcmc::results_maxrate_lq(){
    return RVlowquantile[2];
}

double mcmc::results_timetoendoflag_hq(){
    return RVhighquantile[1];
}

double mcmc::results_max_hq(){
    return RVhighquantile[3];
}

double mcmc::results_maxrate_hq(){
    return RVhighquantile[2];
}




// Functions to compute lag, maximum and maximum rate from modelled data. These are used to sample additional random variables in the MCMC loop, and to compute the posterior mean of these random variables.

double mcmc::results_lag(std::vector<double> data){
    
    
    double intercepts[data.size()-1];
    double slopes[data.size()-1];
    
    for (int i=0; i<data.size()-1; i++){
        slopes[i] = data[i+1] - data[i];
        intercepts[i] = data[i] - i * slopes [i];
    }
    
    int maxslope = 0;
    for (int i=0; i<data.size()-1; i++){
        if (slopes[i] > slopes[maxslope]) maxslope = i;
    }
    
    double y0 = data[0];
    
    double lag = (y0 - intercepts[maxslope]) / slopes[maxslope];
    
    
    
    if (lag <= 0) lag = 0;
    if (lag >= data.size()) lag = data.size();
    if (isnan(lag) || isinf(lag)) lag = 0;
    
    return lag;     
    
}


double mcmc::results_max(std::vector<double> data){
    std::vector<double> q = data;
    double m = q[0];
    for (int i = 0; i<q.size(); i++){
        if (q[i] > m) m=q[i];
    }
    return m-q[0];
}

double mcmc::results_maxrate(std::vector<double> data){
    std::vector<double> q = data;
    std::vector<double> qdash = q;
    qdash[0] = 0;
    for (int i = 1; i<qdash.size(); i++){
        qdash[i] = q[i] - q[i-1];
    }
    double m = qdash[0];
    for (int i = 0; i<qdash.size(); i++){
        if (qdash[i] > m) m=qdash[i];
    }
    return m;
}





// Functions to compute some characteristics of the raw data. Useful for guessing initial parameter vectors to start the chain with, may be used to mcmc_model classes for this purpose.


// Maximum in the raw data.
double mcmc::odmax(){
    double m = 0;
    for (int i=0; i<observed_data.size(); i++){
        if (observed_data[i]>=m) m=observed_data[i];
    }
    return m;
}

// Maximum of smoothed raw data.
double mcmc::odymax(){
    int s=(int)observed_data.size()-10;
    double ymax = (observed_data[s+0]+observed_data[s+1]+observed_data[s+2]+observed_data[s+3]+observed_data[s+4] + observed_data[s+5]+observed_data[s+6]+observed_data[s+7]+observed_data[s+8]+observed_data[s+9])/10;
    if (ymax < ody0()+1) ymax = ody0()+1;
    return ymax;
}


// Best guess for the lag, based on smoothed raw data.
double mcmc::odlag(){
    double* numbers = new double[observed_data.size()];
    double* obsd = new double[observed_data.size()];
    for (int i=0; i<observed_data.size(); i++){
        numbers[i] = i;
        obsd[i] = observed_data[i];
    }
    
    double intercepts[observed_data.size()-21];
    double slopes[observed_data.size()-21];
    
    double cov00, cov01, cov11, sumsq;
    
    for (int i=0; i<observed_data.size()-21; i++){
        gsl_fit_linear(numbers+i, 1, obsd+i, 1, 20, intercepts+i, slopes+i, &cov00, &cov01, &cov11, &sumsq);
    }
    
    int maxslope = 0;
    for (int i=0; i<observed_data.size()-21; i++){
        if (slopes[i] > slopes[maxslope] && slopes[i < 20]) maxslope = i;
    }
    
    logfile << "Max slope at index: " << maxslope+10 << "; intercept = " << intercepts[maxslope]<< "; slope = " << slopes[maxslope] << "\n";
    
    double y0 = ody0();
    
    double lag = (y0 - intercepts[maxslope]) / slopes[maxslope];
    
    logfile << "Assumed y0= " << y0 << "; lag = " << lag << "; \n";
    
    delete[] numbers; delete[] obsd;
    
    if (isnan(lag) || isinf(lag) || lag <= 0 || lag >= observed_data.size()) lag = 1;
    
    return lag;

}


// Minimum in smoothed raw data.
double mcmc::ody0(){
    
    double y0 = (observed_data[0]+observed_data[1]+observed_data[2]+observed_data[3]+observed_data[4] + observed_data[5]+observed_data[6]+observed_data[7]+observed_data[8]+observed_data[9])/10;
    return y0;
}

// Two heuristics to guess rate parameter and any arbitrary parameter. These vary the specified parameter while keeping the others fixed. Thus, the order in which these are called to guess different parameters matters.

// Heuristic best guess for rate parameter. Very simple search that adjusts parameter until the maximum growth rate of the modelled data matches the maximum rate of the smoothed raw data.
double mcmc::odr(std::vector<double> params, const int index, double limitlower, double limitupper){
    double* numbers = new double[observed_data.size()];
    double* obsd = new double[observed_data.size()];
    for (int i=0; i<observed_data.size(); i++){
        numbers[i] = i;
        obsd[i] = observed_data[i];
    }
    
    double intercepts[observed_data.size()-21];
    double slopes[observed_data.size()-21];
    
    double cov00, cov01, cov11, sumsq;
    
    for (int i=0; i<observed_data.size()-21; i++){
        gsl_fit_linear(numbers+i, 1, obsd+i, 1, 20, intercepts+i, slopes+i, &cov00, &cov01, &cov11, &sumsq);
    }
    
    int maxslope = 0;
    for (int i=0; i<observed_data.size()-21; i++){
        if (slopes[i] > slopes[maxslope] && slopes[i < 10]) maxslope = i;
    }
    
    
    std::vector<double> deltaslopes(10);
    if (logspace[index] == 1) {limitlower = log(limitlower); limitupper = log(limitupper);}
    std::vector<double> parameters = params;
    for (int i=0; i<100; i++){
        double deltar = (limitupper-limitlower)/10;
        for (int j=0; j<10; j++){
            parameters[index] = limitlower + j*deltar;
            if (checkparameters(parameters)){
                std::vector<double> data = md(parameters);
                double modelmaxslope = 0;
                for (int k=0; k<observed_data.size()-1; k++){
                    if (data[k+1]-data[k] > modelmaxslope) modelmaxslope = data[k+1]-data[k];
                }
                deltaslopes[j] = modelmaxslope - slopes[maxslope];
                if (deltaslopes[j] < 0) deltaslopes[j] *= -1;
            } else{
                deltaslopes[j] = 999999;
            }
        }
        int m = 0;        
        for (int j=0; j<10; j++){
            if (deltaslopes[j] < deltaslopes[m]) m=j;
        }
        limitlower = limitlower + (m-1)*deltar;
        limitupper = limitlower + (m+1)*deltar;
    }
    
    delete numbers; delete obsd;
    if (logspace[index] == 1) logfile << "Optimisation for rate-parameter " << index << "; Optimal value = " << exp((limitlower+limitupper)/2) << "; \n";
    else logfile << "Optimisation for rate-parameter " << index << "; Optimal value = " << (limitlower+limitupper)/2 << "; \n";
    
    return (limitlower+limitupper)/2;
}


// Generic heuristic to guess any parameter. Adjusts parameter until logprobability is highest.
double mcmc::odoptimiseparameter(std::vector<double> params, const int index, double limitlower, double limitupper){
    std::vector<double> logprobabilities(10);
    if (logspace[index] == 1) {limitlower = log(limitlower); limitupper = log(limitupper);}
    std::vector<double> parameters = params;
    for (int i=0; i<100; i++){
        double deltar = (limitupper-limitlower)/10;
        for (int j=0; j<10; j++){
            parameters[index] = limitlower + j*deltar;
            if (checkparameters(parameters)){
                logprobabilities[j] = logprobability(parameters);
            } else {
                logprobabilities[j] = -1E306;
            }
        }
        int m = 0;        
        for (int j=0; j<10; j++){
            if (logprobabilities[j] > logprobabilities[m]) m=j;
        }
        if (m!=0) limitlower = limitlower + (m-1)*deltar;
        if (m!=10) limitupper = limitlower + (m+1)*deltar;
    }
    if (logspace[index] == 1) logfile << "Optimisation for parameter " << index << "; Optimal value = " << exp((limitlower+limitupper)/2) << "; \n";
    else logfile << "Optimisation for parameter " << index << "; Optimal value = " << (limitlower+limitupper)/2 << "; \n";
    return (limitlower+limitupper)/2;
}





// Functions to calculate DIC, BIC and AICc based on the sampled posterior. Like other functions based on the results of the MCMC chain, we do not check that we've actually run the chain - if we haven't these may behave unexpectedly.


double mcmc::bic(){
    double b = -2 * (log((1/(sqrt(odvariance*2*3.1415926))))+logprobability(maxlikelihood)) + dim * log(observed_data.size());
    if (isnan(b)) b=999999999;
    return b;
    //return dic();
}

double mcmc::dic(){

    double dic = RVresults[4] + 0.5 * RVvariance[4];

    if (acceptance_count > -10 && isnormal(dic)){
        return  dic;
    }
    else return 999999;
}

double mcmc::aic(){
    return 2 * dim - 2 * (log((1/(sqrt(odvariance*2*3.1415926))))+logprobability(maxlikelihood)) + (2*dim*(dim+1))/(observed_data.size() - dim - 1);
}






// Some generic diagnostic output in text form.

std::string mcmc::textoutput(){
    std::stringstream stream;
    stream << "Results for well " <<  label << ": \n ";
    std::vector<double> means = logtolin(results);
    stream << "Accepted: " << acceptance_count << "\n";
    stream << "lag= " << results_timetoendoflag() << "; maxrate= " << results_maxrate() << "; max= " << results_max() << ";";
    stream << "\nDIC score: " << dic() << "; BIC score: " << bic() << ";\n";
    stream << modelname << ":\nEstimated Parameters: ";
    for (int k=0; k<means.size(); k++){
        if(k%2==0) {stream << "\n";}
        stream << parameternames[k] << " = " << means[k] << "(s.d. = " << variance[k] << "), ";
    }
    
    stream << "\nInitial Parameters: ";
    for (int k=0; k<means.size(); k++){
        if(k%3==0) {stream << "\n";}
        stream << parameternames[k] << " = " << logtolin(startparameters)[k] << "; ";
    }
    double auc=0;
    for (int i=0; i<observed_data.size(); i++){
        auc += observed_data[i];
    }
    auc /= observed_data.size();
    stream << "\n od AUC:" << auc << "; " << "; odmax: " << odmax() << "; od size: " << observed_data.size() << ";\n";
    
    std::vector<double> data = md(mean());
    
    double intercepts[data.size()-1];
    double slopes[data.size()-1];
    
    //double cov00, cov01, cov11, sumsq;
    
    stream << "[slopes, intercepts] = ";
    
    for (int i=0; i<data.size()-1; i++){
        //gsl_fit_linear(numbers+i, 1, obsd+i, 1, 20, intercepts+i, slopes+i, &cov00, &cov01, &cov11, &sumsq);
        slopes[i] = data[i+1] - data[i];
        intercepts[i] = data[i] - i * slopes [i];
        stream << " [" << slopes[i] << ", " << intercepts[i] << "];   ";
    }
    
    int maxslope = 0;
    for (int i=0; i<data.size()-1; i++){
        if (slopes[i] > slopes[maxslope]) maxslope = i;
    }
    stream << " \n" << "maxslope = " << maxslope << "\n";
    
    double y0 = data[0];
    
    double lag = (y0 - intercepts[maxslope]) / slopes[maxslope];
    
    stream << "y0 = " << y0 << "\n";
    stream << "intercepts[maxslope] = " << intercepts[maxslope] << "\n";
    stream << "slopes[maxslope] = " << slopes[maxslope] << "\n";
    
    // delete[] numbers; delete[] obsd;
    

    
    
    return stream.str();
}
