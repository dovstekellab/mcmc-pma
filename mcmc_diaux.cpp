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
#include "mcmc_diaux.h"

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



//The diauxic model uses a simple ODE solver. These are the functions used by the solver.
__inline double ydiauxa ( double t, double y, double s1, double s2, double r1, double r2, double k1, double q0, double nu){
    return (q0/(q0+exp(-nu*t)))*((r1)*s1*y) + (((r2)*s2*k1*y)/((k1+s1)));
}
__inline double s1diauxa ( double t, double y, double s1, double s2, double r1, double r2, double k1, double q0, double nu){
    return -(q0/(q0+exp(-nu*t)))*(r1)*s1*y;
}
__inline double s2diauxa ( double t, double y, double s1, double s2, double r1, double r2, double k1, double q0, double nu){
    return -(q0/(q0+exp(-nu*t)))*((r2)*s2*k1*y)/((k1+s1));
}

mcmc_diaux::mcmc_diaux(std::vector<double> od, std::string l) : mcmc(od, l)


{
    dim = 8; 
    parameternames.resize(dim);
    parameternames[0] = "y0";
    parameternames[1] = "s1_0";
    parameternames[2] = "s2_0";
    parameternames[3] = "r1";
    parameternames[4] = "r2";
    parameternames[5] = "k1";
    parameternames[6] = "lag";
    parameternames[7] = "nu";
    modelname = "Diauxic";
    bool logsp[] = {0,1,1,1,1,1,0,1};
    logspace.assign(logsp, logsp+dim);
    double p[] = {0.155,2.3,2.6,0.0011,0.11,0.42,0.99,0.000028};
    double startv[] = {ody0(), ods1(), odymax()-ody0()-ods1(), 0.0002, 0.0003, 0.005, odlag(), 0.9};
    if (startv[2] <= 1){
        startv[2] = 1;
    }
    //double startv[] = {ody0(), 100, 200, 0.0002, 0.0003, 0.005, 60, 0.9};
    std::vector<double> startp(startv, startv+dim);
    startparameters = lintolog(startp);
    startvariances.assign(p, p+dim);
};

double mcmc_diaux::results_lagbyparams(std::vector<double> parameters){
    double t = 0;
    std::vector<double> p = logtolin(parameters);
    t = p[6];
    return t/4;
}

bool mcmc_diaux::do_before_run(){
    
        startparameters[3] = odr(startparameters, 3, 1E-8, 0.0001);
        startparameters[7] = odoptimiseparameter(startparameters, 7, 0.05, 1);
        startparameters[5] = odoptimiseparameter(startparameters, 5, 0.0001, 10);

    return TRUE;
}


/*
 Diauxic model:
 ODEs with two substrates - y, s1, s2.
 Parameters are:
 0 - y at t=0.
 1 - s1 at t=0.
 2 - s2 at t=0.
 3 - r1.
 4 - r2.
 5 - k1.
 6 - t0, end of flat phase.
 7 - q0.
 */

bool mcmc_diaux::checkparameters(std::vector<double> parameters){
    std::vector<double> p = logtolin(parameters);
    
    if (p[0] < 0) return FALSE;
    if (p[1] < 0) return FALSE;
    if (p[2] < 0) return FALSE;
    if (p[3] < 0 || p[3] > 1) return FALSE;
    if (p[4] < 0 || p[4] > 1) return FALSE;
    if (p[5] <= 1E-25 || p[5] > 20) return FALSE;
    if (p[6] < 0 || p[6] > observed_data.size()) return FALSE;
    if (ydiauxa(999, p[0], p[1], p[2], p[3], p[4], p[5], (double)1/(exp(p[7]*p[6])-1), p[7]) > 100) return FALSE;
    if (p[7] < 0.05 || p[7] > 10) return FALSE;
    
    
    return TRUE;
}


// This implements a very (!) simple fourth-order Runge-Kutte ODE solver.
std::vector<double> mcmc_diaux::md(std::vector<double> parameters){
    std::vector<double> model_data(observed_data.size());
    std::vector<double> p = logtolin(parameters);
    
    
    // Set step size (2 steps per time interval)
    int steps=2;
    
    // Initialise variables
    double y = p[0];
    double s1 = p[1];
    double s2 = p[2];
    
    double t = 0;
    
    double f11, f21, f31, f41, f12, f22, f32, f42, f13, f23, f33, f43;
    
    double dt;
    
    double q0 = (double)1/(exp(p[7]*p[6])-1);
    
    
    
    
    dt = (double)1 / steps;
    
    // Iteratively compute values of y, s1, s2 for each step.
    for (int i = 0; i<observed_data.size()*steps; i++){
        
        if (i%steps == 0){
            model_data[i/steps] = y; // Update model_data array at the same time.
        }
        
        
        
        t = (double)(i) / steps;
        
        
        f11 = ydiauxa ( t, y, s1, s2, p[3], p[4], p[5], q0, p[7]  );
        f12 = s1diauxa ( t, y, s1, s2, p[3], p[4], p[5], q0, p[7]  );
        f13 = s2diauxa ( t, y, s1, s2, p[3], p[4], p[5], q0, p[7]  );
        f21 = ydiauxa ( t + dt / 2, y + dt * f11 / 2 , s1 + dt * f12 / 2, s2 + dt * f13 / 2, p[3], p[4], p[5], q0, p[7] );
        f22 = s1diauxa ( t + dt / 2, y + dt * f11 / 2 , s1 + dt * f12 / 2, s2 + dt * f13 / 2, p[3], p[4], p[5], q0, p[7] );
        f23 = s2diauxa ( t + dt / 2, y + dt * f11 / 2 , s1 + dt * f12 / 2, s2 + dt * f13 / 2, p[3], p[4], p[5], q0, p[7] );
        f31 = ydiauxa ( t + dt / 2, y + dt * f21 / 2 , s1 + dt * f22 / 2, s2 + dt * f23 / 2, p[3], p[4], p[5], q0, p[7] );
        f32 = s1diauxa ( t + dt / 2, y + dt * f21 / 2 , s1 + dt * f22 / 2, s2 + dt * f23 / 2, p[3], p[4], p[5], q0, p[7] );
        f33 = s2diauxa ( t + dt / 2, y + dt * f21 / 2 , s1 + dt * f22 / 2, s2 + dt * f23 / 2, p[3], p[4], p[5], q0, p[7] );
        f41 = ydiauxa ( t + dt,     y + dt * f31, s1 + dt * f32 , s2 + dt * f33, p[3], p[4], p[5], q0, p[7] );
        f42 = s1diauxa ( t + dt,     y + dt * f31, s1 + dt * f32 , s2 + dt * f33, p[3], p[4], p[5], q0, p[7] );
        f43 = s2diauxa ( t + dt,     y + dt * f31, s1 + dt * f32 , s2 + dt * f33, p[3], p[4], p[5], q0, p[7] );
        
        y = y + dt * ( f11 + 2.0 * f21 + 2.0 * f31 + f41 ) / 6.0;
        s1 = s1 + dt * ( f12 + 2.0 * f22 + 2.0 * f32 + f42 ) / 6.0;
        s2 = s2 + dt * ( f13 + 2.0 * f23 + 2.0 * f33 + f43 ) / 6.0;
        
        // Due to rounding errors s1, s2 can become negative. Set them to zero if that happens.
        if (s1 <= 0) s1=0;
        if (s2 <= 0) s2=0;
    }
    return model_data;
}







// Some additional functions for initial parameter guesses specific to the diauxic model.


double mcmc_diaux::ods1(){
    
    logfile << "Looking for s_1 in OD\n";
    
    double v=0;
    for (int j=4; j<observed_data.size()-4; j++){
        v += pow(observed_data[j]-(observed_data[j-4]+observed_data[j-3]+observed_data[j-2]+observed_data[j-1]+observed_data[j]+observed_data[j+1]+observed_data[j+2]+observed_data[j+3]+observed_data[j+4])/9, 2);
    }
    v/=observed_data.size()-9;
    if (v<5) v=5;
    double sd = sqrt(v);
    
    
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
        gsl_fit_linear(numbers+i, 1, obsd+i, 1, 20, intercepts+i, slopes + i, &cov00, &cov01, &cov11, &sumsq);
    }
    
    int maxslope = 0;
    for (int i=0; i<observed_data.size()-21; i++){
        if (slopes[i] > slopes[maxslope] && slopes[i < 10]) maxslope = i;
    }
    
    

    int min = 0;
    for (int i=0; i<observed_data.size()-21; i++){
        if(obsd[i+11] > ody0() + 9*sd && obsd[i+11] < odymax() - 12*sd) min=i;
    }
    for (int i=0; i<observed_data.size()-21; i++){
        if(slopes[i] < slopes[min] && obsd[i+11] > ody0() + 12*sd && obsd[i+11] < odymax() - 9*sd) min=i;
    }
    logfile << "minimum in between at: " << min << "\n";
    min += 11;
    logfile << "minimum in between at: " << min << "\n";
    
    double s1 = (observed_data[min-2] + observed_data[min-1] + observed_data[min] + observed_data[min+1] + observed_data[min+2])/5 - ody0();
    
    logfile << "s1 is: " << s1 << "\n";
    
    if (s1 < (odymax()-ody0())/8) s1 = (odymax()-ody0())/4;
    
    if (s1 <= 4) s1 = 4;
    return s1;

    
    
    
    
    
}

double mcmc_diaux::odminslope(){
    
    double v=0;
    for (int j=4; j<observed_data.size()-4; j++){
        v += pow(observed_data[j]-(observed_data[j-4]+observed_data[j-3]+observed_data[j-2]+observed_data[j-1]+observed_data[j]+observed_data[j+1]+observed_data[j+2]+observed_data[j+3]+observed_data[j+4])/9, 2);
    }
    v/=observed_data.size()-9;
    if (v<5) v=5;
    double sd = sqrt(v);
    
    
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
    
    
    
    int min = 0;
    for (int i=0; i<observed_data.size()-21; i++){
        if(obsd[i+11] > ody0() + 9*sd && obsd[i+11] < odymax() - 12*sd) min=i;
    }
    for (int i=0; i<observed_data.size()-21; i++){
        if(slopes[i] < slopes[min] && obsd[i+11] > ody0() + 12*sd && obsd[i+11] < odymax() - 9*sd) min=i;
    }
    min += 11;

    return min;
    
    
    
    
    
    
}

double mcmc_diaux::odlag(){

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
    for (int i=0; i<odminslope()-11; i++){
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





