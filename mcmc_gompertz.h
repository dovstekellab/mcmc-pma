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

#ifndef mcmc_mcmc_gompertz_h
#define mcmc_mcmc_gompertz_h
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
#include <math.h>
#include <time.h>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>







/*
 Class: mcmc_gomp: Specialised class for Gompertz model.
 */
class mcmc_gomp: public mcmc{
public:
    mcmc_gomp(std::vector<double> od, std::string l);
    std::vector<double> md(std::vector<double> parameters);
    double results_lagbyparams(std::vector<double> parameters);
protected:
    bool checkparameters(std::vector<double> parameters);
    virtual bool do_before_run();
};

#endif