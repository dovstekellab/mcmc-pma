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

#ifndef mcmc_oo_csvdata_h
#define mcmc_oo_csvdata_h

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>



class csvdata{
public:
    csvdata(std::string filename);
    std::vector<std::string> labels;
    bool cut(double threshold);
    int size();
    int length();
    std::vector<double> data(int well);
    std::vector<double> uncutdata(int well);
protected:
    std::vector<std::vector<double> > raw_data;
    int original_length;
    std::vector<int> cutoff;
};

#endif
