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






#ifndef ITERATIONS
#define ITERATIONS (500000)
#endif


#define VERBOSITY 1

#define _IC_ dic





#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

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

//Uncomment if you want to use OpenMP
//#include <omp.h>


#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>

#include "csvdata.h"





#include "mcmc.h"

#include "mcmc_baranyi.h"
#include "mcmc_simplebaranyi.h"
#include "mcmc_simplestbaranyi.h"
#include "mcmc_simplebaranyi2.h"

#include "mcmc_diaux.h"

#include "mcmc_nogrowth.h"

#include "mcmc_logistic.h"
#include "mcmc_gompertz.h"




std::string itoa (int i);


std::string itoa( int i ){
    
    std::stringstream s;
    s << std::setw(3) << std::setfill('0') << i;
    return s.str();
    
}


std::string ID;



int main (int argc, const char * argv[])
{
    
    nice(20);
    
    std::cout << "Build: " << __DATE__ << " " << __TIME__ << "\n";

    
    
    //Get filename from CLI parameter:
    std::string fn = " ";
    
    if (argc >=2){
        char temp[256];
        for(int i = 0; (argv[1][i]) != '\0'; i++){
            temp[i] = (argv[1][i]);
            temp[i+1] = '\0';
        }
        fn = temp;
    }
    
    
    //Read data from CSV file into vector:
    csvdata csv(fn);
    
    
    
    
    
    std::stringstream timestamp;
    timestamp << time(NULL);
    std::string filename = fn;
    std::replace(filename.begin(), filename.end(), '/', '-');
    
    std::string basedir = "output/" + filename + "-" + timestamp.str();
    mkdir("output", 0750);
    mkdir("output-archives", 0750);
    mkdir(basedir.c_str(), 0750);
    mkdir((basedir + "/wells").c_str(), 0750);
    
    ID = basedir;
    
    
    //Remove aberrant data at the end of each well. Remove if this is not required for your data.
    csv.cut(0.95);

    

    
    //This gives us a vector of vectors of pointers to be filled with MCMC objects
    std::vector<std::vector<mcmc *> > MCMCs;
    
    
    
    //Now we populate this with our data:
    for (int i=0; i < csv.size(); i++) {
        
        std::vector<mcmc *> temp;
        
        temp.push_back(new mcmc_baranyi(csv.data(i), itoa(i) + "-" + csv.labels[i]));
        temp.push_back(new mcmc_simplestbaranyi(csv.data(i), itoa(i) + "-" + csv.labels[i]));
        temp.push_back(new mcmc_simplebaranyi(csv.data(i), itoa(i) + "-" + csv.labels[i]));
        temp.push_back(new mcmc_simplebaranyi2(csv.data(i), itoa(i) + "-" + csv.labels[i]));
        
        temp.push_back(new mcmc_diaux(csv.data(i), itoa(i) + "-" + csv.labels[i]));
        
        temp.push_back(new mcmc_nogrowth(csv.data(i), itoa(i) + "-" + csv.labels[i]));
        
        //temp.push_back(new mcmc_logistic(csv.data(i), itoa(i) + "-" + csv.labels[i]));
        temp.push_back(new mcmc_gomp(csv.data(i), itoa(i) + "-" + csv.labels[i]));
        
        MCMCs.push_back(temp);
        
    }
    
    
    //Run MCMC for all data:
    
    int m = (int)MCMCs.size();
    int n = (int)MCMCs[0].size();
    
    std::vector<int> finished(m*n);
    
    for (int i=0; i < m*n; i++) {finished[i]=0;}
    int count = 0;
    
    
    std::cout << "Working: " << count << " out of " << m*n << " done. \r" << std::flush;
    
    //Uncomment if you want to use OpenMP.
    // omp_set_num_threads(8);
    
    

    
#pragma omp parallel for schedule(dynamic, 1)   //This should be ignored by non-OpenMP compilers. If not - remove.
    for (int i=0; i < m*n; i++) {
        
        int j = floor(i/n);
        int k = floor(i%n);
        MCMCs[j][k]->run_mcmc(ITERATIONS);
        
        finished[i]=1;
#pragma omp critical
        
        count = 0;
        for (int l=0; l <m*n; l++) {count += finished[l];}
        
        std::cout << "Working: " << count << " out of " << m*n << " done. \r" << std::flush;
        
    }
    
    
    std::cout << "\n";
    
    
    
    // Get an interval of maximum coloration attained at control well to determine growth vs no growth:
    int controlbestmodel = 0;
    for (int j=0; j<MCMCs[2].size(); j++) {
        if (MCMCs[2][j]->_IC_() < MCMCs[2][controlbestmodel]->_IC_() || isnan(MCMCs[2][controlbestmodel]->_IC_())) controlbestmodel=j;
    }
    double controlwell_max_highquantile = MCMCs[2][controlbestmodel]-> results_max_hq();
    
    

    
    
    
    
    
    
    
    
    
    //Now, write the results to disk:
    

    
    
    //Write detailed data for each well into one folder per well:
    for (int i=0; i<MCMCs.size(); i++){
        std::string dir = basedir + "/wells/" + MCMCs[i][0]->label;
        mkdir(dir.c_str(), 0750);
        std::ofstream file((dir + "/data.csv").c_str());
        
        
        for (int t = 0; t<(csv.uncutdata(i)).size(); t++){ file << csv.uncutdata(i)[t] << ","; }
        file << "\n";
        
        for(int j=0; j<MCMCs[i].size(); j++){
            std::vector<double> md;
            md = MCMCs[i][j]->md(MCMCs[i][j]->mean());
            for (int t = 0; t<md.size(); t++){ file << md[t] << ","; }
            file << "\n";
        }
        file.close();
        
        file.open((dir + "/data_initialparameters.csv").c_str());
        for (int t = 0; t<(csv.uncutdata(i)).size(); t++){ file << csv.uncutdata(i)[t] << ","; }
        file << "\n";
        for(int j=0; j<MCMCs[i].size(); j++){
            std::vector<double> md;
            md = MCMCs[i][j]->md(MCMCs[i][j]->startparameters);
            for (int t = 0; t<md.size(); t++){ file << md[t] << ","; }
            file << "\n";
        }
        file.close();
        
        
        
        for(int j=0; j<MCMCs[i].size(); j++){
            file.open((dir + "/parameters_" + itoa(j) + ".txt").c_str());
            file << MCMCs[i][j]->textoutput();
            file.close();
            file.open((dir + "/log_" + itoa(j) + ".txt").c_str());
            file << MCMCs[i][j]->logfile.str();
            file.close();
        }
        
    }
    
    

    
    //Write big CSV file:
    
    std::ofstream file((basedir + "/results_" + filename + ".csv").c_str());
    
    //Header:
    file << "Plate, ";
    file << "Well, ";
    file << "RunTime, ";
    file << "BuildTime, ";
    
    file << "Best-fit, ";
    file << "Best-fit-Max, ";
    file << "Growth/NoGrowth, ";
    file << "Best-fit-BIC, ";
    file << "Lag, ";
    file << "var, ";
    file << "Rate, ";
    file << "var, ";
    file << "Max, ";
    file << "var, ";
    
    for (int j=0; j<MCMCs[0].size(); j++){
        file << ", ";
        file << MCMCs[0][j]->modelname << ", ";
        file << "DIC, ";
        file << "BIC, ";
        file << "accepted, ";
        file << "invalid, ";
        file << "lag, ";
        file << "var, ";
        file << "max rate, ";
        file << "var, ";
        file << "max, ";
        file << "var, ";
        if (MCMCs[0][j]->modelname != "Dummy model"){
            for (int k=0; k<MCMCs[0][j]->dim; k++){
                file << MCMCs[0][j]->parameternames[k] << ", ";
                file << "SD, ";
            }
        }
        file << ", ";
        
    }
    file << ", AUC, \n";
    
    //Body:
    for(int i=0; i<MCMCs.size(); i++){
        
        //Filename
        file << fn << " , ";
        
        //Well
        file << MCMCs[i][0]->label << ", ";
        
        //Run time
        file << timestamp.str() << " , ";
        
        //Build time
        file << __DATE__ << " " << __TIME__ << ", ";
        
        
        
        
        //Best fit (according to DIC):
        int bestmodel = 0;
        for (int j=0; j<MCMCs[i].size(); j++) {
            if (MCMCs[i][j]->_IC_() < MCMCs[i][bestmodel]->_IC_() || isnan(MCMCs[i][bestmodel]->_IC_())) bestmodel=j;
        }
        
        
        if (!isnan(MCMCs[i][bestmodel]->_IC_())){
            file << MCMCs[i][bestmodel]->modelname;
        } else {
            file << "error";
        }
        file << ", ";
        
        
        // Growth vs No-growth:
        
        file << MCMCs[i][bestmodel]-> results_max();
        
        file << ", ";
        
        if (controlwell_max_highquantile >= MCMCs[i][bestmodel]-> results_max()){
            file << "no-growth";
        } else {
            file << "growth";
        }
        
        
        file << ", ";
        
        
        
        
        
        
        //Best fit according to BIC:
        int bestmodel_bic = 0;
        for (int j=0; j<MCMCs[i].size(); j++) {
            if (MCMCs[i][j]->bic() < MCMCs[i][bestmodel_bic]->bic() || isnan(MCMCs[i][bestmodel_bic]->bic())) bestmodel_bic=j;
        }
        
        
        if (!isnan(MCMCs[i][bestmodel_bic]->_IC_())){
            file << MCMCs[i][bestmodel_bic]->modelname;
        } else {
            file << "error";
        }
        file << ", ";
        
        
        
        
        
        
        //Lag:
        file << MCMCs[i][bestmodel]->results_timetoendoflag() << ", ";
        file << MCMCs[i][bestmodel]->results_timetoendoflag_var() << ", ";
        
        //Max Rate:
        file << MCMCs[i][bestmodel]->results_maxrate() << ", ";
        file << MCMCs[i][bestmodel]->results_maxrate_var() << ", ";
        
        //Max:
        file << MCMCs[i][bestmodel]->results_max() << ", ";
        file << MCMCs[i][bestmodel]->results_max_var() << ", ";
        
        
        
        
        
        
        
        
        //Individual models:
        for (int j=0; j<MCMCs[0].size(); j++){
                file << ", ";
                file << MCMCs[i][j]->modelname << ", ";
                file << MCMCs[i][j]->dic() << ", ";
                file << MCMCs[i][j]->bic() << ", ";
                file << MCMCs[i][j]->acceptance_count << ", ";
                file << MCMCs[i][j]->invalid_count << ", ";
                file << MCMCs[i][j]->results_timetoendoflag() << ", ";
                file << MCMCs[i][j]->results_timetoendoflag_var() << ", ";
                file << MCMCs[i][j]->results_maxrate() << ", ";
                file << MCMCs[i][j]->results_maxrate_var() << ", ";
                file << MCMCs[i][j]->results_max() << ", ";
            file << MCMCs[i][j]->results_max_var() << ", ";
            if (MCMCs[0][j]->modelname != "Dummy model"){
                for (int k=0; k<MCMCs[0][j]->dim; k++){
                    file << MCMCs[i][j]->logtolin(MCMCs[i][j]->mean())[k] << ", ";
                    file << MCMCs[i][j]->logvartolinvar(MCMCs[i][j]->variances(), MCMCs[i][j]->mean())[k] << ", ";
                }
            }
                file << ", ";
            
        }
        
        
        
        
        //AUC of the raw data:
        
        file << ", ";
        double auc=0;
        std::vector<double> od = csv.uncutdata(i);
        for (int k=0; k<od.size(); k++){
            auc += od[k];
        }
        auc /= od.size();
        
        file << auc << ", \n";
        
    
        
    }
    
    
    
    
    file.close();
    
    
    
    
    chdir(basedir.c_str());
    
    // Generate some plots in R.
    system("Rscript ~/bin/output.Rscript");
    
    chdir("../..");
    

    
    
    
    
    return 0;
}

