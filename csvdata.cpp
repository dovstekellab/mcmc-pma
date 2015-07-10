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

#include <iostream>
#include "csvdata.h"


csvdata::csvdata(std::string filename){
    std::vector<std::vector<double> > file_data_vectors(0);
    std::ifstream  data(filename.c_str());
    if (!data.good()){
        std::cerr << "\nError: File not found!\n";
        exit(0);
    }
    std::string line;
    std::string        cell;
    
   
    //Get labels:
    std::vector<std::string> l(1);
    while(std::getline(data,line, ‘\n’)) {
        std::vector<std::string> ltemp(0);
        std::stringstream  lineStream(line);
        
        //
        while(std::getline(lineStream,cell,','))
        {
            ltemp.push_back(cell.c_str());
        }
        l = ltemp;
        if (l[0]=="Hour" || l[0]=="Time") break;
    }
             
    int numcols = (int)l.size();
    
    //Get rid of first cell:
    for (int i=0; i<l.size()-1; i++){
        l[i] = l[i+1];
    }
    l.pop_back();
        
    //Get data:
    while(std::getline(data,line, ‘\n’))
    {
        std::stringstream  lineStream(line);
        std::vector<double> current_line(0);
        while(std::getline(lineStream,cell,','))
        {
            current_line.push_back(atof(cell.c_str()));
        }
        if (current_line.size()==numcols){
            file_data_vectors.push_back(current_line);
        }
    }
    std::vector<std::vector<double> > transpose(0);
    for (int i=1; i<file_data_vectors[0].size(); i++){  //starting at i=1 to get rid of Hour/Time row
        std::vector<double> current_row(0);
        for (int j=0; j<file_data_vectors.size(); j++){
            current_row.push_back(file_data_vectors[j][i]);
            
        }
        transpose.push_back(current_row);
    }
    
    
    labels = l;
    raw_data = transpose;
    original_length = (int)raw_data[0].size();
    
    
    std::vector<int> c(raw_data.size());
    for (int i=0; i<raw_data.size(); i++){
        c[i] = -1;
    }

    cutoff = c;
    
}



bool csvdata::cut(double threshold){
    std::vector<int> c(raw_data.size());
    for (int i=0; i<raw_data.size(); i++){
        double v=0;
        for (int j=4; j<raw_data[i].size()-4; j++){
            v += pow(raw_data[i][j]-(raw_data[i][j-4]+raw_data[i][j-3]+raw_data[i][j-2]+raw_data[i][j-1]+raw_data[i][j]+raw_data[i][j+1]+raw_data[i][j+2]+raw_data[i][j+3]+raw_data[i][j+4])/9, 2);
        }
        v/=raw_data[i].size()-9;
        if (v<5) v=5;
        double sd = sqrt(v);
        
        int max = 0;
        int min = 0;
        for (int j=4; j < raw_data[i].size()-4; j++){

            if ((raw_data[i][j-4]+raw_data[i][j-3]+raw_data[i][j-2]+raw_data[i][j-1]+raw_data[i][j]+raw_data[i][j+1]+raw_data[i][j+2]+raw_data[i][j+3]+raw_data[i][j+4])/9 > (raw_data[i][max-4]+raw_data[i][max-3]+raw_data[i][max-2]+raw_data[i][max-1]+raw_data[i][max]+raw_data[i][max+1]+raw_data[i][max+2]+raw_data[i][max+3]+raw_data[i][max+4])/9) max = j;
        }

        
        int k=max;
        bool cut = 0;
        while (k < raw_data[i].size()-4 && cut==0){
            
            int m=0; int n=0;
            for (int j=k; j < raw_data[i].size()-4; j++){
                if (raw_data[i][j] < (raw_data[i][max-4]+raw_data[i][max-3]+raw_data[i][max-2]+raw_data[i][max-1]+raw_data[i][max]+raw_data[i][max+1]+raw_data[i][max+2]+raw_data[i][max+3]+raw_data[i][max+4])/9 - 0.5*sd){
                    m++;
                } else {
                    n++;
                }
            }
            if (n<0.1*(m+n)){
                cut = 1;
                break;
            }
            k++;
        }

        
        if (cut==1 && raw_data[i][max] > raw_data[i][0] + 6*sd) {
            
            if (k<40) k = 40;
        
            c[i] = k;
            
        }
        else {
            c[i] = -1;
        }
    }
    cutoff = c;
    return 1;
}

int csvdata::size(){
    return (int)raw_data.size();
}

int csvdata::length(){
    return original_length;
}


std::vector<double> csvdata::data(int well){
    std::vector<double> d = raw_data[well];
    if (cutoff[well] != -1) d.resize(cutoff[well]);
    return d;
}

std::vector<double> csvdata::uncutdata(int well){
    return raw_data[well];
}

