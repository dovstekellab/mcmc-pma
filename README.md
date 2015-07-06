# mcmc-pma
Bayesian parameter estimation for phenotype microarray data

About: A project for fitting different microbial growth models to Biolog [1] phenotype microarray datasets using Bayesian parameter estimation. It implements an adaptive Metropolis algorithm, together with model choice using DIC. For more information see Gerstgrasser M, S Nicholls, M Stout, K Smart, T Kypraios, C Powell, D J Stekel: A Bayesian approach to analysing phenotype microarray data enables estimation of microbial growth parameters. [2]


Requirements: GNU Scientific Library installed and header files accessible. Tested with GSL version 1.16. Optionally, OpenMP can be used. (See main.cpp comments.) Optionally, with R installed and Rscript accessible, PDF plots will be generated.


To compile: "clang++ *.cpp -lgsl -lc -lm -lgslcblas -lgomp -O3 -o ~/bin/mcmc" For gcc replace clang++ by g++. Tested with clang 3.6.0 and gcc 4.9.2 on Ubuntu Linux and Mac OS X, both on x86-64.


To run: "~/bin/mcmc filename.csv" will run the MCMC analysis for a file containing raw Biolog date in CSV form.


Input file format: The program is designed to work with data in CSV format, obtained from excel files exported by Biolog's software. These contain several lines of metadata followed by the coloration data itself. This includes a header row (with lables for each well / phenotype) and a header column (with labels for each time point). The header cell of the header column always contains the string "Time" or "Hour" in the data files we have worked with. The program will therefore skip all lines until it encounters a line starting with one of these strings. It will read column labels from this line, and ignore the first column containing time labels.

If you are using a different data format, export data to a CSV file with wells in columns, time in rows. Add an extra row and column at the beginning, and put the string "Time" in the first cell in the fist line. Optionally, you can add well labels to the remainder of the first line.

Alternatively, this behaviour is easily modified in the csvdata class constructor in the csvdata.cpp file.


Raw data preprocessing: The data we worked with required removing of artifacts at the end of some phenotypes' data. If this is not required remove the corresponding line in main.cpp.


Code comments: The MCMC algorithm and various helper functions are implemented in the (virtual) class mcmc. Individual models are implemented in derived classes. The csvdata class includes functions to import a CSV file, and basic data preprocessing. main.cpp, mcmc.cpp, and mcmc_baranyi.cpp contain some comments to help explain the code structure. Other models' classes work very similarly to the mcmc_baranyi class.


Notes:

[1] See http://www.biolog.com/ for more information. The authors of this software are not affiliated with Biolog.
[2] Forthcoming.




License:
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

