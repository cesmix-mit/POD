/***************************************************************************
  Compile code: clang++ -std=c++11 -O3 -Wall -llapack -lblas main.cpp -o pod
  /home/linuxbrew/.linuxbrew/opt/llvm@11/bin/clang++ -std=c++11 -stdlib=libc++ -O3 -Wall -llapack -lblas -Wl,-rpath=/home/linuxbrew/.linuxbrew/lib main.cpp -o pod
  Run code: ./../../src/pod pod.txt data.txt  
****************************************************************************/

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <glob.h>

#define _USE_MATH_DEFINES

#include <math.h>

using std::cout;
using std::endl;
using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::ostringstream;

#include "podcommon.h"
#include "cpuArrayOperations.cpp"
#include "readinputfiles.cpp"
#include "podneighborlist.cpp"
#include "snap.cpp"
#include "pod.cpp"
#include "podfitting.cpp"

int main(int argc, char** argv) 
{
    if (argc < 3) {
      printf("Usage: ./pod InputFile1 InputFile2\n");
      return 1;
    }                

    std::string pod_file = std::string(argv[1]);  // pod input file
    std::string data_file = std::string(argv[2]); // data input file           
    std::string coeff_file;                       // coefficient input file           
    
    if (argc > 3)
        coeff_file = std::string(argv[3]); // coefficient input file           
    else
        coeff_file = "";
        
    // data structures defined in podcommon.h
    podstruct pod;        
    snastruct sna;        
    datastruct traindata;    
    datastruct testdata;    
    descriptorstruct desc;            
    neighborstruct nb;
    
    read_input_files(pod, traindata, testdata, pod_file, data_file, coeff_file);                            
        
    if (pod.snaptwojmax > 0) {
        InitSnap(sna, pod.snapelementradius, pod.snapelementweight, pod.rcut, 
            0.0, pod.snaprfac0, pod.snaptwojmax, pod.nelements, pod.snapchemflag);
        //sna.printout();
    }
    
    // allocate memory for data structures
    if ((int) traindata.data_path.size() > 1) 
        allocate_memory(desc, nb, pod, sna, traindata);    
    else if ((int) testdata.data_path.size() > 1)
        allocate_memory(desc, nb, pod, sna, testdata);
        
    if (coeff_file != "") // get POD coefficients from an input file           
        cpuArrayCopy(desc.c, pod.coeff, pod.nd);
    else // compute POD coefficients using least-squares method
        least_squares_fit(desc, nb, pod, sna, traindata);
    
    // calculate errors for the training data set
    if ((traindata.training_analysis) && ((int) traindata.data_path.size() > 1) )
        error_analsysis(desc, nb, pod, sna, traindata, desc.c);    
            
    // calculate errors for the test data set
    if ((testdata.test_analysis) && ((int) testdata.data_path.size() > 1) && (testdata.data_path != traindata.data_path)) 
        error_analsysis(desc, nb, pod, sna, testdata, desc.c);    
    
    // calculate energy and force for the training data set
    if ((traindata.training_calculation) && ((int) traindata.data_path.size() > 1) )
        energyforce_calculation(desc, nb, pod, sna, traindata, desc.c);   
    
    // calculate energy and force for the test data set
    if ((testdata.test_calculation) && ((int) testdata.data_path.size() > 1) && (testdata.data_path != traindata.data_path) )
        energyforce_calculation(desc, nb, pod, sna, testdata, desc.c);   
}





