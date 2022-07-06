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
#include <random>
#include <algorithm>

#define _USE_MATH_DEFINES

#include <math.h>

using std::cout;
using std::endl;
using std::string;
using std::ios;
using std::ofstream;
using std::ifstream;
using std::ostringstream;

void PrintErrorAndExit(const char* errmsg, const char *file, int line ) 
{    
    printf( "%s in %s at line %d\n", errmsg, file, line );
    
#ifdef  USE_MPI       
    MPI_Finalize();    
#endif
    
    exit( 1 );    
}

#define pod_error( errmsg ) (PrintErrorAndExit( errmsg, __FILE__, __LINE__ ))

#include "podcommon.h"
#include "pod.cpp"
#include "podfit.cpp"
#include "pairpod.cpp"

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
            
    if (argc == 3) {
        // create podfit object
        CPODFIT podfit(pod_file, data_file, coeff_file);

        // compute POD coefficients using least-squares method
        podfit.least_squares_fit(podfit.traindata);

        // calculate errors for the training data set
        if ((podfit.traindata.training_analysis) && ((int) podfit.traindata.data_path.size() > 1) )
            podfit.error_analsysis(podfit.traindata, podfit.desc.c);    

        // calculate errors for the test data set
        if ((podfit.testdata.test_analysis) && ((int) podfit.testdata.data_path.size() > 1) && (podfit.testdata.data_path != podfit.traindata.data_path)) 
            podfit.error_analsysis(podfit.testdata, podfit.desc.c);    

        // calculate energy and force for the training data set
        if ((podfit.traindata.training_calculation) && ((int) podfit.traindata.data_path.size() > 1) )
            podfit.energyforce_calculation(podfit.traindata, podfit.desc.c);   

        // calculate energy and force for the test data set
        if ((podfit.testdata.test_calculation) && ((int) podfit.testdata.data_path.size() > 1) && (podfit.testdata.data_path != podfit.traindata.data_path) )
            podfit.energyforce_calculation(podfit.testdata, podfit.desc.c);   
    }
    else {
        CPairPOD pairpod(pod_file, coeff_file, data_file);        
        pairpod.error_analsysis();  
    }    
}





