/***************************************************************************
  Compile code: clang++ -std=c++17 -O3 -Wall -llapack -lblas main.cpp -o pod
  Run code: ./../../src/pod pod.txt data.txt  
****************************************************************************/

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
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
#include "pod.cpp"
#include "podfitting.cpp"

int main(int argc, char** argv) 
{
    if (argc < 3) {
      printf("Usage: ./pod InputFile1 InputFile2\n");
      return 1;
    }                

    std::string pod_file = std::string(argv[1]); // pod input file
    std::string data_file = std::string(argv[2]); // data input file           
    std::string coeff_file = "";
//     if (argc >= 3)
//         coeff_file = std::string(argv[3]); // coefficient input file           
    
    podstruct pod;        
    datastruct traindata;    
    datastruct testdata;    
    descriptorstruct desc;            
    neighborstruct nb;
    
    read_input_files(pod, traindata, testdata, pod_file, data_file, coeff_file);        
                
    allocate_memory(desc, nb, pod, traindata);
        
    //error_analsysis(desc, nb, pod, traindata, pod.coeff);    
    
    least_squares_fit(desc, nb, pod, traindata);
        
    error_analsysis(desc, nb, pod, traindata, desc.c);    
    
    if ((testdata.data_path != "") && (testdata.data_path != traindata.data_path)) {
        error_analsysis(desc, nb, pod, testdata, desc.c);
    }    
}






