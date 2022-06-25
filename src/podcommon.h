#ifndef PODCOMMON_H
#define PODCOMMON_H

using std::string;

#ifdef USE_CUDA
#include "cudamem.h"
#endif

#ifndef CUDA_SYNC
#define CUDA_SYNC
#endif

#ifndef USE_CUDA    
#define cublasHandle_t int
#define cudaEvent_t int
#endif

#define DDOT ddot_
#define DGEMV dgemv_
#define DGEMM dgemm_
#define DGETRF dgetrf_
#define DGETRI dgetri_
#define DSYEV dsyev_
#define DPOSV dposv_

#define NEIGHMASK 0x3FFFFFFF
#define MAX_ATOM_TYPES      30
#define MAX_MOLECULE_TYPES  20
#define MAX_MOLECULE_SIZE   10
#define PODMIN(a,b) ((a) < (b) ? (a) : (b))
#define PODMAX(a,b) ((a) > (b) ? (a) : (b))

//#define _TIMING
#ifdef _TIMING
#define INIT_TIMING auto begin = std::chrono::high_resolution_clock::now(); auto end = std::chrono::high_resolution_clock::now();
#define TIMING_START  begin = std::chrono::high_resolution_clock::now();   
#define TIMING_END    end = std::chrono::high_resolution_clock::now();   
#define TIMING_GET(num) common.timing[num] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;        
#define START_TIMING {CUDA_SYNC; TIMING_START;}       
#define END_TIMING(num) {CUDA_SYNC; TIMING_END; TIMING_GET(num)}   
#define TIMING_GET1(num) CCal.common.timing[num] += std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count()/1e6;        
#define END_TIMING_CCAL(num) {CUDA_SYNC; TIMING_END; TIMING_GET1(num)}   
#else
#define INIT_TIMING
#define START_TIMING
#define END_TIMING(num)
#define END_TIMING_CCAL(num)
#endif

#define CPUFREE(x)                                                           \
{                                                                         \
    if (x != NULL) {                                                      \
        free(x);                                                          \
        x = NULL;                                                         \
    }                                                                     \
}

extern "C" {
    double DNRM2(int*,double*,int*);
    double DDOT(int*,double*,int*,double*,int*);
    void DAXPY(int*,double*,double*,int*,double*,int*);
    void DGEMV(char*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);  
    void DGEMM(char*,char*,int*,int*,int*,double*,double*,int*,
             double*,int*,double*,double*,int*);        
    void DGETRF(int*,int*,double*,int*,int*,int*);
    void DGETRI(int*,double*,int*,int*,double*,int*,int*);
    void DTRSM(char *, char*, char*, char *, int *, int *, double*, double*, int*,
             double*, int*);
    
    void DSYEV( char* jobz, char* uplo, int* n, double* a, int* lda,
        double* w, double* work, int* lwork, int* info );   
    
    void DPOSV( char* uplo, int* n, int* nrhs, double* a, int* lda,
                double* b, int* ldb, int* info );
}

// // global variables for BLAS  
double one = 1.0;
double minusone = -1.0;
double zero = 0.0;
char chn = 'N';
char cht = 'T';
char chl = 'L';
char chu = 'U';
char chv = 'V';
char chr = 'R';
int inc1 = 1;

// global variables for CUBLAS  
double cublasOne[1] = {one};
double cublasMinusone[1] = {minusone};
double cublasZero[1] = {zero};
                
template <typename T> static void TemplateMalloc(T **data, int n, int backend)
{
    if (backend == 0)       // One thread CPU
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));     
    if (backend == 1)  // Open MP
        // allocate the memory on the CPU        
        *data = (T *) malloc(n*sizeof(T));    
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA C                
        // allocate the memory on the GPU            
        CUDA_CHECK( cudaMalloc( (void**)data, n * sizeof(T) ) );
#endif                  
}

template <typename T> static void TemplateCopytoDevice(T *d_data, T *h_data, int n, int backend)
{
    if (backend == 0)       
        //cpuArrayCopy(d_data, h_data, n);
        for (int i=0; i<n; i++) d_data[i] = h_data[i];
    if (backend == 1)         
        for (int i=0; i<n; i++) d_data[i] = h_data[i];
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA          
        cudaCopytoDevice(d_data, h_data, n);
#endif                  
#ifdef USE_HIP            
    if (backend == 3)  // HIP
        hipCopytoDevice(d_data, h_data, n);
#endif                      
}

template <typename T> static void TemplateCopytoHost(T *h_data, T *d_data, int n, int backend)
{
    if (backend == 0)       
        //cpuArrayCopy(h_data, d_data, n);
        for (int i=0; i<n; i++) h_data[i] = d_data[i];
    if (backend == 1)         
        for (int i=0; i<n; i++) h_data[i] = d_data[i];
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA          
        cudaCopytoHost(h_data, d_data, n);
#endif                  
#ifdef USE_HIP            
    if (backend == 3)  // HIP
        hipCopytoHost(h_data, d_data, n);
#endif                      
}

template <typename T> static void TemplateFree(T *data, int backend)
{
    if (backend == 0)       
        CPUFREE(data);
    if (backend == 1)         
        CPUFREE(data);
#ifdef USE_CUDA            
    if (backend == 2)  // CUDA          
        GPUFREE(data);
#endif                  
#ifdef USE_HIP            
    if (backend == 3)  // HIP          
        HIPFREE(data);
#endif                      
}

int checknan(double *a, int n) 
{
    for (int i=0; i<n; i++)
        if (isnan(a[i]) || fabs(a[i])>1e14) 
            return 1;
            
    return 0;
}

void print_matrix(const char* desc, int m, int n, double* a, int lda ) 
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ ) 
    {
        for( j = 0; j < n; j++ ) printf( " %6.12f", a[i+j*lda] );
        printf( "\n" );
    }
}

void print_matrix(const char* desc, int m, int n, int* a, int lda ) 
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ ) 
    {
        for( j = 0; j < n; j++ ) printf( " %d", a[i+j*lda] );
        printf( "\n" );
    }
}


void PrintErrorAndExit(const char* errmsg, const char *file, int line ) 
{    
    printf( "%s in %s at line %d\n", errmsg, file, line );
    
#ifdef  HAVE_MPI       
    MPI_Finalize();    
#endif
    
    exit( 1 );    
}

#define error( errmsg ) (PrintErrorAndExit( errmsg, __FILE__, __LINE__ ))
#define display( errmsg ) (PrintMsg( errmsg, __FILE__, __LINE__ ))

template <typename T> void readarrayfromfile(const char* filename, T **a, int N)
{
    if (N>0) {
        // Open file to read
        ifstream in(filename, ios::in | ios::binary);

        if (!in) {
            error("Unable to open file ");
        }

        if (in) {
            *a = (T*) malloc (sizeof (T)*N);
            in.read( reinterpret_cast<char*>( *a ), sizeof(T)*N );        
        }

        in.close();
    }
}

template <typename T> void writearray2file(const char* filename, T *a, int N, int backend)
{
    if (N>0) {        
        // Open file to read
        ofstream out(filename, ios::out | ios::binary);

        if (!out) {
            error("Unable to open file ");
        }

        if (backend==2) { //GPU
#ifdef  HAVE_CUDA                        
            T *a_host;            
            a_host = (T*) malloc (sizeof (T)*N);            
            
            // transfer data from GPU to CPU to save in a file
            cudaMemcpy(&a_host[0], &a[0], N*sizeof(T), cudaMemcpyDeviceToHost);    
            
            out.write( reinterpret_cast<char*>( &a_host[0] ), sizeof(T) * N );
            
            free(a_host);
#endif            
        }
        else 
            out.write( reinterpret_cast<char*>( &a[0] ), sizeof(T) * N );                    
        
        out.close();
    }
}

struct podstruct {     
    std::vector<std::string> species;    
    int *pbc; //[3] = {1,1,1};
    
    int nelements = 0;
    int onebody = 1;
    int twobody[3] = {5,10,10};
    int threebody[4] = {4,8,8,5}; 
    int fourbody[4] = {0,0,0,0};    
    
    int quadratic22[2] = {0,0};
    int quadratic23[2] = {0,0};
    int quadratic24[2] = {0,0};
    int quadratic33[2] = {0,0};
    int quadratic34[2] = {0,0};
    int quadratic44[2] = {0,0};        
    int cubic333[3] = {0,0,0};
    int cubic444[3] = {0,0,0};
    
    double rin = 0.5;
    double rcut = 4.6;
    double *besselparams; //[3] = {0.0, 2.0, 4.0};        
    double *Phi2, *Phi3, *Phi4, *Lambda2, *Lambda3, *Lambda4;    
    double *coeff;
        
    int nbesselpars = 3;    
    int ns2, ns3, ns4;       // number of snapshots for radial basis functions for linear POD potentials      
    int nc2, nc3, nc4;       // number of chemical  combinations for linear POD potentials      
    int nbf1, nbf2, nbf3, nbf4; // number of basis functions for linear POD potentials      
    int nd1, nd2, nd3, nd4;     // number of descriptors for linear POD potentials 
    int nd22, nd23, nd24, nd33, nd34, nd44; // number of descriptors for quadratic POD potentials    
    int nd333, nd444; // number of descriptors for cubic POD potentials    
    int nrbf3, nabf3, nrbf4, nabf4;    
    int nd;
    
    int snaptwojmax = 0;
    int snapchemflag = 0;
    double snaprfac0 = 0.99363;
    double snapelementradius[10] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    double snapelementweight[10] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    
    void allocatememory(int backend)
    {
        TemplateMalloc(&pbc, 3, backend);
        TemplateMalloc(&besselparams, 3, backend);
    }    
    
    void freememory(int backend)
    {
        TemplateFree(pbc, backend);        
        TemplateFree(besselparams, backend);        
        TemplateFree(Phi2, backend);        
        TemplateFree(Phi3, backend);        
        TemplateFree(Phi4, backend);        
        TemplateFree(Lambda2, backend);        
        TemplateFree(Lambda3, backend);        
        TemplateFree(Lambda4, backend);    
        TemplateFree(coeff, backend);    
    }        
};

struct datastruct {     
    std::string file_format;
    std::string file_extension;    
    std::string data_path;    
    std::vector<std::string> data_files;     
    std::vector<std::string> filenames;     
    
    std::vector<int> num_atom;
    std::vector<int> num_atom_cumsum;
    std::vector<int> num_atom_each_file;
    std::vector<int> num_config;
    std::vector<int> num_config_cumsum;
    int num_atom_sum; 
    int num_atom_min; 
    int num_atom_max; 
    int num_config_sum;    
    
    double *lattice;
    double *energy; 
    double *stress;
    double *position;
    double *force;
    int *atomtype;
    
    int training = 1;
    int normalizeenergy = 1;
    int training_analysis = 1;
    int test_analysis = 1;
    int training_calculation = 0;
    int test_calculation = 0;
    
    double fitting_weights[10] = {0.0, 0.0, 0.0, 1, 1, 0, 0, 0, 0, 0};
    
    void freememory(int backend)
    {
        TemplateFree(lattice, backend);        
        TemplateFree(energy, backend);        
        TemplateFree(stress, backend);        
        TemplateFree(position, backend);        
        TemplateFree(force, backend);        
        TemplateFree(atomtype, backend);        
    }            
};

// struct tempstruct {     
//     double *tmpmem;
//     int *tmpint;        
//     int szd;
//     int szi;
// };

struct neighborstruct {
    int *elemindex;
    int *alist;
    int *pairnum;
    int *pairnum_cumsum;
    int *pairlist;
    double *y;        
    
    int natom;    
    int nalist;    
    int natom_max;
    int sze;
    int sza;
    int szy;    
    int szp;    
    
    void freememory(int backend)
    {
        TemplateFree(elemindex, backend);        
        TemplateFree(alist, backend);        
        TemplateFree(pairnum, backend);        
        TemplateFree(pairnum_cumsum, backend);        
        TemplateFree(pairlist, backend);        
        TemplateFree(y, backend);        
    }                
};
        
struct descriptorstruct {
    double *gd;  // global descriptors    
    double *gdd; // derivatives of global descriptors and peratom descriptors
    double *A;  // least-square matrix for all descriptors    
    double *b;  // least-square vector for all descriptors    
    double *c;  // coefficents of descriptors
    int *tmpint;        
    int szd;
    int szi;
    
    void freememory(int backend)
    {
        TemplateFree(gd, backend);     
        TemplateFree(gdd, backend);     
        TemplateFree(A, backend);        
        TemplateFree(b, backend);       
        TemplateFree(c, backend);       
        TemplateFree(tmpint, backend);        
    }                
};

struct snastruct {        
    int twojmax;
    int ncoeff;
    int idxb_max;
    int idxu_max;
    int idxz_max;
    int idxcg_max;
    int ntypes;
    int nelements;    
    int ndoubles;   // number of multi-element pairs
    int ntriples;   // number of multi-element triplets      
    int bnormflag;
    int chemflag;    
    int switchflag;
    int bzeroflag;
    int wselfallflag;
    
    double wself;
    double rmin0;
    double rfac0;
    double rcutfac;
    double rcutmax;    
        
    int *map=NULL;  // map types to [0,nelements)    
    int *idx_max=NULL; 
    int *idxz=NULL;
    int *idxz_block=NULL;
    int *idxb=NULL;
    int *idxb_block=NULL;
    int *idxu_block=NULL;
    int *idxcg_block=NULL;
    
    double *rcutsq=NULL;    
    double *radelem=NULL;
    double *wjelem=NULL; 
    double *bzero=NULL;
    double *fac=NULL;
    double *rootpqarray=NULL; 
    double *cglist=NULL;
    
    void printout()
    {
        printf("twojmax %d \n", twojmax); 
        printf("ncoeff %d \n", ncoeff);         
        printf("idxb_max %d \n", idxb_max);         
        printf("idxu_max %d \n", idxu_max);         
        printf("idxz_max %d \n", idxz_max); 
        printf("idxcg_max %d \n", idxcg_max);
        printf("ntypes %d \n", ntypes);
        printf("nelements %d \n", nelements);
        printf("ndoubles %d \n", ndoubles);
        printf("ntriples %d \n", ntriples);
        printf("bnormflag %d \n", bnormflag);
        printf("chemflag %d \n", chemflag);
        printf("switchflag %d \n", switchflag);
        printf("bzeroflag %d \n", bzeroflag);
        printf("wselfallflag %d \n", wselfallflag);        
        printf("rfac0 %g \n", rfac0);
        printf("rmin0 %g \n", rmin0);
        printf("rcutfac %g \n", rcutfac);
        printf("rcutmax %g \n", rcutmax);    
        print_matrix( "map:", 1, ntypes+1, map, 1); 
        print_matrix( "radelem:", 1, ntypes+1, radelem, 1); 
        print_matrix( "wjelem:", 1, ntypes+1, wjelem, 1); 
        print_matrix( "rcutsq:", ntypes+1, ntypes+1, rcutsq, ntypes+1); 
        print_matrix( "bzero:", 1, twojmax+1, bzero, 1);
        print_matrix( "fac:", 1, 20, fac, 1);
        print_matrix( "rootpqarray:", twojmax+1, twojmax+1, rootpqarray, (twojmax+1));
        print_matrix( "cglist:", 1, idxcg_max, cglist, 1);            
    }
    
    void freememory(int backend)
    {   
        TemplateFree(map, backend);
        TemplateFree(idx_max, backend);
        TemplateFree(idxz, backend);
        TemplateFree(idxb, backend);
        TemplateFree(idxb_block, backend);
        TemplateFree(idxu_block, backend);
        TemplateFree(idxz_block, backend);
        TemplateFree(idxcg_block, backend);
        
        TemplateFree(rootpqarray, backend);
        TemplateFree(cglist, backend);
        TemplateFree(fac, backend);
        TemplateFree(bzero, backend);
        TemplateFree(wjelem, backend);
        TemplateFree(radelem, backend);
        TemplateFree(rcutsq, backend);
    }                         
};


#endif  
