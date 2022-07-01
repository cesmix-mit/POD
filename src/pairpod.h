
#ifndef PAIR_POD_H
#define PAIR_POD_H

//#include "pod.h"

class CPairPOD  {
private:
    std::vector<std::string> globVector(const std::string& pattern, std::vector<std::string> & files);
        
    bool is_a_number(std::string line);
    
    int latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx);
    
    int podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N, int dim);
    
    void read_data_file(double *inputs, std::string &file_format, std::string &file_extension, 
         std::string &data_path, std::string data_file);

    void get_exyz_files(std::vector<std::string>& files, std::string datapath, std::string extension);
    
    int get_number_atom_exyz(std::vector<int>& num_atom, int& num_atom_sum, std::string file);  
    
    int get_number_atoms(std::vector<int>& num_atom, std::vector<int> &num_atom_sum, std::vector<int>& num_config, std::vector<std::string> training_files);  
    
    void read_exyz_file(double *lattice, double *stress, double *energy, double *pos, double *vel, 
            double *forces, int *atomtype, std::string file, std::vector<std::string> species);
    
    void get_data(std::vector<std::string> species);

    void read_data_files(std::string data_file, std::vector<std::string> species);    
public:
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

        double *lattice=NULL;
        double *energy=NULL; 
        double *stress=NULL;
        double *position=NULL;
        double *velocity=NULL;
        double *force=NULL;
        int *atomtype=NULL;

        int analysis = 0;
        int runMD = 0;
        int savecalculation = 0;
        int savefrequency = 0;

        void copydatainfo(datastruct &data) {
            data.data_path = data_path;        
            data.file_format = file_format;
            data.file_extension = file_extension;              
            data.data_files = data_files;
            data.filenames = filenames;            
            data.analysis = analysis;
            data.runMD = runMD;
            data.savecalculation = savecalculation;
        }                

        void freememory(int backend)
        {
            TemplateFree(lattice, backend);        
            TemplateFree(energy, backend);        
            TemplateFree(stress, backend);        
            TemplateFree(position, backend);  
            TemplateFree(velocity, backend);  
            TemplateFree(force, backend);        
            TemplateFree(atomtype, backend);        
        }            
    };

    datastruct data;
    class CPOD *podptr;
     
    CPairPOD();
    ~CPairPOD();
    
    // constructor 
    CPairPOD(std::string pod_file, std::string coeff_file, std::string data_file); 
    
    int backend=1;
    int dim = 3;
    int atommemory = 0;
    int podpairlist=0;
    int lammpspairlist=0;
    
    int *atomtype=NULL;
    double *pos;
    double *vel;
    double energy;        
    double *force;
    double *stress;
    
    double *gd=NULL;         // global linear descriptors 
    double *podcoeff=NULL;   // POD coefficients  
    double *effcoeff=NULL;   // effective coefficients
    int ncoeff;              // number of coefficients
    int ndesc;               // number of linear descriptors 
    
    void compute(int, int);
    void settings(int, char **);
    void coeff(int, char **);
    void init_style();
    double init_one(int, int);
    double memory_usage();  
    void *extract(const char *, int &); 
    
    void free_memory();    
    void allocate_memory();
    void allocate_memory(int Ni, int Ng, int Nij);
    void allocate_memory(datastruct data);
    
    
    int podfullneighborlist(double *y, int *alist, int *pairlist, int *pairnum, int *pairnumsum, 
        double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx);

    //void podNeighPairs(double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int natom);
        
    void podNeighPairs(int istart, int iend);
        
    void lammpsNeighPairs(double **x, int **firstneigh, int *atomtype, int *numneigh, 
            int *ilist, int istart, int iend);
            
    double podenergyforce();   
    double lammpsenergyforce(double **f, double **x, int **firstneigh, int *atomtype, 
            int *numneigh, int *ilist);     
    
protected:       
    double *rij=NULL;     // xj - xi    
    int *idxi=NULL;
    int *ai=NULL;
    int *aj=NULL;
    int *ti=NULL;
    int *tj=NULL;
            
    double *y=NULL;     // positions of own and ghost atoms    
    int *pairlist; //     
    int *pairnum;
    int *pairnumsum=NULL;    
    int *atomlist;    
    
    double *tmpmem;
    int *tmpint;
    int szd;    // size of tmpmem
    int szi;    // size of tmpint

    int nboxatom;    // number of atom in the simulation box    
    int nlocalatom;   // number of owned atoms
    int nghostatom; // number of ghost atoms
    int ntotalatom; // number of owned + ghost atoms
    int nmaxatom;  // maximum number of atoms (nmaxatom >= ntotalatom)   
    int nblock; // number of atoms per computation block    
    int nblockmax; // maximum number of atoms per computation block    
    int nij;    //  number of atom pairs per computation block    
    int nijmax;  // maximum number of atom pairs per computation block      
    
    int compblocks; // number of computation blocks    
    int *istart;    
    int pairstart;
    
    double **scale;         // for thermodynamic integration  
};

#endif
