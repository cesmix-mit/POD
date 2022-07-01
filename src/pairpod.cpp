
#include "pairpod.h"

// #include "atom.h"
// #include "comm.h"
// #include "error.h"
// #include "force.h"
// #include "memory.h"
// #include "neigh_list.h"
// #include "neighbor.h"
// #include "tokenizer.h"

#include <cmath>
#include <cstring>

//using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
// CPairPOD::CPairPOD(LAMMPS *lmp) : Pair(lmp)
// {
//   single_enable = 0;
//   restartinfo = 0;
//   one_coeff = 1;
//   manybody_flag = 1;
//   centroidstressflag = CENTROID_NOTAVAIL;
// 
//   radelem = nullptr;
//   wjelem = nullptr;
//   coeffelem = nullptr;
//   sinnerelem = nullptr;
//   dinnerelem = nullptr;
// 
//   beta_max = 0;
//   beta = nullptr;
//   bispectrum = nullptr;
//   podptr = nullptr;
// }


/* ----------------------------------------------------------------------
   This version is a straightforward implementation
   ---------------------------------------------------------------------- */
void CPairPOD::compute(int eflag, int vflag)
{

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void CPairPOD::settings(int narg, char ** /* arg */)
{
//   if (narg > 0)
//     error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void CPairPOD::coeff(int narg, char **arg)
{
//   if (!allocated) allocate();
//   if (narg != 4 + atom->ntypes) error->all(FLERR,"Incorrect args for pair coefficients");
// 
//   map_element2type(narg-4,arg+4);
// 
//   // read podcoeff and podparam files
// 
//   read_files(arg[2],arg[3]);

  // set default scaling
  int n = podptr->pod.nelements;
  for (int ii = 0; ii < n+1; ii++)
    for (int jj = 0; jj < n+1; jj++)
      scale[ii][jj] = 1.0;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void CPairPOD::init_style()
{
//   if (force->newton_pair == 0)
//     error->all(FLERR,"Pair style POD requires newton pair on");
// 
//   // need a full neighbor list
// 
//   neighbor->add_request(this, NeighConst::REQ_FULL);
// 
//   podptr->init();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double CPairPOD::init_one(int i, int j)
{
//   if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");
    scale[j][i] = scale[i][j];
    
    return 0.0;
}

/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double CPairPOD::memory_usage()
{
  //double bytes = Pair::memory_usage();

  double bytes = 0.0;
  
  return bytes;
}

void *CPairPOD::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  return nullptr;
}

// constructor 
CPairPOD::CPairPOD(std::string pod_file, std::string coeff_file, std::string data_file)
{
    podptr = new CPOD(pod_file, coeff_file);
    
    this->read_data_files(data_file, podptr->pod.species);       
    
    if ((int) data.data_path.size() > 1) 
        this->allocate_memory(data);    
                
    // get POD coefficients from an input file           
    if (coeff_file != "") {        
        TemplateMalloc(&podcoeff, podptr->pod.nd, backend);
        podptr->podArrayCopy(podcoeff, podptr->pod.coeff, podptr->pod.nd);    
    }
}

// destructor        
CPairPOD::~CPairPOD()
{
    data.freememory(1);
}

bool CPairPOD::is_a_number(std::string line)
{    
    return isdigit(line.at(0));
}

void CPairPOD::read_data_file(double *inputs, std::string &file_format, std::string &file_extension, 
        std::string &data_path, std::string data_file)
{
    std::ifstream file_in(data_file);
    if (!file_in) {pod_error("Error: Data input file is not found");}
    
    std::string line;        
    while (std::getline(file_in, line)) // Read next line to `line`, stop if no more lines.
    {                                            
        if (line != "") {
            std::string s, s2;
            double d;
            
            std::istringstream ss_line(line);                                    
            ss_line >> s;
                                    
            if (s == "error_analysis_for_data_set") {
                ss_line >> d;           
                inputs[0] = d;               
            }
            else if (s == "save_calculation_in_binary_files") {
                ss_line >> d;           
                inputs[1] = d;               
            }
            else if (s == "save_frequency") {
                ss_line >> d;           
                inputs[2] = d;               
            }
            else if (s == "run_molecular_dynamics_simulation") {
                ss_line >> d;           
                inputs[3] = d;               
            }
            else if (s != "#") {
                ss_line >> s2;                
                if (s == "file_format") file_format = s2;
                if (s == "file_extension") file_extension = s2;
                if (s == "path_to_data_set") {
                    data_path = s2;                    
                    while(ss_line >> s2){    
                        data_path = data_path + " " + s2;   
                    }                        
                    //data_path.erase(std::remove(data_path.begin(), data_path.end(), '"'), data_path.end());                                        
                    data_path.erase(data_path.begin());                    
                    data_path.erase(data_path.end()-1);                                        
                }                
            }
        }
    }        
    file_in.close();
    
    std::cout<<"**************** Begin of Data File ****************"<<std::endl;
    std::cout<<"file format: "<<file_format<<std::endl;
    std::cout<<"file extension: "<<file_extension<<std::endl;
    std::cout<<"path to data set: "<<data_path<<std::endl;
    std::cout<<"error analysis for data set: "<<inputs[0]<<std::endl;
    std::cout<<"save calculation fin binary files: "<<inputs[1]<<std::endl;
    std::cout<<"save frequency: "<<inputs[2]<<std::endl;
    std::cout<<"run MD simulation: "<<inputs[3]<<std::endl;
    std::cout<<"**************** End of Data File ****************"<<std::endl<<std::endl;
}

std::vector<std::string> CPairPOD::globVector(const std::string& pattern, std::vector<std::string> & files)
{
    glob_t glob_result;
    glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
      std::string s = string(glob_result.gl_pathv[i]);
      //std::cout << s << "\n";
      files.push_back(s);
    }
    globfree(&glob_result);
    return files;
}

void CPairPOD::get_exyz_files(std::vector<std::string>& files, std::string datapath, std::string extension)  
{
    std::vector<std::string> res = this->globVector(datapath + "/*." + extension, files);
} 

int CPairPOD::get_number_atom_exyz(std::vector<int>& num_atom, int& num_atom_sum, std::string file)  
{
    std::ifstream datafile(file);
    if (!datafile) {/*error*/}    
               
    std::string line; 
    int num_configs = 0;
    num_atom_sum = 0;
    while (std::getline(datafile, line)) // Read next line to `line`, stop if no more lines.
    {
        int d;
        if (this->is_a_number(line)) {            
            d = std::stoi(line);
            num_atom.push_back(d);                                        
            num_configs += 1;
            num_atom_sum += d;
        }        
    }    
    datafile.close();
    
    return num_configs;
}

int CPairPOD::get_number_atoms(std::vector<int>& num_atom, std::vector<int> &num_atom_sum, std::vector<int>& num_config, std::vector<std::string> training_files)  
{
    int nfiles = training_files.size(); // number of files
    int d, n;
    
    for (int i=0; i<nfiles; i++) {
        d = this->get_number_atom_exyz(num_atom, n, training_files[i]);  
        num_config.push_back(d);  
        num_atom_sum.push_back(n);  
    }
        
    int num_atom_all = 0;
    for (int i=0; i< (int) num_atom.size(); i++)
        num_atom_all += num_atom[i];
    
    return num_atom_all;
}

void CPairPOD::read_exyz_file(double *lattice, double *stress, double *energy, double *pos, 
        double *vel, double *forces, int *atomtype, std::string file, std::vector<std::string> species)  
{
    std::ifstream datafile(file);
    if (!datafile) {/*error*/}    
               
    std::string substr1 = "nergy";
    //std::string substr2 = "ttice";
    std::string substr3 = "tress";
        
    int cfi = 0, nat=0, k = 0, ns = species.size();
    double d; 
    
    std::string line;           
    while (std::getline(datafile, line)) // Read next line to `line`, stop if no more lines.
    {
        if (line.substr(1,6) == "attice") {
            int index1 = line.find(substr1);
            //int index2 = line.find(substr2);
            int index3 = line.find(substr3);
                                    
            if (line.substr(index1,6) == "nergy=") {      
                std::string s1 = line.substr(index1+6,-1);
                int ind = s1.find(" ");
                energy[cfi] = stod(line.substr(index1+6,ind));
            }
            
            //std::cout<<line<<std::endl;
            if (line.substr(1,6) == "attice") {      
                //std::string s1 = line.substr(0,-1);
                int ind1 = line.find("\"");
                int ind2 = line.find("\"",ind1+1);
                std::istringstream ss(line.substr(ind1+1,ind2-ind1));    
                //std::cout<<line.substr(ind1+1,ind2-ind1-1)<<std::endl;
                k = 0;
                while(ss >> d){            
                    lattice[k + 9*cfi] = d;
                    k += 1;
                }                
            }
            
            if (line.substr(index3,6) == "tress=") {      
                std::string s1 = line.substr(index3+7,-1);
                int ind = s1.find("\"");
                std::istringstream ss(line.substr(index3+7,ind));    
                k = 0;
                while(ss >> d){            
                    stress[k + 9*cfi] = d;
                    k += 1;
                }                
            }
                        
            cfi += 1;            
        }
        else if (!this->is_a_number(line)) {            
            std::string s0;
            std::istringstream ss(line);                        
            ss >> s0;
                    
            for (int ii=0; ii<ns; ii++)
                if (species[ii] == s0)
                    atomtype[nat] = ii+1;
                
            k = 0;
            while(ss >> d){            
                if (k <=2 ) pos[k + 3*nat] = d;
                if (k > 2 && k <5) vel[k-3 + 3*nat] = d;
                if (k > 5 ) forces[k-6 + 3*nat] = d;
                k += 1;
            }          
            nat += 1;             
        }                
    }    
    datafile.close();    
}

void CPairPOD::get_data(std::vector<std::string> species)
{    
    this->get_exyz_files(data.data_files, data.data_path, data.file_extension);
    data.num_atom_sum = this->get_number_atoms(data.num_atom, data.num_atom_each_file, data.num_config, data.data_files);            
    data.num_config_sum = data.num_atom.size();
    
    std::cout<<"data file     |    number of configurations     |     number of atoms "<<std::endl;
    for (int i=0; i< (int) data.data_files.size(); i++) {
        string filename = data.data_files[i].substr(data.data_path.size()+1,data.data_files[i].size());
        data.filenames.push_back(filename.c_str());                                        
        std::cout<<data.filenames[i]<<"   |   "<<data.num_config[i]<<"   |   "<<data.num_atom_each_file[i]<< std::endl;  
        //std::cout<<data.data_files[i].substr(data.data_path.size()+1,data.data_files[i].size())<<std::endl;
    }    
    std::cout << "number of files: " <<data.data_files.size() << std::endl;   
    std::cout << "number of configurations in all files: " <<data.num_config_sum << std::endl;   
    std::cout << "number of atoms in all files: " <<data.num_atom_sum << std::endl;       
    
    int n = data.num_config_sum;
    data.lattice = (double *) malloc(9*n*sizeof(double));
    data.stress = (double *) malloc(9*n*sizeof(double));
    data.energy = (double *) malloc(n*sizeof(double));    
    n = data.num_atom_sum;
    data.position = (double *) malloc(3*n*sizeof(double));
    data.force = (double *) malloc(3*n*sizeof(double));
    data.atomtype = (int *) malloc(n*sizeof(int));    
    
    int nfiles = data.data_files.size(); // number of files    
    int nconfigs = 0;
    int natoms = 0;
    for (int i=0; i<nfiles; i++) {        
        this->read_exyz_file(&data.lattice[9*nconfigs], &data.stress[9*nconfigs], &data.energy[nconfigs], 
                &data.position[3*natoms], &data.velocity[3*natoms], &data.force[3*natoms], &data.atomtype[natoms], 
                data.data_files[i], species);          
        nconfigs += data.num_config[i];
        natoms += data.num_atom_each_file[i];
    }    
    
    int len = data.num_atom.size();
    data.num_atom_min = podptr->podArrayMin(&data.num_atom[0], len);    
    data.num_atom_max = podptr->podArrayMax(&data.num_atom[0], len);    
    data.num_atom_cumsum.resize(len+1);
    podptr->podCumsum(&data.num_atom_cumsum[0], &data.num_atom[0], len+1); 
    
    data.num_config_cumsum.resize(nfiles+1);
    podptr->podCumsum(&data.num_config_cumsum[0], &data.num_config[0], nfiles+1); 
    
    std::cout << "minimum number of atoms: " <<data.num_atom_min << std::endl;   
    std::cout << "maximum number of atoms: " <<data.num_atom_max << std::endl;   
}

void CPairPOD::read_data_files(std::string data_file, std::vector<std::string> species)
{        
    double inputs[100];
    
    // read data input file to datastruct
    this->read_data_file(inputs, data.file_format, data.file_extension, data.data_path, data_file);
        
    data.analysis = (int) inputs[0];
    data.savecalculation = (int) inputs[1];
    data.savefrequency = (int) inputs[2];
    data.runMD = (int) inputs[3];    

    if ((int) data.data_path.size() > 1) 
        this->get_data(species);
    else 
        pod_error("data set is not found");                 
}

int CPairPOD::latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
{
    int m=0, n=0, p=0;
    if (pbc[0] == 1) m = (int) ceil(rcut/a1[0]);    
    if (pbc[1] == 1) n = (int) ceil(rcut/a2[1]);    
    if (pbc[2] == 1) p = (int) ceil(rcut/a3[2]);                
        
    // index for the center lattice
    int ind = m + (2*m+1)*(n) + (2*m+1)*(2*n+1)*(p);
    
    // number of lattices
    int nl = (2*m+1)*(2*n+1)*(2*p+1);            
        
    //y = zeros(3, nx*nl)
    for (int j=0; j<3*nx; j++)
        y[j] = x[j];
    int q = nx;
        
    for (int i = 0; i < (2*p+1); i++)
        for (int j = 0; j < (2*n+1); j++)
            for (int k = 0; k < (2*m+1); k++) {
                int ii = k + (2*m+1)*j + (2*m+1)*(2*n+1)*i;                
                if (ii != ind) {                    
                    double x0 = a1[0]*(k - m) + a2[0]*(j - n) + a3[0]*(i - p);
                    double x1 = a1[1]*(k - m) + a2[1]*(j - n) + a3[1]*(i - p);
                    double x2 = a1[2]*(k - m) + a2[2]*(j - n) + a3[2]*(i - p);       
                    for (int jj=0; jj<nx; jj++) {
                        y[0+3*q] = x0 + x[0+3*jj]; 
                        y[1+3*q] = x1 + x[1+3*jj]; 
                        y[2+3*q] = x2 + x[2+3*jj]; 
                        q = q + 1;                               
                    }
                }
            }                
    
    //alist = zeros(Int32,nx*nl)
    for (int i=0; i <nl; i++)
        for (int j=0; j<nx; j++) 
            alist[j + nx*i] = j;            
    
    return nl;
}

int CPairPOD::podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N, int dim)
{
    int k = 0;
    for (int i = 0; i<nx; i++) {
        double *ri = &r[i*dim];
        int inc = 0;
        for (int j=0; j<N; j++) {
            double *rj = &r[dim*j];                        
            double rijsq = (ri[0]-rj[0])*(ri[0]-rj[0]) + (ri[1]-rj[1])*(ri[1]-rj[1]) + (ri[2]-rj[2])*((ri[2]-rj[2]));
            if  ((rijsq > 1e-12) && (rijsq <= rcutsq))  { 
                inc += 1;                                
                neighlist[k] = j;          
                k += 1;                                                  
            }
        }
        numneigh[i] = inc; 
    }
    return k; 
}

int CPairPOD::podfullneighborlist(double *xy, int *alist, int *neighlist, int *numneigh, int *numneighsum, 
        double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
{
    double rcutsq = rcut*rcut;    
    int dim = 3, nl = 0, nn = 0;
    
    // number of lattices
    nl = this->latticecoords(xy, alist, x, a1, a2, a3, rcut, pbc, nx);        
    int N = nx*nl;
            
    // total number of neighbors
   nn = this->podneighborlist(neighlist, numneigh, xy, rcutsq, nx, N, dim);
    
   podptr->podCumsum(numneighsum, numneigh, nx+1); 
       
   return nn;
}   

void CPairPOD::allocate_memory(datastruct data)
{
    //int nd = podptr->pod.nd;    
    int dim = 3;
    int natom_max = data.num_atom_max;    
    int nd1 = podptr->pod.nd1;
    int nd2 = podptr->pod.nd2;
    int nd3 = podptr->pod.nd3;    
    int nd4 = podptr->pod.nd4;
    //int nelements = podptr->pod.nelements;
    int nbesselpars = podptr->pod.nbesselpars;
    int nrbf2 = podptr->pod.nbf2;
    int nabf3 = podptr->pod.nabf3;
    int nrbf3 = podptr->pod.nrbf3;
    int *pdegree2 = podptr->pod.twobody;
    int *pdegree3 = podptr->pod.threebody;
    int *pbc = podptr->pod.pbc;
    double rcut = podptr->pod.rcut;
                
    int Nij=0;
    int m=0, n=0, p=0, nl=0, ny=0, na=0, np=0;
    
    for (int ci=0; ci<(int) data.num_atom.size(); ci++)
    {
        int natom = data.num_atom[ci];    
        double *lattice = &data.lattice[9*ci];
        double *a1 = &lattice[0];
        double *a2 = &lattice[3];
        double *a3 = &lattice[6];
        if (pbc[0] == 1) m = (int) ceil(rcut/a1[0]);    
        if (pbc[1] == 1) n = (int) ceil(rcut/a2[1]);    
        if (pbc[2] == 1) p = (int) ceil(rcut/a3[2]);         
    
        // number of lattices
        nl = (2*m+1)*(2*n+1)*(2*p+1);              
        ny = PODMAX(ny,dim*natom*nl);
        na = PODMAX(na, natom*nl);            
        np = PODMAX(np, natom*natom*nl);        
    }
        
    y = (double*) malloc (sizeof (double)*(ny));
    atomlist = (int*) malloc (sizeof (int)*(na));    
    pairnum = (int*) malloc (sizeof (int)*(natom_max));
    pairnumsum = (int*) malloc (sizeof (int)*(natom_max+1));
    pairlist = (int*) malloc (sizeof (int)*(np));                   
        
    szd = 0, szi=0;
    int szsnap=0;
    for (int ci=0; ci<(int) data.num_atom.size(); ci++)
    {
        int natom = data.num_atom[ci];    
        int natom_cumsum = data.num_atom_cumsum[ci];    
        double *x = &data.position[dim*natom_cumsum];
        double *lattice = &data.lattice[9*ci];
        double *a1 = &lattice[0];
        double *a2 = &lattice[3];
        double *a3 = &lattice[6];
            
        // neighbor list
        Nij = podfullneighborlist(y, atomlist, pairlist, pairnum, pairnumsum, x, a1, a2, a3, rcut, pbc, natom);
    
        int ns2 = pdegree2[0]*nbesselpars + pdegree2[1];
        int ns3 = pdegree3[0]*nbesselpars + pdegree3[1];

        int szd1 = 3*Nij+ (1+dim)*Nij*PODMAX(nrbf2+ns2,nrbf3+ns3) + (nabf3+1)*7;
        int szi1 = 6*Nij + 2*natom+1;                
        szd = PODMAX(szd, szd1);   
        szi = PODMAX(szi, szi1);   
        
        if (podptr->sna.twojmax>0) {
            szd1 = 0;
            szd1 += Nij*dim; // rij
            szd1 += PODMAX(2*podptr->sna.idxu_max*Nij, 2*podptr->sna.idxz_max*podptr->sna.ndoubles*natom); // (Ur, Ui) and (Zr, Zi)             
            szd1 += 2*podptr->sna.idxu_max*dim*Nij; // dUr, dUi            
            szd1 += PODMAX(podptr->sna.idxb_max*podptr->sna.ntriples*dim*Nij, 2*podptr->sna.idxu_max*podptr->sna.nelements*natom); // dblist and (Utotr, Utoti)                                  
            szsnap = PODMAX(szsnap, szd1);   
        }        
    }            
    szd = PODMAX(szsnap, szd);   
    szd = natom_max*(nd1+nd2+nd3+nd4) + szd;  
    
    tmpmem = (double*) malloc (sizeof (double)*(szd)); // [ldd qdd]
    tmpint = (int*) malloc (sizeof (int)*(szi));        
}

void CPairPOD::free_memory()
{
    if (atommemory) {
        TemplateFree(atomtype, backend);
        TemplateFree(pos, backend);
        TemplateFree(vel, backend);
        TemplateFree(force, backend);
        TemplateFree(stress, backend);
    }
    
    TemplateFree(rij, backend);
    TemplateFree(idxi, backend);
    TemplateFree(ai, backend);
    TemplateFree(aj, backend);
    TemplateFree(ti, backend);
    TemplateFree(tj, backend);
    
    if (podpairlist) {
        TemplateFree(y, backend);
        TemplateFree(pairlist, backend);
        TemplateFree(pairnum, backend);
        TemplateFree(pairnumsum, backend);
        TemplateFree(atomlist, backend);                
    }
    
    TemplateFree(tmpmem, backend);
    TemplateFree(tmpint, backend);    
}

void CPairPOD::allocate_memory()
{
    if (atommemory) {
        TemplateMalloc(&atomtype, nmaxatom, backend);
        TemplateMalloc(&pos, dim*nmaxatom, backend);
        TemplateMalloc(&vel, dim*nmaxatom, backend);
        TemplateMalloc(&force, dim*nmaxatom, backend);
        TemplateMalloc(&stress, 9, backend);
    }
    
    TemplateMalloc(&rij, dim*nijmax, backend);
    TemplateMalloc(&idxi, nijmax, backend);
    TemplateMalloc(&ai, nijmax, backend);
    TemplateMalloc(&aj, nijmax, backend);
    TemplateMalloc(&ti, nijmax, backend);
    TemplateMalloc(&tj, nijmax, backend);

    if (podpairlist) {
        TemplateMalloc(&y, dim*nmaxatom, backend);        
        TemplateMalloc(&pairnum, nmaxatom, backend);
        TemplateMalloc(&pairnumsum, nmaxatom+1, backend);
        TemplateMalloc(&pairlist, nijmax, backend);
        TemplateMalloc(&atomlist, nijmax, backend);                
    }
    
    TemplateMalloc(&tmpmem, szd, backend);
    TemplateMalloc(&tmpint, szi, backend);
}


void CPairPOD::allocate_memory(int Ni, int Ng, int Nij)
{
    nlocalatom = Ni;
    nghostatom = Ng;
    ntotalatom = Ni + Ng;
    nij = Nij;    
    if ((nij > nijmax) || (ntotalatom > nmaxatom)) {
        nijmax = nij;
        nmaxatom = ntotalatom;   
        this->free_memory();
        this->allocate_memory();
    }
}

void CPairPOD::podNeighPairs(int istart, int iend)
{
    nblock = iend - istart;
    
    nij = 0;
    for (int ii=0; ii<nblock; ii++) {  // for each atom i in the simulation box     
        int gi = istart + ii;       // atom i
        int itype = atomtype[gi];
        int start = pairnumsum[gi];   
        int m = pairnumsum[gi+1] - start;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = start + l;
            int gj = pairlist[k];  // atom j                                                               
            idxi[k]      = ii;
            ai[k]        = atomlist[gi];
            aj[k]        = atomlist[gj];          
            ti[k]        = itype;       
            tj[k]        = atomtype[aj[k]];              
            rij[k*3+0]   = y[gj*3+0] -  y[gi*3+0];  // xj - xi            
            rij[k*3+1]   = y[gj*3+1] -  y[gi*3+1];  // xj - xi            
            rij[k*3+2]   = y[gj*3+2] -  y[gi*3+2];  // xj - xi            
            nij += 1;
        }
    }                       
    
    pairstart = istart;
}

void CPairPOD::lammpsNeighPairs(double **x, int **firstneigh, int *atype, int *numneigh, 
        int *ilist, int istart, int iend)
{
    nblock = iend - istart;    
    double rcutsq = podptr->pod.rcut*podptr->pod.rcut;
    
    nij = 0;
    for (int ii=0; ii<nblock; ii++) {  // for each atom i in the simulation box             
        int gi = ilist[istart+ii];       // atom i
        int itype = atype[gi];
        int m = numneigh[gi];
        pairnumsum[ii+1] = 0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int gj = firstneigh[gi][l];  // atom j     
            double delx   = x[gj][0] -  x[gi][0];  // xj - xi            
            double dely   = x[gj][1] -  x[gi][1];  // xj - xi            
            double delz   = x[gj][2] -  x[gi][2];  // xj - xi            
            double rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < rcutsq && rsq > 1e-20) {
                rij[nij*3 + 0] = delx;
                rij[nij*3 + 1] = dely;
                rij[nij*3 + 2] = delz;
                idxi[nij]      = ii;
                ai[nij]        = gi;
                aj[nij]        = gj;          
                ti[nij]        = itype;       
                tj[nij]        = atype[gj];       
                nij++;
                pairnumsum[ii+1] += 1;
            }                        
        }
    }    
    
    pairnumsum[0] = 0;
    for (int ii=0; ii<nblock; ii++)
        pairnumsum[ii+1] = pairnumsum[ii+1] + pairnumsum[ii];
              
    pairstart = 0;
}

double CPairPOD::podenergyforce()
{    
    energy = 0.0;
    for (int i = 0; i< compblocks; i++) {
        podNeighPairs(istart[i], istart[i+1]);
    
        energy += podptr->energyforce_calculation(force, podcoeff, effcoeff, gd, rij, 
            tmpmem, &pairnumsum[pairstart], atomtype, idxi, ai, aj, ti, tj, nblock, nij);
    }
    
    return energy;
}

double CPairPOD::lammpsenergyforce(double **f, double **x, int **firstneigh, int *atype, 
            int *numneigh, int *ilist)
{
    energy = 0.0;
    for (int i = 0; i< compblocks; i++) {            
        lammpsNeighPairs(x, firstneigh, atype, numneigh, ilist, istart[i], istart[i+1]);
        
        energy += podptr->energyforce_calculation(force, podcoeff, effcoeff, gd, rij, 
            tmpmem, &pairnumsum[pairstart], atype, idxi, ai, aj, ti, tj, nblock, nij);
    }
    
    for (int i = 0; i<ntotalatom; i++) {
        f[i][0] = force[0+3*i];
        f[i][1] = force[1+3*i];
        f[i][2] = force[2+3*i];
    }
    
    return energy;    
}
