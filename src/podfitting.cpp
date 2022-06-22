#ifndef __PODFITTING
#define __PODFITTING
        
void allocate_memory(descriptorstruct &desc, neighborstruct &nb, podstruct pod, datastruct data)
{
    int nd = pod.nd;
    desc.gd = (double *) malloc(nd*sizeof(double));
    desc.A = (double *) malloc(nd*nd*sizeof(double));    
    desc.b = (double *) malloc(nd*sizeof(double));    
    desc.c = (double *) malloc(nd*sizeof(double));    
    cpuArraySetValue(desc.A, 0.0, nd*nd);
    cpuArraySetValue(desc.b, 0.0, nd);    
    cpuArraySetValue(desc.c, 0.0, nd);    
    
    int dim = 3;
    int natom_max = data.num_atom_max;    
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;    
    int nd4 = pod.nd4;
    int nelements = pod.nelements;
    int nbesselpars = pod.nbesselpars;
    int nabf3 = pod.nabf3;
    int nrbf3 = pod.nrbf3;
    int *pdegree3 = pod.threebody;
    int *pbc = pod.pbc;
    double rcut = pod.rcut;
                
    int Nj=0, Nij=0, Nijk=0;
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
        
    nb.y = (double*) malloc (sizeof (double)*(ny));
    nb.alist = (int*) malloc (sizeof (int)*(na));    
    nb.pairnum = (int*) malloc (sizeof (int)*(natom_max));
    nb.pairnum_cumsum = (int*) malloc (sizeof (int)*(natom_max+1));
    nb.pairlist = (int*) malloc (sizeof (int)*(np));           
    nb.elemindex = (int*) malloc (sizeof (int)*(nelements*nelements));     
        
    int k = 1;
    for (int i=0; i < nelements; i++) 
        for (int j=0; j < nelements; j++) {
            nb.elemindex[i + nelements*j] = k;
            nb.elemindex[j + nelements*i] = k;
            k += 1;
        }            
    
    nb.natom_max = natom_max;
    nb.sze = nelements*nelements;
    nb.sza = na;
    nb.szy = ny;    
    nb.szp = np;    
    
    std::cout<<"**************** Begin of Memory Allocation ****************"<<std::endl;
    
    int szd = 0, szi=0;
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
        Nij = podfullneighborlist(nb.y, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum, x, a1, a2, a3, rcut, pbc, natom);
    
        int ns3 = pdegree3[0]*nbesselpars + pdegree3[1];

        Nj=0, Nijk=0;
        for (int i=0; i < natom; i++) {
            Nj = (Nj > nb.pairnum[i]) ? Nj : nb.pairnum[i];
            Nijk +=  (nb.pairnum[i]-1)*nb.pairnum[i]/2;
        }

        int szd1 = 3*Nij+ (1+dim)*Nij*(nrbf3+ns3) + 2*(1+dim)*Nijk*nrbf3 + 4*Nj*nrbf3 + 2*dim*Nijk + 7*Nijk*nabf3;
        int szi1 = 4*Nij + 2*natom+1 + 2*Nijk + (Nj-1)*Nj + 6*Nijk;                
        szd = PODMAX(szd, szd1);   
        szi = PODMAX(szi, szi1);   
        //std::cout<<"Ni, Nj, Nij, Nijk: "<<natom<<", "<<Nj<<", "<<Nij<<", "<<Nijk<<std::endl;
    }        
    
    szd = PODMAX(natom_max*(nd1+nd2+nd3+nd4) + szd, dim*natom_max*(nd-nd1-nd2-nd3-nd4));
    szd = dim*natom_max*(nd1+nd2+nd3+nd4) + szd;        
    
    // gdd includes linear descriptors derivatives, quadratic descriptors derivatives and temporary memory
    desc.gdd = (double*) malloc (sizeof (double)*(szd)); // [ldd qdd]
    desc.tmpint = (int*) malloc (sizeof (int)*(szi));    
    desc.szd = szd;
    desc.szi = szi;        

    std::cout<<"maximum number of atoms in periodic domain: "<<natom_max<<std::endl;
    std::cout<<"maximum number of atoms in extended domain: "<<nb.sza<<std::endl;
    std::cout<<"maximum number of neighbors in extended domain: "<<nb.szp<<std::endl;
    std::cout<<"size of double memory: "<<szd<<std::endl;
    std::cout<<"size of int memory: "<<szi<<std::endl;
    std::cout<<"size of descriptor matrix: "<<nd<<" x "<<nd<<std::endl;
    std::cout<<"**************** End of Memory Allocation ****************"<<std::endl<<std::endl;
}

void linear_descriptors(descriptorstruct &desc, neighborstruct &nb, podstruct pod, datastruct data, int ci)
{
    int dim = 3;    
    int nelements = pod.nelements;
    int nbesselpars = pod.nbesselpars;
    int nrbf2 = pod.nbf2;
    int nabf3 = pod.nabf3;
    int nrbf3 = pod.nrbf3;
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int nd4 = pod.nd4;
    int *pdegree2 = pod.twobody;
    int *pdegree3 = pod.threebody;
    int *pbc = pod.pbc;
    double rin = pod.rin;
    double rcut = pod.rcut;
    double *Phi2 = pod.Phi2;
    double *Phi3 = pod.Phi3;
    double *besselparams = pod.besselparams;
        
    int natom = data.num_atom[ci];    
    int natom_cumsum = data.num_atom_cumsum[ci];    
    int *atomtype = &data.atomtype[natom_cumsum];
    double *position = &data.position[dim*natom_cumsum];
    double *lattice = &data.lattice[9*ci];
    double *a1 = &lattice[0];
    double *a2 = &lattice[3];
    double *a3 = &lattice[6];
                
    // neighbor list
    int Nij = podfullneighborlist(nb.y, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum, 
                position, a1, a2, a3, rcut, pbc, natom);
    
    double *fatom1 = &desc.gdd[0];
    double *fatom2 = &desc.gdd[dim*natom*(nd1)];
    double *fatom3 = &desc.gdd[dim*natom*(nd1+nd2)];    
    double *eatom1 = &desc.gdd[dim*natom*(nd1+nd2+nd3)];
    double *eatom2 = &desc.gdd[dim*natom*(nd1+nd2+nd3)+natom*nd1];
    double *eatom3 = &desc.gdd[dim*natom*(nd1+nd2+nd3)+natom*(nd1+nd2)];
    double *tmpmem = &desc.gdd[dim*natom*(nd1+nd2+nd3+nd4)+natom*(nd1+nd2+nd3+nd4)];
    int *tmpint = &desc.tmpint[0];   
    
//     std::cout<<dim*natom*(nd1+nd2+nd3)+natom*(nd1+nd2)<<std::endl;
//     std::cout<<dim*natom*(nd1+nd2+nd3+nd4)+natom*(nd1+nd2+nd3+nd4)<<std::endl;
//     std::cout<<dim<<" "<<natom<<" "<<nd1<<" "<<nd2<<" "<<nd3<<" "<<nd4<<std::endl;
//     error("here");
    
    cpuArraySetValue(eatom1, 0.0, natom*(nd1+nd2+nd3+nd4));
    cpuArraySetValue(fatom1, 0.0, dim*natom*(nd1+nd2+nd3+nd4));    
    
    // peratom descriptors for one-body, two-body, and three-body linear potentials
    poddesc(eatom1, fatom1, eatom2, fatom2, eatom3, fatom3, nb.y, Phi2, Phi3, besselparams, 
            tmpmem, rin, rcut, atomtype, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum, 
            nb.elemindex, pdegree2, pdegree3, tmpint, nbesselpars, nrbf2, nrbf3, nabf3, 
            nelements, Nij, natom);            
        
    // global descriptors for one-body, two-body, and three-body linear potentials
    int nd123 = nd1+nd2+nd3;    
    cpuArraySetValue(tmpmem, 1.0, natom);
    DGEMV(&cht, &natom, &nd123, &one, eatom1, &natom, tmpmem, &inc1, &zero, desc.gd, &inc1);    
    
//     writearray2file("Phi2.bin", Phi2, pod.ns2*pod.ns2, 1);
//     writearray2file("Phi3.bin", Phi3, pod.ns3*pod.ns3, 1);
//     writearray2file("gd123.bin", desc.gd, nd123, 1);
//     writearray2file("eatom.bin", eatom1, natom*nd123, 1);
//     writearray2file("fatom.bin", fatom1, dim*natom*nd123, 1);

    //print_matrix( "Two-body eigevectors:", pod.ns2, nrbf2, Phi2, pod.ns2);        
    //print_matrix( "Three-body eigevectors:", pod.ns3, nrbf3, Phi3, pod.ns3);        
//     print_matrix( "Lattice vector 1:", 1, 3, a1, 1);        
//     print_matrix( "Lattice vector 2:", 1, 3, a2, 1);        
//     print_matrix( "Lattice vector 3:", 1, 3, a3, 1); 
//     print_matrix( "Atom positions:", 3, natom, position, 3);        
//     print_matrix( "Atom types:", 1, natom, atomtype, 1);       
//     print_matrix( "Periodic boundary conditions:", 1, 3, pbc, 1); 
//     print_matrix( "Bessel parameters:", 1, 3, besselparams, 1); 
//     print_matrix( "Two-body potential:", 1, 3, pdegree2, 1); 
//     print_matrix( "Three-body potential:", 1, 4, pdegree3, 1); 
//     std::cout<<std::endl<<"Inner cut-off radius: "<<rin<<std::endl;
//     std::cout<<std::endl<<"Outer cut-off radius: "<<rcut<<std::endl;
//     std::cout<<std::endl<<"Number of bessel parameters: "<<nbesselpars<<std::endl;
//     std::cout<<std::endl<<"Number of two-body basis functions: "<<nrbf2<<std::endl;
//     std::cout<<std::endl<<"Number of three-body radial basis functions: "<<nrbf3<<std::endl;
//     std::cout<<std::endl<<"Number of three-body angular basis functions: "<<nabf3<<std::endl;
//     std::cout<<std::endl<<"Number of elements: "<<nelements<<std::endl;
//     std::cout<<std::endl<<"Number of atoms: "<<natom<<std::endl;
//     
//     print_matrix( "POD descriptors:", natom, nd1+nd2+nd3, eatom1, natom); 
    
//     print_matrix( "One-body descriptors:", natom, nd1, eatom1, natom); 
//     print_matrix( "One-body descriptors derivarives:", 3*natom, nd1, fatom1, 3*natom); 
//     print_matrix( "Two-body descriptors:", natom, nd2, eatom2, natom); 
//     print_matrix( "Two-body descriptors derivarives:", 3*natom, nd2, fatom2, 3*natom); 
//     print_matrix( "Three-body descriptors:", natom, nd3, eatom3, natom); 
//     print_matrix( "Three-body descriptors derivarives:", 3*natom, nd3, fatom3, 3*natom); 
//     error("here");
    
}

void quadratic_descriptors(descriptorstruct &desc, podstruct pod, datastruct data, int ci)
{    
    int dim = 3;
    int natom = data.num_atom[ci];    
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int nd4 = pod.nd4;
    int nd22 = pod.nd22;
    int nd23 = pod.nd23;
    int nd24 = pod.nd24;
    int nd33 = pod.nd33;
    int nd34 = pod.nd34;
    int nd44 = pod.nd44;    
    int nd123 = nd1+nd2+nd3;    
    int nd1234 = nd1+nd2+nd3+nd4;    
    
    double *fatom2 = &desc.gdd[dim*natom*(nd1)];
    double *fatom3 = &desc.gdd[dim*natom*(nd1+nd2)];    
    double *fatom4 = &desc.gdd[dim*natom*(nd123)];    
    
    // global descriptors for four-body quadratic22 potential
    if (nd22>0) {
        int nq2 = pod.quadratic22[0]*pod.nc2;        
        quadratic_descriptors(&desc.gd[nd1234], &desc.gdd[dim*natom*nd1234], 
                &desc.gd[nd1], &desc.gd[nd1], fatom2, fatom2, nq2, nq2, dim*natom);
    }
    
    // global descriptors for four-body quadratic23 potential
    if (nd23>0) {
        int nq2 = pod.quadratic23[0]*pod.nc2;        
        int nq3 = pod.quadratic23[1]*pod.nc3;              
        quadratic_descriptors(&desc.gd[nd1234+nd22], &desc.gdd[dim*natom*(nd1234+nd22)], 
                &desc.gd[nd1], &desc.gd[nd1+nd2], fatom2, fatom3, nq2, nq3, dim*natom);        
    }
    
    // global descriptors for five-body quadratic24 potential
    if (nd24>0) {
        int nq2 = pod.quadratic24[0]*pod.nc2;        
        int nq4 = pod.quadratic24[1]*pod.nc4;              
        quadratic_descriptors(&desc.gd[nd1234+nd22+nd23], &desc.gdd[dim*natom*(nd1234+nd22+nd23)], 
                &desc.gd[nd1], &desc.gd[nd1+nd2+nd3], fatom2, fatom4, nq2, nq4, dim*natom);        
    }
    
    // global descriptors for five-body quadratic33 potential
    if (nd33>0) {
        int nq3 = pod.quadratic33[0]*pod.nc3;        
        quadratic_descriptors(&desc.gd[nd1234+nd22+nd23+nd24], &desc.gdd[dim*natom*(nd1234+nd22+nd23+nd24)], 
                &desc.gd[nd1+nd2], &desc.gd[nd1+nd2], fatom3, fatom3, nq3, nq3, dim*natom);        
    }

    // global descriptors for six-body quadratic34 potential
    if (nd34>0) {
        int nq3 = pod.quadratic34[0]*pod.nc3; 
        int nq4 = pod.quadratic34[1]*pod.nc4; 
        quadratic_descriptors(&desc.gd[nd1234+nd22+nd23+nd24+nd33], &desc.gdd[dim*natom*(nd1234+nd22+nd23+nd24+nd33)], 
                &desc.gd[nd1+nd2], &desc.gd[nd1+nd2+nd3], fatom3, fatom4, nq3, nq4, dim*natom);        
    }

    // global descriptors for seven-body quadratic44 potential
    if (nd44>0) {
        int nq4 = pod.quadratic44[1]*pod.nc4; 
        quadratic_descriptors(&desc.gd[nd1234+nd22+nd23+nd24+nd33+nd34], &desc.gdd[dim*natom*(nd1234+nd22+nd23+nd24+nd33+nd34)], 
                &desc.gd[nd1+nd2+nd3], &desc.gd[nd1+nd2+nd3], fatom4, fatom4, nq4, nq4, dim*natom);        
    }        
    
//     writearray2file("fatom2.bin", fatom2, dim*natom*nd2, 1);
//     writearray2file("fatom3.bin", fatom3, dim*natom*nd3, 1);
//     writearray2file("gdesc.bin", desc.gd, pod.nd, 1);
//     writearray2file("gddesc.bin", desc.gdd, dim*natom*pod.nd, 1);
}

void least_squares_matrix(descriptorstruct &desc, podstruct pod, datastruct data, int ci)
{        
    int dim = 3;
    int natom = data.num_atom[ci];    
    int natom_cumsum = data.num_atom_cumsum[ci];    
    int nd = pod.nd;
    int nforce = dim*natom;
    
    // compute energy weight and force weight
    double normconst = 1.0;
    if (data.normalizeenergy==1) normconst = 1.0/natom;        
    double we = data.fitting_weights[0];
    double wf = data.fitting_weights[1];
    double we2 = (we*we)*(normconst*normconst);
    double wf2 = (wf*wf);
            
    // get energy and force from the training data set
    double energy = data.energy[ci];    
    double *force = &data.force[dim*natom_cumsum];
        
    // least-square matrix for all descriptors: A = A + (we*we)*(gd^T * gd)      
    cpuKron(desc.A, desc.gd, desc.gd, we2, nd, nd);
                
    // least-square matrix for all descriptors derivatives: A =  A + (wf*wf) * (gdd^T * gdd)
    DGEMM(&cht, &chn, &nd, &nd, &nforce, &wf2, desc.gdd, &nforce, desc.gdd, &nforce, &one, desc.A, &nd);    
        
    // least-square vector for all descriptors: b = b + (we*we*energy)*gd    
    double wee = we2*energy;
    //DAXPY(&nd, &wee, desc.gd, &inc1, desc.b, &inc1);
    for (int i = 0; i< nd; i++)
        desc.b[i] += wee*desc.gd[i];    
    
    // least-square vector for all descriptors derivatives: b = b + (wf*wf) * (gdd^T * f)
    DGEMV(&cht, &nforce, &nd, &wf2, desc.gdd, &nforce, force, &inc1, &one, desc.b, &inc1);    
    
//     writearray2file("A.bin", desc.A, pod.nd*pod.nd, 1);
//     writearray2file("b.bin", desc.b, pod.nd, 1);
//     error("here");
}

void InverseMatrix(double* A, double* work, int* ipiv, int n)
{
    int lwork = n*n;
    int info;
    DGETRF(&n,&n,A,&n,ipiv,&info);
    DGETRI(&n,A,&n,ipiv,work,&lwork,&info);
}

void least_squares_fit(descriptorstruct &desc, neighborstruct &nb, podstruct pod, datastruct data)
{        
    std::cout<<"**************** Begin of Least-Squares Fitting ****************"<<std::endl;
    
    // loop over each configuration in the training data set
    for (int ci=0; ci < (int) data.num_atom.size(); ci++) {
        if ((ci % 100)==0) std::cout<<"Configuration: # "<<ci+1<<std::endl;
        
        // compute linear POD descriptors
        linear_descriptors(desc, nb, pod, data, ci);
        
        // compute quadratic POD descriptors
        quadratic_descriptors(desc, pod, data, ci);        
        
        // assemble the least-squares linear system
        least_squares_matrix(desc, pod, data, ci);          
    }
    
    int nd = pod.nd;
    
//     print_matrix( "Least-squares matrix:", nd, nd, desc.A, nd); 
//     print_matrix( "Least-squares vector:", 1, nd, desc.b, 1); 
    
    for (int i = 0; i<nd; i++) {       
        desc.c[i] = desc.b[i];
        desc.A[i + nd*i] = desc.A[i + nd*i]*(1.0 + 1e-12);
    }
    
    // solving the linear system A * c = b
    int nrhs=1, info;    
    DPOSV(&chu, &nd, &nrhs, desc.A, &nd, desc.c, &nd, &info);
    
//     double *work = &desc.gdd[0];  
//     int *ipiv = &desc.tmpint[0];
            
//     // compute A^{-1}  
//     InverseMatrix(desc.A, work, ipiv, nd);    
//     
//     // c = A^{-1} * b
//     DGEMV(&chn, &nd, &nd, &one, desc.A, &nd, desc.b, &inc1, &zero, desc.c, &inc1);    
    
    print_matrix( "Least-squares coefficient vector:", 1, nd, desc.c, 1); 
//     error("here");
    
    // save coefficients into a text file
    std::ofstream myfile ("coefficients.txt");
    if (myfile.is_open())
    {
        myfile << "POD_coefficients: " + std::to_string(nd) + "\n";
        for(int count = 0; count < nd; count ++)
            myfile << std::setprecision(20)<< desc.c[count] << std::endl;        
        myfile.close();
    }
    else std::cout << "Unable to open file";
    
    std::cout<<"**************** End of Least-Squares Fitting ****************"<<std::endl<<std::endl;
}

void print_analysis(datastruct data, double *outarray, double *errors)
{                
    string s = "All files";
    int nfiles = data.data_files.size();    // number of files    
    int lm = s.size();    
    for (int i = 0; i < nfiles; i++) 
        lm = PODMAX(lm, (int) data.filenames[i].size());                
    lm = lm + 2;
    
    std::string filename = data.training ? "training_errors.txt" : "test_errors.txt";   
    std::ofstream myfile (filename);
    if (!myfile.is_open()) std::cout << "Unable to open file";
    
    filename = data.training ? "training_analysis.txt" : "test_analysis.txt";
    std::ofstream mfile (filename);
    if (!mfile.is_open()) std::cout << "Unable to open file";
        
    std::string sa = "**************** Begin of Error Analysis for the Training Data Set ****************";
    std::string sb = "**************** Begin of Error Analysis for the Test Data Set ****************";            
    std::string mystr = (data.training) ? sa : sb;    
    std::cout<<mystr<<std::endl;            
    myfile <<mystr+"\n";
    
    sa = "----------------------------------------------------------------------------------------\n";
    sb = "  File      | # configs | # atoms  | MAE energy | RMSE energy | MAE force | RMSE force |\n";
    std::cout<<sa; myfile <<sa;
    std::cout<<sb; myfile <<sb;
    std::cout<<sa; myfile <<sa;   
    
    int ci=0, m=8, nc=0, nf=0;    
    for (int file = 0; file < nfiles; file++) {        
        mfile<<"# " + data.filenames[file] + "\n";
        sb = "|  config  |  # atoms  |  energy  | DFT energy | energy error |  force  | DFT force | force error |\n";
        mfile <<sb;
        
        int nforceall = 0;
        int nconfigs = data.num_config[file];
        nc += nconfigs;
        for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file
            
            mfile <<"     "; 
            for(int count = 0; count < m; count ++)
                mfile << outarray[count + m*ci] << "      ";  
            mfile<<std::endl;
            
            nforceall += 3*data.num_atom[ci];            
            ci += 1;
        }
        nf += nforceall;
        
        int q = file+1;
        string s = data.filenames[file];
        s = s + std::string(lm-s.size(), ' ');       
        string s1 = std::to_string(nconfigs);
        s1 = s1 + std::string(PODMAX(6- (int) s1.size(),1), ' ');    
        s = s + "   " + s1;
        s1 = std::to_string(nforceall/3);
        s1 = s1 + std::string(PODMAX(7 - (int) s1.size(),1), ' ');    
        s = s + "   " + s1;
        s1 = std::to_string(errors[0 + 4*q]);
        s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');    
        s = s + "   " + s1;
        s1 = std::to_string(errors[1 + 4*q]);
        s1 = s1 + std::string(PODMAX(10 - (int)  s1.size(),1), ' ');    
        s = s + "   " + s1;
        s1 = std::to_string(errors[2 + 4*q]);
        s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');    
        s = s + "   " + s1;
        s1 = std::to_string(errors[3 + 4*q]);
        s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');    
        s = s + "   " + s1 + "\n";        
        std::cout<<s; 
        myfile <<s;        
    }
    std::cout<<sa; myfile <<sa;
        
    s = s + std::string(PODMAX(lm - (int) s.size(),1), ' ');       
    string s1 = std::to_string(nc);
    s1 = s1 + std::string(PODMAX(6- (int) s1.size(),1), ' ');    
    s = s + "   " + s1;
    s1 = std::to_string(nf/3);
    s1 = s1 + std::string(PODMAX(7 - (int) s1.size(),1), ' ');    
    s = s + "   " + s1;
    s1 = std::to_string(errors[0]);
    s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');    
    s = s + "   " + s1;
    s1 = std::to_string(errors[1]);
    s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');    
    s = s + "   " + s1;
    s1 = std::to_string(errors[2]);
    s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');    
    s = s + "   " + s1;
    s1 = std::to_string(errors[3]);
    s1 = s1 + std::string(PODMAX(10 - (int) s1.size(),1), ' ');    
    s = s + "   " + s1 + "\n";        
    std::cout<<s; myfile <<s;
    std::cout<<sa; myfile <<sa;
    
    sa = "**************** End of Error Analysis for the Training Data Set ****************";        
    sb = "**************** End of Error Analysis for the Test Data Set ****************";     
    mystr = (data.training) ? sa : sb;    
    std::cout<<mystr<<std::endl;            
    myfile <<mystr+"\n";    
    myfile.close();
    mfile.close();    
}

void error_analsysis(descriptorstruct &desc, neighborstruct &nb, podstruct pod, datastruct data, double *coeff)
{                
    int dim = 3;
    double energy;
    double force[dim*data.num_atom_max];
    
    int nfiles = data.data_files.size();    // number of files    
    int num_configs = data.num_atom.size(); // number of configurations in all files    
    int nd1234 = pod.nd1 + pod.nd2 + pod.nd3 + pod.nd4; 
            
    int m = 8; 
    double outarray[m*num_configs];            
    double errors[4*(nfiles+1)];        
    for (int i=0; i<4*(nfiles+1); i++)
        errors[i] = 0.0;
                    
    std::cout<<"**************** Begin of Error Calculation ****************"<<std::endl;
    
    int ci = 0; // configuration counter    
    int nc = 0, nf = 0;    
    for (int file = 0; file < nfiles; file++) { // loop over each file in the training data set
        double emae=0.0, essr=0.0, fmae=0.0, fssr=0.0;
        int nforceall = 0;
                
        int nconfigs = data.num_config[file];
        nc += nconfigs;
        for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file
            if ((ci % 100)==0) std::cout<<"Configuration: # "<<ci+1<<std::endl;
            
            int natom = data.num_atom[ci];
            int nforce = dim*natom;

            // compute linear POD descriptors
            linear_descriptors(desc, nb, pod, data, ci);

//             // compute quadratic POD descriptors
//             quadratic_descriptors(desc, pod, data, ci);        
//             
//             // calculate energy and force
//             energy = calculate_energyforce(force, desc.gd, desc.gdd, coeff, pod.nd, natom);
            
            // calculate energy and force
            energy = calculate_energyforce(force, desc.gd, desc.gdd, coeff, &desc.gdd[nforce*nd1234], 
                        pod.quadratic22, pod.quadratic23, pod.quadratic24, pod.quadratic33, 
                        pod.quadratic34, pod.quadratic44, pod.nd1, pod.nd2, pod.nd3, pod.nd4, 
                        pod.nelements, pod.nc2, pod.nc3, pod.nc4, natom);
            
//            print_matrix( "global descriptors:", pod.nd, 1, desc.gd, pod.nd); 
            
            double DFTenergy = data.energy[ci];   
            int natom_cumsum = data.num_atom_cumsum[ci];    
            double *DFTforce = &data.force[dim*natom_cumsum];     

//             print_matrix( "predicted force:", 3, natom, force, 3); 
//             print_matrix( "DFT force:", 3, natom, DFTforce, 3); 
            
            outarray[0 + m*ci] = ci+1;
            outarray[1 + m*ci] = natom;
            outarray[2 + m*ci] = energy;
            outarray[3 + m*ci] = DFTenergy;        
            outarray[4 + m*ci] = fabs(DFTenergy-energy)/natom;        
            outarray[5 + m*ci] = cpuArrayNorm(force, nforce);
            outarray[6 + m*ci] = cpuArrayNorm(DFTforce, nforce);

            double diff, sum = 0.0, ssr = 0.0;
            for (int j=0; j<dim*natom; j++) {
                diff = DFTforce[j] - force[j]; 
                sum += fabs(diff);
                ssr += diff*diff;
            }
            outarray[7 + m*ci] = sum/nforce;
            //outarray[8 + m*ci] = sqrt(ssr/nforce);        
                                                
            emae += outarray[4 + m*ci];
            essr += outarray[4 + m*ci]*outarray[4 + m*ci];                    
            fmae += sum;
            fssr += ssr;            
            nforceall += nforce;                        
            ci += 1;             
            //std::cout<<"configuration: "<<ci<<",   predicted energy: "<<energy<<",   DFT energy: "<<DFTenergy<<std::endl;
            //error("here");
        }
        int q = file + 1;
        errors[0 + 4*q] = emae/nconfigs; 
        errors[1 + 4*q] = sqrt(essr/nconfigs);
        errors[2 + 4*q] = fmae/nforceall; 
        errors[3 + 4*q] = sqrt(fssr/nforceall); 
                
        nf += nforceall;
        errors[0] += emae; 
        errors[1] += essr; 
        errors[2] += fmae; 
        errors[3] += fssr;         
    }   
    errors[0] = errors[0]/nc;
    errors[1] = sqrt(errors[1]/nc);
    errors[2] = errors[2]/nf;
    errors[3] = sqrt(errors[3]/nf);
    
    std::cout<<"**************** End of Error Calculation ****************"<<std::endl<<std::endl;
    
    print_analysis(data, outarray, errors);    
}

#endif

