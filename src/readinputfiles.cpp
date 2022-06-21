#ifndef __READINPUTFILES
#define __READINPUTFILES

void snapshots(double *rbf, double *xij, double *besselparams, double rin, double rcut, 
        int besseldegree, int inversedegree, int nbesselpars, int N)
{
    double rmax = rcut-rin;
    for (int n=0; n<N; n++) {            
        double dij = xij[n];    

        double r = dij - rin;        
        double y = r/rmax;            
        double y2 = y*y;
        double y3 = 1.0 - y2*y;
        double y4 = y3*y3 + 1e-6;
        double y5 = pow(y4, 0.5);
        double y6 = exp(-1.0/y5);
        double fcut = y6/exp(-1.0);
        
        for (int j=0; j<nbesselpars; j++) {            
            double alpha = besselparams[j];    
            if (fabs(alpha) <= 1.0e-6) alpha = 1e-3;                        
            double x =  (1.0 - exp(-alpha*r/rmax))/(1.0-exp(-alpha));
            //if (n==0) std::cout<<x<<" "<<r<<std::endl;
            
            for (int i=0; i<besseldegree; i++) {
                double a = (i+1)*M_PI;
                double b = (sqrt(2.0/(rmax))/(i+1));
                int nij = n + N*i + N*besseldegree*j;            
                rbf[nij] = b*fcut*sin(a*x)/r;
            }
        }

        for (int i=0; i<inversedegree; i++) {
            int p = besseldegree*nbesselpars + i;
            int nij = n + N*p;     
            double a = pow(dij, (double) (i+1.0));
            rbf[nij] = fcut/a;
        }
        
        //if (n==0) std::cout<<dij<<" "<<rin<<" "<<rcut<<" "<<fcut<<" "<<rbf[0]<<std::endl;
    }
}

void eigenvaluedecomposition(double *Phi, double *Lambda, double *besselparams, double rin, double rcut, 
        int besseldegree, int inversedegree, int nbesselpars, int N)
{
    int ns = besseldegree*nbesselpars + inversedegree;
    
    double *xij = (double *) malloc(N*sizeof(double));
    double *S = (double *) malloc(N*ns*sizeof(double));
    double *Q = (double *) malloc(N*ns*sizeof(double));
    double *A = (double *) malloc(ns*ns*sizeof(double));
    double *b = (double *) malloc(ns*sizeof(double));
    
    for (int i=0; i<N; i++)
        xij[i] = (rin+1e-6) + (rcut-rin-1e-6)*(i*1.0/(N-1));        
    
    snapshots(S, xij, besselparams, rin, rcut, besseldegree, inversedegree, nbesselpars, N);
    
    //print_matrix( "xij", 1, 20, xij, 1);
    //print_matrix( "Snapshot", 10, ns, S, N);
    
//     char chn = 'N';
//     char cht = 'T';
    double alpha = 1.0, beta = 0.0;    
    DGEMM(&cht, &chn, &ns, &ns, &N, &alpha, S, &N, S, &N, &beta, A, &ns);    
        
    for (int i=0; i<ns*ns; i++)
        A[i] = A[i]*(1.0/N);        
    
    // Declaring Function Input for DSYEV
//     char jobz = 'V';    // 'V':  Compute eigenvalues and eigenvectors.
//     char uplo = 'U';    // 'U':  Upper triangle of A is stored.
    int lwork = ns * ns;  // The length of the array work, lwork >= max(1,3*N-1).
    int info = 1;       // = 0:  successful exit
    double work[ns*ns];         
    DSYEV(&chv, &chu, &ns, A, &ns, b, work, &lwork, &info);
    
    // order eigenvalues and eigenvectors from largest to smallest
    for (int j=0; j<ns; j++)
        for (int i=0; i<ns; i++)        
            Phi[i + ns*(ns-j-1)] = A[i + ns*j];        

    for (int i=0; i<ns; i++)        
        Lambda[(ns-i-1)] = b[i];        
        
//     /* Print eigenvalues */
//     print_matrix( "Eigenvalues", 1, ns, Lambda, 1 );
//     
//     /* Print eigenvectors */
//     print_matrix( "Eigenvectors (stored columnwise)", ns, ns, Phi, ns);
    
    DGEMM(&chn, &chn, &N, &ns, &ns, &alpha, S, &N, Phi, &ns, &beta, Q, &N);        
    for (int i=0; i<(N-1); i++)
        xij[i] = xij[i+1] - xij[i];
    double area;
    for (int m=0; m<ns; m++) {
        area = 0.0;
        for (int i=0; i<(N-1); i++)            
            area += 0.5*xij[i]*(Q[i + N*m]*Q[i + N*m] + Q[i+1 + N*m]*Q[i+1 + N*m]);
        for (int i=0; i<ns; i++)            
            Phi[i + ns*m] = Phi[i + ns*m]/sqrt(area);
    }            
        
    free(xij); free(S); free(A); free(b); free(Q);
}

void read_pod(podstruct &pod, std::string pod_file)
{
    pod.nbesselpars = 3;
    pod.besselparams = (double *) malloc(3*sizeof(double));
    pod.pbc = (int *) malloc(3*sizeof(int));
    
    std::ifstream file_in(pod_file);
    if (!file_in) {error("Error: POD input file is not found");}
        
    std::string line;        
    while (std::getline(file_in, line)) // Read next line to `line`, stop if no more lines.
    {                                            
        if (line != "") {
            std::string s;
            double d;
            
            std::istringstream ss_line(line);                                    
            ss_line >> s;
            if (s == "species") {
                std::string element;
                while(ss_line >> element){                                        
                    pod.species.push_back(element);                    
                    pod.nelements += 1;
                }                        
            }
            
            if (s == "pbc") {
                int i, j = 0;                
                while(ss_line >> i) {                                        
                    pod.pbc[j++] = i;                    
                }                        
            }
            
            if ((s != "#") && (s != "species") && (s != "pbc")) {
                ss_line >> d;
                
                if (s == "rin") pod.rin = d;
                if (s == "rcut") pod.rcut = d;
                if (s == "bessel_scaling_parameter1") pod.besselparams[0] = d;
                if (s == "bessel_scaling_parameter2") pod.besselparams[1] = d;
                if (s == "bessel_scaling_parameter3") pod.besselparams[2] = d;                
                if (s == "onebody") pod.onebody = (int) d;
                if (s == "twobody_bessel_polynomial_degree") pod.twobody[0] = (int) d;
                if (s == "twobody_inverse_polynomial_degree") pod.twobody[1] = (int) d;
                if (s == "twobody_number_radial_basis_functions") pod.twobody[2] = (int) d;
                if (s == "threebody_bessel_polynomial_degree") pod.threebody[0] = (int) d;
                if (s == "threebody_inverse_polynomial_degree") pod.threebody[1] = (int) d;
                if (s == "threebody_number_radial_basis_functions") pod.threebody[2] = (int) d;
                if (s == "threebody_number_angular_basis_functions") pod.threebody[3] = (int) d;
                if (s == "fourbody_bessel_polynomial_degree") pod.fourbody[0] = (int) d;
                if (s == "fourbody_inverse_polynomial_degree") pod.fourbody[1] = (int) d;
                if (s == "fourbody_number_radial_basis_functions") pod.fourbody[2] = (int) d;
                if (s == "fourbody_number_spherical_harmonic_basis_functions") pod.fourbody[3] = (int) d;
                if (s == "quadratic22_number_twobody_basis_functions") pod.quadratic22[0] = (int) d;
                if (s == "quadratic22_number_twobody_basis_functions") pod.quadratic22[1] = (int) d;
                if (s == "quadratic23_number_twobody_basis_functions") pod.quadratic23[0] = (int) d;
                if (s == "quadratic23_number_threebody_basis_functions") pod.quadratic23[1] = (int) d;
                if (s == "quadratic24_number_twobody_basis_functions") pod.quadratic24[0] = (int) d;
                if (s == "quadratic24_number_fourbody_basis_functions") pod.quadratic24[1] = (int) d;
                if (s == "quadratic33_number_threebody_basis_functions") pod.quadratic33[0] = (int) d;
                if (s == "quadratic33_number_threebody_basis_functions") pod.quadratic33[1] = (int) d;
                if (s == "quadratic34_number_threebody_basis_functions") pod.quadratic34[0] = (int) d;
                if (s == "quadratic34_number_fourbody_basis_functions") pod.quadratic34[1] = (int) d;
                if (s == "quadratic44_number_fourbody_basis_functions") pod.quadratic44[0] = (int) d;
                if (s == "quadratic44_number_fourbody_basis_functions") pod.quadratic44[1] = (int) d;
            }
        }        
    }          
    file_in.close();
    
    // number of snapshots
    pod.ns2 = pod.nbesselpars*pod.twobody[0] + pod.twobody[1];
    pod.ns3 = pod.nbesselpars*pod.threebody[0] + pod.threebody[1];
    pod.ns4 = pod.nbesselpars*pod.fourbody[0] + pod.fourbody[1];
    
    // allocate memory for eigenvectors and eigenvalues
    pod.Phi2 = (double *) malloc(pod.ns2*pod.ns2*sizeof(double));
    pod.Lambda2 = (double *) malloc(pod.ns2*sizeof(double));
    pod.Phi3 = (double *) malloc(pod.ns3*pod.ns3*sizeof(double));
    pod.Lambda3 = (double *) malloc(pod.ns3*sizeof(double));
    pod.Phi4 = (double *) malloc(pod.ns4*pod.ns4*sizeof(double));
    pod.Lambda4 = (double *) malloc(pod.ns4*sizeof(double));    
    
    if (pod.ns2>0) {
        eigenvaluedecomposition(pod.Phi2, pod.Lambda2, pod.besselparams, pod.rin, pod.rcut, 
            pod.twobody[0], pod.twobody[1], pod.nbesselpars, 2000);                
            
//         /* Print eigenvalues */
//         print_matrix( "Eigenvalues for two-body potential:", 1, pod.ns2, pod.Lambda2, 1 );
// 
//         /* Print eigenvectors */
//         print_matrix( "Eigenvectors for two-body potential:", pod.ns2, pod.ns2, pod.Phi2, pod.ns2);        
    }
    if (pod.ns3>0) {
        eigenvaluedecomposition(pod.Phi3, pod.Lambda3, pod.besselparams, pod.rin, pod.rcut, 
            pod.threebody[0], pod.threebody[1], pod.nbesselpars, 2000);
        
//         /* Print eigenvalues */
//         print_matrix( "Eigenvalues for three-body potential:", 1, pod.ns3, pod.Lambda3, 1 );
// 
//         /* Print eigenvectors */
//         print_matrix( "Eigenvectors for three-body potential:", pod.ns3, pod.ns3, pod.Phi3, pod.ns3);        
    }
    if (pod.ns4>0) {
        eigenvaluedecomposition(pod.Phi4, pod.Lambda4, pod.besselparams, pod.rin, pod.rcut, 
            pod.fourbody[0], pod.fourbody[1], pod.nbesselpars, 2000);        
//         /* Print eigenvalues */
//         print_matrix( "Eigenvalues for four-body potential:", 1, pod.ns4, pod.Lambda4, 1 );
// 
//         /* Print eigenvectors */
//         print_matrix( "Eigenvectors for four-body potential:", pod.ns4, pod.ns4, pod.Phi4, pod.ns4);        
    }
    
    // number of chemical combinations
    pod.nc2 = pod.nelements*(pod.nelements+1)/2;
    pod.nc3 = pod.nelements*pod.nelements*(pod.nelements+1)/2;
    pod.nc4 = pod.nelements*(pod.nelements+1)/2;
        
    // number of basis functions and descriptors for one-body potential
    if (pod.onebody==1) {
        pod.nbf1 = 1;
        pod.nd1 = pod.nelements;
    }
    else {
        pod.nbf1 = 0;
        pod.nd1 = 0;        
    }    
    
    // number of basis functions and descriptors for two-body potential
    pod.nbf2 = pod.twobody[2];
    pod.nd2 = pod.nbf2*pod.nc2;
    
    // number of basis functions and descriptors for three-body potential
    pod.nrbf3 = pod.threebody[2];
    pod.nabf3 = pod.threebody[3];
    pod.nbf3 = pod.nrbf3*(1 + pod.nabf3);
    pod.nd3 = pod.nbf3*pod.nc3;

    // number of basis functions and descriptors for four-body potential
    pod.nrbf4 = pod.fourbody[2];
    pod.nabf4 = pod.fourbody[3];
    pod.nbf4 = pod.nrbf4*(1 + pod.nabf4);
    pod.nd4 = pod.nbf4*pod.nc4;
        
    pod.quadratic23[0] = PODMIN(pod.quadratic23[0], pod.nbf2);
    pod.quadratic23[1] = PODMIN(pod.quadratic23[1], pod.nbf3);
    pod.quadratic33[0] = PODMIN(pod.quadratic33[0], pod.nbf3);
    pod.quadratic33[1] = PODMIN(pod.quadratic33[1], pod.nbf3);
    pod.quadratic34[0] = PODMIN(pod.quadratic34[0], pod.nbf3);
    pod.quadratic34[1] = PODMIN(pod.quadratic34[1], pod.nbf4);
    pod.quadratic44[0] = PODMIN(pod.quadratic44[0], pod.nbf4);
    pod.quadratic44[1] = PODMIN(pod.quadratic44[1], pod.nbf4);
    
    // number of descriptors for quadratic potentials
    pod.nd22 = pod.quadratic22[0]*pod.quadratic22[1]*pod.nc2*pod.nc2;
    pod.nd23 = pod.quadratic23[0]*pod.quadratic23[1]*pod.nc2*pod.nc3;
    pod.nd24 = pod.quadratic24[0]*pod.quadratic24[1]*pod.nc2*pod.nc4;
    pod.nd33 = pod.quadratic33[0]*pod.quadratic33[1]*pod.nc3*pod.nc3;
    pod.nd34 = pod.quadratic34[0]*pod.quadratic34[1]*pod.nc3*pod.nc4;
    pod.nd44 = pod.quadratic44[0]*pod.quadratic44[1]*pod.nc4*pod.nc4;
        
    // total number of descriptors for all POD potentials
    pod.nd = pod.nd1 + pod.nd2 + pod.nd3 + pod.nd4 + pod.nd22 + pod.nd23 + pod.nd24 + pod.nd33 + pod.nd34 + pod.nd44; 
            
    std::cout<<"**************** Begin of POD Potentials ****************"<<std::endl;
    std::cout<<"species: ";
    for (int i=0; i<pod.nelements; i++)
        std::cout<<pod.species[i]<<" ";
    std::cout<<std::endl;      
    std::cout<<"periodic boundary conditions:  "<<pod.pbc[0]<<"  "<<pod.pbc[1]<<"  "<<pod.pbc[2]<<std::endl;
    std::cout<<"inner cut-off radius: "<<pod.rin<<std::endl;
    std::cout<<"outer cut-off radius: "<<pod.rcut<<std::endl;
    std::cout<<"bessel parameters: "<<pod.besselparams[0]<<"  "<<pod.besselparams[1]<<"  "<<pod.besselparams[2]<<std::endl;
    std::cout<<"one-body potential: "<<pod.onebody<<std::endl;
    std::cout<<"two-body potential: "<<pod.twobody[0]<<"  "<<pod.twobody[1]<<"  "<<pod.twobody[2]<<std::endl;
    std::cout<<"three-body potential: "<<pod.threebody[0]<<"  "<<pod.threebody[1]<<"  "<<pod.threebody[2]<<"  "<<pod.threebody[3]<<std::endl;
    std::cout<<"four-body potential: "<<pod.fourbody[0]<<"  "<<pod.fourbody[1]<<"  "<<pod.fourbody[2]<<"  "<<pod.fourbody[3]<<std::endl;
    std::cout<<"three-body quadratic22 potential: "<<pod.quadratic22[0]<<"  "<<pod.quadratic22[1]<<std::endl;
    std::cout<<"four-body quadratic23 potential: "<<pod.quadratic23[0]<<"  "<<pod.quadratic23[1]<<std::endl;
    std::cout<<"five-body quadratic24 potential: "<<pod.quadratic24[0]<<"  "<<pod.quadratic24[1]<<std::endl;
    std::cout<<"five-body quadratic33 potential: "<<pod.quadratic33[0]<<"  "<<pod.quadratic33[1]<<std::endl;
    std::cout<<"six-body quadratic34 potential: "<<pod.quadratic34[0]<<"  "<<pod.quadratic34[1]<<std::endl;
    std::cout<<"seven-body quadratic44 potential: "<<pod.quadratic44[0]<<"  "<<pod.quadratic44[1]<<std::endl;    
    std::cout<<"number of snapshots for two-body potential: "<<pod.ns2<<std::endl;
    std::cout<<"number of snapshots for three-body potential: "<<pod.ns3<<std::endl;
    std::cout<<"number of snapshots for four-body potential: "<<pod.ns4<<std::endl;    
    std::cout<<"number of basis functions for one-body potential: "<<pod.nbf1<<std::endl;
    std::cout<<"number of basis functions for two-body potential: "<<pod.nbf2<<std::endl;
    std::cout<<"number of basis functions for three-body potential: "<<pod.nbf3<<std::endl;
    std::cout<<"number of basis functions for four-body potential: "<<pod.nbf4<<std::endl;
    std::cout<<"number of descriptors for one-body potential: "<<pod.nd1<<std::endl;
    std::cout<<"number of descriptors for two-body potential: "<<pod.nd2<<std::endl;
    std::cout<<"number of descriptors for three-body potential: "<<pod.nd3<<std::endl;
    std::cout<<"number of descriptors for four-body potential: "<<pod.nd4<<std::endl;
    std::cout<<"number of descriptors for three-body quadratic22 potential: "<<pod.nd22<<std::endl;
    std::cout<<"number of descriptors for four-body quadratic23 potential: "<<pod.nd23<<std::endl;
    std::cout<<"number of descriptors for five-body quadratic24 potential: "<<pod.nd24<<std::endl;
    std::cout<<"number of descriptors for five-body quadratic33 potential: "<<pod.nd33<<std::endl;
    std::cout<<"number of descriptors for six-body quadratic34 potential: "<<pod.nd34<<std::endl;
    std::cout<<"number of descriptors for seven-body quadratic44 potential: "<<pod.nd44<<std::endl;
    std::cout<<"total number of descriptors for all POD potentials: "<<pod.nd<<std::endl;    
    std::cout<<"**************** End of POD Potentials ****************"<<std::endl<<std::endl;
}

bool is_a_number(std::string line)
{    
    return isdigit(line.at(0));
}

void read_data_file(double *fitting_weights, std::string &file_format, std::string &file_extension, 
        std::string &test_path, std::string &training_path, std::string data_file)
{
    std::ifstream file_in(data_file);
    if (!file_in) {error("Error: Data input file is not found");}
    
    std::string line;        
    while (std::getline(file_in, line)) // Read next line to `line`, stop if no more lines.
    {                                            
        if (line != "") {
            std::string s, s2;
            double d;
            
            std::istringstream ss_line(line);                                    
            ss_line >> s;
                        
            if (s == "fitting_weight_energy") {
                ss_line >> d;           
                fitting_weights[0] = d;               
            }
            else if (s == "fitting_weight_force") {
                ss_line >> d;           
                fitting_weights[1] = d;               
            }
            else if (s == "fitting_weight_stress") {
                ss_line >> d;           
                fitting_weights[2] = d;               
            }
            else if (s != "#") {
                ss_line >> s2;                
                if (s == "file_format") file_format = s2;
                if (s == "file_extension") file_extension = s2;
                if (s == "path_to_training_data_set") {
                    training_path = s2;                    
                    while(ss_line >> s2){    
                        training_path = training_path + " " + s2;   
                    }                        
                    training_path.erase(remove(training_path.begin(), training_path.end(), '"'), training_path.end());                    
                }                
                if (s == "path_to_test_data_set") {
                    test_path = s2;    
                    while (ss_line >> s2) {             
                        test_path = test_path + " " + s2;
                    }
                    test_path.erase(remove(test_path.begin(), test_path.end(), '"'), test_path.end());
                }
            }
        }
    }        
    file_in.close();
    
    std::cout<<"**************** Begin of Data File ****************"<<std::endl;
    std::cout<<"file format: "<<file_format<<std::endl;
    std::cout<<"file extension: "<<file_extension<<std::endl;
    std::cout<<"path to training data set: "<<training_path<<std::endl;
    std::cout<<"path to test data set: "<<test_path<<std::endl;    
    std::cout<<"fitting weight for energy: "<<fitting_weights[0]<<std::endl;
    std::cout<<"fitting weight for force: "<<fitting_weights[1]<<std::endl;
    std::cout<<"fitting weight for stress: "<<fitting_weights[2]<<std::endl;
    std::cout<<"**************** End of Data File ****************"<<std::endl<<std::endl;
}

void read_coeff_file(podstruct &pod, std::string coeff_file)
{
    std::ifstream file_in(coeff_file);
    if (!file_in) {error("Error: Coefficient input file is not found");}
    
    int ncoeff=0;
    std::string line;        
    while (std::getline(file_in, line)) // Read next line to `line`, stop if no more lines.
    {                                            
        if (line != "") {
            std::string s;
            int n;
            
            std::istringstream ss_line(line);                                    
            ss_line >> s;
                        
            if (s == "POD_coefficients:") {
                ss_line >> n;           
                ncoeff = n;   
                break;
            }
        }
    }        
    
    pod.coeff = (double *) malloc(ncoeff*sizeof(double));
    
    int k = 0;
    while (std::getline(file_in, line)) // Read next line to `line`, stop if no more lines.
    {                                            
        if (line != "") {
            double d;            
            std::istringstream ss_line(line);                                    
            ss_line >> d;                                    
            pod.coeff[k] = d;
            k += 1;
        }
    }        
    
    file_in.close();
    
    std::cout<<"**************** Begin of Coefficient File ****************"<<std::endl;    
    std::cout<<"number of POD coefficients: "<<ncoeff<<std::endl;
    print_matrix( "POD coefficient vector:", 1, ncoeff, pod.coeff, 1); 
    std::cout<<"**************** End of Coefficient File ****************"<<std::endl<<std::endl;
}

void get_exyz_files(std::vector<std::string>& files, std::string datapath, std::string extension)  
{
    int m = extension.length();
    for (const auto & entry : std::filesystem::directory_iterator(datapath.c_str())) {
        std::string filename = entry.path();        
        int n = filename.length();           
        std::string ext = filename.substr(n-m,n);   
        if (ext == extension) files.push_back(filename.c_str());                                        
    }        
}

int get_number_atom_exyz(std::vector<int>& num_atom, int& num_atom_sum, std::string file)  
{
    std::ifstream datafile(file);
    if (!datafile) {/*error*/}    
               
    std::string line; 
    int num_configs = 0;
    num_atom_sum = 0;
    while (std::getline(datafile, line)) // Read next line to `line`, stop if no more lines.
    {
        int d;
        if (is_a_number(line)) {            
            d = std::stoi(line);
            num_atom.push_back(d);                                        
            num_configs += 1;
            num_atom_sum += d;
        }        
    }    
    datafile.close();
    
    return num_configs;
}

int get_number_atoms(std::vector<int>& num_atom, std::vector<int> &num_atom_sum, std::vector<int>& num_config, std::vector<std::string> training_files)  
{
    int nfiles = training_files.size(); // number of files
    int d, n;
    
    for (int i=0; i<nfiles; i++) {
        d = get_number_atom_exyz(num_atom, n, training_files[i]);  
        num_config.push_back(d);  
        num_atom_sum.push_back(n);  
    }
        
    int num_atom_all = 0;
    for (int i=0; i<num_atom.size(); i++)
        num_atom_all += num_atom[i];
    
    return num_atom_all;
}

void read_exyz_file(double *lattice, double *stress, double *energy, double *pos, double *forces, 
        int *atomtype, std::string file, std::vector<std::string> species)  
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
        else if (!is_a_number(line)) {            
            std::string s0;
            std::istringstream ss(line);                        
            ss >> s0;
                    
            for (int ii=0; ii<ns; ii++)
                if (species[ii] == s0)
                    atomtype[nat] = ii+1;
                
            k = 0;
            while(ss >> d){            
                if (k <=2 ) pos[k + 3*nat] = d;
                if (k > 2 ) forces[k-3 + 3*nat] = d;
                k += 1;
            }          
            nat += 1;             
        }                
    }    
    datafile.close();    
}

void get_data(datastruct &data, std::vector<std::string> species)
{    
    get_exyz_files(data.data_files, data.data_path, data.file_extension);
    data.num_atom_sum = get_number_atoms(data.num_atom, data.num_atom_each_file, data.num_config, data.data_files);            
    data.num_config_sum = data.num_atom.size();
    
    std::cout<<"data file     |    number of configurations     |     number of atoms "<<std::endl;
    for (int i=0; i<data.data_files.size(); i++) {
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
        read_exyz_file(&data.lattice[9*nconfigs], &data.stress[9*nconfigs], &data.energy[nconfigs], 
                &data.position[3*natoms], &data.force[3*natoms], &data.atomtype[natoms], 
                data.data_files[i], species);          
        nconfigs += data.num_config[i];
        natoms += data.num_atom_each_file[i];
    }    
    
    int len = data.num_atom.size();
    data.num_atom_min = cpuArrayMin(&data.num_atom[0], len);    
    data.num_atom_max = cpuArrayMax(&data.num_atom[0], len);    
    data.num_atom_cumsum.resize(len+1);
    cpuCumsum(&data.num_atom_cumsum[0], &data.num_atom[0], len+1); 
    
    data.num_config_cumsum.resize(nfiles+1);
    cpuCumsum(&data.num_config_cumsum[0], &data.num_config[0], nfiles+1); 
    
    std::cout << "minimum number of atoms: " <<data.num_atom_min << std::endl;   
    std::cout << "maximum number of atoms: " <<data.num_atom_max << std::endl;   
}

void read_input_files(podstruct &pod, datastruct &traindata, datastruct &testdata, 
        std::string pod_file, std::string data_file, std::string coeff_file)
{
    // read pod input file to podstruct
    read_pod(pod, pod_file);    
        
    // read pod coefficient file to podstruct
    if (coeff_file != "")
        read_coeff_file(pod, coeff_file);    
    
    // read data input file to datastruct
    read_data_file(traindata.fitting_weights, traindata.file_format, traindata.file_extension, 
            testdata.data_path, traindata.data_path, data_file);
    
    std::cout<<"**************** Begin of Training Data Set ****************"<<std::endl;
    get_data(traindata, pod.species);
    std::cout<<"**************** End of Training Data Set ****************"<<std::endl<<std::endl;
                
    if ((testdata.data_path != "") && (testdata.data_path != traindata.data_path)) {
        testdata.training = 0;
        testdata.file_format = traindata.file_format;
        testdata.file_extension = traindata.file_extension;         
        std::cout<<"**************** Begin of Test Data Set ****************"<<std::endl;
        get_data(testdata, pod.species);
        std::cout<<"**************** End of Test Data Set ****************"<<std::endl<<std::endl;    
    }
    else {
        testdata.data_path = traindata.data_path;
    }        
}

void get_position(double *x, datastruct &data, int ci)
{
    int inum = data.num_atom[ci];    
    int start = data.num_atom_cumsum[ci];          
    cpuArrayCopy(x, &data.position[3*start], 3*inum);        
}

void get_force(double *x, datastruct &data, int ci)
{
    int inum = data.num_atom[ci];    
    int start = data.num_atom_cumsum[ci];          
    cpuArrayCopy(x, &data.force[3*start], 3*inum);        
}

#endif


