void podsnapshots(double *rbf, double *xij, double *besselparams, double rin, double rcut, 
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
    }
}

void podeigenvaluedecomposition(double *Phi, double *Lambda, double *besselparams, double rin, double rcut, 
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
    
    podsnapshots(S, xij, besselparams, rin, rcut, besseldegree, inversedegree, nbesselpars, N);
        
    char chn = 'N';
    char cht = 'T';
    char chv = 'V';
    char chu = 'U';
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


void CPOD::read_pod(std::string pod_file)
{
    pod.nbesselpars = 3;
    pod.besselparams = (double *) malloc(3*sizeof(double));
    pod.pbc = (int *) malloc(3*sizeof(int));
    
    std::ifstream file_in(pod_file);
    if (!file_in) {pod_error("Error: POD input file is not found");}
        
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
                if (s == "bessel_polynomial_degree") pod.besseldegree = (int) d;
                if (s == "inverse_polynomial_degree") pod.inversedegree = (int) d;
                if (s == "onebody") pod.onebody = (int) d;
                if (s == "twobody_bessel_polynomial_degree") pod.twobody[0] = (int) d;
                if (s == "twobody_inverse_polynomial_degree") pod.twobody[1] = (int) d;
                if (s == "twobody_number_radial_basis_functions") pod.twobody[2] = (int) d;
                if (s == "threebody_bessel_polynomial_degree") pod.threebody[0] = (int) d;
                if (s == "threebody_inverse_polynomial_degree") pod.threebody[1] = (int) d;
                if (s == "threebody_number_radial_basis_functions") pod.threebody[2] = (int) d;
                if (s == "threebody_number_angular_basis_functions") pod.threebody[3] = (int) (d-1);
                if (s == "fourbody_bessel_polynomial_degree") pod.fourbody[0] = (int) d;
                if (s == "fourbody_inverse_polynomial_degree") pod.fourbody[1] = (int) d;
                if (s == "fourbody_number_radial_basis_functions") pod.fourbody[2] = (int) d;
                if (s == "fourbody_snap_twojmax") pod.snaptwojmax = (int) d;
                if (s == "fourbody_snap_chemflag") pod.snapchemflag = (int) d;
                if (s == "fourbody_snap_rfac0") pod.snaprfac0 = d;
                if (s == "fourbody_snap_neighbor_weight1") pod.snapelementweight[0] = d;
                if (s == "fourbody_snap_neighbor_weight2") pod.snapelementweight[1] = d;
                if (s == "fourbody_snap_neighbor_weight3") pod.snapelementweight[2] = d;
                if (s == "fourbody_snap_neighbor_weight4") pod.snapelementweight[3] = d;
                if (s == "fourbody_snap_neighbor_weight5") pod.snapelementweight[4] = d;
                //if (s == "fourbody_number_spherical_harmonic_basis_functions") pod.fourbody[3] = (int) d;
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
                if (s == "cubic234_number_twobody_basis_functions") pod.cubic234[0] = (int) d;
                if (s == "cubic234_number_threebody_basis_functions") pod.cubic234[1] = (int) d;
                if (s == "cubic234_number_fourbody_basis_functions") pod.cubic234[2] = (int) d;
                if (s == "cubic333_number_threebody_basis_functions") pod.cubic333[0] = (int) d;
                if (s == "cubic444_number_fourbody_basis_functions") pod.cubic444[0] = (int) d;
            }
        }        
    }          
    file_in.close();
        
    pod.twobody[0] = pod.besseldegree;
    pod.twobody[1] = pod.inversedegree;
    pod.threebody[0] = pod.besseldegree;
    pod.threebody[1] = pod.inversedegree;
    
    // number of snapshots
    pod.ns2 = pod.nbesselpars*pod.twobody[0] + pod.twobody[1];
    pod.ns3 = pod.nbesselpars*pod.threebody[0] + pod.threebody[1];
    pod.ns4 = pod.nbesselpars*pod.fourbody[0] + pod.fourbody[1];
    
    for (int i = 0; i < pod.nbesselpars; i++)
        if (fabs(pod.besselparams[i]) < 1e-3) pod.besselparams[i] = 1e-3;
            
    // allocate memory for eigenvectors and eigenvalues
    pod.Phi2 = (double *) malloc(pod.ns2*pod.ns2*sizeof(double));
    pod.Lambda2 = (double *) malloc(pod.ns2*sizeof(double));
    pod.Phi3 = (double *) malloc(pod.ns3*pod.ns3*sizeof(double));
    pod.Lambda3 = (double *) malloc(pod.ns3*sizeof(double));
    pod.Phi4 = (double *) malloc(pod.ns4*pod.ns4*sizeof(double));
    pod.Lambda4 = (double *) malloc(pod.ns4*sizeof(double));    
    
    if (pod.ns2>0) {
        podeigenvaluedecomposition(pod.Phi2, pod.Lambda2, pod.besselparams, pod.rin, pod.rcut, 
            pod.twobody[0], pod.twobody[1], pod.nbesselpars, 2000);                
            
//         /* Print eigenvalues */
//         print_matrix( "Eigenvalues for two-body potential:", 1, pod.ns2, pod.Lambda2, 1 );
// 
//         /* Print eigenvectors */
//         print_matrix( "Eigenvectors for two-body potential:", pod.ns2, pod.ns2, pod.Phi2, pod.ns2);        
    }
    if (pod.ns3>0) {
        podeigenvaluedecomposition(pod.Phi3, pod.Lambda3, pod.besselparams, pod.rin, pod.rcut, 
            pod.threebody[0], pod.threebody[1], pod.nbesselpars, 2000);        
    }
    if (pod.ns4>0) {
        podeigenvaluedecomposition(pod.Phi4, pod.Lambda4, pod.besselparams, pod.rin, pod.rcut, 
            pod.fourbody[0], pod.fourbody[1], pod.nbesselpars, 2000);        
    }
    
    // number of chemical combinations
    pod.nc2 = pod.nelements*(pod.nelements+1)/2;
    pod.nc3 = pod.nelements*pod.nelements*(pod.nelements+1)/2;            
    pod.nc4 = pod.snapchemflag ? pod.nelements*pod.nelements*pod.nelements*pod.nelements : pod.nelements;
            
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
    int twojmax = pod.snaptwojmax;    
    int idxb_count = 0;    
    if (twojmax > 0) {
        for(int j1 = 0; j1 <= twojmax; j1++)
            for(int j2 = 0; j2 <= j1; j2++)
                for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
                    if (j >= j1) idxb_count++;
    }
    pod.nbf4 = idxb_count;
    pod.nd4 = pod.nbf4*pod.nc4;
        
    pod.quadratic22[0] = PODMIN(pod.quadratic22[0], pod.nbf2);
    pod.quadratic22[1] = PODMIN(pod.quadratic22[1], pod.nbf2);
    pod.quadratic23[0] = PODMIN(pod.quadratic23[0], pod.nbf2);
    pod.quadratic23[1] = PODMIN(pod.quadratic23[1], pod.nbf3);
    pod.quadratic24[0] = PODMIN(pod.quadratic24[0], pod.nbf2);
    pod.quadratic24[1] = PODMIN(pod.quadratic24[1], pod.nbf4);
    pod.quadratic33[0] = PODMIN(pod.quadratic33[0], pod.nbf3);
    pod.quadratic33[1] = PODMIN(pod.quadratic33[1], pod.nbf3);
    pod.quadratic34[0] = PODMIN(pod.quadratic34[0], pod.nbf3);
    pod.quadratic34[1] = PODMIN(pod.quadratic34[1], pod.nbf4);
    pod.quadratic44[0] = PODMIN(pod.quadratic44[0], pod.nbf4);
    pod.quadratic44[1] = PODMIN(pod.quadratic44[1], pod.nbf4);
    
    pod.cubic234[0] = PODMIN(pod.cubic234[0], pod.nbf2);
    pod.cubic234[1] = PODMIN(pod.cubic234[1], pod.nbf3);
    pod.cubic234[2] = PODMIN(pod.cubic234[2], pod.nbf4);    
    pod.cubic333[0] = PODMIN(pod.cubic333[0], pod.nbf3);
    pod.cubic333[1] = PODMIN(pod.cubic333[0], pod.nbf3);
    pod.cubic333[2] = PODMIN(pod.cubic333[0], pod.nbf3);
    pod.cubic444[0] = PODMIN(pod.cubic444[0], pod.nbf4);
    pod.cubic444[1] = PODMIN(pod.cubic444[0], pod.nbf4);
    pod.cubic444[2] = PODMIN(pod.cubic444[0], pod.nbf4);
    
    // number of descriptors for quadratic POD potentials        
    pod.nd22 = pod.quadratic22[0]*pod.quadratic22[1]*pod.nc2*pod.nc2;
    pod.nd23 = pod.quadratic23[0]*pod.quadratic23[1]*pod.nc2*pod.nc3;
    pod.nd24 = pod.quadratic24[0]*pod.quadratic24[1]*pod.nc2*pod.nc4;    
    pod.nd33 = pod.quadratic33[0]*pod.quadratic33[1]*pod.nc3*pod.nc3;
    pod.nd34 = pod.quadratic34[0]*pod.quadratic34[1]*pod.nc3*pod.nc4;
    pod.nd44 = pod.quadratic44[0]*pod.quadratic44[1]*pod.nc4*pod.nc4;
    
    int nq;
    nq = pod.quadratic22[0]*pod.nc2; pod.nd22 = nq*(nq+1)/2;
    nq = pod.quadratic33[0]*pod.nc3; pod.nd33 = nq*(nq+1)/2;
    nq = pod.quadratic44[0]*pod.nc4; pod.nd44 = nq*(nq+1)/2;
    
    // number of descriptors for cubic POD potentials        
    pod.nd234 = pod.cubic234[0]*pod.cubic234[1]*pod.cubic234[2]*pod.nc2*pod.nc3*pod.nc4;
    nq = pod.cubic333[0]*pod.nc3; pod.nd333 = nq*(nq+1)*(nq+2)/6;    
    nq = pod.cubic444[0]*pod.nc4; pod.nd444 = nq*(nq+1)*(nq+2)/6;    
    
    // total number of descriptors for all POD potentials
    pod.nd = pod.nd1 + pod.nd2 + pod.nd3 + pod.nd4 + pod.nd22 + pod.nd23 + pod.nd24 + 
             pod.nd33 + pod.nd34 + pod.nd44 + pod.nd234 + pod.nd333 + pod.nd444; 
            
    int nelements = pod.nelements;
    pod.elemindex = (int*) malloc (sizeof (int)*(nelements*nelements));     
        
    int k = 1;
    for (int i=0; i < nelements; i++) 
        for (int j=i; j < nelements; j++) {
            pod.elemindex[i + nelements*j] = k;
            pod.elemindex[j + nelements*i] = k;
            k += 1;
        }          
    
    std::cout<<"**************** Begin of POD Potentials ****************"<<std::endl;
    std::cout<<"species: ";
    for (int i=0; i<pod.nelements; i++)
        std::cout<<pod.species[i]<<" ";
    std::cout<<std::endl;      
    std::cout<<"periodic boundary conditions:  "<<pod.pbc[0]<<"  "<<pod.pbc[1]<<"  "<<pod.pbc[2]<<std::endl;
    std::cout<<"inner cut-off radius: "<<pod.rin<<std::endl;
    std::cout<<"outer cut-off radius: "<<pod.rcut<<std::endl;
    std::cout<<"bessel parameters: "<<pod.besselparams[0]<<"  "<<pod.besselparams[1]<<"  "<<pod.besselparams[2]<<std::endl;
    std::cout<<"bessel polynomial degree: "<<pod.besseldegree<<std::endl;
    std::cout<<"inverse polynomial degree: "<<pod.inversedegree<<std::endl;
    std::cout<<"one-body potential: "<<pod.onebody<<std::endl;
    std::cout<<"two-body potential: "<<pod.twobody[0]<<"  "<<pod.twobody[1]<<"  "<<pod.twobody[2]<<std::endl;
    std::cout<<"three-body potential: "<<pod.threebody[0]<<"  "<<pod.threebody[1]<<"  "<<pod.threebody[2]<<"  "<<pod.threebody[3]+1<<std::endl;
    std::cout<<"four-body SNAP potential: "<<pod.snaptwojmax<<"  "<<pod.snapchemflag<<std::endl;
    std::cout<<"three-body quadratic22 potential: "<<pod.quadratic22[0]<<"  "<<pod.quadratic22[1]<<std::endl;
    std::cout<<"four-body quadratic23 potential: "<<pod.quadratic23[0]<<"  "<<pod.quadratic23[1]<<std::endl;
    std::cout<<"five-body quadratic24 potential: "<<pod.quadratic24[0]<<"  "<<pod.quadratic24[1]<<std::endl;
    std::cout<<"five-body quadratic33 potential: "<<pod.quadratic33[0]<<"  "<<pod.quadratic33[1]<<std::endl;
    std::cout<<"six-body quadratic34 potential: "<<pod.quadratic34[0]<<"  "<<pod.quadratic34[1]<<std::endl;
    std::cout<<"seven-body quadratic44 potential: "<<pod.quadratic44[0]<<"  "<<pod.quadratic44[1]<<std::endl;    
    std::cout<<"seven-body cubic234 potential: "<<pod.cubic234[0]<<"  "<<pod.cubic234[1]<<"  "<<pod.cubic234[2]<<std::endl;
    std::cout<<"seven-body cubic333 potential: "<<pod.cubic333[0]<<"  "<<pod.cubic333[1]<<"  "<<pod.cubic333[2]<<std::endl;    
    std::cout<<"ten-body cubic444 potential: "<<pod.cubic444[0]<<"  "<<pod.cubic444[1]<<"  "<<pod.cubic444[2]<<std::endl;    
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
    std::cout<<"number of descriptors for seven-body cubic234 potential: "<<pod.nd333<<std::endl;
    std::cout<<"number of descriptors for seven-body cubic333 potential: "<<pod.nd234<<std::endl;
    std::cout<<"number of descriptors for ten-body cubic444 potential: "<<pod.nd444<<std::endl;
    std::cout<<"total number of descriptors for all POD potentials: "<<pod.nd<<std::endl;    
    std::cout<<"**************** End of POD Potentials ****************"<<std::endl<<std::endl;
}

void CPOD::read_coeff_file(std::string coeff_file)
{
    std::ifstream file_in(coeff_file);
    if (!file_in) {pod_error("Error: Coefficient input file is not found");}
    
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

