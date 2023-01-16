#ifndef __PODNEIGHBORLIST
#define __PODNEIGHBORLIST

int latticecoords(double *y, int *alist, double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx)
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


int neighborlist(int *ai, int *aj, int *numneigh, double *r, double rcutsq, int nx, int N, int dim)
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
                ai[k] = i+1;
                aj[k] = j+1;          
                k += 1;                                                  
            }
        }
        numneigh[i] = inc; 
    }
    return k; 
}

int podneighborlist(int *neighlist, int *numneigh, double *r, double rcutsq, int nx, int N, int dim)
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

int podfullneighborlist(double *y, int *alist, int *neighlist, int *numneigh, int *numneighsum, 
        double *x, double *a1, double *a2, double *a3, double rcut, int *pbc, int nx, int *nll)
{
    double rcutsq = rcut*rcut;    
    int dim = 3, nl = 0, nn = 0;
    
    // number of lattices
    nl = latticecoords(y, alist, x, a1, a2, a3, rcut, pbc, nx);        
    *nll = nl;
    int N = nx*nl;
            
    // total number of neighbors
   nn = podneighborlist(neighlist, numneigh, y, rcutsq, nx, N, dim);
    
   cpuCumsum(numneighsum, numneigh, nx+1); 
       
   return nn;
}   

int podNeighPairList(int *pairnum, int *pairnumsum, int *pairlist, double *x, double rcutsq, 
        int *neighlist, int *neighnumsum, int inum,  int dim)
{    
    int sumcount = 0;
    for (int ii=0; ii<inum; ii++) {
        int i = ii;       // atom i
        int n1 = neighnumsum[i];    
        int m = neighnumsum[i+1] - n1;
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i {
            int j = neighlist[n1 + l];         
            // distance between atom i and atom j                                    
            double xij0 = x[j*dim] - x[i*dim];  // xj - xi
            double xij1 = x[j*dim+1] - x[i*dim+1]; // xj - xi               
            double xij2 = x[j*dim+2] - x[i*dim+2]; // xj - xi               
            double dij = xij0*xij0 + xij1*xij1 + xij2*xij2;                        
            if (dij < rcutsq && dij>1e-20) {
                pairlist[count + sumcount] = j;  // atom j     
                count += 1;
            }
        }        
        pairnum[ii] = count;       
        sumcount += count; 
    }    

    cpuCumsum(pairnumsum, pairnum, inum+1); 

    return sumcount; 
};

void podNeighPairs(double *xij, double *x, int *ai, int *aj,  int *ti, int *tj, 
        int *pairlist, int *pairnumsum, int *atomtype, int *alist, int inum, int dim)
{        
    // std::cout << "(k, d) grid is\n";
    std::cout << "Intermediate values are:" << endl;
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ii;       // atom i
        int itype = atomtype[i];        
        int start = pairnumsum[ii];   
        int m = pairnumsum[ii+1] - start; // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
	    int k = start + l;                        
            int j = pairlist[k];  // atom j
            ai[k]        = i;
            aj[k]        = alist[j];          
            ti[k]        = itype;       
            tj[k]        = atomtype[alist[j]];
            std::cout << k << " " << j << " " << i << " ";
	    std::cout << "ai: " << ai[k] << ", aj: " << aj[k] << ", ti: " << ti[k] << ", tj: " << tj[k] << endl;
            for (int d=0; d<dim; d++) {
                xij[k*dim+d]   = x[j*dim+d] -  x[i*dim+d];  // xj - xi            
		        // std::cout << xij[k*dim+d] << " ";
		std::cout << "k: " << k << " - " << xij[k*dim+d] << " j: " << j << " - " << x[j*dim+d] << " i: " << i << " - " << x[i*dim+d] << "\n";
            }
	    // std::cout << "\n";
        }
    }    
};

void podTripletnum(int* tripletnum, int* pairnum, int length) 
{	
    // pairnum = Int32.(pairnumsum[2:end]-pairnumsum[1:end-1]);
    // tripletnum  = Int32.((pairnum.-1).*pairnum/2);
	for (int ii = 0; ii < length; ++ii)	
		tripletnum[ii] = (pairnum[ii]-1)*pairnum[ii]/2;       
}

void podNeighTripletList(int *tripletlist, int *tripletnum, int *tripletnumsum, int *pairlist, 
     int *pairnum, int *pairnumsum, int inum)
{        
	for (int ii = 0; ii < inum; ++ii)	
		tripletnum[ii] = (pairnum[ii]-1)*pairnum[ii]/2;       

    cpuCumsum(tripletnumsum, tripletnum, inum+1); 

    int Nijk = tripletnumsum[inum];

    for (int ii=0; ii<inum; ii++) {
        //int i = ii;       // atom i
        int p = pairnum[ii];      // number of pairs (i,j) around i         
        int q = tripletnumsum[ii];
        int s = pairnumsum[ii];
        int count = 0;
        for (int lj=0; lj<p ; lj++) {   // loop over each atom j around atom i
            int gj = pairlist[lj + s];  // atom j           
            for (int lk=lj+1; lk<p; lk++) { // loop over each atom k around atom i (k > j)
                int gk = pairlist[lk + s];  // atom k
                tripletlist[(count + q)] = gj;
                tripletlist[(count + q) + Nijk] = gk;                
                count += 1;
            }
        }                        
    }
}

void podNeighTriplets(double *xij, double *xik, double *x, int *ai, int *aj, int *ak,  
      int *ti, int *tj, int *tk, int *tripletlist, int *tripletnumsum, 
      int *alist,  int *atomtype, int inum, int dim)
{        
    int Nijk = tripletnumsum[inum];
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ii;       // atom i
        int itype = atomtype[i];        
        int start = tripletnumsum[ii];   
        int m = tripletnumsum[ii+1]-tripletnumsum[ii];        // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom pair (j,k) around atom i
            // int gj = tripletlist[0 + 2*(l + start)];  // ghost index of atom j  
            // int gk = tripletlist[1 + 2*(l + start)];  // ghost index of atom k  
            int gj = tripletlist[(l + start)];  // ghost index of atom j  
            int gk = tripletlist[(l + start) + Nijk];  // ghost index of atom k  
            int j = alist[gj];  // atom j
            int k = alist[gk];  // atom k
            int n = start + l;       
            ai[n]        = i;
            aj[n]        = j;    
            ak[n]        = k;    
            ti[n]        = itype;       
            tj[n]        = atomtype[j];     
            tk[n]        = atomtype[k];     
            for (int d=0; d<dim; d++) {
                xij[n*dim+d]   = x[gj*dim+d] - x[i*dim+d];  // xj - xi  
                xik[n*dim+d]   = x[gk*dim+d] - x[i*dim+d];  // xk - xi  
            }
        }
    }    
}


#endif

