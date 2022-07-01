void CPOD::linear_descriptors(double *gd, double *efatom, double *y, double *tmpmem, int *atomtype, 
            int *alist, int *pairlist, int *pairnum, int *pairnumsum, int *tmpint, int natom, int Nij)
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
    int nd1234 = nd1+nd2+nd3+nd4;    
    int *pdegree2 = pod.twobody;
    int *elemindex = pod.elemindex;
    double rin = pod.rin;
    double rcut = pod.rcut;
    double *Phi2 = pod.Phi2;
    double *besselparams = pod.besselparams;        
    
    double *fatom1 = &efatom[0];
    double *fatom2 = &efatom[dim*natom*(nd1)];
    double *fatom3 = &efatom[dim*natom*(nd1+nd2)];    
    double *fatom4 = &efatom[dim*natom*(nd1+nd2+nd3)];    
    double *eatom1 = &efatom[dim*natom*(nd1+nd2+nd3+nd4)];
    double *eatom2 = &efatom[dim*natom*(nd1+nd2+nd3+nd4)+natom*nd1];
    double *eatom3 = &efatom[dim*natom*(nd1+nd2+nd3+nd4)+natom*(nd1+nd2)];
    double *eatom4 = &efatom[dim*natom*(nd1+nd2+nd3+nd4)+natom*(nd1+nd2+nd3)];    
        
    podArraySetValue(fatom1, 0.0, (1+dim)*natom*(nd1+nd2+nd3+nd4));    
    
//     // peratom descriptors for one-body, two-body, and three-body linear potentials
//     this->poddesc(eatom1, fatom1, eatom2, fatom2, eatom3, fatom3, y, Phi2, besselparams, 
//             tmpmem, rin, rcut, atomtype, alist, pairlist, pairnum, pairnumsum, 
//             elemindex, pdegree2, tmpint, nbesselpars, nrbf2, nrbf3, nabf3, 
//             nelements, Nij, natom);                    
//     
//     if (pod.snaptwojmax>0) 
//         this->snapdesc(eatom4, fatom4, y, tmpmem, atomtype, alist, 
//                 pairlist, pairnumsum, tmpint, natom, Nij);            

    double *rij = &tmpmem[0]; // 3*Nij
    int *ai = &tmpint[0];     // Nij
    int *aj = &tmpint[Nij];   // Nij 
    int *ti = &tmpint[2*Nij]; // Nij
    int *tj = &tmpint[3*Nij]; // Nij
    this->podNeighPairs(rij, y, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, 
                alist, natom, dim);
    
    // peratom descriptors for one-body, two-body, and three-body linear potentials
    this->poddesc(eatom1, fatom1, eatom2, fatom2, eatom3, fatom3, rij, Phi2, besselparams, 
            &tmpmem[3*Nij], rin, rcut, pairnumsum, atomtype, ai, aj, ti, tj, elemindex, pdegree2, 
            nbesselpars, nrbf2, nrbf3, nabf3, nelements, Nij, natom);                    
    
    if (pod.snaptwojmax>0) 
        this->snapdesc(eatom4, fatom4, rij, &tmpmem[3*Nij], atomtype, ai, aj, ti, tj, natom, Nij);            
    
    // global descriptors for one-body, two-body, three-body, and four-bodt linear potentials    
    podArraySetValue(tmpmem, 1.0, natom);
    
    char cht = 'T';
    double one = 1.0, zero = 0.0;    
    int inc1 = 1;
    DGEMV(&cht, &natom, &nd1234, &one, eatom1, &natom, tmpmem, &inc1, &zero, gd, &inc1);            
}

void CPOD::quadratic_descriptors(double* d23, double *dd23, double* d2, double *d3, double* dd2, double *dd3, 
        int M2, int M3, int N)
{
    for (int m3 = 0; m3<M3; m3++)
        for (int m2 = 0; m2<M2; m2++)
        {
            int m = m2 + M2*m3;
            d23[m] = d2[m2]*d3[m3];                
            for (int n=0; n<N; n++)
                dd23[n + N*m] = d2[m2]*dd3[n + N*m3] + dd2[n + N*m2]*d3[m3];
        }
}

void CPOD::quadratic_descriptors(double* d33, double *dd33, double *d3, double *dd3, int M3, int N)
{
    int m = 0;
    for (int m3 = 0; m3<M3; m3++)
        for (int m2 = m3; m2<M3; m2++)
        {            
            d33[m] = d3[m2]*d3[m3];                
            for (int n=0; n<N; n++)
                dd33[n + N*m] = d3[m2]*dd3[n + N*m3] + dd3[n + N*m2]*d3[m3];
            m += 1;
        }
}

void CPOD::cubic_descriptors(double* d234, double *dd234, double* d2, double *d3, double *d4, 
        double* dd2, double *dd3, double *dd4, int M2, int M3, int M4, int N)
{
    for (int m4 = 0; m4<M4; m4++)
        for (int m3 = 0; m3<M3; m3++)
            for (int m2 = 0; m2<M2; m2++)
            {
                int m = m2 + M2*m3 + M2*M3*m4;
                d234[m] = d2[m2]*d3[m3]*d4[m4];                
                for (int n=0; n<N; n++)
                    dd234[n + N*m] = d2[m2]*d3[m3]*dd4[n + N*m4] + 
                                     d2[m2]*dd3[n + N*m3]*d4[m4] + 
                                     dd2[n + N*m2]*d3[m3]*d4[m4];
            }
}

void CPOD::cubic_descriptors(double* d333, double *Dd333, double *d3, double *Dd3, int M3, int N)
{
    int m = 0;
    for (int m3 = 0; m3<M3; m3++)
        for (int m2 = m3; m2<M3; m2++)
            for (int m1 = m2; m1<M3; m1++)
            {            
                d333[m] = d3[m1]*d3[m2]*d3[m3];                
                for (int n=0; n<N; n++)
                    Dd333[n + N*m] = d3[m1]*d3[m2]*Dd3[n + N*m3] + d3[m1]*Dd3[n + N*m2]*d3[m3] + Dd3[n + N*m1]*d3[m2]*d3[m3];
                m += 1;
            }
}

double CPOD::quadratic_coefficients(double *c2, double *c3, double *d2, double *d3, 
        double *coeff23, int *quadratic, int nc2, int nc3)
{        
    int nd2 = quadratic[0]*nc2;
    int nd3 = quadratic[1]*nc3;
        
    double energy = 0.0;    
    int m = 0;
    for (int j=0; j< nd3; j++)
        for (int k=0; k< nd2; k++) {
            energy += coeff23[m]*d3[j]*d2[k];  
            c2[k] += coeff23[m]*d3[j];
            c3[j] += coeff23[m]*d2[k];
            m += 1;        
        }
            
    return energy;
}

double CPOD::quadratic_coefficients(double *c3, double *d3, double *coeff33, 
        int *quadratic, int nc3)
{        
    int nd3 = quadratic[0]*nc3;
        
    double energy = 0.0;    
    int m = 0;
    for (int j=0; j< nd3; j++)
        for (int k=j; k< nd3; k++) {
            energy += coeff33[m]*d3[j]*d3[k];  
            c3[k] += coeff33[m]*d3[j];
            c3[j] += coeff33[m]*d3[k];
            m += 1;        
        }
            
    return energy;
}

double CPOD::cubic_coefficients(double *c2, double *c3, double *c4, double *d2, double *d3, double *d4, 
        double *coeff234, int *cubic, int nc2, int nc3, int nc4)
{        
    int nd2 = cubic[0]*nc2;
    int nd3 = cubic[1]*nc3;
    int nd4 = cubic[2]*nc4;
    
    double energy = 0.0;    
    int m = 0;
    for (int i=0; i< nd4; i++)
        for (int j=0; j< nd3; j++)
            for (int k=0; k< nd2; k++) {
                energy += coeff234[m]*d4[i]*d3[j]*d2[k];  
                c2[k] += coeff234[m]*d4[i]*d3[j];
                c3[j] += coeff234[m]*d4[i]*d2[k];
                c4[i] += coeff234[m]*d3[j]*d2[k];                
                m += 1;        
            }
            
    return energy;
}

double CPOD::cubic_coefficients(double *c3, double *d3, double *coeff333, int *cubic, int nc3)
{        
    int nd3 = cubic[0]*nc3;
        
    double energy = 0.0;
    
    int m = 0;
    for (int i=0; i< nd3; i++)
        for (int j=i; j< nd3; j++)
            for (int k=j; k< nd3; k++) {
                energy += coeff333[m]*d3[i]*d3[j]*d3[k];  
                c3[k] += coeff333[m]*d3[i]*d3[j];
                c3[j] += coeff333[m]*d3[i]*d3[k];
                c3[i] += coeff333[m]*d3[j]*d3[k];                
                m += 1;        
            }
                                
    return energy;
}

double CPOD::calculate_energyforce(double *force, double *gd, double *gdd, double *coeff, double *tmp, int natom)
{        
    int dim = 3;    
    int nforce = dim*natom;
    int nd1 = pod.nd1;
    int nd2 = pod.nd2;
    int nd3 = pod.nd3;
    int nd4 = pod.nd4;
    int nd1234 = nd1+nd2+nd3+nd4;        
    int nd22 = pod.nd22;
    int nd23 = pod.nd23;
    int nd24 = pod.nd24;
    int nd33 = pod.nd33;
    int nd34 = pod.nd34;
    int nd44 = pod.nd44;
    int nd234 = pod.nd234;
    int nd333 = pod.nd333;
    int nd444 = pod.nd444;    
    int nc2 = pod.nc2;
    int nc3 = pod.nc3;
    int nc4 = pod.nc4;
        
    // two-body, three-body, and four-body descriptors
    double *d2 = &gd[nd1];
    double *d3 = &gd[nd1+nd2];
    double *d4 = &gd[nd1+nd2+nd3];
        
    // quadratic and cubic POD coefficients
    double *coeff22 = &coeff[nd1234];
    double *coeff23 = &coeff[nd1234+nd22];
    double *coeff24 = &coeff[nd1234+nd22+nd23];
    double *coeff33 = &coeff[nd1234+nd22+nd23+nd24];
    double *coeff34 = &coeff[nd1234+nd22+nd23+nd24+nd33];
    double *coeff44 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34];  
    double *coeff234 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44];  
    double *coeff333 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234];  
    double *coeff444 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34+nd44+nd234+nd333];  
            
    // effective POD coefficients for calculating force 
    double *c1 = &tmp[0];
    double *c2 = &tmp[nd1];  
    double *c3 = &tmp[nd1+nd2];
    double *c4 = &tmp[nd1+nd2+nd3];
        
    // calculate energy for linear potentials
    double energy = 0.0;
    for (int i=0; i< nd1234; i++) {
        c1[i] = 0.0;    
        energy += coeff[i]*gd[i];    
    }
    
    // calculate energy for quadratic22 potential
    if (nd22>0) energy += this->quadratic_coefficients(c2, d2, coeff22, pod.quadratic22, nc2);
    
    // calculate energy for quadratic23 potential
    if (nd23>0) energy += this->quadratic_coefficients(c2, c3, d2, d3, coeff23, pod.quadratic23, nc2, nc3);
    
    // calculate energy for quadratic24 potential
    if (nd24>0) energy += this->quadratic_coefficients(c2, c4, d2, d4, coeff24, pod.quadratic24, nc2, nc4);
    
    // calculate energy for quadratic33 potential
    if (nd33>0) energy += this->quadratic_coefficients(c3, d3, coeff33, pod.quadratic33, nc3);
    
    // calculate energy for quadratic34 potential
    if (nd34>0) energy += this->quadratic_coefficients(c3, c4, d3, d4, coeff34, pod.quadratic34, nc3, nc4);
    
    // calculate energy for quadratic44 potential
    if (nd44>0) energy += this->quadratic_coefficients(c4, d4, coeff44, pod.quadratic44, nc4);
    
    // calculate energy for cubic234 potential
    if (nd234>0) energy += this->cubic_coefficients(c2, c3, c4, d2, d3, d4, coeff234, pod.cubic234, nc2, nc3, nc4);
    
    // calculate energy for cubic333 potential
    if (nd333>0) energy += this->cubic_coefficients(c3, d3, coeff333, pod.cubic333, nc3);

    // calculate energy for cubic444 potential
    if (nd444>0) energy += this->cubic_coefficients(c4, d4, coeff444, pod.cubic444, nc4);
    
    // calculate effective POD coefficients
    for (int i=0; i< nd1234; i++) c1[i] += coeff[i];    
    
    // calculate force = gdd * c1
    char chn = 'N';
    double one = 1.0, zero = 0.0;    
    int inc1 = 1;
    DGEMV(&chn, &nforce, &nd1234, &one, gdd, &nforce, c1, &inc1, &zero, force, &inc1);        
        
    return energy;
}

double CPOD::energyforce_calculation(double *force, double *gd, double *gdd, double *coeff, double *y, 
    int *atomtype, int *alist, int *pairlist, int *pairnum, int *pairnumsum, int *tmpint, int natom, int Nij)         
{
    int dim = 3;
    int nd1234 = pod.nd1+pod.nd2+pod.nd3+pod.nd4;        
    double *tmpmem = &gdd[dim*natom*nd1234+natom*nd1234];
    
    // calculate POD and SNAP descriptors and their derivatives
    this->linear_descriptors(gd, gdd, y, tmpmem, atomtype, alist, 
            pairlist, pairnum, pairnumsum, tmpint, natom, Nij);                
    
    // calculate energy and force
    double energy = 0.0;
    energy = this->calculate_energyforce(force, gd, gdd, coeff, &gdd[dim*natom*nd1234], natom);
            
    return energy;
}

void CPOD::podNeighPairs(double *xij, double *x, int *ai, int *aj,  int *ti, int *tj, 
        int *pairlist, int *pairnumsum, int *atomtype, int *alist, int inum, int dim)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ii;       // atom i
        int itype = atomtype[i];        
        int start = pairnumsum[ii];   
        int m = pairnumsum[ii+1] - start; // number of neighbors around i             
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = pairlist[l + start];  // atom j              
            int k = start + l;                                     
            ai[k]        = i;
            aj[k]        = alist[j];          
            ti[k]        = itype;       
            tj[k]        = atomtype[alist[j]];        
            for (int d=0; d<dim; d++) 
                xij[k*dim+d]   = x[j*dim+d] -  x[i*dim+d];  // xj - xi            
        }
    }    
};

void CPOD::podradialbasis(double *rbf, double *drbf, double *xij, double *besselparams, double rin, 
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N)
{
    for (int n=0; n<N; n++) {    
        double xij1 = xij[0+3*n];
        double xij2 = xij[1+3*n];
        double xij3 = xij[2+3*n];

        double dij = pow(xij1*xij1 + xij2*xij2 + xij3*xij3, 0.5);    
        double dr1 = xij1/dij;    
        double dr2 = xij2/dij;    
        double dr3 = xij3/dij;    

        double r = dij - rin;        
        double y = r/rmax;    
        double y2 = y*y;
        double y3 = 1.0 - y2*y;
        double y4 = y3*y3 + 1e-6;
        double y5 = pow(y4, 0.5);
        double y6 = exp(-1.0/y5);
        double y7 = pow(y4, 1.5);
        double fcut = y6/exp(-1.0);
        double dfcut = ((3.0/(rmax*exp(-1.0)))*(y2)*y6*(y*y2 - 1.0))/y7;

        for (int j=0; j<nbesselpars; j++) {            
            double alpha = besselparams[j];    
            if (fabs(alpha) <= 1.0e-6) alpha = 1e-3;                        
            double x =  (1.0 - exp(-alpha*r/rmax))/(1.0-exp(-alpha));
            double dx = (alpha/rmax)*exp(-(alpha*r/rmax))/(1.0 - exp(-alpha));        

            for (int i=0; i<besseldegree; i++) {
                double a = (i+1)*M_PI;
                double b = (sqrt(2.0/(rmax))/(i+1));
                int nij = n + N*i + N*besseldegree*j;            
                rbf[nij] = b*fcut*sin(a*x)/r;
                double drbfdr = b*(dfcut*sin(a*x)/r - fcut*sin(a*x)/(r*r) + a*cos(a*x)*fcut*dx/r);
                drbf[0 + 3*nij] = drbfdr*dr1;
                drbf[1 + 3*nij] = drbfdr*dr2;
                drbf[2 + 3*nij] = drbfdr*dr3;
            }
        }

        for (int i=0; i<inversedegree; i++) {
            int p = besseldegree*nbesselpars + i;
            int nij = n + N*p;     
            double a = pow(dij, (double) (i+1.0));
            rbf[nij] = fcut/a;
            double drbfdr = dfcut/a - (i+1.0)*fcut/(a*dij);  
            drbf[0 + 3*nij] = drbfdr*dr1;
            drbf[1 + 3*nij] = drbfdr*dr2;
            drbf[2 + 3*nij] = drbfdr*dr3;
        }
    }       
}

void CPOD::podtally2b(double *eatom, double *fatom, double *eij, double *fij, int *ai, int *aj, 
        int *ti, int *tj, int *elemindex, int nelements, int nbf, int natom, int N)
{
    int nelements2 = nelements*(nelements+1)/2;
    for (int n=0; n<N; n++) {
        int i1 = ai[n];
        int j1 = aj[n];
        int typei = ti[n]-1;
        int typej = tj[n]-1;
        //int mij = (elemindex[typei + typej*nelements]-1)*nbf;
        for (int m=0; m<nbf; m++) {               
            //int im = i1 + natom*(m + mij);
            //int jm = j1 + natom*(m + mij);
            int im =  i1 + natom*((elemindex[typei + typej*nelements] - 1) + nelements2*m);
            int jm =  j1 + natom*((elemindex[typei + typej*nelements] - 1) + nelements2*m);
            int nm = n + N*m;
            eatom[im] += eij[nm];
            fatom[0 + 3*im] += fij[0 + 3*nm];
            fatom[1 + 3*im] += fij[1 + 3*nm];
            fatom[2 + 3*im] += fij[2 + 3*nm];
            fatom[0 + 3*jm] -= fij[0 + 3*nm];
            fatom[1 + 3*jm] -= fij[1 + 3*nm];
            fatom[2 + 3*jm] -= fij[2 + 3*nm];          
        }
    }
}

void CPOD::pod1body(double *eatom, double *fatom, int *atomtype, int nelements, int natom)
{
    for (int m=1; m<=nelements; m++)       
        for (int i=0; i<natom; i++)         
            eatom[i + natom*(m-1)] = (atomtype[i] == m) ? 1.0 : 0.0;
        
    for (int i=0; i<3*natom*nelements; i++)  
        fatom[i] = 0.0;
}

void CPOD::pod3body(double *eatom, double *fatom, double *yij, double *e2ij, double *f2ij, double *tmpmem, 
             int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti, int *tj, int nrbf, int nabf, 
             int nelements, int natom, int Nij)
{   
    int dim = 3, nabf1 = nabf + 1;
    int nelements2 = nelements*(nelements+1)/2;
    int n, nijk, nijk3;
    
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, sinthe, theta, dtheta; 
    double tm, tm1, tm2, dct1, dct2, dct3, dct4, dct5, dct6;
    double uj, uk, rbf, drbf1, drbf2, drbf3, drbf4, drbf5, drbf6;
    double eijk, fj1, fj2, fj3, fk1, fk2, fk3;
            
    double *abf = &tmpmem[0];
    double *dabf1 = &tmpmem[nabf1];
    double *dabf2 = &tmpmem[2*nabf1];
    double *dabf3 = &tmpmem[3*nabf1];
    double *dabf4 = &tmpmem[4*nabf1];
    double *dabf5 = &tmpmem[5*nabf1];
    double *dabf6 = &tmpmem[6*nabf1];
    
    for (int ii=0; ii<natom; ii++) {
        int numneigh = pairnumsum[ii+1] - pairnumsum[ii];      // number of pairs (i,j) around i         
        int s = pairnumsum[ii];        
        for (int lj=0; lj<numneigh ; lj++) {   // loop over each atom j around atom i            
            int ij = lj + s;
            int i = ai[ij];  // atom i                        
            int j = aj[ij];  // atom j    
            int typei = ti[ij] - 1;           
            int typej = tj[ij] - 1;                   
            xij1 = yij[0+dim*ij];  // xj - xi           
            xij2 = yij[1+dim*ij];  // xj - xi           
            xij3 = yij[2+dim*ij];  // xj - xi           
            rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
            rij = pow(rijsq, 0.5);                         
            for (int lk=lj+1; lk<numneigh; lk++) { // loop over each atom k around atom i (k > j)
                int ik = lk + s;
                int k = aj[ik];  // atom k                       
                int typek = tj[ik] - 1;                         
                xik1 = yij[0+dim*ik];  // xk - xi           
                xik2 = yij[1+dim*ik];  // xk - xi           
                xik3 = yij[2+dim*ik];  // xk - xi           s
                riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;                    
                rik = pow(riksq, 0.5); 

                xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;                
                costhe = xdot/(rij*rik);    
                costhe = costhe > 1.0 ? 1.0 : costhe;
                costhe = costhe < -1.0 ? -1.0 : costhe;
                xdot = costhe*(rij*rik);

                sinthe = pow(1.0 - costhe*costhe, 0.5);
                sinthe = sinthe > 1e-12 ? sinthe : 1e-12;    
                theta = acos(costhe);            
                dtheta = -1.0/sinthe; 

                tm1 = pow(rijsq,1.5)*rik;
                tm2 = rij*pow(riksq,1.5);
                tm1 = 1.0/tm1;
                tm2 = 1.0/tm2;
                dct1 = (xik1*rijsq - xij1*xdot)*tm1; 
                dct2 = (xik2*rijsq - xij2*xdot)*tm1;
                dct3 = (xik3*rijsq - xij3*xdot)*tm1;
                dct4 = (xij1*riksq - xik1*xdot)*tm2;
                dct5 = (xij2*riksq - xik2*xdot)*tm2;
                dct6 = (xij3*riksq - xik3*xdot)*tm2;

                for (int p=0; p <nabf1; p++) {
                    abf[p] = cos(p*theta);                
                    tm = -p*sin(p*theta)*dtheta;                    
                    dabf1[p] = tm*dct1;
                    dabf2[p] = tm*dct2;
                    dabf3[p] = tm*dct3;
                    dabf4[p] = tm*dct4;
                    dabf5[p] = tm*dct5;
                    dabf6[p] = tm*dct6;        
                }

                for (int m=0; m<nrbf; m++) {
                    uj = e2ij[lj + s + Nij*m];
                    uk = e2ij[lk + s + Nij*m];
                    rbf = uj*uk;
                    drbf1 = f2ij[0 + dim*(lj + s) + dim*Nij*m]*uk;
                    drbf2 = f2ij[1 + dim*(lj + s) + dim*Nij*m]*uk;
                    drbf3 = f2ij[2 + dim*(lj + s) + dim*Nij*m]*uk;                                                        
                    drbf4 = f2ij[0 + dim*(lk + s) + dim*Nij*m]*uj;
                    drbf5 = f2ij[1 + dim*(lk + s) + dim*Nij*m]*uj;
                    drbf6 = f2ij[2 + dim*(lk + s) + dim*Nij*m]*uj;     

                    for (int p=0; p <nabf1; p++) {
                        eijk = rbf*abf[p];
                        fj1 = drbf1*abf[p] + rbf*dabf1[p];
                        fj2 = drbf2*abf[p] + rbf*dabf2[p];
                        fj3 = drbf3*abf[p] + rbf*dabf3[p];
                        fk1 = drbf4*abf[p] + rbf*dabf4[p];
                        fk2 = drbf5*abf[p] + rbf*dabf5[p];
                        fk3 = drbf6*abf[p] + rbf*dabf6[p];

                        n = p + (nabf1)*m;
                        nijk = natom*((elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*n);
                        eatom[i + nijk] += eijk;

                        nijk3 = 3*i + 3*nijk;                        
                        fatom[0 + nijk3] += fj1 + fk1;
                        fatom[1 + nijk3] += fj2 + fk2;
                        fatom[2 + nijk3] += fj3 + fk3;

                        nijk3 = 3*j + 3*nijk;    
                        fatom[0 + nijk3] -= fj1;
                        fatom[1 + nijk3] -= fj2;
                        fatom[2 + nijk3] -= fj3;

                        nijk3 = 3*k + 3*nijk;    
                        fatom[0 + nijk3] -= fk1;   
                        fatom[1 + nijk3] -= fk2;   
                        fatom[2 + nijk3] -= fk3;                           

                    }                    
                }
            }
        }
    }
}

void CPOD::poddesc(double *eatom1, double *fatom1, double *eatom2, double *fatom2, double *eatom3, 
            double *fatom3, double *rij, double *Phi, double *besselparams, double *tmpmem, double rin, 
            double rcut, int *pairnumsum, int *atomtype, int *ai, int *aj, int *ti, int *tj, int *elemindex, 
            int *pdegree, int nbesselpars, int nrbf2, int nrbf3, int nabf, int nelements, int Nij, int natom)
{       
    int nrbf = PODMAX(nrbf2, nrbf3);
    int ns = pdegree[0]*nbesselpars + pdegree[1];
    
    double *e2ij = &tmpmem[0]; // Nij*nrbf
    double *f2ij = &tmpmem[Nij*nrbf]; // dim*Nij*nrbf
    double *e2ijt = &tmpmem[4*Nij*nrbf]; // Nij*ns
    double *f2ijt = &tmpmem[4*Nij*nrbf+Nij*ns]; // dim*Nij*ns    

    // orthogonal radial basis functions
    this->podradialbasis(e2ijt, f2ijt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nij);
    this->podMatMul(e2ij, e2ijt, Phi, Nij, ns, nrbf);
    this->podMatMul(f2ij, f2ijt, Phi, 3*Nij, ns, nrbf);

    // one-body descriptors
    this->pod1body(eatom1, fatom1, atomtype, nelements, natom);

    // two-body descriptors
    this->podtally2b(eatom2, fatom2, e2ij, f2ij, ai, aj, ti, tj, elemindex, nelements, nrbf2, natom, Nij);   
    
    // three-body descriptors
    this->pod3body(eatom3, fatom3, rij, e2ij, f2ij, &tmpmem[4*Nij*nrbf], elemindex, pairnumsum, 
             ai, aj, ti, tj, nrbf3, nabf, nelements, natom, Nij);            
}
