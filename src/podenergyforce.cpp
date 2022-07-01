void CPOD::podNeighPairs(double *rij, double *x, int *idxi, int *ai, int *aj,  int *ti, int *tj, 
        int *pairnumsum, int *atomtype, int *jlist, int *alist, int inum)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        //int gi = ilist[ii];       // atom i
        int gi = ii;       // atom i
        int itype = atomtype[gi];
        int start = pairnumsum[ii];   
        int m = pairnumsum[ii+1] - start;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int k = start + l;
            int gj = jlist[k];  // atom j                                                               
            idxi[k]      = ii;
            ai[k]        = alist[gi];
            aj[k]        = alist[gj];          
            ti[k]        = itype;       
            tj[k]        = atomtype[aj[k]];              
            rij[k*3+0]   = x[gj*3+0] -  x[gi*3+0];  // xj - xi            
            rij[k*3+1]   = x[gj*3+1] -  x[gi*3+1];  // xj - xi            
            rij[k*3+2]   = x[gj*3+2] -  x[gi*3+2];  // xj - xi            
        }
    }            
};

int CPOD::lammpsNeighPairs(double *rij, double **x, double rcutsq, int *idxi, int *ai, int *aj,  int *ti, int *tj, 
        int *pairnumsum, int *atomtype, int *numneigh, int *ilist, int **jlist, int inum)
{  
    
    int ninside = 0;
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int gi = ilist[ii];       // atom i
        int itype = atomtype[gi];
        int m = numneigh[gi];
        pairnumsum[ii+1] = 0;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int gj = jlist[gi][l];  // atom j     
            double delx   = x[gj][0] -  x[gi][0];  // xj - xi            
            double dely   = x[gj][1] -  x[gi][1];  // xj - xi            
            double delz   = x[gj][2] -  x[gi][2];  // xj - xi            
            double rsq = delx*delx + dely*dely + delz*delz;
            if (rsq < rcutsq && rsq > 1e-20) {
                rij[ninside*3 + 0] = delx;
                rij[ninside*3 + 1] = dely;
                rij[ninside*3 + 2] = delz;
                idxi[ninside]      = ii;
                ai[ninside]        = gi;
                aj[ninside]        = gj;          
                ti[ninside]        = itype;       
                tj[ninside]        = atomtype[gj];       
                ninside++;
                pairnumsum[ii+1] += 1;
            }                        
        }
    }    
    
    pairnumsum[0] = 0;
    for (int ii=0; ii<inum; ii++)
        pairnumsum[ii+1] = pairnumsum[ii+1] + pairnumsum[ii];
    
    
    return ninside;
};

void CPOD::podradialbasis(double *rbf, double *xij, double *besselparams, double rin, 
        double rmax, int besseldegree, int inversedegree, int nbesselpars, int N)
{
    for (int n=0; n<N; n++) {    
        double xij1 = xij[0+3*n];
        double xij2 = xij[1+3*n];
        double xij3 = xij[2+3*n];

        double dij = pow(xij1*xij1 + xij2*xij2 + xij3*xij3, 0.5);    
        double r = dij - rin;        
        double y = r/rmax;    
        double y2 = y*y;
        double y3 = 1.0 - y2*y;
        double y4 = y3*y3 + 1e-6;
        double y5 = pow(y4, 0.5);
        double y6 = exp(-1.0/y5);
        double fcut = y6/exp(-1.0);

        for (int j=0; j<nbesselpars; j++) {            
            double x =  (1.0 - exp(-besselparams[j]*r/rmax))/(1.0-exp(-besselparams[j]));
            for (int i=0; i<besseldegree; i++)                 
                rbf[n + N*i + N*besseldegree*j] = ((sqrt(2.0/(rmax))/(i+1)))*fcut*sin((i+1)*M_PI*x)/r;            
        }

        for (int i=0; i<inversedegree; i++) {
            int p = besseldegree*nbesselpars + i;
            double a = pow(dij, (double) (i+1.0));
            rbf[n + N*p] = fcut/a;
        }
    }                    
}

void CPOD::podtally2b(double *eatom, double *eij, int *idxi, int *ti, int *tj, int *elemindex, 
        int nelements, int nbf, int natom, int N)
{
    int nelements2 = nelements*(nelements+1)/2;
    for (int n=0; n<N; n++) {
        int i1 = idxi[n];
        int typei = ti[n]-1;
        int typej = tj[n]-1;
        for (int m=0; m<nbf; m++) {               
            int im =  i1 + natom*((elemindex[typei + typej*nelements] - 1) + nelements2*m);
            int nm = n + N*m;
            eatom[im] += eij[nm];
        }
    }
}

void CPOD::pod1body(double *eatom, int *atomtype, int nelements, int natom)
{
    for (int m=1; m<=nelements; m++)       
        for (int i=0; i<natom; i++)         
            eatom[i + natom*(m-1)] = (atomtype[i] == m) ? 1.0 : 0.0;        
}

void CPOD::pod3body(double *eatom, double *yij, double *e2ij, double *tmpmem, int *elemindex, int *pairnumsum, 
        int *idxi, int *ti, int *tj, int nrbf, int nabf, int nelements, int natom, int Nij)
{   
    int dim = 3, nabf1 = nabf + 1;
    int nelements2 = nelements*(nelements+1)/2;
    int n, nijk;
    
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, theta; 
    double uj, uk, rbf;
            
    double *abf = &tmpmem[0];
    
    for (int ii=0; ii<natom; ii++) {
        int numneigh = pairnumsum[ii+1] - pairnumsum[ii];      // number of pairs (i,j) around i         
        int s = pairnumsum[ii];        
        for (int lj=0; lj<numneigh ; lj++) {   // loop over each atom j around atom i            
            int ij = lj + s;
            int i = idxi[ij];  // atom i                        
            int typei = ti[ij] - 1;           
            int typej = tj[ij] - 1;                   
            xij1 = yij[0+dim*ij];  // xj - xi           
            xij2 = yij[1+dim*ij];  // xj - xi           
            xij3 = yij[2+dim*ij];  // xj - xi           
            rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
            rij = pow(rijsq, 0.5);                         
            for (int lk=lj+1; lk<numneigh; lk++) { // loop over each atom k around atom i (k > j)
                int ik = lk + s;
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
                theta = acos(costhe);            

                for (int p=0; p <nabf1; p++) 
                    abf[p] = cos(p*theta);                

                for (int m=0; m<nrbf; m++) {
                    uj = e2ij[lj + s + Nij*m];
                    uk = e2ij[lk + s + Nij*m];
                    rbf = uj*uk;
                    for (int p=0; p <nabf1; p++) {                        
                        n = p + (nabf1)*m;
                        nijk = natom*((elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*n);
                        eatom[i + nijk] += rbf*abf[p];
                    }                    
                }
            }
        }
    }
}


void CPOD::poddesc_ij(double *eatom1, double *eatom2, double *eatom3, double *rij, double *Phi, double *besselparams, 
            double *tmpmem, double rin, double rcut, int *pairnumsum, int *atomtype, int *idxi, int *ti, int *tj, 
            int *elemindex, int *pdegree, int nbesselpars, int nrbf2, int nrbf3, int nabf, int nelements, int Nij, int natom)
{       
    int nrbf = PODMAX(nrbf2, nrbf3);
    int ns = pdegree[0]*nbesselpars + pdegree[1];
    
    double *e2ij = &tmpmem[0]; // Nij*nrbf
    double *e2ijt = &tmpmem[Nij*nrbf]; // Nij*ns

    // orthogonal radial basis functions
    this->podradialbasis(e2ijt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nij);
    this->podMatMul(e2ij, e2ijt, Phi, Nij, ns, nrbf);
  
    
    // one-body descriptors
    this->pod1body(eatom1, atomtype, nelements, natom);

    this->podtally2b(eatom2, e2ij, idxi, ti, tj, elemindex, nelements, nrbf2, natom, Nij);   

    // three-body descriptors
    this->pod3body(eatom3, rij, e2ij, &tmpmem[Nij*nrbf], elemindex, pairnumsum, 
             idxi, ti, tj, nrbf3, nabf, nelements, natom, Nij);     
}

void CPOD::snapComputeUij(double *Sr, double *Si, double *rootpqarray, double *rij, 
        double *wjelem, double *radelem, double rmin0, double rfac0, double rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)                
{    
  for(int ij=0; ij<ijnum; ij++) {        
    double x = rij[ij*3+0];
    double y = rij[ij*3+1];
    double z = rij[ij*3+2];
    double rsq = x * x + y * y + z * z;
    double r = sqrt(rsq);

    double rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    double rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    double theta0 = (r - rmin0) * rscale0;
    double z0 = r / tan(theta0);                
            
    double sfac = 0.0; 
    if (switch_flag == 0) {
        sfac = 1.0;
    }
    else if (switch_flag == 1) {
        if (r <= rmin0) {
            sfac = 1.0;
        }
        else if(r > rcutij) {
            sfac = 0.0;
        }
        else {
            double rcutfac0 = M_PI / (rcutij - rmin0);
            sfac =  0.5 * (cos((r - rmin0) * rcutfac0) + 1.0);   
        }
    } 
    sfac *= wjelem[tj[ij]];
    
    double r0inv;
    double a_r, a_i, b_r, b_i;
    double rootpq;
    int jdim = twojmax + 1;
  
    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;
    
    Sr[ij+0*ijnum] = 1.0;
    Si[ij+0*ijnum] = 0.0;
    for (int j = 1; j <= twojmax; j++) {
        int jju = idxu_block[j];
        int jjup = idxu_block[j-1];
        
        // fill in left side of matrix layer from previous layer
        for (int mb = 0; 2*mb <= j; mb++) {
            Sr[ij+jju*ijnum] = 0.0;
            Si[ij+jju*ijnum] = 0.0;
            for (int ma = 0; ma < j; ma++) {
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                int njju = ij+jju*ijnum;
                int njju1 = ij+(jju+1)*ijnum;
                int njjup = ij+jjup*ijnum;
                double u_r = Sr[njjup];
                double u_i = Si[njjup];

                Sr[njju] += rootpq * (a_r * u_r + a_i * u_i);
                Si[njju] += rootpq * (a_r * u_i - a_i * u_r);

                rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
                Sr[njju1] = -rootpq * (b_r * u_r + b_i * u_i);
                Si[njju1] = -rootpq * (b_r * u_i - b_i * u_r);
                jju++;
                jjup++;
            }
            jju++;
        }
                   
        jju = idxu_block[j];
        jjup = jju+(j+1)*(j+1)-1;
        int mbpar = 1;
        for (int mb = 0; 2*mb <= j; mb++) {
            int mapar = mbpar;
            for (int ma = 0; ma <= j; ma++) {
                int njju = ij+jju*ijnum;
                int njjup = ij+jjup*ijnum;
                if (mapar == 1) {
                    Sr[njjup] = Sr[njju];
                    Si[njjup] = -Si[njju];
                } else {
                    Sr[njjup] = -Sr[njju];
                    Si[njjup] =  Si[njju];
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }        
    }        
    
    for (int k=0; k<idxu_max; k++) {
        int ijk = ij + ijnum*k;
        Sr[ijk] = sfac*Sr[ijk];
        Si[ijk] = sfac*Si[ijk];
    }            
  }
};

void CPOD::snapdesc_ij(double *blist, double *rij, double *tmpmem, int *atomtype, int *idxi, 
        int *ti, int *tj, int natom, int Nij)            
{    
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int switchflag = sna.switchflag;    
    //int bzeroflag = sna.bzeroflag;
    int wselfallflag = sna.wselfallflag;
    int nelem = (chemflag) ? nelements : 1;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    //int *idxb_block = sna.idxb_block;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    double wself = sna.wself;
    double rmin0 = sna.rmin0;
    double rfac0 = sna.rfac0;
    double rcutfac = sna.rcutfac;
    double *rootpqarray = sna.rootpqarray;
    double *cglist = sna.cglist;
    double *radelem = sna.radelem;
    double *wjelem = sna.wjelem; 
            
    int ne = 0;
    double *Ur = &tmpmem[ne]; 
    double *Zr = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *Ui = &tmpmem[ne]; 
    double *Zi = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *Utotr = &tmpmem[ne];
    ne += idxu_max*nelements*natom;
    double *Utoti = &tmpmem[ne];        
                    
    this->snapComputeUij(Ur, Ui, rootpqarray, rij, wjelem, radelem, rmin0, 
         rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, Nij, switchflag);
        
    this->snapZeroUarraytot2(Utotr, Utoti, wself, idxu_block, atomtype, map, idxi, wselfallflag, 
            chemflag, idxu_max, nelem, twojmax, natom);

    this->snapAddUarraytot(Utotr, Utoti, Ur, Ui, map, idxi, tj, idxu_max, natom, Nij, chemflag);    
        
    this->snapComputeZi2(Zr, Zi, Utotr, Utoti, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, natom);

    this->snapComputeBi1(blist, Zr, Zi, Utotr, Utoti, idxb, idxu_block, idxz_block, twojmax, idxb_max, 
            idxu_max, idxz_max, nelem, natom);                        
}

void CPOD::linear_descriptors_ij(double *gd, double *eatom, double *rij, double *tmpmem, int *pairnumsum,
        int *atomtype, int *idxi, int *ti, int *tj, int natom, int Nij)
{
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
    
    double *eatom1 = &eatom[0];
    double *eatom2 = &eatom[0+natom*nd1];
    double *eatom3 = &eatom[0+natom*(nd1+nd2)];
    double *eatom4 = &eatom[0+natom*(nd1+nd2+nd3)];    
        
    podArraySetValue(eatom1, 0.0, natom*nd1234);    
    
    // peratom descriptors for one-body, two-body, and three-body linear potentials
    this->poddesc_ij(eatom1, eatom2, eatom3, rij, Phi2, besselparams, 
            tmpmem, rin, rcut, pairnumsum, atomtype, idxi, ti, tj, elemindex, pdegree2, 
            nbesselpars, nrbf2, nrbf3, nabf3, nelements, Nij, natom);                    
    
    // peratom snap descriptors
    if (pod.snaptwojmax>0) 
        this->snapdesc_ij(eatom4, rij, tmpmem, atomtype, idxi, ti, tj, natom, Nij);                
    
    // global descriptors for one-body, two-body, three-body, and four-bodt linear potentials    
    podArraySetValue(tmpmem, 1.0, natom);
    
    char cht = 'T';
    double one = 1.0;    
    int inc1 = 1;
    DGEMV(&cht, &natom, &nd1234, &one, eatom1, &natom, tmpmem, &inc1, &one, gd, &inc1);                        
}

double CPOD::calculate_energy(double *effectivecoeff, double *gd, double *coeff, int natom)
{        
    //int dim = 3;    
    //int nforce = dim*natom;
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
                    
    // calculate energy for linear potentials
    double energy = 0.0;
    for (int i=0; i< nd1234; i++) {
        effectivecoeff[i] = 0.0;    
        energy += coeff[i]*gd[i];    
    }
    
    // effective POD coefficients for calculating force 
    //double *c1 = &effectivecoeff[0];
    double *c2 = &effectivecoeff[nd1];  
    double *c3 = &effectivecoeff[nd1+nd2];
    double *c4 = &effectivecoeff[nd1+nd2+nd3];
    
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
    for (int i=0; i< nd1234; i++) effectivecoeff[i] += coeff[i];    
                
    return energy;
}

void CPOD::pod2body_force(double *force, double *fij, double *coeff2, int *ai, int *aj, 
        int *ti, int *tj, int *elemindex, int nelements, int nbf, int natom, int N)
{
    int nelements2 = nelements*(nelements+1)/2;
    for (int n=0; n<N; n++) {
        int i1 = ai[n];
        int j1 = aj[n];
        int typei = ti[n]-1;
        int typej = tj[n]-1;
        for (int m=0; m<nbf; m++) {               
            int im =  3*i1;
            int jm =  3*j1;
            int nm = n + N*m;
            int km = (elemindex[typei + typej*nelements] - 1) + nelements2*m;
            double ce = coeff2[km];
            force[0 + im] += fij[0 + 3*nm]*ce;
            force[1 + im] += fij[1 + 3*nm]*ce;
            force[2 + im] += fij[2 + 3*nm]*ce;
            force[0 + jm] -= fij[0 + 3*nm]*ce;
            force[1 + jm] -= fij[1 + 3*nm]*ce;
            force[2 + jm] -= fij[2 + 3*nm]*ce;          
        }
    }
}

void CPOD::pod3body_force(double *force, double *yij, double *e2ij, double *f2ij, double *coeff3, double *tmpmem, 
             int *elemindex, int *pairnumsum, int *ai, int *aj, int *ti, int *tj, int nrbf, int nabf, 
             int nelements, int natom, int Nij)
{   
    int dim = 3, nabf1 = nabf + 1;
    int nelements2 = nelements*(nelements+1)/2;
    int n, c, nijk3;
    
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, sinthe, theta, dtheta; 
    double tm, tm1, tm2, dct1, dct2, dct3, dct4, dct5, dct6;
    double uj, uk, rbf, drbf1, drbf2, drbf3, drbf4, drbf5, drbf6;
    double fj1, fj2, fj3, fk1, fk2, fk3;
            
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
                        tm = abf[p];
                        fj1 = drbf1*tm + rbf*dabf1[p];
                        fj2 = drbf2*tm + rbf*dabf2[p];
                        fj3 = drbf3*tm + rbf*dabf3[p];
                        fk1 = drbf4*tm + rbf*dabf4[p];
                        fk2 = drbf5*tm + rbf*dabf5[p];
                        fk3 = drbf6*tm + rbf*dabf6[p];

                        n = p + (nabf1)*m;
                        c = (elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*n;
                        tm = coeff3[c];
                        
                        nijk3 = 3*i;                        
                        force[0 + nijk3] += (fj1 + fk1)*tm;
                        force[1 + nijk3] += (fj2 + fk2)*tm;
                        force[2 + nijk3] += (fj3 + fk3)*tm;

                        nijk3 = 3*j;    
                        force[0 + nijk3] -= fj1*tm;
                        force[1 + nijk3] -= fj2*tm;
                        force[2 + nijk3] -= fj3*tm;

                        nijk3 = 3*k;    
                        force[0 + nijk3] -= fk1*tm;   
                        force[1 + nijk3] -= fk2*tm;   
                        force[2 + nijk3] -= fk3*tm;                           
                    }                    
                }
            }
        }
    }
}

void CPOD::snapTallyForce(double *force, double *dbdr, double *coeff4,
        int *ai, int *aj, int *ti, int ijnum, int ncoeff, int ntype)
{           
    int N2 = ijnum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ij = idx%ijnum;
        int icoeff = (idx-ij)/ijnum;        
        int i = ai[ij]; // index of atom i
        int j = aj[ij]; // index of atom i
        int itype = ti[ij]; // element type of atom i       
        int n = ncoeff*(itype-1);        
        int nij = ijnum*3*icoeff;
        
        double bix = dbdr[ij + ijnum*0 + nij];
        double biy = dbdr[ij + ijnum*1 + nij];
        double biz = dbdr[ij + ijnum*2 + nij];      
        double ce = coeff4[icoeff + n];
        
        force[0 + 3*i] += bix*ce; 
        force[1 + 3*i] += biy*ce;
        force[2 + 3*i] += biz*ce;
        force[0 + 3*j] -= bix*ce;
        force[1 + 3*j] -= biy*ce;
        force[2 + 3*j] -= biz*ce;        
    }
}

void CPOD::pod4body_force(double *force, double *rij, double *coeff4, double *tmpmem, int *atomtype, 
        int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij)            
{    
    int dim = 3;    
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int switchflag = sna.switchflag;    
    int wselfallflag = sna.wselfallflag;
    int nelem = (chemflag) ? nelements : 1;
    
    int *map = sna.map;
    int *idxz = sna.idxz;
    int *idxz_block = sna.idxz_block;
    int *idxb = sna.idxb;
    int *idxu_block = sna.idxu_block;
    int *idxcg_block = sna.idxcg_block;   
    
    double wself = sna.wself;
    double rmin0 = sna.rmin0;
    double rfac0 = sna.rfac0;
    double rcutfac = sna.rcutfac;
    double *rootpqarray = sna.rootpqarray;
    double *cglist = sna.cglist;
    double *radelem = sna.radelem;
    double *wjelem = sna.wjelem; 
            
    int ne = 0;
    double *Ur = &tmpmem[ne]; 
    double *Zr = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *Ui = &tmpmem[ne]; 
    double *Zi = &tmpmem[ne]; 
    ne += PODMAX(idxu_max*Nij, idxz_max*ndoubles*natom); 
    double *dUr = &tmpmem[ne];
    ne += idxu_max*dim*Nij;
    double *dUi = &tmpmem[ne];
    ne += idxu_max*dim*Nij;    
    double *dblist = &tmpmem[ne]; // idxb_max*ntriples*dim*Nij          
    double *Utotr = &tmpmem[ne];
    ne += idxu_max*nelements*natom;
    double *Utoti = &tmpmem[ne];        
                    
    this->snapComputeUlist(Ur, Ui, dUr, dUi, rootpqarray, rij, wjelem, radelem, rmin0, 
         rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, Nij, switchflag);
    
    this->snapZeroUarraytot2(Utotr, Utoti, wself, idxu_block, atomtype, map, idxi, wselfallflag, 
            chemflag, idxu_max, nelem, twojmax, natom);

    this->snapAddUarraytot(Utotr, Utoti, Ur, Ui, map, idxi, tj, idxu_max, natom, Nij, chemflag);    
        
    this->snapComputeZi2(Zr, Zi, Utotr, Utoti, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, natom);
            
    this->snapComputeDbidrj(dblist, Zr, Zi, dUr, dUi, idxb, idxu_block, idxz_block, map, idxi, tj, 
            twojmax, idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, natom, Nij);
    
    this->snapTallyForce(force, dblist, coeff4, ai, aj, ti, Nij, ncoeff, ntypes);           
}

void CPOD::calculate_force(double *force, double *effectivecoeff, double *rij, double *tmpmem, int *pairnumsum,
        int *atomtype, int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij)
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
    int *pdegree = pod.twobody;
    int *elemindex = pod.elemindex;
    double rin = pod.rin;
    double rcut = pod.rcut;
    double *Phi = pod.Phi2;
    double *besselparams = pod.besselparams;        
    
    // effective POD coefficients for calculating force 
    double *coeff2 = &effectivecoeff[nd1];  
    double *coeff3 = &effectivecoeff[nd1+nd2];
    double *coeff4 = &effectivecoeff[nd1+nd2+nd3];    
    
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
    
    podArraySetValue(force, 0.0, dim*natom);    
    
    this->pod2body_force(force, f2ij, coeff2, ai, aj, ti, tj, elemindex, nelements, nrbf2, natom, Nij);
    
    this->pod3body_force(force, rij, e2ij, f2ij, coeff3, &tmpmem[4*Nij*nrbf], elemindex, pairnumsum, ai, aj, 
            ti, tj, nrbf3, nabf3, nelements, natom, Nij);
        
    if (pod.snaptwojmax>0) 
        this->pod4body_force(force, rij, coeff4, tmpmem, atomtype, idxi, ai, aj, ti, tj, natom, Nij);       
            
//     double eatom[natom*nd3];
//     double fatom[3*natom*nd3];    
//     podArraySetValue(eatom, 0.0, natom*nd3);    
//     podArraySetValue(fatom, 0.0, dim*natom*nd3);    
//     this->pod3body(eatom, fatom, rij, e2ij, f2ij,  &tmpmem[4*Nij*nrbf], elemindex, pairnumsum, ai, aj, 
//             ti, tj, nrbf3, nabf3, nelements, natom, Nij);
// 
//     
//     // calculate force = gdd * c1
//     char chn = 'N';
//     double one = 1.0, zero = 0.0;    
//     int inc1 = 1;
//     int nforce = 3*natom;
//     cout<<nd3<<endl;
//     DGEMV(&chn, &nforce, &nd3, &one, fatom, &nforce, coeff3, &inc1, &one, force, &inc1);        
//     
//     print_matrix( "fatom3", nforce, nd3, fatom, nforce);
//     
// //     print_matrix( "coeff2", 1, nd2, coeff2, 1);
//     print_matrix( "coeff3", 1, nd3, coeff3, 1);
    
//     if (pod.snaptwojmax>0) {
//         this->pod4body_force(force, coeff4, rij, tmpmem, atomtype, idxi, ai, aj, ti, tj, natom, Nij);       
//         int nd4 = pod.nd4;
//         double blist[natom*nd4]; 
//         double bd[dim*natom*nd4]; 
//         this->snapdesc_ij(blist, bd, rij, tmpmem, atomtype, ai, aj, ti, tj, natom, Nij);            
//      
//         char chn = 'N';
//         double one = 1.0, zero = 0.0;    
//         int inc1 = 1;
//         int nforce = 3*natom;
//         DGEMV(&chn, &nforce, &nd4, &one, bd, &nforce, coeff4, &inc1, &one, force, &inc1);                
//    }
}

double CPOD::energyforce_calculation(double *force, double *podcoeff, double *effectivecoeff, double *gd, double *rij, 
        double *tmpmem, int *pairnumsum, int *atomtype, int *idxi, int *ai, int *aj, int *ti, int *tj, int natom, int Nij)
{       
    int nd1234 = pod.nd1+pod.nd2+pod.nd3+pod.nd4;
    double *eatom = &tmpmem[0];
                
    podArraySetValue(gd, 0.0, nd1234);    
    this->linear_descriptors_ij(gd, eatom, rij, &tmpmem[natom*nd1234], pairnumsum, atomtype, idxi, ti, tj, natom, Nij);
    
    // Need to do MPI_Allreduce on gd for paralell
    
    double energy = this->calculate_energy(effectivecoeff, gd, podcoeff, natom);    
        
    this->calculate_force(force, effectivecoeff, rij, tmpmem, pairnumsum, atomtype, idxi, ai, aj, ti, tj, natom, Nij);
    
    return energy;
}



// void snapComputeUi(double *Utotr, double *Utoti, double *rootpqarray, double *rij, double *wjelem, double *radelem, 
//         double rmin0, double rfac0, double rcutfac, int *idxu_block, int *map, int *aii, int *ti, int *tj, 
//         int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag)                
// {    
//   for(int ij=0; ij<ijnum; ij++) {        
//     double x = rij[ij*3+0];
//     double y = rij[ij*3+1];
//     double z = rij[ij*3+2];    
//     double r = sqrt(x * x + y * y + z * z);
// 
//     double rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
//     //double rscale0 = rfac0 * M_PI / (rcutij - rmin0);
//     //double theta0 = (r - rmin0) * rscale0;
//     double z0 = r / tan((r - rmin0) * rfac0 * M_PI / (rcutij - rmin0));                
//             
//     double sfac = 0.0;
//     if (switchflag == 0) 
//         sfac = 1.0;    
//     else if (switchflag == 1) {
//         if (r <= rmin0) {
//             sfac = 1.0;
//         }
//         else if(r > rcutij) {
//             sfac = 1.0;
//         }
//         else {
//             double rcutfac = M_PI / (rcutij - rmin0);
//             sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
//         }
//     } 
//     sfac *= wjelem[tj[ij]];
//     
//     double a_r, a_i, b_r, b_i;
//     double rootpq;
// 
//     //double r0inv;    
//     rcutij = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = rcutij * z0;
//     a_i = -rcutij * z;
//     b_r = rcutij * y;
//     b_i = -rcutij * x;
// 
//     // 2Jmax = 10
//     double Pr[11], Pi[11], Qr[9], Qi[9];
//     Pr[0] = 1.0;
//     Pi[0] = 0.0;    
//     
//     int jdim = twojmax + 1;
//     int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
//     int i = aii[ij] + njelem;                        
//     Utotr[i] += sfac; // atomic add   
//     
//     int mb = 0;    
//     for (int j = 1; j <= twojmax; j++) {        
//         // fill in left side of matrix layer from previous layer
//         int ma = 0;
//         // x y z z0 
//         // double p_r, p_i; // -> x, y
//         // double u_r = Pr[ma]; // -> z
//         // double u_i = Pi[ma]; // -> z0
//         z = Pr[ma];
//         z0 = Pi[ma];
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//         Pr[ma] = rootpq * (a_r * z + a_i * z0);
//         Pi[ma] = rootpq * (a_r * z0 - a_i * z);            
//         for (ma = 1; ma < j; ma++) {
//             x = Pr[ma];
//             y = Pi[ma];
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             rcutij = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma] = rootpq * (a_r * x + a_i * y) -rcutij * (b_r * z + b_i * z0);
//             Pi[ma] = rootpq * (a_r * y - a_i * x) -rcutij * (b_r * z0 - b_i * z);
//             z = x;
//             z0 = y;
//         }
//         ma = j;
//         rcutij = rootpqarray[ma*jdim + (j - mb)];
//         Pr[ma] = -rcutij * (b_r * z + b_i * z0);
//         Pi[ma] = -rcutij * (b_r * z0 - b_i * z);                        
//                                 
//         if (j==(2*mb+1)) { // store Qr, Qi, for the next mb level
//             int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//             for (int ma = 0; ma <= j; ma++) {     
//                 if (mapar == 1) {                    
//                     Qr[j-ma] = Pr[ma];
//                     Qi[j-ma] = -Pi[ma];
//                 } else {
//                     Qr[j-ma] = -Pr[ma];
//                     Qi[j-ma] =  Pi[ma];
//                 }
//                 mapar = -mapar;
//             }                                                
//         }
//         
//         int k =  1 + (j+1)*mb;
//         for (int ma = 2; ma <= j; ma++)
//             k += ma*ma;                    
//         for (int ma = 0; ma <= j; ma++) {
//             int in = i + inum*k;                
//             Utotr[in] += sfac*Pr[ma]; // atomic add   
//             Utoti[in] += sfac*Pi[ma]; // atomic add                       
//             k += 1;
//         }                   
//     }
//     
//     for (mb = 1; 2*mb <= twojmax; mb++) {     
//         for (int ma = 0; ma < 2*mb; ma++) {                      
//             Pr[ma] = Qr[ma];
//             Pi[ma] = Qi[ma];
//         }                
//         for (int j = 2*mb; j <= twojmax; j++) { 
//             int ma = 0;
//             //double p_r, p_i;
//             //double u_r = Pr[ma];
//             //double u_i = Pi[ma];
//             z = Pr[ma];
//             z0 = Pi[ma];
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             Pr[ma] = rootpq * (a_r * z + a_i * z0);
//             Pi[ma] = rootpq * (a_r * z0 - a_i * z);            
//             for (ma = 1; ma < j; ma++) {
//                 x = Pr[ma];
//                 y = Pi[ma];
//                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//                 rcutij = rootpqarray[ma*jdim + (j - mb)];
//                 Pr[ma] = rootpq * (a_r * x + a_i * y) -rcutij * (b_r * z + b_i * z0);
//                 Pi[ma] = rootpq * (a_r * y - a_i * x) -rcutij * (b_r * z0 - b_i * z);
//                 z = x;
//                 z0 = y;
//             }
//             ma = j;
//             rcutij = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma] = -rcutij * (b_r * z + b_i * z0);
//             Pi[ma] = -rcutij * (b_r * z0 - b_i * z);       
//             
//             if (j==(2*mb)) {
//                 int mapar = 1;
//                 for (int ma = 0; ma <= j/2; ma++) {
//                     if (mapar == 1) {                    
//                         Pr[j/2+ma] = Pr[j/2-ma];
//                         Pi[j/2+ma] = -Pi[j/2-ma];
//                     } else {
//                         Pr[j/2+ma] = -Pr[j/2-ma];
//                         Pi[j/2+ma] = Pi[j/2-ma];
//                     }
//                     mapar = -mapar;        
//                 }                                                        
//             }
//             
//             // store Qr, Qi, for the next mb level
//             if (j==(2*mb+1)) {
//                 int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//                 for (int ma = 0; ma <= j; ma++) {     
//                     if (mapar == 1) {                    
//                         Qr[j-ma] = Pr[ma];
//                         Qi[j-ma] = -Pi[ma];
//                     } else {
//                         Qr[j-ma] = -Pr[ma];
//                         Qi[j-ma] =  Pi[ma];
//                     }
//                     mapar = -mapar;
//                 }                                                
//             }
//             
//             int k =  1 + (j+1)*mb;
//             for (int ma = 2; ma <= j; ma++)
//                 k += ma*ma;            
//             for (int ma = 0; ma <= j; ma++) {
//                 int in = i + inum*k;                
//                 Utotr[in] += sfac*Pr[ma]; // atomic add   
//                 Utoti[in] += sfac*Pi[ma]; // atomic add                       
//                 k += 1; 
//             }                                                           
//         }
//     }        
//   }
// };
// 
// void snapAddWself2Ui(double *Utotr, double *Utoti, double wself, int *idxu_block, int *type, int *map, 
//         int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, int twojmax, int inum)
// {
//     int N1 = inum;
//     int N2 = N1*(twojmax+1);
//     int N3 = N2*nelements;                                
//     
//     for (int idx=0; idx < N3; idx++) {
//         int l = idx%N2;  // inum*(twojmax+1)
//         int ii = l%N1;    // inum
//         int j = (l-ii)/N1; // (twojmax+1)
//         int jelem = (idx-l)/N2; // nelements   
//         int ielem = (chemflag) ? map[type[ii]]: 0;                
//         int nmax = ii + inum*idxu_max*jelem;
//         
//         int jju = idxu_block[j];                
//         for (int mb = 0; mb <= j; mb++) {
//             for (int ma = 0; ma <= j; ma++) {                
//                 if (jelem == ielem || wselfall_flag)
//                     if (ma==mb)                        
//                         Utotr[inum*jju + nmax] += wself;                                     
//                 jju++;                                
//             }
//         }
//         
//         // copy left side to right side with inversion symmetry VMK 4.4(2)
//         // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
//         
//         jju = idxu_block[j];        
//         int jjup = jju+(j+1)*(j+1)-1;
//         int mbpar = 1;
//         for (int mb = 0; 2*mb < j; mb++) {
//             int mapar = mbpar;
//             for (int ma = 0; ma <= j; ma++) {
//                 int njju =  inum*jju + nmax;
//                 int njjup = inum*jjup + nmax;
//                 if (mapar == 1) {
//                     Utotr[njjup] = Utotr[njju];
//                     Utoti[njjup] = -Utoti[njju];
//                 } else {
//                     Utotr[njjup] = -Utotr[njju];
//                     Utoti[njjup] =  Utoti[njju];
//                 }
//                 mapar = -mapar;
//                 jju++;
//                 jjup--;
//             }
//             mbpar = -mbpar;
//         }        
//     }                    
// };
// 
// void snapComputeYi(double *ylist_r, double *ylist_i, double *Utotr, double *Utoti, double *cglist, double* beta, 
//         int *map, int *type, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
//         int idxb_max, int idxu_max, int idxz_max, int nelements, int ncoeffall, int bnorm_flag, int inum)
// {    
//     
//     int N1 = idxu_max*nelements*inum;
//     //printf("%i %i %i %i %i %i\n", inum, idxu_max, nelements, N1, idxz_max, idxb_max);
//     podArraySetValue(ylist_r, (double) 0.0, N1);
//     podArraySetValue(ylist_i, (double) 0.0, N1);
//     //printf("%i %i %i %i %i %i\n", inum, idxu_max, nelements, N1, idxz_max, idxb_max);
//     
//     int jdim = twojmax + 1;         
//     int N2 = idxz_max*inum;                          
//     for (int idx=0; idx < N2; idx++) {
//       int ii = idx%inum;              
//       int jjz = (idx-ii)/inum;         
//       int jjz10 = jjz*10;
//       const int j1 = idxz[jjz10+0];
//       const int j2 = idxz[jjz10+1];
//       const int j = idxz[jjz10+2];
//       const int ma1min = idxz[jjz10+3];
//       const int ma2max = idxz[jjz10+4];
//       const int na = idxz[jjz10+5];
//       const int mb1min = idxz[jjz10+6];
//       const int mb2max = idxz[jjz10+7];
//       const int nb = idxz[jjz10+8];          
//       const int jju = idxz[jjz10+9];
//       const double *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
//       
//       int itype = type[ii]; // element type of atom i
//       int ielem = map[itype];  // index of that element type                        
// 
//       for(int elem1 = 0; elem1 < nelements; elem1++)
//         for (int elem2 = 0; elem2 < nelements; elem2++) {        
//           
//             double ztmp_r = 0.0;
//             double ztmp_i = 0.0;
//             
//             int ii1 = ii + inum*idxu_max*elem1;
//             int ii2 = ii + inum*idxu_max*elem2;
//             int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
//             int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
//             int icgb = mb1min * (j2 + 1) + mb2max;
//             for (int ib = 0; ib < nb; ib++) {
//                 int ma1 = ma1min;
//                 int ma2 = ma2max;
//                 int icga = ma1min * (j2 + 1) + ma2max;
//                 for (int ia = 0; ia < na; ia++) {
//                     ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
//                     ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
//                   ma1++;
//                   ma2--;
//                   icga += j2;
//                 } // end loop over ia
// 
//                 jju1 += j1 + 1;
//                 jju2 -= j2 + 1;
//                 icgb += j2;
//             } // end loop over ib
// 
//             if (bnorm_flag) {
//               ztmp_i /= j+1;
//               ztmp_r /= j+1;
//             }            
//                                
//             for(int elem3 = 0; elem3 < nelements; elem3++) {
//               int itriple;  
//               double betaj;
//               if (j >= j1) {
//                 const int jjb = idxb_block[j + j2*jdim + j1*jdim*jdim]; 
//                 itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb + 1 + ielem*ncoeffall;
//                 if (j1 == j) {
//                   if (j2 == j) betaj = 3*beta[itriple];
//                   else betaj = 2*beta[itriple];
//                 } else betaj = beta[itriple];          
//               } else if (j >= j2) {
//                 const int jjb = idxb_block[j1 + j2*jdim + j*jdim*jdim];
//                 itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb + 1 + ielem*ncoeffall;
//                 if (j2 == j) betaj = 2*beta[itriple];
//                 else betaj = beta[itriple];
//               } else {
//                 const int jjb = idxb_block[j1 + j*jdim + j2*jdim*jdim];
//                 itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb + 1 + ielem*ncoeffall;
//                 betaj = beta[itriple];
//               }
//               
//               if (!bnorm_flag && j1 > j)
//                 betaj *= (j1 + 1) / (j + 1.0);
//                          
//               //printf("%i %i %i %i %i %i\n", inum, ii, jju, elem1, elem2, elem3);
//               ylist_r[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_r;
//               ylist_i[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_i;        
//            }
//         }         
//     }  
// }
// 
// void snapComputeFi(double *fatom, double *ylist_r, double *ylist_i, double *rootpqarray, double *rij, 
//         double *wjelem, double *radelem,  double rmin0, double rfac0, double rcutfac, int *idxu_block, int *map, int *aii, int *ai, int *aj, 
//         int *ti, int *tj, int twojmax, int idxu_max, int inum, int anum, int ijnum, int switchflag, int chemflag) 
// {                 
//   for(int ij=0; ij<ijnum; ij++) {        
//     double x = rij[ij*3+0];
//     double y = rij[ij*3+1];
//     double z = rij[ij*3+2];    
//     double rsq = x * x + y * y + z * z;
//     double r = sqrt(rsq);
//     double rinv = 1.0 / r;
//     double ux = x * rinv;
//     double uy = y * rinv;
//     double uz = z * rinv;
// 
//     double rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
//     double z0 = r / tan((r - rmin0) * rfac0 * M_PI / (rcutij - rmin0));                
//     double dz0dr = z0 / r - (r*rfac0 * M_PI / (rcutij - rmin0)) * (rsq + z0 * z0) / rsq;
// 
//     double sfac = 0.0, dsfac = 0.0;        
//     if (switchflag == 0) {
//         sfac = 1.0;
//         dsfac = 0.0;
//     }
//     else if (switchflag == 1) {
//         if (r <= rmin0) {
//             sfac = 1.0;
//             dsfac = 0.0;
//         }
//         else if(r > rcutij) {
//             sfac = 1.0;
//             dsfac = 0.0;
//         }
//         else {
//             double rcutfac = M_PI / (rcutij - rmin0);
//             sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
//             dsfac = -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
//         }
//     } 
//     sfac *= wjelem[tj[ij]];
//     dsfac *= wjelem[tj[ij]];
//     
//     double a_r, a_i, b_r, b_i, rootpq;        
//     rcutij = 1.0 / sqrt(r * r + z0 * z0);
//     a_r = rcutij * z0;
//     a_i = -rcutij * z;
//     b_r = rcutij * y;
//     b_i = -rcutij * x;
// 
//     double u_r, u_i, ux_r, ux_i, uy_r, uy_i, uz_r, uz_i;
//     double w_r, w_i, wx_r, wx_i, wy_r, wy_i, wz_r, wz_i;
//     u_r = -pow(rcutij, 3.0) * (r + z0 * dz0dr);
//     wx_r = u_r * ux;
//     wy_r = u_r * uy;
//     wz_r = u_r * uz;
//     ux_r = dz0dr * ux;
//     uy_r = dz0dr * uy;
//     uz_r = dz0dr * uz;
// 
//     double dardx, daidx, dardy, daidy, dardz, daidz;
//     dardx = ux_r * rcutij + z0 * wx_r;
//     daidx = -z * wx_r;
//     dardy = uy_r * rcutij + z0 * wy_r;
//     daidy = -z * wy_r;
//     dardz = uz_r * rcutij + z0 * wz_r;
//     daidz = -z * wz_r;    
//     daidz += -rcutij;
// 
//     double dbrdx, dbidx, dbrdy, dbidy, dbrdz, dbidz;
//     dbrdx = y * wx_r;
//     dbidx = -x * wx_r;    
//     dbrdy = y * wy_r;
//     dbidy = -x * wy_r;    
//     dbrdz = y * wz_r;
//     dbidz = -x * wz_r;        
//     dbidx += -rcutij;
//     dbrdy += rcutij;
//     
//     // 2Jmax = 10    
//     double Pr[11], Pi[11], Qr[9], Qi[9];
//     double Prx[11], Pix[11], Qrx[9], Qix[9];
//     double Pry[11], Piy[11], Qry[9], Qiy[9];        
//     double Prz[11], Piz[11], Qrz[9], Qiz[9];
//     Pr[0] = 1.0; Pi[0] = 0.0;    
//     Prx[0] = 0.0; Pix[0] = 0.0;    
//     Pry[0] = 0.0; Piy[0] = 0.0;    
//     Prz[0] = 0.0; Piz[0] = 0.0;        
//     
//     int jdim = twojmax + 1;
//     int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
//     int i = aii[ij] + njelem;                        
//                     
//     double dedx, dedy, dedz;       
//     u_r = 0.5*ylist_r[i];  
//     //e    =  sfac*u_r;
//     dedx = (dsfac * ux) * u_r;
//     dedy = (dsfac * uy) * u_r;
//     dedz = (dsfac * uz) * u_r; 
//             
//     int j, k, ma, mb, mapar;    
//     mb = 0;
//     for (j = 1; j <= twojmax; j++) {        
//         // fill in left side of matrix layer from previous layer
//         ma = 0;
//         u_r = Pr[ma];
//         u_i = Pi[ma];
//         ux_r = Prx[ma];
//         ux_i = Pix[ma];            
//         uy_r = Pry[ma];
//         uy_i = Piy[ma];            
//         uz_r = Prz[ma];
//         uz_i = Piz[ma];                    
//         rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//         Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
//         Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);        
//         Prx[ma] = rootpq * (dardx * u_r + daidx * u_i + a_r * ux_r + a_i * ux_i);
//         Pix[ma] = rootpq * (dardx * u_i - daidx * u_r + a_r * ux_i - a_i * ux_r);
//         Pry[ma] = rootpq * (dardy * u_r + daidy * u_i + a_r * uy_r + a_i * uy_i);
//         Piy[ma] = rootpq * (dardy * u_i - daidy * u_r + a_r * uy_i - a_i * uy_r);
//         Prz[ma] = rootpq * (dardz * u_r + daidz * u_i + a_r * uz_r + a_i * uz_i);
//         Piz[ma] = rootpq * (dardz * u_i - daidz * u_r + a_r * uz_i - a_i * uz_r);                    
//         for (ma = 1; ma < j; ma++) {
//             w_r = Pr[ma];
//             w_i = Pi[ma];
//             wx_r = Prx[ma];
//             wx_i = Pix[ma];            
//             wy_r = Pry[ma];
//             wy_i = Piy[ma];            
//             wz_r = Prz[ma];
//             wz_i = Piz[ma];                                        
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             rcutij = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
//             Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
//             Prx[ma] = rootpq * (dardx * w_r + daidx * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
//             Pix[ma] = rootpq * (dardx * w_i - daidx * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
//             Pry[ma] = rootpq * (dardy * w_r + daidy * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
//             Piy[ma] = rootpq * (dardy * w_i - daidy * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
//             Prz[ma] = rootpq * (dardz * w_r + daidz * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
//             Piz[ma] = rootpq * (dardz * w_i - daidz * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);            
//             u_r = w_r;
//             u_i = w_i;
//             ux_r = wx_r;
//             ux_i = wx_i;
//             uy_r = wy_r;
//             uy_i = wy_i;
//             uz_r = wz_r;
//             uz_i = wz_i;            
//         }
//         ma = j;
//         rcutij = rootpqarray[ma*jdim + (j - mb)];
//         Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
//         Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);                        
//         Prx[ma] =-rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
//         Pix[ma] =-rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
//         Pry[ma] =-rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
//         Piy[ma] =-rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
//         Prz[ma] =-rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
//         Piz[ma] =-rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);
//                                 
//         if (j==(2*mb+1)) { // store Qr, Qi, for the next mb level
//             mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//             for (ma = 0; ma <= j; ma++) {     
//                 if (mapar == 1) {                    
//                     Qr[j-ma] = Pr[ma];
//                     Qi[j-ma] = -Pi[ma];
//                     Qrx[j-ma] =  Prx[ma];
//                     Qix[j-ma] = -Pix[ma];
//                     Qry[j-ma] =  Pry[ma];
//                     Qiy[j-ma] = -Piy[ma];
//                     Qrz[j-ma] =  Prz[ma];
//                     Qiz[j-ma] = -Piz[ma];                    
//                 } else {
//                     Qr[j-ma] = -Pr[ma];
//                     Qi[j-ma] =  Pi[ma];
//                     Qrx[j-ma] = -Prx[ma];
//                     Qix[j-ma] =  Pix[ma];
//                     Qry[j-ma] = -Pry[ma];
//                     Qiy[j-ma] =  Piy[ma];
//                     Qrz[j-ma] = -Prz[ma];
//                     Qiz[j-ma] =  Piz[ma];                    
//                 }
//                 mapar = -mapar;
//             }                              
//         }
//         
//         k =  1 + (j+1)*mb;
//         for (ma = 2; ma <= j; ma++)
//             k += ma*ma;                    
//         for (ma = 0; ma <= j; ma++) {                            
//             rsq = ylist_r[i + inum*k]; 
//             rinv = ylist_i[i + inum*k];
//             //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
//             dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
//             dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
//             dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;            
//             k += 1;
//         }                   
//     }
//     
//     for (mb = 1; 2*mb <= twojmax; mb++) {     
//         for (ma = 0; ma < 2*mb; ma++) {                      
//             Pr[ma] = Qr[ma];
//             Pi[ma] = Qi[ma];
//             Prx[ma] = Qrx[ma];
//             Pix[ma] = Qix[ma];
//             Pry[ma] = Qry[ma];
//             Piy[ma] = Qiy[ma];
//             Prz[ma] = Qrz[ma];
//             Piz[ma] = Qiz[ma];            
//         }                
//         for (j = 2*mb; j <= twojmax; j++) { 
//             ma = 0;
//             u_r = Pr[ma];
//             u_i = Pi[ma];
//             ux_r = Prx[ma];
//             ux_i = Pix[ma];            
//             uy_r = Pry[ma];
//             uy_i = Piy[ma];            
//             uz_r = Prz[ma];
//             uz_i = Piz[ma];                                
//             rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//             Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
//             Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);            
//             Prx[ma] = rootpq * (dardx * u_r + daidx * u_i + a_r * ux_r + a_i * ux_i);
//             Pix[ma] = rootpq * (dardx * u_i - daidx * u_r + a_r * ux_i - a_i * ux_r);
//             Pry[ma] = rootpq * (dardy * u_r + daidy * u_i + a_r * uy_r + a_i * uy_i);
//             Piy[ma] = rootpq * (dardy * u_i - daidy * u_r + a_r * uy_i - a_i * uy_r);
//             Prz[ma] = rootpq * (dardz * u_r + daidz * u_i + a_r * uz_r + a_i * uz_i);
//             Piz[ma] = rootpq * (dardz * u_i - daidz * u_r + a_r * uz_i - a_i * uz_r);                                
//             for (ma = 1; ma < j; ma++) {
//                 w_r = Pr[ma];
//                 w_i = Pi[ma];
//                 wx_r = Prx[ma];
//                 wx_i = Pix[ma];            
//                 wy_r = Pry[ma];
//                 wy_i = Piy[ma];            
//                 wz_r = Prz[ma];
//                 wz_i = Piz[ma];                                            
//                 rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
//                 rcutij = rootpqarray[ma*jdim + (j - mb)];
//                 Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
//                 Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
//                 Prx[ma] = rootpq * (dardx * w_r + daidx * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
//                 Pix[ma] = rootpq * (dardx * w_i - daidx * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
//                 Pry[ma] = rootpq * (dardy * w_r + daidy * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
//                 Piy[ma] = rootpq * (dardy * w_i - daidy * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
//                 Prz[ma] = rootpq * (dardz * w_r + daidz * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
//                 Piz[ma] = rootpq * (dardz * w_i - daidz * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);                            
//                 u_r = w_r;
//                 u_i = w_i;
//                 ux_r = wx_r;
//                 ux_i = wx_i;
//                 uy_r = wy_r;
//                 uy_i = wy_i;
//                 uz_r = wz_r;
//                 uz_i = wz_i;                
//             }
//             ma = j;
//             rcutij = rootpqarray[ma*jdim + (j - mb)];
//             Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
//             Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);       
//             Prx[ma] =-rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
//             Pix[ma] =-rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
//             Pry[ma] =-rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
//             Piy[ma] =-rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
//             Prz[ma] =-rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
//             Piz[ma] =-rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);
//             
//             if (j==(2*mb)) {
//                 mapar = 1;
//                 for (ma = 0; ma <= j/2; ma++) {
//                     if (mapar == 1) {                    
//                         Pr[j/2+ma] = Pr[j/2-ma];
//                         Pi[j/2+ma] = -Pi[j/2-ma];
//                         Prx[j/2+ma] = Prx[j/2-ma];
//                         Pix[j/2+ma] = -Pix[j/2-ma];
//                         Pry[j/2+ma] = Pry[j/2-ma];
//                         Piy[j/2+ma] = -Piy[j/2-ma];
//                         Prz[j/2+ma] = Prz[j/2-ma];
//                         Piz[j/2+ma] = -Piz[j/2-ma];                        
//                     } else {
//                         Pr[j/2+ma] = -Pr[j/2-ma];
//                         Pi[j/2+ma] = Pi[j/2-ma];
//                         Prx[j/2+ma] = -Prx[j/2-ma];
//                         Pix[j/2+ma] =  Pix[j/2-ma];
//                         Pry[j/2+ma] = -Pry[j/2-ma];
//                         Piy[j/2+ma] =  Piy[j/2-ma];
//                         Prz[j/2+ma] = -Prz[j/2-ma];
//                         Piz[j/2+ma] =  Piz[j/2-ma];                        
//                     }
//                     mapar = -mapar;        
//                 }                                                        
//             }
//             
//             if (j==(2*mb)) {
//                 mapar = 1;
//                 for (ma = 0; ma <= j; ma++) {
//                     if (mapar == 1) {                    
//                         Prx[j/2+ma] = Prx[j/2-ma];
//                         Pix[j/2+ma] = -Pix[j/2-ma];
//                         Pry[j/2+ma] = Pry[j/2-ma];
//                         Piy[j/2+ma] = -Piy[j/2-ma];
//                         Prz[j/2+ma] = Prz[j/2-ma];
//                         Piz[j/2+ma] = -Piz[j/2-ma];                        
//                     } else {
//                         Prx[j/2+ma] = -Prx[j/2-ma];
//                         Pix[j/2+ma] =  Pix[j/2-ma];
//                         Pry[j/2+ma] = -Pry[j/2-ma];
//                         Piy[j/2+ma] =  Piy[j/2-ma];
//                         Prz[j/2+ma] = -Prz[j/2-ma];
//                         Piz[j/2+ma] =  Piz[j/2-ma];                        
//                     }
//                     mapar = -mapar;        
//                 }                                                        
//             }
//                         
//             // store Qr, Qi, for the next mb level
//             if (j==(2*mb+1)) {
//                 mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
//                 for (ma = 0; ma <= j; ma++) {     
//                     if (mapar == 1) {                    
//                         Qr[j-ma] = Pr[ma];
//                         Qi[j-ma] = -Pi[ma];  
//                         Qrx[j-ma] =  Prx[ma];
//                         Qix[j-ma] = -Pix[ma];
//                         Qry[j-ma] =  Pry[ma];
//                         Qiy[j-ma] = -Piy[ma];
//                         Qrz[j-ma] =  Prz[ma];
//                         Qiz[j-ma] = -Piz[ma];                                            
//                     } else {
//                         Qr[j-ma] = -Pr[ma];
//                         Qi[j-ma] =  Pi[ma];
//                         Qrx[j-ma] = -Prx[ma];
//                         Qix[j-ma] =  Pix[ma];
//                         Qry[j-ma] = -Pry[ma];
//                         Qiy[j-ma] =  Piy[ma];
//                         Qrz[j-ma] = -Prz[ma];
//                         Qiz[j-ma] =  Piz[ma];                                            
//                     }
//                     mapar = -mapar;
//                 }                                                
//             }
//             
//             k =  1 + (j+1)*mb;
//             for (ma = 2; ma <= j; ma++)
//                 k += ma*ma;                            
//             if (j==(2*mb)) {
//                 for (ma = 0; ma < mb; ma++) {
//                     rsq = ylist_r[i + inum*k]; 
//                     rinv = ylist_i[i + inum*k];         
//                     //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
//                     dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
//                     dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
//                     dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
//                     k += 1; 
//                 }     
//                 ma = mb;
//                 rsq = 0.5*ylist_r[i + inum*k]; 
//                 rinv = 0.5*ylist_i[i + inum*k];        
//                 //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
//                 dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
//                 dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
//                 dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
//             }
//             else {
//                 for (ma = 0; ma <= j; ma++) {
//                     rsq = ylist_r[i + inum*k]; 
//                     rinv = ylist_i[i + inum*k];             
//                     //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
//                     dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
//                     dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
//                     dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
//                     k += 1; 
//                 }                
//             }            
//         }
//     }      
// 
//     ma = ai[ij];        
//     mb = aj[ij];                    
//     
//     dedx = 2.0*dedx;      
//     dedy = 2.0*dedy;      
//     dedz = 2.0*dedz;             
//     fatom[0+3*ma] += dedx;
//     fatom[1+3*ma] += dedy;
//     fatom[2+3*ma] += dedz;
//     fatom[0+3*mb] -= dedx;
//     fatom[1+3*mb] -= dedy;
//     fatom[2+3*mb] -= dedz;        
//   }   
// }
// 
// void CPOD::pod4body_force(double *force, double *coeffelem, double *rij, double *tmpmem, int *atomtype, 
//         int *aii, int *ai, int *aj, int *ti, int *tj, int natom, int Nij)            
// {
//     int inum = natom;
//     
//     //int dim = 3;    
//     //int idxcg_max = sna.idxcg_max;
//     int idxu_max = sna.idxu_max;
//     int idxb_max = sna.idxb_max;
//     int idxz_max = sna.idxz_max;    
//     int twojmax = sna.twojmax;
//     int ncoeff = sna.ncoeff;
//     //int ncoeffall = sna.ncoeffall;
//     //int ntypes = sna.ntypes;
//     int nelements = sna.nelements;    
//     //int ndoubles = sna.ndoubles;   
//     //int ntriples = sna.ntriples;   
//     int bnormflag = sna.bnormflag;
//     int chemflag = sna.chemflag;    
//     int switchflag = sna.switchflag;    
//     //int bzeroflag = sna.bzeroflag;
//     int wselfallflag = sna.wselfallflag;
//     int nelem = (chemflag) ? nelements : 1;
//     
//     int *map = sna.map;
//     int *idxz = sna.idxz;
//     //int *idxz_block = sna.idxz_block;
//     //int *idxb = sna.idxb;
//     int *idxb_block = sna.idxb_block;
//     int *idxu_block = sna.idxu_block;
//     int *idxcg_block = sna.idxcg_block;   
//     
//     double wself = sna.wself;
//     double rmin0 = sna.rmin0;
//     double rfac0 = sna.rfac0;
//     double rcutfac = sna.rcutfac;
//     //double rcutmax = sna.rcutmax;        
//     //double *bzero = sna.bzero;
//     double *rootpqarray = sna.rootpqarray;
//     double *cglist = sna.cglist;
//     //double *rcutsq = sna.rcutsq;    
//     double *radelem = sna.radelem;
//     double *wjelem = sna.wjelem; 
//     
//     int ne = 0;
//     double *ulisttot_r = &tmpmem[ne];
//     ne += idxu_max*nelements*natom;
//     double *ulisttot_i = &tmpmem[ne];    
//     ne += idxu_max*nelements*natom;        
//     double *ylist_r = &tmpmem[ne]; 
//     ne += idxu_max*nelements*natom;
//     double *ylist_i = &tmpmem[ne];     
//     
//     podArraySetValue(ulisttot_r, 0.0, natom*idxu_max*nelem);
//     podArraySetValue(ulisttot_i, 0.0, natom*idxu_max*nelem);  
// 
//     snapComputeUi(ulisttot_r, ulisttot_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
//           idxu_block, map, aii, ti, tj, twojmax, idxu_max, natom, Nij, switchflag, chemflag);                
// 
//     snapAddWself2Ui(ulisttot_r, ulisttot_i, wself, idxu_block, atomtype,
//             map, ai, wselfallflag, chemflag, idxu_max, nelem, twojmax, natom);        
// 
//     snapComputeYi(ylist_r, ylist_i, ulisttot_r, ulisttot_i, cglist, coeffelem, map, atomtype, 
//           idxz, idxb_block, idxu_block, idxcg_block, twojmax, idxb_max, idxu_max, idxz_max, 
//           nelem, ncoeff, bnormflag, natom);                
// 
//     snapComputeFi(force, ylist_r, ylist_i, rootpqarray, rij, wjelem, radelem, rmin0, rfac0, rcutfac, 
//         idxu_block, map, aii, ai, aj, ti, tj, twojmax, idxu_max, natom, inum, Nij, switchflag, chemflag);        
// }
