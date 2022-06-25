#ifndef __POD
#define __POD

double calculate_energyforce(double *force, double *gd, double *gdd, double *coeff, int nd, int natom)
{        
    int dim = 3;    
    int nforce = dim*natom;

    // calculate energy = gd * coeff
    double energy = 0.0;
    for (int i=0; i< nd; i++)
        energy += coeff[i]*gd[i];    
    
    // calculate force = gdd * coeff
    DGEMV(&chn, &nforce, &nd, &one, gdd, &nforce, coeff, &inc1, &zero, force, &inc1);        
    
    return energy;
}

double quadratic_energyforce(double *force, double *d2, double *d3, double *dd2, 
        double *dd3, double *coeff23, double *tmp, int *quadratic, int nc2, int nc3, int natom)
{        
    int dim = 3;    
    int nforce = dim*natom;
    int nd2 = quadratic[0]*nc2;
    int nd3 = quadratic[1]*nc3;
    
    double *c2 = &tmp[0];  
    double *c3 = &tmp[nd2];
    
    // calculate c2 = coeff23 * d3
    DGEMV(&chn, &nd2, &nd3, &one, coeff23, &nd2, d3, &inc1, &zero, c2, &inc1);            
    
    // calculate c3 = coeff23^T * d2
    DGEMV(&cht, &nd2, &nd3, &one, coeff23, &nd2, d2, &inc1, &zero, c3, &inc1);            
    
    // calculate energy = c2 * d2 
    double energy = 0.0;
    for (int i=0; i< nd2; i++)
        energy += c2[i]*d2[i];    
    
    // calculate force = force + dd2 * c2
    DGEMV(&chn, &nforce, &nd2, &one, dd2, &nforce, c2, &inc1, &one, force, &inc1);        

    // calculate force = force + dd3 * c3
    DGEMV(&chn, &nforce, &nd3, &one, dd3, &nforce, c3, &inc1, &one, force, &inc1);        
    
    return energy;
}

double quadratic_energyforce(double *force, double *d3, double *dd3, double *coeff33, 
        double *tmp, int *quadratic, int nc3, int natom)
{        
    int dim = 3;    
    int nforce = dim*natom;
    int nd3 = quadratic[0]*nc3;
        
    // calculate energy 
    double energy = 0.0;
    double *c33 = &tmp[0];  
    double *c3 = &tmp[nd3*nd3];  
    
    int m = 0;
    for (int i=0; i< nd3; i++)
        for (int j=i; j< nd3; j++) {
            energy += coeff33[m]*d3[i]*d3[j];  
            double cij = coeff33[m];
            if (i==j) c33[i + nd3*j] = 2*cij;
            else {
                c33[i + nd3*j] = cij;
                c33[j + nd3*i] = cij;
            }                
            m += 1;        
        }
                    
    // calculate c3 = c33 * d3
    DGEMV(&chn, &nd3, &nd3, &one, c33, &nd3, d3, &inc1, &zero, c3, &inc1);            
    
    // calculate force = force + dd3 * c3
    DGEMV(&chn, &nforce, &nd3, &one, dd3, &nforce, c3, &inc1, &one, force, &inc1);        
    
    return energy;
}

double cubic_energyforce(double *force, double *d3, double *dd3, double *coeff333, 
        double *tmp, int *cubic, int nc3, int natom)
{        
    int dim = 3;    
    int nforce = dim*natom;
    int nd3 = cubic[0]*nc3;
        
    // calculate energy 
    double energy = 0.0;
    double *c3 = &tmp[0];  
    for (int i=0; i< nd3; i++)
        c3[i] = 0.0;
    
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
                            
    // calculate force = force + dd3 * c3
    DGEMV(&chn, &nforce, &nd3, &one, dd3, &nforce, c3, &inc1, &one, force, &inc1);        
    
    return energy;
}

double calculate_energyforce(double *force, double *gd, double *gdd, double *coeff, double *tmp, 
        int *quadratic22, int *quadratic23, int *quadratic24, int *quadratic33, int *quadratic34,
        int *quadratic44, int nd1, int nd2, int nd3, int nd4, int nc1, int nc2, int nc3, int nc4, 
        int natom)
{        
    int dim = 3;    
    int nforce = dim*natom;
    int nd1234 = nd1+nd2+nd3+nd4;
    
    int nd22 = quadratic22[0]*quadratic22[1]*nc2*nc2;
    int nd23 = quadratic23[0]*quadratic23[1]*nc2*nc3;
    int nd24 = quadratic24[0]*quadratic24[1]*nc2*nc4;
    int nd33 = quadratic33[0]*quadratic33[1]*nc3*nc3;
    int nd34 = quadratic34[0]*quadratic34[1]*nc3*nc4;
    int nd44 = quadratic44[0]*quadratic44[1]*nc4*nc4;

    int nq;
    nq = quadratic22[0]*nc2; nd22 = nq*(nq+1)/2;
    nq = quadratic33[0]*nc3; nd33 = nq*(nq+1)/2;
    nq = quadratic44[0]*nc4; nd44 = nq*(nq+1)/2;
    
    // two-body, three-body, and four-body descriptors
    double *d2 = &gd[nd1];
    double *d3 = &gd[nd1+nd2];
    double *d4 = &gd[nd1+nd2+nd3];
    
    // two-body, three-body, and four-body descriptors derivatives
    double *dd2 = &gdd[nforce*nd1];
    double *dd3 = &gdd[nforce*(nd1+nd2)];
    double *dd4 = &gdd[nforce*(nd1+nd2+nd3)];
    
    // quadratic POD coefficients
    double *coeff22 = &coeff[nd1234];
    double *coeff23 = &coeff[nd1234+nd22];
    double *coeff24 = &coeff[nd1234+nd22+nd23];
    double *coeff33 = &coeff[nd1234+nd22+nd23+nd24];
    double *coeff34 = &coeff[nd1234+nd22+nd23+nd24+nd33];
    double *coeff44 = &coeff[nd1234+nd22+nd23+nd24+nd33+nd34];    
        
    // calculate energy and force for linear potentials
    double energy = calculate_energyforce(force, gd, gdd, coeff, nd1234, natom);    

    // calculate energy and force for quadratic22 potential
//     if (nd22>0) energy += quadratic_energyforce(force, d2, d2, dd2, dd2, 
//                             coeff22, tmp, quadratic22, nc2, nc2, natom);
    if (nd22>0) energy += quadratic_energyforce(force, d2, dd2, 
                            coeff22, tmp, quadratic22, nc2, natom);
    
    // calculate energy and force for quadratic23 potential
    if (nd23>0) energy += quadratic_energyforce(force, d2, d3, dd2, dd3, 
                            coeff23, tmp, quadratic23, nc2, nc3, natom);

    // calculate energy and force for quadratic24 potential
    if (nd24>0) energy += quadratic_energyforce(force, d2, d4, dd2, dd4, 
                            coeff24, tmp, quadratic24, nc2, nc4, natom);
    
    // calculate energy and force for quadratic33 potential
//     if (nd33>0) energy += quadratic_energyforce(force, d3, d3, dd3, dd3, 
//                             coeff33, tmp, quadratic33, nc3, nc3, natom);    
    if (nd33>0) energy += quadratic_energyforce(force, d3, dd3, 
                            coeff33, tmp, quadratic33, nc3, natom);    
    
    // calculate energy and force for quadratic34 potential
    if (nd34>0) energy += quadratic_energyforce(force, d3, d4, dd3, dd4, 
                            coeff34, tmp, quadratic34, nc3, nc4, natom);    
    
    // calculate energy and force for quadratic44 potential
//     if (nd44>0) energy += quadratic_energyforce(force, d4, d4, dd4, dd4, 
//                             coeff44, tmp, quadratic44, nc4, nc4, natom);    
    if (nd44>0) energy += quadratic_energyforce(force, d4, dd4, 
                            coeff44, tmp, quadratic44, nc4, natom);    
    
    return energy;
}

void makeindjk(int *indj, int *indk, int n)
{    
    int k1 = 1;
    for (int i=1; i<=(n-1); i++) {
        for (int j=k1; j<=(k1+n-i-1); j++) {        
            indj[j-1] = i-1;
            indk[j-1] = i+1+j-k1-1;
        }
        k1 = k1 + (n-i);
    }         
}

void makejk(double *uij, double *uik, double *wij, double *wik, double *e2ij, double *f2ij, 
        double *ei, double *f1, double *f2, double *f3, int *pairnum, int *tripletnum, 
        int *indj, int *indk, int Nj, int M, int inum, int Nij, int Nijk)
{
    for (int ii =0; ii<inum; ii++) {    
        int i = ii;     
        int n2 = pairnum[i+1] - pairnum[i]; 
        int n3 = tripletnum[i+1] - tripletnum[i];         

        for (int i1=0; i1<M; i1++)
            for (int i2=0; i2<n2; i2++) {        
                int ind1 = i2 + Nj*i1;    
                int ind2 = pairnum[i] + i2 + Nij*i1;
                ei[ind1] = e2ij[ind2];
                f1[ind1] = f2ij[0 + 3*ind2];
                f2[ind1] = f2ij[1 + 3*ind2];
                f3[ind1] = f2ij[2 + 3*ind2];
            }
            
        makeindjk(indj, indk, n2);

        for (int i1=0; i1<M; i1++)
            for (int i2=0; i2<n3; i2++) {        
                int jj = indj[i2] + Nj*i1;             
                int kk = indk[i2] + Nj*i1;                    
                int ind3 = tripletnum[i] + i2 + Nijk*i1;                
                uij[ind3] = ei[jj];
                uik[ind3] = ei[kk];
                wij[0 + 3*ind3] = f1[jj];
                wij[1 + 3*ind3] = f2[jj];
                wij[2 + 3*ind3] = f3[jj];
                wik[0 + 3*ind3] = f1[kk];                       
                wik[1 + 3*ind3] = f2[kk];                       
                wik[2 + 3*ind3] = f3[kk];         
            }
    }
}

void radialbasis(double *rbf, double *drbf, double *xij, double *besselparams, double rin, 
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
    //&rbf[0]
    //&rbf[N*besseldegree*nbesselpars]
}

void cosinbasis(double *abf, double *dabf, double *xij, double *xik, int nabf, int N)
{
    double xij1, xij2, xij3, xik1, xik2, xik3;
    double xdot, rijsq, riksq, rij, rik;
    double costhe, sinthe, theta, dtheta; 
    double tm, tm1, tm2, dct1, dct2, dct3, dct4, dct5, dct6;
    int np6;

    for (int n=0; n<N; n++) {
        xij1 = xij[0+3*n];
        xij2 = xij[1+3*n];
        xij3 = xij[2+3*n];
        xik1 = xik[0+3*n];
        xik2 = xik[1+3*n];
        xik3 = xik[2+3*n];

        xdot  = xij1*xik1 + xij2*xik2 + xij3*xik3;
        rijsq = xij1*xij1 + xij2*xij2 + xij3*xij3;    
        riksq = xik1*xik1 + xik2*xik2 + xik3*xik3;    
        rij = pow(rijsq, 0.5); 
        rik = pow(riksq, 0.5); 

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

        for (int p=0; p <nabf; p++) {
            abf[n + N*p] = cos((p+1)*theta);                
            tm = -(p+1)*sin((p+1)*theta)*dtheta;
            np6 = 6*n + 6*N*p;
            dabf[0 + np6] = tm*dct1;
            dabf[1 + np6] = tm*dct2;
            dabf[2 + np6] = tm*dct3;
            dabf[3 + np6] = tm*dct4;
            dabf[4 + np6] = tm*dct5;
            dabf[5 + np6] = tm*dct6;            
        }
    }
}

void podtally2b(double *eatom, double *fatom, double *eij, double *fij, int *ai, int *aj, 
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

void podtally3b(double *eatom, double *fatom, double *uij, double *uik, double *uijk, double *wij, 
        double *wik, double *wijk, int *ai, int *aj, int *ak, int *ti, int *tj, int *tk, int *elemindex, 
        int nelements, int nrbf, int nabf, int natom, int N)
{
    double fij0, fij1, fij2; // xij0, xij1, xij2;
    double fik0, fik1, fik2; // xik0, xik1, xik2;
    double eijk, uj, uk, ujk, wij0, wij1, wij2, wik0, wik1, wik2;
            
    int nelements2 = nelements*(nelements+1)/2;
    //int M = nrbf*(1+nabf);

    int K = -1;
    for (int m =0; m<nrbf; m++) {                    
        K += 1;
        for (int n=0; n<N; n++) {        
            int nm = n + N*m;
            uj = uij[nm];
            uk = uik[nm];
            ujk = uj*uk;
            double *wij3 = &wij[3*nm];
            double *wik3 = &wik[3*nm];
            wij0 = wij3[0];
            wij1 = wij3[1];
            wij2 = wij3[2];
            wik0 = wik3[0];
            wik1 = wik3[1];
            wik2 = wik3[2];

            eijk = uj*uk;
            fij0 = wij0*uk;
            fij1 = wij1*uk;
            fij2 = wij2*uk;
            fik0 = uj*wik0;
            fik1 = uj*wik1;
            fik2 = uj*wik2;                        

            int typei = ti[n]-1;
            int typej = tj[n]-1;
            int typek = tk[n]-1;
            //int nijk =  natom*(K + (elemindex[typej + typek*nelements] - 1)*M + (typei)*M*nelements2);
            int nijk =  natom*((elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*K);
            int nK3 = 3*nijk;

            int i1 = ai[n];
            int k = i1 + nijk;            
            int k0 = 0 + 3*i1 + nK3;
            int k1 = 1 + 3*i1 + nK3;
            int k2 = 2 + 3*i1 + nK3;            

            eatom[k] += eijk;
            fatom[k0] += fij0 + fik0;
            fatom[k1] += fij1 + fik1;
            fatom[k2] += fij2 + fik2;              
            
            i1 = aj[n];
            k0 = 0 + 3*i1 + nK3;
            k1 = 1 + 3*i1 + nK3;
            k2 = 2 + 3*i1 + nK3;            
            fatom[k0] -= fij0;
            fatom[k1] -= fij1;
            fatom[k2] -= fij2;

            i1 = ak[n];
            k0 = 0 + 3*i1 + nK3;
            k1 = 1 + 3*i1 + nK3;
            k2 = 2 + 3*i1 + nK3;
            fatom[k0] -= fik0;   
            fatom[k1] -= fik1;   
            fatom[k2] -= fik2;   
        }

        for (int p=0; p<nabf; p++) {               
            K = K + 1;                
            for (int n=0; n<N; n++) {       
                int nm = n + N*m;
                int np = n + N*p;
                int np0 = 0 + 6*n + 6*N*p;
                int np1 = 1 + 6*n + 6*N*p;
                int np2 = 2 + 6*n + 6*N*p;
                int np3 = 3 + 6*n + 6*N*p;
                int np4 = 4 + 6*n + 6*N*p;
                int np5 = 5 + 6*n + 6*N*p;      

                uj = uij[nm];
                uk = uik[nm];
                ujk = uj*uk;
                double *wij3 = &wij[3*nm];
                double *wik3 = &wik[3*nm];
                wij0 = wij3[0];
                wij1 = wij3[1];
                wij2 = wij3[2];
                wik0 = wik3[0];
                wik1 = wik3[1];
                wik2 = wik3[2];

                double u = uijk[np];   
                eijk = ujk*u;                
                fij0 = wij0*uk*u + ujk*wijk[np0];
                fij1 = wij1*uk*u + ujk*wijk[np1];
                fij2 = wij2*uk*u + ujk*wijk[np2];
                fik0 = uj*wik0*u + ujk*wijk[np3];
                fik1 = uj*wik1*u + ujk*wijk[np4];
                fik2 = uj*wik2*u + ujk*wijk[np5];
           
                int typei = ti[n]-1;
                int typej = tj[n]-1;
                int typek = tk[n]-1;
                //int nijk = natom*(K + (elemindex[typej + typek*nelements] - 1)*M + (typei)*M*nelements2);
                int nijk =  natom*((elemindex[typej + typek*nelements] - 1) + nelements2*typei + nelements2*nelements*K);
                int nK3 = 3*nijk;

                int i1 = ai[n];
                int k = i1 + nijk;
                int k0 = 0 + 3*i1 + nK3;
                int k1 = 1 + 3*i1 + nK3;
                int k2 = 2 + 3*i1 + nK3;

                eatom[k] += eijk;
                fatom[k0] += fij0 + fik0;
                fatom[k1] += fij1 + fik1;
                fatom[k2] += fij2 + fik2;
                
                i1 = aj[n];
                k0 = 0 + 3*i1 + nK3;
                k1 = 1 + 3*i1 + nK3;
                k2 = 2 + 3*i1 + nK3;
                fatom[k0] -= fij0;
                fatom[k1] -= fij1;
                fatom[k2] -= fij2;

                i1 = ak[n];
                k0 = 0 + 3*i1 + nK3;
                k1 = 1 + 3*i1 + nK3;
                k2 = 2 + 3*i1 + nK3;
                fatom[k0] -= fik0;   
                fatom[k1] -= fik1;   
                fatom[k2] -= fik2;   
           }
        }
    }
}

void matmul(double *c, double *a, double *b, int r1, int c1, int c2) 
{
    int i, j, k;

    for(j = 0; j < c2; j++)        
        for(i = 0; i < r1; i++)
            c[i + r1*j] = 0.0;        
    
    for(j = 0; j < c2; j++)
        for(i = 0; i < r1; i++)        
            for(k = 0; k < c1; k++)            
                c[i + r1*j] += a[i + r1*k] * b[k + c1*j];            
}

void pod1body(double *eatom, double *fatom, int *atomtype, int nelements, int natom)
{
    for (int m=1; m<=nelements; m++)       
        for (int i=0; i<natom; i++)         
            eatom[i + natom*(m-1)] = (atomtype[i] == m) ? 1.0 : 0.0;
        
    for (int i=0; i<3*natom*nelements; i++)  
        fatom[i] = 0.0;
}

void pod2body(double *eatom, double *fatom, double *y, double *Phi, double *besselparams, double *tmpmem, 
             double rin, double rcut, int *tmpint, int *elemindex, int *pairlist, int *pairnumsum, int *atomtype, 
             int *alist, int *pdegree, int nbesselpars, int ns, int nrbf, int nelements, int natom, int Nij)
{
    int dim = 3;
    
    double *rij = &tmpmem[0]; // 3*Nij
    int *ai = &tmpint[0];     // Nij
    int *aj = &tmpint[Nij];   // Nij 
    int *ti = &tmpint[2*Nij]; // Nij
    int *tj = &tmpint[3*Nij]; // Nij
    podNeighPairs(rij, y, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, 
                alist, natom, dim);

    double *e2ij = &tmpmem[3*Nij]; // Nij*nrbf
    double *f2ij = &tmpmem[3*Nij+Nij*nrbf]; // dim*Nij*nrbf
    double *e2ijt = &tmpmem[3*Nij+4*Nij*nrbf]; // Nij*ns
    double *f2ijt = &tmpmem[3*Nij+4*Nij*nrbf+Nij*ns]; // dim*Nij*ns    
    radialbasis(e2ijt, f2ijt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nij);
    matmul(e2ij, e2ijt, Phi, Nij, ns, nrbf);
    matmul(f2ij, f2ijt, Phi, 3*Nij, ns, nrbf);
        
    podtally2b(eatom, fatom, e2ij, f2ij, ai, aj, ti, tj, elemindex, nelements, nrbf, natom, Nij);   
}

void pod3body(double *eatom, double *fatom, double *y, double *Phi, double *besselparams, double *tmpmem, 
             double rin, double rcut, int *tmpint, int *elemindex, int *pairlist, int *pairnum, 
             int *pairnumsum, int *atomtype, int *alist, int *pdegree, int nbesselpars, int ns, 
             int nrbf, int nabf, int nelements, int natom, int Nj, int Nij, int Nijk)
{
    int dim = 3;
    
    double *rij = &tmpmem[0]; // 3*Nij
    int *ai = &tmpint[0];     // Nij
    int *aj = &tmpint[Nij];   // Nij 
    int *ti = &tmpint[2*Nij]; // Nij
    int *tj = &tmpint[3*Nij]; // Nij
    podNeighPairs(rij, y, ai, aj, ti, tj, pairlist, pairnumsum, atomtype, 
                alist, natom, dim);

    int *tripletnum = &tmpint[4*Nij];  // natom
    int *tripletnumsum = &tmpint[4*Nij+natom]; // natom+1
    int *tripletlist = &tmpint[4*Nij+2*natom+1]; // 2*Nijk
    podNeighTripletList(tripletlist, tripletnum, tripletnumsum, pairlist, 
            pairnum, pairnumsum, natom);        

    double *e2ij = &tmpmem[3*Nij]; // Nij*nrbf
    double *f2ij = &tmpmem[3*Nij+Nij*nrbf]; // dim*Nij*nrbf
    double *e2ijt = &tmpmem[3*Nij+4*Nij*nrbf]; // Nij*ns
    double *f2ijt = &tmpmem[3*Nij+4*Nij*nrbf+Nij*ns]; // dim*Nij*ns    

    radialbasis(e2ijt, f2ijt, rij, besselparams, rin, rcut-rin, pdegree[0], pdegree[1], nbesselpars, Nij);
    matmul(e2ij, e2ijt, Phi, Nij, ns, nrbf);
    matmul(f2ij, f2ijt, Phi, 3*Nij, ns, nrbf);

    int n = 3*Nij+ (1+dim)*Nij*nrbf;
    double *uij = &tmpmem[n]; // Nijk*nrbf
    double *uik = &tmpmem[n+Nijk*nrbf]; // Nijk*nrbf
    double *wij = &tmpmem[n+2*Nijk*nrbf]; // dim*Nijk*nrbf
    double *wik = &tmpmem[n+(2+dim)*Nijk*nrbf]; // dim*Nijk*nrbf

    n = 3*Nij+ (1+dim)*Nij*nrbf + 2*(1+dim)*Nijk*nrbf;
    double *ei = &tmpmem[n]; // Nj*nrbf 
    double *f1 = &tmpmem[n + Nj*nrbf]; // Nj*nrbf
    double *f2 = &tmpmem[n + 2*Nj*nrbf]; // Nj*nrbf
    double *f3 = &tmpmem[n + 3*Nj*nrbf]; // Nj*nrbf
    
    int m = 4*Nij+2*natom+1+2*Nijk;
    int *indj = &tmpint[m]; //(Nj-1)*Nj/2 
    int *indk = &tmpint[m+(Nj-1)*Nj/2]; //(Nj-1)*Nj/2
    makejk(uij, uik, wij, wik, e2ij, f2ij, ei, f1, f2, f3, pairnumsum, tripletnumsum, 
            indj, indk, Nj, nrbf, natom, Nij, Nijk);

    n = 3*Nij+ (1+dim)*Nij*nrbf + 2*(1+dim)*Nijk*nrbf + 4*Nj*nrbf;
    double *xij = &tmpmem[n]; // dim*Nijk
    double *xik = &tmpmem[n+dim*Nijk]; // dim*Nijk

    m = 4*Nij + 2*natom+1 + 2*Nijk + (Nj-1)*Nj;
    ai = &tmpint[m];     // Nijk
    aj = &tmpint[m+Nijk];   // Nijk 
    ti = &tmpint[m+2*Nijk]; // Nijk
    tj = &tmpint[m+3*Nijk]; // Nijk
    int *ak = &tmpint[m+4*Nijk]; // Nijk
    int *tk = &tmpint[m+5*Nijk]; // Nijk
    podNeighTriplets(xij, xik, y, ai, aj, ak, ti, tj, tk, tripletlist, tripletnumsum, 
        alist, atomtype,  natom, dim);                    
    
    n = 3*Nij+ (1+dim)*Nij*nrbf + 2*(1+dim)*Nijk*nrbf + 4*Nj*nrbf + 2*dim*Nijk;
    double *uijk = &tmpmem[n]; // Nijk*nabf
    double *wijk = &tmpmem[n+Nijk*nabf]; // 6*Nijk*nabf
    cosinbasis(uijk, wijk, xij, xik, nabf, Nijk);

    podtally3b(eatom, fatom, uij, uik, uijk, wij, wik, wijk, ai, aj, ak, 
        ti, tj, tk, elemindex, nelements, nrbf, nabf, natom, Nijk);   
}

void poddesc(double *eatom1, double *fatom1, double *eatom2, double *fatom2, double *eatom3, 
            double *fatom3, double *y, double *Phi2, double *Phi3, double *besselparams, 
            double *tmpmem, double rin, double rcut, int *atomtype, int *alist, int *pairlist, 
            int *pairnum, int *pairnumsum, int *elemindex, int *pdegree2, int *pdegree3, int *tmpint, 
            int nbesselpars, int nrbf2, int nrbf3, int nabf, int nelements, int Nij, int natom)
{    
    int Nj=0, Nijk=0;
    
    // one-body descriptors
    pod1body(eatom1, fatom1, atomtype, nelements, natom);
        
    //print_matrix( "One-body descriptors:", natom, nelements, eatom1, natom); 
    //print_matrix( "One-body descriptors derivarives:", 3*natom, nelements, fatom1, 3*natom); 
    
    int ns2 = pdegree2[0]*nbesselpars + pdegree2[1];
    int ns3 = pdegree3[0]*nbesselpars + pdegree3[1];

    for (int i=0; i < natom; i++) {
        Nj = (Nj > pairnum[i]) ? Nj : pairnum[i];
        Nijk +=  (pairnum[i]-1)*pairnum[i]/2;
    }

    int k = 1;
    for (int i=0; i < nelements; i++) 
        for (int j=i; j < nelements; j++) {
            elemindex[i + nelements*j] = k;
            elemindex[j + nelements*i] = k;
            k += 1;
        }            

    // two-body descriptors
    pod2body(eatom2, fatom2, y, Phi2, besselparams, tmpmem, rin, rcut, tmpint, 
             elemindex, pairlist, pairnumsum, atomtype, alist, pdegree2, 
             nbesselpars, ns2, nrbf2, nelements, natom, Nij); 

//     print_matrix( "One-body descriptors:", natom, nelements, eatom1, natom); 
//     print_matrix( "Two-body descriptors:", natom, nrbf2*3, eatom2, natom); 
    
    // three-body descriptors
    pod3body(eatom3, fatom3, y, Phi3, besselparams, tmpmem, rin, rcut, tmpint, 
             elemindex, pairlist, pairnum, pairnumsum, atomtype, alist, pdegree3, 
             nbesselpars, ns3, nrbf3, nabf, nelements, natom, Nj, Nij, Nijk);    
        
//     print_matrix( "One-body descriptors:", natom, nelements, eatom1, natom); 
//     print_matrix( "Two-body descriptors:", natom, nrbf2*3, eatom2, natom); 
//     
//     print_matrix( "element indices:", nelements, nelements, elemindex, nelements); 
//     
//     error("here");
}

void linear_descriptors(descriptorstruct &desc, neighborstruct &nb, podstruct pod, 
        snastruct sna, double *lattice, double *position, int *atomtype, int natom)
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
    
    cpuArraySetValue(eatom1, 0.0, natom*(nd1+nd2+nd3+nd4));
    cpuArraySetValue(fatom1, 0.0, dim*natom*(nd1+nd2+nd3+nd4));    
    
    // peratom descriptors for one-body, two-body, and three-body linear potentials
    poddesc(eatom1, fatom1, eatom2, fatom2, eatom3, fatom3, nb.y, Phi2, Phi3, besselparams, 
            tmpmem, rin, rcut, atomtype, nb.alist, nb.pairlist, nb.pairnum, nb.pairnum_cumsum, 
            nb.elemindex, pdegree2, pdegree3, tmpint, nbesselpars, nrbf2, nrbf3, nabf3, 
            nelements, Nij, natom);            
        
//     snapCompute(double *bi, double *bd, snastruct &sna, neighborstruct &nb, 
//         double *y, double *tmpmem, int *atomtype, int *tmpint, int natom, int Nij)            

    // global descriptors for one-body, two-body, three-body, and four-bodt linear potentials
    int nd1234 = nd1+nd2+nd3+nd4;    
    cpuArraySetValue(tmpmem, 1.0, natom);
    DGEMV(&cht, &natom, &nd1234, &one, eatom1, &natom, tmpmem, &inc1, &zero, desc.gd, &inc1);            
}

void quadratic_descriptors(double* d23, double *dd23, double* d2, double *d3, double* dd2, double *dd3, 
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

void quadratic_descriptors(double* d33, double *dd33, double *d3, double *dd3, int M3, int N)
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

void cubic_descriptors(double* d333, double *Dd333, double *d3, double *Dd3, int M3, int N)
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

void energyforce_calculation(descriptorstruct &desc, neighborstruct &nb, podstruct pod, 
        snastruct sna, datastruct data, double *coeff)
{                
    int dim = 3;
    double energy;
    double force[1+dim*data.num_atom_max];
    
    int nfiles = data.data_files.size();    // number of files    
    int nd1234 = pod.nd1 + pod.nd2 + pod.nd3 + pod.nd4; 
                                
    std::cout<<"**************** Begin of Energy/Force Calculation ****************"<<std::endl;
    
    int ci = 0; // configuration counter    
    for (int file = 0; file < nfiles; file++) { // loop over each file in the data set
        
        int nconfigs = data.num_config[file];
        for (int ii=0; ii < nconfigs; ii++) { // loop over each configuration in a file
            if ((ci % 100)==0) std::cout<<"Configuration: # "<<ci+1<<std::endl;
            
            int natom = data.num_atom[ci];
            int nforce = dim*natom;
            int natom_cumsum = data.num_atom_cumsum[ci];    
            int *atomtype = &data.atomtype[natom_cumsum];
            double *position = &data.position[dim*natom_cumsum];
            double *lattice = &data.lattice[9*ci];
                        
            // compute linear POD descriptors
            linear_descriptors(desc, nb, pod, sna, lattice, position, atomtype, natom);
                    
            // calculate energy and force
            energy = calculate_energyforce(&force[1], desc.gd, desc.gdd, coeff, &desc.gdd[nforce*nd1234], 
                        pod.quadratic22, pod.quadratic23, pod.quadratic24, pod.quadratic33, 
                        pod.quadratic34, pod.quadratic44, pod.nd1, pod.nd2, pod.nd3, pod.nd4, 
                        pod.nelements, pod.nc2, pod.nc3, pod.nc4, natom);
                        
            if (pod.nd333>0) {
                energy += cubic_energyforce(force, &desc.gd[pod.nd1+pod.nd2], &desc.gdd[nforce*(pod.nd1+pod.nd2)], 
                            &coeff[nd1234+pod.nd22+pod.nd23+pod.nd24+pod.nd33+pod.nd34+pod.nd44],
                            &desc.gdd[nforce*nd1234], pod.cubic333, pod.nc3, natom);    
            }
            
            if (pod.nd444>0) {
                energy += cubic_energyforce(force, &desc.gd[pod.nd1+pod.nd2+pod.nd3], &desc.gdd[nforce*(pod.nd1+pod.nd2+pod.nd3)], 
                            &coeff[nd1234+pod.nd22+pod.nd23+pod.nd24+pod.nd33+pod.nd34+pod.nd44+pod.nd333],
                            &desc.gdd[nforce*nd1234], pod.cubic444, pod.nc4, natom);    
            }
            
            ci += 1;             
            
            // save energy and force into a binary file
            force[0] = energy;
            string filename = "energyforce_config" + std::to_string(ci) + ".bin";
            writearray2file(filename.c_str(), force, 1 + dim*natom, 1);
        }
    }       
    std::cout<<"**************** End of Energy/Force Calculation ****************"<<std::endl<<std::endl;    
}




// void poddesc(double *eatom1, double *fatom1, double *eatom2, double *fatom2, double *eatom3, 
//             double *fatom3, double *x, double *y, double *a1, double *a2, double *a3, double *Phi2, 
//             double *Phi3, double *besselparams, double *tmpmem, double rin, double rcut, int *atomtype, 
//             int *pbc, int *alist, int *pairlist, int *pairnum, int *pairnumsum, int *elemindex, 
//             int *pdegree2, int *pdegree3, int *tmpint, int nbesselpars, int nrbf2, int nrbf3, int nabf,
//             int nelements, int natom)
// {    
//     int Nj=0, Nij=0, Nijk=0;
//     
//     // neighbor list
//     Nij = podfullneighborlist(y, alist, pairlist, pairnum, pairnumsum, x, a1, a2, a3, rcut, pbc, natom);
//     
//     // one-body descriptors
//     pod1body(eatom1, fatom1, atomtype, nelements, natom);
//     
//     if (Nij>0) {
//         int ns2 = pdegree2[0]*nbesselpars + pdegree2[1];
//         int ns3 = pdegree3[0]*nbesselpars + pdegree3[1];
//         
//         for (int i=0; i < natom; i++) {
//             Nj = (Nj > pairnum[i]) ? Nj : pairnum[i];
//             Nijk +=  (pairnum[i]-1)*pairnum[i]/2;
//         }
//         
//         int k = 1;
//         for (int i=0; i < nelements; i++) 
//             for (int j=i; j < nelements; j++) {
//                 elemindex[i + nelements*j] = k;
//                 elemindex[j + nelements*i] = k;
//                 k += 1;
//             }            
//         
//         // two-body descriptors
//         pod2body(eatom2, fatom2, y, Phi2, besselparams, tmpmem, rin, rcut, tmpint, 
//                  elemindex, pairlist, pairnumsum, atomtype, alist, pdegree2, 
//                  nbesselpars, ns2, nrbf2, nelements, natom, Nij); 
//                 
//         // three-body descriptors
//         pod3body(eatom3, fatom3, y, Phi3, besselparams, tmpmem, rin, rcut, tmpint, 
//                  elemindex, pairlist, pairnum, pairnumsum, atomtype, alist, pdegree3, 
//                  nbesselpars, ns3, nrbf3, nabf, nelements, natom, Nj, Nij, Nijk);         
//     }        
// }
// 
// void poddesc(double *eatom1, double *fatom1, double *eatom2, double *fatom2, double *eatom3, 
//             double *fatom3, double *x, double *a1, double *a2, double *a3, double *Phi2, double *Phi3, 
//             double *besselparams, double rin, double rcut, int *atomtype, int *pbc, int *pdegree2, 
//             int *pdegree3, int nbesselpars, int nrbf2, int nrbf3, int nabf, int nelements, int natom)
// {
//     int Nj=0, Nij=0, Nijk=0;
//     int m=0, n=0, p=0, dim=3;
//     if (pbc[0] == 1) m = (int) ceil(rcut/a1[0]);    
//     if (pbc[1] == 1) n = (int) ceil(rcut/a2[1]);    
//     if (pbc[2] == 1) p = (int) ceil(rcut/a3[2]);                
//     // number of lattices
//     int nl = (2*m+1)*(2*n+1)*(2*p+1);      
//     
//     double *y, *tmpmem;    
//     int *alist, *pairlist, *pairnum, *pairnumsum, *elemindex, *tmpint;
//     y = (double*) malloc (sizeof (double)*(dim*natom*nl));
//     alist = (int*) malloc (sizeof (int)*(natom*nl));    
//     pairnum = (int*) malloc (sizeof (int)*(natom));
//     pairnumsum = (int*) malloc (sizeof (int)*(natom+1));
//     pairlist = (int*) malloc (sizeof (int)*(natom*PODMIN(natom*nl,1000)));
//     elemindex = (int*) malloc (sizeof (int)*(nelements*nelements));
//     
//     // neighbor list
//     Nij = podfullneighborlist(y, alist, pairlist, pairnum, pairnumsum, x, a1, a2, a3, rcut, pbc, natom);
//     
//     // one-body descriptors
//     pod1body(eatom1, fatom1, atomtype, nelements, natom);
//     
//     if (Nij>0) {
//         int ns2 = pdegree2[0]*nbesselpars + pdegree2[1];
//         int ns3 = pdegree3[0]*nbesselpars + pdegree3[1];
//         
//         for (int i=0; i < natom; i++) {
//             Nj = (Nj > pairnum[i]) ? Nj : pairnum[i];
//             Nijk +=  (pairnum[i]-1)*pairnum[i]/2;
//         }
// 
//         int szd = 3*Nij+ (1+dim)*Nij*(nrbf3+ns3) + 2*(1+dim)*Nijk*nrbf3 + 4*Nj*nrbf3 + 2*dim*Nijk + 7*Nijk*nabf;
//         int szi = 4*Nij + 2*natom+1 + 2*Nijk + (Nj-1)*Nj + 6*Nijk;
//         tmpmem = (double*) malloc (sizeof (double)*(szd));
//         tmpint = (int*) malloc (sizeof (int)*(szi));
//         
//         int k = 1;
//         for (int i=0; i < nelements; i++) 
//             for (int j=i; j < nelements; j++) {
//                 elemindex[i + nelements*j] = k;
//                 elemindex[j + nelements*i] = k;
//                 k += 1;
//             }            
//         
//         // two-body descriptors
//         pod2body(eatom2, fatom2, y, Phi2, besselparams, tmpmem, rin, rcut, tmpint, 
//                  elemindex, pairlist, pairnumsum, atomtype, alist, pdegree2, 
//                  nbesselpars, ns2, nrbf2, nelements, natom, Nij); 
//                 
//         // three-body descriptors
//         pod3body(eatom3, fatom3, y, Phi3, besselparams, tmpmem, rin, rcut, tmpint, 
//                  elemindex, pairlist, pairnum, pairnumsum, atomtype, alist, pdegree3, 
//                  nbesselpars, ns3, nrbf3, nabf, nelements, natom, Nj, Nij, Nijk); 
//         
//         free(tmpint); free(tmpmem);
//     }
//     
//     free(y); free(alist); free(pairlist); free(pairnum); free(pairnumsum); free(elemindex);    
// }
// 
#endif

