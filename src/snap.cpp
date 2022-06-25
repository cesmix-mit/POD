#ifndef __SNAP
#define __SNAP

void snapBuildIndexList(int *idx_max, int *idxz, int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, 
        int *idxcg_block, int twojmax)
{
  // index list for cglist

  int jdim = twojmax + 1;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
        idxcg_block[j + j2*jdim + j1*jdim*jdim] = idxcg_count;  
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idx_max[0] = idxcg_count;
          
  int idxu_count = 0;

  for(int j = 0; j <= twojmax; j++) {
    idxu_block[j] = idxu_count;
    for(int mb = 0; mb <= j; mb++)
      for(int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  //idxu_max = idxu_count;
  idx_max[1] = idxu_count;
  
  // index list for beta and B

  int idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  int idxb_max = idxb_count;
  idx_max[2] = idxb_max;

  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          idxb[idxb_count*3 + 0] = j1;
          idxb[idxb_count*3 + 1] = j2;
          idxb[idxb_count*3 + 2] = j;  
          idxb_count++;
        }
  
  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          idxb_block[j + j2*jdim + j1*jdim*jdim] = idxb_count;    
          idxb_count++;
        }
      }

  // index list for zlist

  int idxz_count = 0;

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  int idxz_max = idxz_count;
  idx_max[3] = idxz_max;

  idxz_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
        idxz_block[j + j2*jdim + j1*jdim*jdim] = idxz_count;    

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {

            idxz[idxz_count*10 + 0] = j1;
            idxz[idxz_count*10 + 1] = j2;
            idxz[idxz_count*10 + 2] = j;
            idxz[idxz_count*10 + 3] = PODMAX(0, (2 * ma - j - j2 + j1) / 2);
            idxz[idxz_count*10 + 4] = (2 * ma - j - (2 * idxz[idxz_count*10 + 3] - j1) + j2) / 2;
            idxz[idxz_count*10 + 5] = PODMIN(j1, (2 * ma - j + j2 + j1) / 2) - idxz[idxz_count*10 + 3] + 1;
            idxz[idxz_count*10 + 6] = PODMAX(0, (2 * mb - j - j2 + j1) / 2);
            idxz[idxz_count*10 + 7] = (2 * mb - j - (2 * idxz[idxz_count*10 + 6] - j1) + j2) / 2;
            idxz[idxz_count*10 + 8] = PODMIN(j1, (2 * mb - j + j2 + j1) / 2) - idxz[idxz_count*10 + 6] + 1;
            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = idxu_block[j] + (j+1)*mb + ma;
            idxz[idxz_count*10 + 9] = jju;
              
            idxz_count++;
          }
      }
};

void snapInitRootpqArray(double *rootpqarray, int twojmax)
{
  int jdim = twojmax + 1;
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      rootpqarray[p*jdim + q] = sqrt(((double) p)/q);
};

double snapDeltacg(double *factorial, int j1, int j2, int j)
{
  double sfaccg = factorial[(j1 + j2 + j) / 2 + 1];
  return sqrt(factorial[(j1 + j2 - j) / 2] *
              factorial[(j1 - j2 + j) / 2] *
              factorial[(-j1 + j2 + j) / 2] / sfaccg);
};

void snapInitClebschGordan(double *cglist, double *factorial, int twojmax)
{
  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) {
              cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = PODMAX(0, PODMAX(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= PODMIN((j1 + j2 - j) / 2,
                          PODMIN((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (factorial[z] *
                 factorial[(j1 + j2 - j) / 2 - z] *
                 factorial[(j1 - aa2) / 2 - z] *
                 factorial[(j2 + bb2) / 2 - z] *
                 factorial[(j - j2 + aa2) / 2 + z] *
                 factorial[(j - j1 - bb2) / 2 + z]);
            }

            cc2 = 2 * m - j;
            dcg = snapDeltacg(factorial, j1, j2, j);
            sfaccg = sqrt(factorial[(j1 + aa2) / 2] *
                          factorial[(j1 - aa2) / 2] *
                          factorial[(j2 + bb2) / 2] *
                          factorial[(j2 - bb2) / 2] *
                          factorial[(j  + cc2) / 2] *
                          factorial[(j  - cc2) / 2] *
                          (j + 1));

            cglist[idxcg_count] = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
}

void snapInitSna(double *rootpqarray, double *cglist, double *factorial, int *idx_max, int *idxz, 
      int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax)
{
    snapBuildIndexList(idx_max, idxz, idxz_block, idxb, 
            idxb_block, idxu_block, idxcg_block, twojmax);
    
    snapInitRootpqArray(rootpqarray, twojmax);
    snapInitClebschGordan(cglist, factorial, twojmax);        
}

void snapComputeUij(double *Sr, double *Si, double *dSr, double *dSi, double *rootpqarray, double *rij, 
        double *wjelem, double *radelem, double rmin0, double rfac0, double rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)                
{    
  double *Srx = &dSr[0];  
  double *Sry = &dSr[idxu_max*ijnum];  
  double *Srz = &dSr[2*idxu_max*ijnum];  
  double *Six = &dSi[0];  
  double *Siy = &dSi[idxu_max*ijnum];  
  double *Siz = &dSi[2*idxu_max*ijnum];  

  for(int ij=0; ij<ijnum; ij++) {        
    double x = rij[ij*3+0];
    double y = rij[ij*3+1];
    double z = rij[ij*3+2];
    double rsq = x * x + y * y + z * z;
    double r = sqrt(rsq);
    double rinv = 1.0 / r;
    double ux = x * rinv;
    double uy = y * rinv;
    double uz = z * rinv;

    double rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    double rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    double theta0 = (r - rmin0) * rscale0;
    double z0 = r / tan(theta0);                
    double dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
            
    double sfac = 0.0, dsfac = 0.0;        
    if (switch_flag == 0) {
        sfac = 1.0;
        dsfac = 0.0;
    }
    else if (switch_flag == 1) {
        if (r <= rmin0) {
            sfac = 1.0;
            dsfac = 0.0;
        }
        else if(r > rcutij) {
            sfac = 0.0;
            dsfac = 0.0;
        }
        else {
            double rcutfac0 = M_PI / (rcutij - rmin0);
            sfac =  0.5 * (cos((r - rmin0) * rcutfac0) + 1.0);   
            dsfac = -0.5 * sin((r - rmin0) * rcutfac0) * rcutfac0;                    
//             double y = (r - rmin0)/(rcutij - rmin0);    
//             double y2 = y*y;
//             double y3 = 1.0 - y2*y;
//             double y4 = y3*y3 + 1e-6;
//             double y5 = pow(y4, 0.5);
//             double y6 = exp(-1.0/y5);
//             double y7 = pow(y4, 1.5);
//             sfac = y6/exp(-1.0);
//             dsfac = ((3.0/((rcutij - rmin0)*exp(-1.0)))*(y2)*y6*(y*y2 - 1.0))/y7;            
        }
    } 
    sfac *= wjelem[tj[ij]];
    dsfac *= wjelem[tj[ij]];

    //sfac = 1.0; 
    //dsfac = 0.0;
    
    double r0inv, dr0invdr;
    double a_r, a_i, b_r, b_i;
    double da_r[3], da_i[3], db_r[3], db_i[3];
    double dz0[3], dr0inv[3];
    double rootpq;
    int jdim = twojmax + 1;
  
    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;

    dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

    dr0inv[0] = dr0invdr * ux;
    dr0inv[1] = dr0invdr * uy;
    dr0inv[2] = dr0invdr * uz;

    dz0[0] = dz0dr * ux;
    dz0[1] = dz0dr * uy;
    dz0[2] = dz0dr * uz;

    for (int k = 0; k < 3; k++) {
        da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
        da_i[k] = -z * dr0inv[k];
    }
    da_i[2] += -r0inv;

    for (int k = 0; k < 3; k++) {
        db_r[k] = y * dr0inv[k];
        db_i[k] = -x * dr0inv[k];
    }
    db_i[0] += -r0inv;
    db_r[1] += r0inv;
    
    Sr[ij+0*ijnum] = 1.0;
    Si[ij+0*ijnum] = 0.0;
    Srx[ij+0*ijnum] = 0.0;
    Six[ij+0*ijnum] = 0.0;
    Sry[ij+0*ijnum] = 0.0;
    Siy[ij+0*ijnum] = 0.0;
    Srz[ij+0*ijnum] = 0.0;
    Siz[ij+0*ijnum] = 0.0;
    for (int j = 1; j <= twojmax; j++) {
        int jju = idxu_block[j];
        int jjup = idxu_block[j-1];
        
        // fill in left side of matrix layer from previous layer
        for (int mb = 0; 2*mb <= j; mb++) {
            Sr[ij+jju*ijnum] = 0.0;
            Si[ij+jju*ijnum] = 0.0;
            Srx[ij+jju*ijnum] = 0.0;
            Six[ij+jju*ijnum] = 0.0;
            Sry[ij+jju*ijnum] = 0.0;
            Siy[ij+jju*ijnum] = 0.0;
            Srz[ij+jju*ijnum] = 0.0;
            Siz[ij+jju*ijnum] = 0.0;
            for (int ma = 0; ma < j; ma++) {
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                int njju = ij+jju*ijnum;
                int njju1 = ij+(jju+1)*ijnum;
                int njjup = ij+jjup*ijnum;
                double u_r = Sr[njjup];
                double u_i = Si[njjup];
                double ux_r = Srx[njjup];
                double ux_i = Six[njjup];
                double uy_r = Sry[njjup];
                double uy_i = Siy[njjup];
                double uz_r = Srz[njjup];
                double uz_i = Siz[njjup];

                Sr[njju] += rootpq * (a_r * u_r + a_i * u_i);
                Si[njju] += rootpq * (a_r * u_i - a_i * u_r);
                Srx[njju] += rootpq * (da_r[0] * u_r + da_i[0] * u_i + a_r * ux_r + a_i * ux_i);
                Six[njju] += rootpq * (da_r[0] * u_i - da_i[0] * u_r + a_r * ux_i - a_i * ux_r);
                Sry[njju] += rootpq * (da_r[1] * u_r + da_i[1] * u_i + a_r * uy_r + a_i * uy_i);
                Siy[njju] += rootpq * (da_r[1] * u_i - da_i[1] * u_r + a_r * uy_i - a_i * uy_r);
                Srz[njju] += rootpq * (da_r[2] * u_r + da_i[2] * u_i + a_r * uz_r + a_i * uz_i);
                Siz[njju] += rootpq * (da_r[2] * u_i - da_i[2] * u_r + a_r * uz_i - a_i * uz_r);

                rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
                Sr[njju1] = -rootpq * (b_r * u_r + b_i * u_i);
                Si[njju1] = -rootpq * (b_r * u_i - b_i * u_r);
                Srx[njju1] = -rootpq * (db_r[0] * u_r + db_i[0] * u_i + b_r * ux_r + b_i * ux_i);
                Six[njju1] = -rootpq * (db_r[0] * u_i - db_i[0] * u_r + b_r * ux_i - b_i * ux_r);
                Sry[njju1] = -rootpq * (db_r[1] * u_r + db_i[1] * u_i + b_r * uy_r + b_i * uy_i);
                Siy[njju1] = -rootpq * (db_r[1] * u_i - db_i[1] * u_r + b_r * uy_i - b_i * uy_r);
                Srz[njju1] = -rootpq * (db_r[2] * u_r + db_i[2] * u_i + b_r * uz_r + b_i * uz_i);
                Siz[njju1] = -rootpq * (db_r[2] * u_i - db_i[2] * u_r + b_r * uz_i - b_i * uz_r);
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
                    if (j%2==1 && mb==(j/2)) {                    
                    Srx[njjup] =  Srx[njju];
                    Six[njjup] = -Six[njju];
                    Sry[njjup] =  Sry[njju];
                    Siy[njjup] = -Siy[njju];
                    Srz[njjup] =  Srz[njju];
                    Siz[njjup] = -Siz[njju];
                    }
                } else {
                    Sr[njjup] = -Sr[njju];
                    Si[njjup] =  Si[njju];
                    if (j%2==1 && mb==(j/2)) {
                    Srx[njjup] = -Srx[njju];
                    Six[njjup] =  Six[njju];
                    Sry[njjup] = -Sry[njju];
                    Siy[njjup] =  Siy[njju];
                    Srz[njjup] = -Srz[njju];
                    Siz[njjup] =  Siz[njju];                    
                    }
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }        
    }        

    for (int j = 0; j <= twojmax; j++) {
        int jju = idxu_block[j];
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            int ijk = ij+jju*ijnum;               
            Srx[ijk] = dsfac * Sr[ijk] * ux + sfac * Srx[ijk]; 
            Six[ijk] = dsfac * Si[ijk] * ux + sfac * Six[ijk]; 
            Sry[ijk] = dsfac * Sr[ijk] * uy + sfac * Sry[ijk]; 
            Siy[ijk] = dsfac * Si[ijk] * uy + sfac * Siy[ijk]; 
            Srz[ijk] = dsfac * Sr[ijk] * uz + sfac * Srz[ijk]; 
            Siz[ijk] = dsfac * Si[ijk] * uz + sfac * Siz[ijk];                  
            jju++;
          }
    }
    
    for (int k=0; k<idxu_max; k++) {
        int ijk = ij + ijnum*k;
        Sr[ijk] = sfac*Sr[ijk];
        Si[ijk] = sfac*Si[ijk];
    }            
  }
};

void snapZeroUarraytot2(double *Stotr, double *Stoti, double wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
         int twojmax, int inum)
{
    int N1 = inum;
    int N2 = N1*(twojmax+1);
    int N3 = N2*nelements;                                
    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;  // inum*(twojmax+1)
        int ii = l%N1;    // inum
        int j = (l-ii)/N1; // (twojmax+1)
        int jelem = (idx-l)/N2; // nelements   
        int ielem = (chemflag) ? map[type[ai[ii]]]: 0;                
        int jju = idxu_block[j];        
        for (int mb = 0; mb <= j; mb++) {
            for (int ma = 0; ma <= j; ma++) {
                int n = ii + inum*jju + inum*idxu_max*jelem;        
                Stotr[n] = 0.0;
                Stoti[n] = 0.0;
                if (jelem == ielem || wselfall_flag)
                    if (ma==mb)
                        Stotr[n] = wself; ///// double check this
                jju++;
            }
        }
        
    }                    
};

 void snapKernelAddUarraytot(double *Stotr, double *Stoti, double *Sr, double *Si, 
        int *map, int *ai, int *tj, int inum, int ijnum, int N1, int N2, int chemflag)
{    
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;  // ijnum
        int jju = (idx-ij)/ijnum;    // idxu_max   
        int jelem = (chemflag) ? map[tj[ij]] : 0;     
        int i = ai[ij] + inum*jju + N1*jelem;                
        Stotr[i] += Sr[idx];
        Stoti[i] += Si[idx];                
    }
};
 void snapKernelAddUarraytot(double *Stotr, double *Stoti, double *Sr, double *Si, 
        int *ai, int inum, int ijnum, int N2)
{    
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;  // ijnum
        int jju = (idx-ij)/ijnum;    // idxu_max        
        int i = ai[ij] + inum*jju;                
        Stotr[i] += Sr[idx];
        Stoti[i] += Si[idx];                
    }
};
void snapAddUarraytot(double *Stotr, double *Stoti, double *Sr, 
        double *Si, int *map, int *ai, int *tj, int idxu_max, int inum, int ijnum, int chemflag)
{   
    int N1 = inum*idxu_max;    
    int N2 = ijnum*idxu_max;    
    if (chemflag==0) {
        snapKernelAddUarraytot(Stotr, Stoti, Sr, Si, ai, inum, ijnum, N2);          
    } else
        snapKernelAddUarraytot(Stotr, Stoti, Sr, Si, map, ai, tj, inum, ijnum, N1, N2, chemflag);  
};

void snapComputeZi2(double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        double *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int inum)
{
    int jdim = twojmax + 1;    
    int N1 = inum;    
    int N2 = N1*idxz_max;
    int N3 = N2*nelements*nelements;                                
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;   //  inum*idxz_max
        int ii = l%inum;    // inum
        int jjz = (l-ii)/inum; // idxz_max
        int ielem = (idx-l)/N2;  // nelements*nelements  
        int elem2 = ielem%nelements; // nelements
        int elem1 = (ielem-elem2)/nelements; // nelements
              
        const int j1 = idxz[jjz*10+0];
        const int j2 = idxz[jjz*10+1];
        const int j = idxz[jjz*10+2];
        const int ma1min = idxz[jjz*10+3];
        const int ma2max = idxz[jjz*10+4];
        const int na = idxz[jjz*10+5];
        const int mb1min = idxz[jjz*10+6];
        const int mb2max = idxz[jjz*10+7];
        const int nb = idxz[jjz*10+8];
        int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
        int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
        int icgb = mb1min * (j2 + 1) + mb2max;

        const double *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
        double qr = 0.0;
        double qi = 0.0;          
        for (int ib = 0; ib < nb; ib++) {
            double suma1_r = 0.0;
            double suma1_i = 0.0;

            // Stotr: inum*idxu_max*nelements  
            const double *u1_r = &Stotr[ii + inum*jju1 + inum*idxu_max*elem1];
            const double *u1_i = &Stoti[ii + inum*jju1 + inum*idxu_max*elem1];
            const double *u2_r = &Stotr[ii + inum*jju2 + inum*idxu_max*elem2];
            const double *u2_i = &Stoti[ii + inum*jju2 + inum*idxu_max*elem2];

            int ma1 = ma1min;
            int ma2 = ma2max;
            int icga = ma1min * (j2 + 1) + ma2max;

            for (int ia = 0; ia < na; ia++) {
                suma1_r += cgblock[icga] * (u1_r[inum*ma1] * u2_r[inum*ma2] - u1_i[inum*ma1] * u2_i[inum*ma2]);
                suma1_i += cgblock[icga] * (u1_r[inum*ma1] * u2_i[inum*ma2] + u1_i[inum*ma1] * u2_r[inum*ma2]);
                ma1++;
                ma2--;
                icga += j2;
            } // end loop over ia

            qr += cgblock[icgb] * suma1_r;
            qi += cgblock[icgb] * suma1_i;

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
        } // end loop over ib
        
        if (bnorm_flag) {
            qr /= (j+1);
            qi /= (j+1);
        }        
        
        zlist_r[idx] = qr;
        zlist_i[idx] = qi;          
    }
};

void snapKernelComputeBi1(double *blist, double *zlist_r, double *zlist_i, 
        double *Stotr, double *Stoti, int *idxb, int *idxu_block, int *idxz_block, int jdim,         
        int nelements, int nelemsq, int nz_max, int nu_max, int inum, int N2, int N3)
{    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;        
        int ii = l%inum;
        int jjb = (l-ii)/inum;
        int jelem = (idx-l)/N2;                    
        int k = jelem%nelemsq;      
        int elem3 = k%nelements;
        int elem2 = (k-elem3)/nelements;
        int elem1 = (jelem-k)/nelemsq;    
        //int itriple = elem3 + nelements*elem2 + nelemsq*elem1;    
        int idouble = elem2 + nelements*elem1;  
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  

        int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
        int jju = idxu_block[j];
        int idu;
        int idz;
        double sumzu = 0.0;
        for (int mb = 0; 2 * mb < j; mb++)
            for (int ma = 0; ma <= j; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            } // end loop over ma, mb

        // For j even, handle middle column
        if (j % 2 == 0) {
            int mb = j / 2;
            for (int ma = 0; ma < mb; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            }
            idu = ii + inum*jju + nu_max*elem3;        
            idz = ii + inum*jjz + nz_max*idouble;        
            sumzu += 0.5 * (Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz]);
        } // end if jeven

        blist[idx] = 2.0 * sumzu;                
    }
}
void snapKernelComputeBi2(double *blist, double *bzero,int *ilist, int *type,
       int *map, int *idxb, int nelements, int nb_max, int inum, int N2, int chemflag)
{        
    for (int idx=0; idx < N2; idx++) {        
        int ii = idx%inum;        
        int jjb = (idx-ii)/inum;    
        
        int ielem = (chemflag) ? map[type[ilist[ii]]]: 0;                
        int itriple = (ielem*nelements+ielem)*nelements+ielem;

        const int j = idxb[jjb*3 + 2];  
        blist[ii + inum*jjb + nb_max*itriple] -= bzero[j];                
    }
}
void snapKernelComputeBi4(double *blist, double *bzero,
       int *idxb, int inum, int N2, int N3)
{        
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;
        int ii = l%inum;        
        int jjb = (l-ii)/inum;    
        int j = idxb[jjb*3 + 2];  
        blist[idx] -= bzero[j];        
    }
}
void snapComputeBi2(double *blist, double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        double *bzero, int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum)
{                
    int nelemsq = nelements*nelements;
    int nz_max = idxz_max*inum;
    int nu_max = idxu_max*inum;
    int nb_max = idxb_max*inum;
    int N2 = inum*idxb_max;
    int N3 = N2*nelements*nelemsq;
    int jdim = twojmax+1;

    snapKernelComputeBi1(blist, zlist_r, zlist_i, Stotr, Stoti, 
            idxb, idxu_block, idxz_block, jdim, nelements, nelemsq, nz_max, nu_max, inum, N2, N3);

    if (bzero_flag) {
        if (!wselfall_flag) {
            snapKernelComputeBi2(blist, bzero, ilist, type, map, 
                    idxb, nelements, nb_max, inum, N2, chemflag);
        }
        else {
            snapKernelComputeBi4(blist, bzero, idxb, inum, N2, N3);            
        }
    }
};

void snapComputeDbidrj(double *dblist, double *zlist_r, double *zlist_i, 
        double *dulist_r, double *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int inum, int ijnum)
{                    
    int nz_max = idxz_max*inum;
    int nb_max = idxb_max*ijnum;    
    int nu_max = idxu_max*ijnum;    
    int N2 = ijnum*idxb_max;
    int jdim = twojmax+1;

    for (int i=0; i<nb_max*3*nelements*nelements*nelements; i++)
        dblist[i] = 0.0;
    //snapArraySetValue(dblist, (double) 0.0, nb_max*3*nelements*nelements*nelements);
        
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;              
        int jjb = (idx-ij)/ijnum;                              
        int elem3 = (chemflag) ? map[tj[ij]] : 0;//(chemflag) ? map[type[alist[aj[ij]]]] : 0;
        int i = ai[ij]; // atom i
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  

       // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)
        for(int elem1 = 0; elem1 < nelements; elem1++)
            for(int elem2 = 0; elem2 < nelements; elem2++) {

            int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
            int jju = idxu_block[j];
            int idouble = elem1*nelements+elem2;
            int itriple = (elem1*nelements+elem2)*nelements+elem3;
            int nimax = nz_max*idouble;                      

            double *dbdr = &dblist[nb_max*3*itriple];
            double sumzdu_r[3];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j; mb++)
              for (int ma = 0; ma <= j; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j even, handle middle column

            if (j % 2 == 0) {
              int mb = j / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;
                jjz++;
                jju++;
              }

              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              double z_r = zlist_r[n1];
              double z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if jeven
                                
            for (int k = 0; k < 3; k++)
              dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
            
            // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
            double j1fac = (j + 1) / (j1 + 1.0);
            idouble = elem1*nelements+elem2;
            itriple = (elem3*nelements+elem2)*nelements+elem1;            
            //jjz = idxz_block[j][j2][j1];
            jjz = idxz_block[j1 + j2*jdim + j*jdim*jdim];
            jju = idxu_block[j1];

            //dbdr = &dblist[nb_max*3*itriple];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j1; mb++)
              for (int ma = 0; ma <= j1; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                       
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j1 even, handle middle column

            if (j1 % 2 == 0) {
              int mb = j1 / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                       
                jjz++;
                jju++;
              }
              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              double z_r = zlist_r[n1];
              double z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if j1even

            for (int k = 0; k < 3; k++)
              if (bnorm_flag)
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
              else
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k] * j1fac;

            // Sum over Conj(dudr(j2,ma2,mb2))*z(j,j1,j2,ma2,mb2)
            double j2fac = (j + 1) / (j2 + 1.0);
            idouble = elem2*nelements+elem1;
            itriple = (elem1*nelements+elem3)*nelements+elem2;
            //jjz = idxz_block[j][j1][j2];
            jjz = idxz_block[j2 + j1*jdim + j*jdim*jdim];        
            jju = idxu_block[j2];
            
            //dbdr = &dblist[nb_max*3*itriple];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j2; mb++)
              for (int ma = 0; ma <= j2; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                                            
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j2 even, handle middle column

            if (j2 % 2 == 0) {
              int mb = j2 / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                                                                 
                jjz++;
                jju++;
              }
              
              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              double z_r = zlist_r[n1];
              double z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if j2even

            for (int k = 0; k < 3; k++)
              if (bnorm_flag)
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
              else
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k] * j2fac;
          }        
    }
}

void snapTallyBispectrum(double *bi, double *bispectrum, int *ilist, 
        int *type, int inum, int ncoeff, int ntype)
{      
    //snapArraySetValue(bi, (double) 0.0, nperdim*ntype);
    for (int i=0; i<ncoeff*ntype; i++)
        bi[i] = 0.0;
    
    int N2 = inum*ncoeff;
    for (int idx=0; idx<N2; idx++) {        
        int ii = idx%inum;
        int icoeff = (idx-ii)/inum;
        int i = ilist[ii]; // index of atom i
        int itype = type[i]; // element type of atom i        
        int n = ncoeff*(itype-1);
        bi[icoeff+n] += bispectrum[ii + inum*icoeff];                  
    }
}

void snapTallyBispectrumDeriv(double *db, double *dbdr, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int ntype)
{   
    //snapArraySetValue(db, (double) 0.0, inum*3*nperdim*ntype);
    for (int i=0; i<inum*3*ncoeff*ntype; i++)
        db[i] = 0.0;
        
    int N2 = ijnum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ij = idx%ijnum;
        int icoeff = (idx-ij)/ijnum;        
        int i = ai[ij]; // index of atom i
        int j = aj[ij]; // index of atom i
        int itype = ti[ij]; // element type of atom i       
        int n = ncoeff*(itype-1);        
        int nii = inum*3*(icoeff + n);  
        int nij = ijnum*3*icoeff;
        
        //printf("%i %i %i %i %i %i %i %i %i %i %i \n", idx, ij, icoeff, ii, i, j, itype, n, nii, nij, quadraticflag);

        double bix = dbdr[ij + ijnum*0 + nij];
        double biy = dbdr[ij + ijnum*1 + nij];
        double biz = dbdr[ij + ijnum*2 + nij];        
        db[0 + 3*i + nii] += bix; 
        db[1 + 3*i + nii] += biy;
        db[2 + 3*i + nii] += biz;
        db[0 + 3*j + nii] -= bix;
        db[1 + 3*j + nii] -= biy;
        db[2 + 3*j + nii] -= biz;        
    }
}

void snapSetup(snastruct &sna, int twojmax, int ntypes)
{
    int backend = 1;
    sna.twojmax = twojmax;
    sna.ntypes = ntypes;
    
    int jdim = twojmax + 1;    
    int jdimpq = twojmax + 2;  
    
    TemplateMalloc(&sna.map, ntypes+1, backend);
    TemplateMalloc(&sna.idxcg_block, jdim*jdim*jdim, backend);
    TemplateMalloc(&sna.idxz_block, jdim*jdim*jdim, backend);
    TemplateMalloc(&sna.idxb_block, jdim*jdim*jdim, backend);
    TemplateMalloc(&sna.idxu_block, jdim, backend);   
    TemplateMalloc(&sna.idx_max, 5, backend);     
    
    int idxb_count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
          if (j >= j1) idxb_count++;
    
    int idxz_count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
      for(int j2 = 0; j2 <= j1; j2++)
        for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2)
          for (int mb = 0; 2*mb <= j; mb++)
            for (int ma = 0; ma <= j; ma++)
              idxz_count++;
    
    int idxcg_count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
        for(int j2 = 0; j2 <= j1; j2++)
            for(int j = j1 - j2; j <= PODMIN(twojmax, j1 + j2); j += 2) {
                for (int m1 = 0; m1 <= j1; m1++)
                  for (int m2 = 0; m2 <= j2; m2++)
                    idxcg_count++;
            }
   
    TemplateMalloc(&sna.idxz, idxz_count*10, backend);        
    TemplateMalloc(&sna.idxb, idxb_count*3, backend);        
    
    TemplateMalloc(&sna.rcutsq, (ntypes+1)*(ntypes+1), backend);
    TemplateMalloc(&sna.radelem, ntypes+1, backend);
    TemplateMalloc(&sna.wjelem, ntypes+1, backend);
        
    TemplateMalloc(&sna.rootpqarray, jdimpq*jdimpq, backend);      
    TemplateMalloc(&sna.cglist, idxcg_count, backend);        
    TemplateMalloc(&sna.bzero, jdim, backend); 
    TemplateMalloc(&sna.fac, 168, backend); 
        
    for (int i=0; i<jdimpq*jdimpq; i++)
        sna.rootpqarray[i] = 0;
    
    double facn = 1.0;
    for (int i=0; i<168; i++) {
        sna.fac[i] = facn;
        facn = facn*(i+1);
    }
    
    snapInitSna(sna.rootpqarray, sna.cglist, sna.fac, sna.idx_max, sna.idxz, 
        sna.idxz_block, sna.idxb, sna.idxb_block, sna.idxu_block, sna.idxcg_block, sna.twojmax);    

    sna.idxcg_max = sna.idx_max[0];
    sna.idxu_max = sna.idx_max[1];
    sna.idxb_max = sna.idx_max[2];
    sna.idxz_max = sna.idx_max[3];        
}

void InitSnap(snastruct &sna, double *elemradius, double *elemweight, double rcutfac, 
        double rmin0, double rfac0, int twojmax, int ntypes, int chemflag)
{
  int backend=1;  
  int bzeroflag = 0;  
  int switchflag = 1;
  int wselfallflag = 0; 
  int bnormflag = chemflag; 
  double wself=1.0;  
              
  // Calculate maximum cutoff for all elements
  double rcutmax = 0.0;
  for (int ielem = 0; ielem < ntypes; ielem++)
    rcutmax = PODMAX(2.0*elemradius[ielem]*rcutfac,rcutmax);
      
  snapSetup(sna, twojmax, ntypes);  
  TemplateCopytoDevice(&sna.radelem[1], elemradius, ntypes, backend);
  TemplateCopytoDevice(&sna.wjelem[1], elemweight, ntypes, backend);  
  cpuArrayFill(&sna.map[1], (int) 0, ntypes);
  
  double cutsq[100];
  for (int i=0; i<ntypes; i++)
      for (int j=0; j<ntypes; j++) {
          double cut = (elemradius[i] + elemradius[j])*rcutfac;
          cutsq[j+1 + (i+1)*(ntypes+1)] = cut*cut;
      }
  TemplateCopytoDevice(sna.rcutsq, cutsq, (ntypes+1)*(ntypes+1), backend);  
  
  if (bzeroflag) {
    double www = wself*wself*wself;
    double bzero[100];
    for (int j = 0; j <= twojmax; j++)
      if (bnormflag)
        bzero[j] = www;
      else
        bzero[j] = www*(j+1);
    TemplateCopytoDevice(sna.bzero, bzero, twojmax+1, backend);
  }
  
  int nelements = ntypes;
  if (!chemflag)    
    nelements = 1;
  
  sna.nelements = nelements;    
  sna.ndoubles = nelements*nelements;   // number of multi-element pairs
  sna.ntriples = nelements*nelements*nelements;   // number of multi-element triplets      
  sna.bnormflag = bnormflag;
  sna.chemflag = chemflag;    
  sna.switchflag = switchflag;
  sna.bzeroflag = bzeroflag;
  sna.wselfallflag = wselfallflag;
  sna.wself = wself;
  sna.rmin0 = rmin0;
  sna.rfac0 = rfac0;
  sna.rcutfac = rcutfac;
  sna.rcutmax = rcutmax;      
  sna.ncoeff = sna.idxb_max*sna.ntriples;
}

void snapCompute(double *blist, double *bd, snastruct &sna, neighborstruct &nb, 
        double *y, double *tmpmem, int *atomtype, int *tmpint, int natom, int Nij)            
{    
    int dim = 3;    
    //int idxcg_max = sna.idxcg_max;
    int idxu_max = sna.idxu_max;
    int idxb_max = sna.idxb_max;
    int idxz_max = sna.idxz_max;    
    int twojmax = sna.twojmax;
    int ncoeff = sna.ncoeff;
    //int ncoeffall = sna.ncoeffall;
    int ntypes = sna.ntypes;
    int nelements = sna.nelements;    
    int ndoubles = sna.ndoubles;   
    //int ntriples = sna.ntriples;   
    int bnormflag = sna.bnormflag;
    int chemflag = sna.chemflag;    
    int switchflag = sna.switchflag;    
    int bzeroflag = sna.bzeroflag;
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
    //double rcutmax = sna.rcutmax;        
    double *bzero = sna.bzero;
    double *rootpqarray = sna.rootpqarray;
    double *cglist = sna.cglist;
    //double *rcutsq = sna.rcutsq;    
    double *radelem = sna.radelem;
    double *wjelem = sna.wjelem; 
        
    int *ai = &tmpint[0];     // Nij
    int *aj = &tmpint[Nij];   // Nij 
    int *ti = &tmpint[2*Nij]; // Nij
    int *tj = &tmpint[3*Nij]; // Nij
    
    int ne = 0;
    double *rij = &tmpmem[ne]; // Nij*dim    
    ne += Nij*dim; 
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
    //double *blist = &tmpmem[ne]; // idxb_max*ntriples*natom          
    //ne += idxb_max*ntriples*natom;    
    double *dblist = &tmpmem[ne]; // idxb_max*ntriples*dim*Nij          
    double *Utotr = &tmpmem[ne];
    ne += idxu_max*nelements*natom;
    double *Utoti = &tmpmem[ne];        
    
    podNeighPairs(rij, y, ai, aj, ti, tj, nb.pairlist, nb.pairnum_cumsum, atomtype, 
                nb.alist, natom, dim);        
                
    snapComputeUij(Ur, Ui, dUr, dUi, rootpqarray, rij, wjelem, radelem, rmin0, 
         rfac0, rcutfac, idxu_block, ti, tj, twojmax, idxu_max, Nij, switchflag);
    
    snapZeroUarraytot2(Utotr, Utoti, wself, idxu_block, atomtype, map, ai, wselfallflag, 
            chemflag, idxu_max, nelem, twojmax, natom);

    snapAddUarraytot(Utotr, Utoti, Ur, Ui, map, ai, tj, idxu_max, natom, Nij, chemflag);    
        
    snapComputeZi2(Zr, Zi, Utotr, Utoti, cglist, idxz, idxu_block, 
          idxcg_block, twojmax, idxu_max, idxz_max, nelem, bnormflag, natom);

    snapComputeBi2(blist, Zr, Zi, Utotr, Utoti, bzero, nb.alist, atomtype, 
          map, idxb, idxu_block, idxz_block, twojmax, idxb_max, idxu_max, idxz_max, 
          nelem, bzeroflag,  wselfallflag, chemflag, natom);           
            
    snapComputeDbidrj(dblist, Zr, Zi, dUr, dUi, idxb, idxu_block, idxz_block, map, ai, tj, 
            twojmax, idxb_max, idxu_max, idxz_max, nelements, bnormflag, chemflag, natom, Nij);
    
    //snapTallyBispectrum(bi, blist, nb.alist, atomtype, natom, ncoeff, ntypes);    
          
    snapTallyBispectrumDeriv(bd, dblist, ai, aj, ti, natom, Nij, ncoeff, ntypes);    
                    
//    print_matrix( "SNAP descriptors:", natom, ncoeff, blist, natom); 
    
//     print_matrix( "One-body descriptors:", natom, nelements, eatom1, natom); 
//     print_matrix( "Two-body descriptors:", natom, nrbf2*3, eatom2, natom); 
//     
//     print_matrix( "element indices:", nelements, nelements, elemindex, nelements); 
//     
//     error("here");
}

#endif

