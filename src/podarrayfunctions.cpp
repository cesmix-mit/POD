void CPOD::print_matrix(const char* desc, int m, int n, double **a, int lda ) 
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ ) 
    {
        for( j = 0; j < n; j++ ) printf( " %6.12f", a[j][i] );
        printf( "\n" );
    }
}

void CPOD::print_matrix(const char* desc, int m, int n, double* a, int lda ) 
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ ) 
    {
        for( j = 0; j < n; j++ ) printf( " %6.12f", a[i+j*lda] );
        printf( "\n" );
    }
}

void CPOD::print_matrix(const char* desc, int m, int n, int* a, int lda ) 
{
    int i, j;
    printf( "\n %s\n", desc );

    for( i = 0; i < m; i++ ) 
    {
        for( j = 0; j < n; j++ ) printf( " %d", a[i+j*lda] );
        printf( "\n" );
    }
}

void CPOD::podMatMul(double *c, double *a, double *b, int r1, int c1, int c2) 
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

void CPOD::podArrayFill(int* output, int start, int length) 
{	
	for (int j = 0; j < length; ++j)	
		output[j] = start + j;
}

double CPOD::podArrayMin(double *a, int n)
{
    double b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]<b)
            b = a[i];    
    return b;
}

double CPOD::podArrayMax(double *a, int n)
{
    double b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]>b)
            b = a[i];    
    return b;
}

int CPOD::podArrayMin(int *a, int n)
{
    int b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]<b)
            b = a[i];    
    return b;
}

int CPOD::podArrayMax(int *a, int n)
{
    int b = a[0];
    for (int i=1; i<n; i++)
        if (a[i]>b)
            b = a[i];    
    return b;
}

double CPOD::podArraySum(double *a, int n)
{
    double e = a[0];
    for (int i=1; i<n; i++)        
        e += a[i];    
    return e;
}

void CPOD::podKron(double *C, double *A, double *B, double alpha, int M1, int M2)
{            
    int M = M1*M2;        
    for (int idx=0; idx<M; idx++)     
    {
        int ib = idx%M2;
        int ia = (idx-ib)/M2;        
        C[idx] += alpha*A[ia]*B[ib];        
    }
}

void CPOD::podCumsum(int* output, int* input, int length) 
{
	output[0] = 0; 
	for (int j = 1; j < length; ++j)	
		output[j] = input[j - 1] + output[j - 1];	
}

double CPOD::podArrayNorm(double *a, int n)
{
    double e = a[0]*a[0];
    for (int i=1; i<n; i++)        
        e += a[i]*a[i];    
    return sqrt(e);
}

double CPOD::podArrayErrorNorm(double *a, double *b, int n)
{
    double e = (a[0]-b[0])*(a[0]-b[0]);
    for (int i=1; i<n; i++)        
        e += (a[i]-b[i])*(a[i]-b[i]);    
    return sqrt(e);
}

void CPOD::podArraySetValue(double *y, double a, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = a;        
}

void CPOD::podArrayCopy(double *y, double *x, int n)
{    
    for (int i=0; i<n; i++) 
        y[i] = x[i];        
}

void CPOD::rotation_matrix(double *Rmat, double alpha, double beta, double gamma)
{    
    double ca = cos(alpha);
    double cb = cos(beta);
    double cg = cos(gamma);
    double sa = sin(alpha);
    double sb = sin(beta);
    double sg = sin(gamma);

    Rmat[0] = ca*cg*cb - sa*sg;
    Rmat[3] = -ca*cb*sg - sa*cg;
    Rmat[6] = ca*sb;
    
    Rmat[1] = sa*cg*cb + ca*sg;
    Rmat[4] = -sa*cb*sg + ca*cg;
    Rmat[7] = sa*sb;
    
    Rmat[2] = -sb*cg;
    Rmat[5] = sb*sg;
    Rmat[8] = cb;    
}

void CPOD::matrix33_multiplication(double *xrot, double *Rmat, double *x, int natom)
{
    double x1, x2, x3;
    for (int i=0; i < natom; i++) {
        x1 = x[0 + 3*i];
        x2 = x[1 + 3*i];
        x3 = x[2 + 3*i];
        xrot[0 + 3*i] = Rmat[0]*x1 + Rmat[3]*x2 + Rmat[6]*x3;
        xrot[1 + 3*i] = Rmat[1]*x1 + Rmat[4]*x2 + Rmat[7]*x3;
        xrot[2 + 3*i] = Rmat[2]*x1 + Rmat[5]*x2 + Rmat[8]*x3;
    }
}

void CPOD::matrix33_inverse(double *invA, double *A1, double *A2, double *A3)
{                 
    double a11 = A1[0];
    double a21 = A1[1];
    double a31 = A1[2];
    double a12 = A2[0];
    double a22 = A2[1];
    double a32 = A2[2];
    double a13 = A3[0];
    double a23 = A3[1];
    double a33 = A3[2];        
    double detA = (a11*a22*a33 - a11*a23*a32 - a12*a21*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31);

    invA[0] = (a22*a33 - a23*a32)/detA;
    invA[1] = (a23*a31 - a21*a33)/detA;
    invA[2] = (a21*a32 - a22*a31)/detA;
    invA[3] = (a13*a32 - a12*a33)/detA;
    invA[4] = (a11*a33 - a13*a31)/detA;
    invA[5] = (a12*a31 - a11*a32)/detA;
    invA[6] = (a12*a23 - a13*a22)/detA;
    invA[7] = (a13*a21 - a11*a23)/detA;
    invA[8] = (a11*a22 - a12*a21)/detA;            
}

void CPOD::triclinic_lattice_conversion(double *a, double *b, double *c, double *A, double *B, double *C)
{
    double Anorm = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
    double Bnorm = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);
    double Cnorm = sqrt(C[0]*C[0] + C[1]*C[1] + C[2]*C[2]);
    
    double Ahat[3];
    Ahat[0] = A[0]/Anorm; Ahat[1] = A[1]/Anorm; Ahat[2] = A[2]/Anorm;
 
    double ax = Anorm;
    double bx = B[0]*Ahat[0] + B[1]*Ahat[1] + B[2]*Ahat[2]; //dot(B,Ahat);
    double by = sqrt(Bnorm*Bnorm - bx*bx); //sqrt(Bnorm^2 - bx^2);// #norm(cross(Ahat,B));
    double cx = C[0]*Ahat[0] + C[1]*Ahat[1] + C[2]*Ahat[2]; // dot(C,Ahat);
    double cy = (B[0]*C[0] + B[1]*C[1] + B[2]*C[2] - bx*cx)/by; // (dot(B, C) - bx*cx)/by;
    double cz = sqrt(Cnorm*Cnorm - cx*cx - cy*cy); // sqrt(Cnorm^2 - cx^2 - cy^2);    
    
    a[0] = ax; a[1] = 0.0; a[2] = 0.0;
    b[0] = bx; b[1] = by;  b[2] = 0.0;
    c[0] = cx; c[1] = cy;  c[2] = cz;
}


