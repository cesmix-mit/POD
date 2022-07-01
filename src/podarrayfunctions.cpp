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
