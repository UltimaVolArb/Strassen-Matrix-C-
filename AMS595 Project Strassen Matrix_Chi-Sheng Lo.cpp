//AMS595 Final Project: Strassen matrix in C++
//Chi-Sheng Lo 
//Student ID: 114031563
//email:  Chi-Sheng.Lo@stonybrook.edu

#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;
int n, r;
double *x, *y;
//Function for addition of matrices A and B begins
double *add(double *A, double *B, int m, int col, int str1 = 0, int stl1 = 0, int str2 = 0, int stl2 = 0)
{
    
    double *R = new double[m * m]; 
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            R[i * m + j] = A[(i + str1) * col + j + stl1] + B[(i + str2) * col + j + stl2];
    return R;
}
//Subtraction for one of the A,B matrices
double *sub(double *A, double *B, int m, int col, int str1 = 0, int stl1 = 0, int str2 = 0, int stl2 = 0)
{
    
    double *R = new double[m * m]; 
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            R[i * m + j] = A[(i + str1) * col + j + stl1] - B[(i + str2) * col + j + stl2];
    return R;
}

//Cut function: one of the matrices from A is being copied and then return one new matrix
double *cut(double *A, int m, int col, int str, int stl)
{
   //saving result in mx*m 
    double *R = new double[m * m]; 
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            R[i * m + j] = A[(i + str) * col + j + stl];
    return R;
}
//multiplication: for time complexityO(n^{3})
//Multiplication as m*m matrix
double *mul(double *A, double *B, int m)
{
    
    double *R = new double[m * m]; 
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
        {
            R[i * m + j] = 0;
            for (int k = 0; k < m; k++)
                R[i * m + j] += A[i * m + k] * B[k * m + j];
        }
    return R;
}
//generate 2 n*n matrices and random initialize range in [-10,10]
void create() 
{
    x = new double[n * n];
    y = new double[n * n];
    for (int i = 0; i < n * n; i++)
    {
        x[i] = (rand() % 20 - 10);
        y[i] = (rand() % 20 - 10);
    }
}
//Combines four (m/2)*(m/2) matrices intor one big matrix
double *merge(double *C11, double *C12, double *C21, double *C22, int m)
{
    /*
    R = [ C11   C12  ]
        [ C21   C22  ]
    */
    //R size is m*m
    double *R = new double[m * m]; 
    for (int i = 0; i < m / 2; i++)
        for (int j = 0; j < m / 2; j++)
        {
            R[i * m + j] = C11[i * m / 2 + j];
            R[i * m + j + m / 2] = C12[i * m / 2 + j];
            R[(i + m / 2) * m + j] = C21[i * m / 2 + j];
            R[(i + m / 2) * m + j + m / 2] = C22[i * m / 2 + j];
        }
    return R;
}

//Strassen matrix algo function; aka divide and conquer
//If current size is smaller than the threshold, then use matrix multiplication
//otherwise, use divide and conquer to conduct recursive calculation
double *strassen(double *A, double *B, int m)
{
    //If < then threshold, then switch to naive algo
	if (m <= r || m <= 2) 
        return mul(A, B, m);
        
    double *S1 = sub(B, B, m / 2, m, 0, m / 2, m / 2, m / 2); //S1=B12-B22
    double *S2 = add(A, A, m / 2, m, 0, 0, 0, m / 2);         //S2=A11+A12
    double *S3 = add(A, A, m / 2, m, m / 2, 0, m / 2, m / 2); //S3=A21+A22
    double *S4 = sub(B, B, m / 2, m, m / 2, 0, 0, 0);         //S4=B21-B11
    double *S5 = add(A, A, m / 2, m, 0, 0, m / 2, m / 2);     //S5=A11+A22
    double *S6 = add(B, B, m / 2, m, 0, 0, m / 2, m / 2);     //S6=B11+B22
    double *S7 = sub(A, A, m / 2, m, 0, m / 2, m / 2, m / 2); //S7=A12-A22
    double *S8 = add(B, B, m / 2, m, m / 2, 0, m / 2, m / 2); //S8=B21+B22
    double *S9 = sub(A, A, m / 2, m, 0, 0, m / 2, 0);         //S9=A11-A21
    double *S10 = add(B, B, m / 2, m, 0, 0, 0, m / 2);        //S10=B11+B12

    double *A11 = cut(A, m / 2, m, 0, 0);
    double *A22 = cut(A, m / 2, m, m / 2, m / 2);
    double *B11 = cut(B, m / 2, m, 0, 0);
    double *B22 = cut(B, m / 2, m, m / 2, m / 2);
  
    double *P1 = strassen(A11, S1, m / 2); //P1=A11*S1
    double *P2 = strassen(S2, B22, m / 2); //P2=S2*B22
    double *P3 = strassen(S3, B11, m / 2); //P3=S3*B11
    double *P4 = strassen(A22, S4, m / 2); //P4=A22*S4
    double *P5 = strassen(S5, S6, m / 2);  //P5=S5*S6
    double *P6 = strassen(S7, S8, m / 2);  //P6=S7*S8
    double *P7 = strassen(S9, S10, m / 2); //P7=S9*S10 and then calcualte the last four matrices for final result

    double *C11 = add(sub(add(P5, P4, m / 2, m / 2), P2, m / 2, m / 2), P6, m / 2, m / 2); //C11=P5+P4-P2+P6
    double *C12 = add(P1, P2, m / 2, m / 2);                                               //C12=P1+P2
    double *C21 = add(P3, P4, m / 2, m / 2);                                               //C21=P3+P4
    double *C22 = sub(sub(add(P5, P1, m / 2, m / 2), P3, m / 2, m / 2), P7, m / 2, m / 2); //C33=P5+P1-P3-P7
   
    double *R = merge(C11, C12, C21, C22, m);
    return R;
}
void show(double *a)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cout << a[i * n + j] << " ";
        cout << endl;
    }
}

//Determine whether n matches the rule, n must >2 and must be 2's exponential power
bool is2n(int n)
{
    if (n < 2)
        return false;
    while (n != 1)
    {
        if (n % 2 == 1)
            return false;
        n /= 2;
    }
    return true;
}
int main()
{
    cin >> n; //n*n matrix
	cin >> r; //Threshold
    if (!is2n(n))
    {
        cout << "n must be greater than or equal to 2, and is the exponent of 2." << endl;
        return 0;
    }
	cout<<"waiting for matrix multiplication...."<<endl;
    create(); //Create a timer
    clock_t start, finish;

    start = clock();
    double *r = strassen(x, y, n);
    finish = clock();
    //show(r);

    double duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << duration << endl;
    return 0;
}

//For the inputs I experimented: 1024 versus 2, 4, 6, 8, 16, 32, 64, 128, 256, 512, 1024, 2048. 
//In the execution pop-up exe file, just type, for example: 1024 64, then enter
