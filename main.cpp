using namespace std;
#include <cmath>
#include <iostream>
#include <armadillo>
#include <time.h>
using namespace arma;


double offdiag(int *k , int *l , int n ); // Finner maxverdi for ikke-diagonal verdi
void Jacobi_rotate(int k , int l , int n );
void initialize1();

int k, l, n = 200;
double rho_min = 0.0, rho_max = 5.0; // Dette gir best verdi for å finne egenverdiene
double h = (rho_max - rho_min)/(n+1); // Må huske på å indeksere med N+1 ikke N
double omega = 5.0; // Vi tester også med 0.01 0.5 1.0 og 5.0

mat A(n, n);
vec rho_i(n);
vec rho2_i(n);

FILE * f;
int main()
{
    double tolerance = 1.0E-8;
    unsigned int iterations = 0, i;
    double max_iter = n*n*n;
    clock_t t;

    for( i = 0; i < n; i++)
         rho_i(i) = rho_min + (i+1)*h; // Husk på (i+1) ikke bare i

    initialize1();


    t = clock();
    double max_offdiag = offdiag(&k, &l, n);

    while(max_offdiag > tolerance && (double)iterations <=max_iter )
    {
        Jacobi_rotate (k, l, n);
        iterations++;
        max_offdiag = offdiag (&k , &l , n ) ;
    }

    rho2_i = sort(A.diag());
    t = clock() -t;

    f = fopen("log.txt", "a+");
    fprintf(f, "Execution time %.5f for N=%d w = %f\n", (float) t/CLOCKS_PER_SEC, n, omega);
    for(i = 0; i < 3; i++)
             fprintf(f, "%f\n", rho2_i[i]);

    fprintf(f, "Iterations: %d\n", iterations);
    fclose(f);
    return 0;
}

double offdiag(int *k , int *l , int n )
{
    double max = 0.0;
    for (int i = 0; i < n ; i++ )
    {
        for ( int j = i +1; j < n ; j++ )
        {
            double aij = fabs(A(i, j)) ;
            if (aij > max)
            {
                max = aij ; *k= i ; *l = j ;
             }
        }
    }

    return max ;
}

void Jacobi_rotate(int k , int l , int n )
{
    double s , c ;
    double a_kk, a_ll, a_ik, a_il;

    if (A(k, l) != 0.0 )
    {
        double t , tau ;
        tau = (A(l,l) - A(k,k))/(2*A(k, l));
        if (tau >= 0 )
         {
                t = 1.0/(tau + sqrt(1.0 + tau *tau)) ;
         }
        else
        {
            t = -1.0/(-tau + sqrt(1.0 + (tau *tau))) ;
        }
         // fprintf(stdout, "%f\n", atan(t) *180/3.1415);

        c = 1.0/sqrt(1+(t*t)) ;
        s = c*t ;
    }
    else
    {
        c = 1.0 ;
        s = 0.0 ;
    }

    a_kk = A(k, k);
    a_ll = A(l,l);
    A(k, k) = c*c*a_kk - 2.0* c*s*A(k, l)+ s*s*a_ll;
    A(l, l) = s*s*a_kk + 2.0* c*s*A(k, l)+ c*c*a_ll;
    A(k, l) = 0.0;
    A(l, k) = 0.0;

    for (int i = 0; i < n ; i ++ )
    {

         if (i != k && i != l )
         {
               a_ik = A(i, k) ;
               a_il = A(i, l) ;
               A(i, k) = c*a_ik - s*a_il ;
               A(k, i) = A(i, k) ;
               A(i, l) = c*a_il + s*a_ik ;
               A(l, i) = A(i, l);
          }
      }

    return;
}


void initialize1()
{
    int i , j;

    for( i = 0; i < n; i++)
         for(j = 0; j < n; j++)
        {
               if(i == j)
                {
                   A(i, j) = 2.0/(h*h) + omega*omega*rho_i(i)*rho_i(i) + 1.0/rho_i(i);
                }
                 else if((j -1 == i) || (i-1 == j))
                {
                    A(i, j) = -1.0/(h*h);

                 }
                else A(i, j) = 0.0;
        }

}



