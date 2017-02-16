#include <adolc/adolc.h>
#include <iostream>
#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;
#define XEPS 1e-4

//////////////////////////////////////////////////////////////
// Function to compute log-likelihood for bivariate copula
// Input:
// n         sample size
// u         first variable of data set
// v         second variable of data set
// theta     parameter for bivariate copula
//////////////////////////////////////////////////////////////
adouble LL(int n, adouble* u, adouble* v, adouble theta)
{
    int j;
    adouble dat[2], ll=0.0, f, ff;

    //Compute log-likelihood:
    f = 0.0;
    for(j=0;j<n;j++)
    {
        dat[0] = u[j]; dat[1] = v[j];
        //f=log1p(theta))-(1.0+theta)*log(dat[0]*dat[1])-(2.0+1.0/(theta))*log(pow(dat[0],-theta)+pow(dat[1],-theta)-1.0);
        f+=log(1+theta)-(1.0+theta)*(log(dat[0])+log(dat[1]))-(2.0+1.0/(theta))*log(pow(dat[0],-theta)+pow(dat[1],-theta)-1.0);
    }
    // ADOL-C does not allow constant 0 in condassign in its newest verion.
    // Question remains on accurate evaluation when theta is very small. Careful analysis is needed.
    //condassign(ff, theta-1e-10, f, 0); 
    ff = f;

    ll += ff;

    //if(theta <= 1e-10) ll = 0;
    //else
    //{
    //    for(j=0;j<n;j++)
    //    {
    //        dat[0] = u[j]; dat[1] = v[j];
    //        //f=log1p(theta))-(1.0+theta)*log(dat[0]*dat[1])-(2.0+1.0/(theta))*log(pow(dat[0],-theta)+pow(dat[1],-theta)-1.0);
    //        f=log(1+theta)-(1.0+theta)*log(dat[0]*dat[1])-(2.0+1.0/(theta))*log(pow(dat[0],-theta)+pow(dat[1],-theta)-1.0);
    //        ll += f;
    //    }
    //}
    return ll;
}

//[[Rcpp::export]]
double DVineAD_LL(NumericVector u, NumericVector v, double theta) {
    adouble sumloglik = 0;
    int n = u.size();
    adouble* pu = new adouble[n];
    adouble* pv = new adouble[n]; 
    for (int i=0; i<n; i++) {
        pu[i] = u[i];
        pv[i] = v[i];
    }
    sumloglik = LL(n, pu, pv, theta); 
    delete[] pu;
    delete[] pv;
    return sumloglik.value();
}

//////////////////////////////////////////////////////////////
// Function to compute h-function for vine simulation and estimation
// Input:
// n        Sample size 
// u, v     Two arrays for samples in two dimensions    
// theta    parameter for the bivariate copula
// out      output
//////////////////////////////////////////////////////////////

void Hfunc(int n, adouble* u, adouble* v, adouble theta, adouble* out)
{
    int j;
    adouble x;

    for(j=0;j<n;j++)
    {
            x = pow(u[j],-theta)+pow(v[j],-theta)-1.0 ;
            x = pow(v[j],-theta-1.0)*pow(x,-1.0-1.0/(theta));
            condassign(out[j], theta-XEPS, x, u[j]);
    }
}

// Compute the negative log-likelihood, also used for automatic differentiation
adouble computeLogLik(int n, int d, adouble*** H1, adouble*** H2, adouble** Theta, int TruncLevel)
{
    adouble loglik(0);
    for (int i=0; i<=(d-2); i++) //Level 0
    {
       loglik -= LL(n, H1[0][i], H2[0][i], Theta[0][i]); 
    }

   //for (int i=1; i<=(d-2); i++)
   for (int i=1; i<=(TruncLevel-1); i++)
   {
       for (int j=0; j<=(d-i-2); j++)
       { 
           Hfunc(n,H1[i-1][j],H2[i-1][j], Theta[i-1][j], H1[i][j]);
           Hfunc(n,H2[i-1][j+1],H1[i-1][j+1], Theta[i-1][j+1], H2[i][j]);
           loglik -= LL(n, H1[i][j], H2[i][j], Theta[i][j]);
       }
   }

    return loglik;
}

//[[Rcpp::export]]
NumericVector DVineVec2RVineVec(int d, NumericVector dvinegrad, int TruncLevel)
{ 
    int len = dvinegrad.size();
    NumericVector rvinegrad(len);
    int k;

    double** mat = new double*[TruncLevel];

    k=0;
    for (int i=0; i<=(TruncLevel-1); i++) {
        mat[i] = new double[d-i-1];
        for (int j=0; j<d-i-1; j++)
        {
            mat[i][j] = dvinegrad[k];
            k++;
        }
    }

    k=0;
    for (int j=(d-2); j>=0; j--) {
        for (int i=std::max(0,j-TruncLevel+1); i<=j; i++) {
            rvinegrad[k] = mat[j-i][i];
            k++;
        }
    }
    for (int i=0; i<=(TruncLevel-1); i++) {
        delete[] mat[i];
    }
    delete[] mat;
    return rvinegrad;
}

//[[Rcpp::export]]
void DVineAD_EstablishTape(NumericMatrix X, NumericVector theta, int TruncLevel)
{
    int d = X.ncol(); //theta.size();
    //int dd = d*(d-1)/2;
    int dd = (d-1+d-TruncLevel)*TruncLevel/2;
    if (dd!= theta.size()){
        Rcout<<"Warning: dimensions do not match!"<<std::endl;
    }
    int n = X.nrow();
    double numloglik = 0;

    trace_on(0);
    //Store theta values to a triangle matrix Theta
    //adouble** Theta = new adouble*[d-1]; 
    adouble** Theta = new adouble*[TruncLevel]; 

    int k = 0; 
    //for (int i=0; i<=(d-2); i++) {
    for (int i=0; i<=(TruncLevel-1); i++) {
        Theta[i] = new adouble[d-i-1];
        for (int j=0; j<=(d-i-2); j++) {
            Theta[i][j] <<= theta[k];
            k++;
        }
    }

    // Claim vector arrays storing pseudo-samples
    //adouble*** H1 = new adouble**[d-1];
    //adouble*** H2 = new adouble**[d-1]; 
    adouble*** H1 = new adouble**[TruncLevel];
    adouble*** H2 = new adouble**[TruncLevel]; 
    for (int i=0; i<=(TruncLevel-1); i++) {
        H1[i] = new adouble*[d-i-1];
        H2[i] = new adouble*[d-i-1];
        for (int j=0; j<=(d-i-2); j++) {
            H1[i][j] = new adouble[n];
            H2[i][j] = new adouble[n];
        }
    } 

   //First level pseudo-samples are just original samples 
   for (int j=0; j<=(d-2); j++) {
       for (int k=0; k<=(n-1); k++) {
            H1[0][j][k] = X(k,j);
            H2[0][j][k] = X(k,j+1);
       }
   }

   adouble loglik = computeLogLik(n, d, H1, H2, Theta, TruncLevel);
   loglik >>= numloglik;

    //Free H1 and H2
    for (int i=0; i<=(TruncLevel-1); i++) {
        for (int j=0; j<=(d-i-2); j++) {
            delete[] H1[i][j];
            delete[] H2[i][j];
        }
        delete[] H1[i];
        delete[] H2[i]; 
    }

    delete[] H1;
    delete[] H2;

    for (int i=0; i<=(TruncLevel-1); i++) {
        delete[] Theta[i];
    }
    delete[] Theta;

    //trace_off(1);
    trace_off();
}

//[[Rcpp::export]]
double DVineAD_LogLik(NumericVector theta)
{
    double loglik = 0.0;
    double* ptheta;
    ptheta = new double[theta.length()];
    for (int i=0; i<theta.length(); i++) {
        ptheta[i] = theta[i];
    }
    function(0,1,theta.length(), ptheta, &loglik);
    delete[] ptheta;
    return loglik;
}


//[[Rcpp::export]]
NumericVector DVineAD_LogLikGrad(NumericVector theta)
{
    int dd = theta.length();

    double* ptheta = new double[dd];
    for (int i=0; i<dd; i++) 
        ptheta[i] = theta[i];
    double* grad = new double[dd];
    gradient(0, dd, ptheta, grad);
    NumericVector outgrad(dd);
    
    for (int i=0; i<dd; i++)
        outgrad[i] = grad[i];

    delete[] ptheta;
    delete[] grad;
    
    return outgrad;
}

//[[Rcpp::export]]
NumericMatrix DVineAD_LogLikHess(NumericVector theta)
{
    int dd = theta.length();

    double* ptheta = new double[dd];
    for (int i=0; i<dd; i++) 
        ptheta[i] = theta[i];
    double** hess = new double*[dd];
    for (int i=0; i<dd; i++) { 
        hess[i] = new double[dd];
    }
    for (int i=0; i<dd; i++){
        for (int j=0; j<dd; j++) {
            hess[i][j] = 0;
        }
    }

    hessian(0, dd, ptheta, hess);
    NumericMatrix outhess(dd, dd);

    for (int i=0; i<dd; i++){
        for (int j=0; j<dd; j++) {
            outhess(i,j) = hess[i][j];
        }
    }

    delete[] ptheta;
    for (int i=0; i<dd; i++) { 
        delete[] hess[i];
    } 
    delete[] hess;
    return outhess;
}
