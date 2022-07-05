// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include<math.h>
#include <algorithm>//min based on this headfile

using namespace Rcpp;
using namespace arma;
using namespace std;//min based on this namespace

vec Epsilon(mat Xn,vec Yn,mat Wn,double rho,vec beta);
vec Expectile(mat Xn,vec Yn,mat Wn,double tau,vec beta,double rho);
field<mat> MHS(mat Xn,vec Yn,mat Wn,double tau,int m,vec beta0,
               double S,double C,double D,double sigma0,double rho0,mat B0);

// [[Rcpp::export]]
vec Epsilon(mat Xn,vec Yn,mat Wn,double rho,vec beta) {
  int n=Yn.n_elem;
  vec temp(n);
  temp = Yn-rho*Wn*Yn-Xn*beta;
  return temp;
}

// [[Rcpp::export]]
vec Expectile(mat Xn,vec Yn,mat Wn,double tau,vec beta,double rho){
  int n=Yn.n_elem;
  //vec temp = compare(Xn,Yn,Wn,rho,beta);
  mat I(n, n, fill::eye);
  uvec temp1=(I-rho*Wn)*Yn>Xn*beta;
  vec temp2 = conv_to<vec>::from(temp1);
  vec vecpow = temp2.for_each([](vec::elem_type& val) { val=pow(-1.0,val); });
  vec comoperator = tau*vecpow+conv_to<vec>::from(temp1);
  return comoperator;
}

// [[Rcpp::export]]   
field<mat> MHS(mat Xn,vec Yn,mat Wn,double tau,int m,vec beta0,
               double S,double C,double D,double sigma0,double rho0,mat B0){
  int n=Yn.n_elem;
  int i;
  int d=Xn.n_cols + 2;
  int ks=0, kbeta=0, krho=0; 
  mat BETA(m,d,fill::zeros); 
  field <mat> Result(6,1);
  vec betaold(d-2),betanew(d-2),t0(n),g1(n),g2(n),t1(n),t2(n);
  double Te,U,sigmaold,sigmanew,rhoold,rhonew;
  double te1,te2;
  mat A1(n,n),A2(n,n);
  mat I(n, n, fill::eye);
  tau = 1-tau;
  vec temp(d, fill::none); temp.fill(0.5);
  BETA.row(0) = temp.t();
  for(i=1;i<m;i++){
    //arma_rng::set_seed_random();//Ensure that the random number generated each time is different
    // but  when called from R, the RNG seed has to be set at the R level via set.seed()
    betaold = BETA(i-1,span(2,d-1)).t();//vector is stored by column,while here needs a row
    rhoold = BETA(i-1,0);
    sigmaold = BETA(i-1,1);
    double rnorm = randn();

    sigmanew = sigmaold + S*rnorm;
    t0 = Epsilon(Xn,Yn,Wn,rhoold,betaold);
    te1 = -log(sigmaold)*n/2.0 -det(t0.t()*diagmat(Expectile(Xn,Yn,Wn,tau,betaold,rhoold))*t0)/sigmaold -0.5*pow((sigmaold-sigma0),2)/10000;//not I
    te2 = -log(sigmanew)*n/2.0 -det(t0.t()*diagmat(Expectile(Xn,Yn,Wn,tau,betaold,rhoold))*t0)/sigmanew -0.5*pow((sigmanew-sigma0),2)/10000;
    if(te2>=te1){Te=2;}
    else{Te=exp(te2-te1);}
    Te = min(1.0,Te);
    U = randu();
    if(U<Te){
      BETA(i,1)=sigmanew;
      ks=ks+1;
    }
    else{BETA(i,1)=sigmaold;}

    vec tempv(d-2,fill::randn);
    betanew = betaold + C*tempv;
    t1 = Epsilon(Xn,Yn,Wn,rhoold,betaold);
    t2 = Epsilon(Xn,Yn,Wn,rhoold,betanew);
    te1 =  -det(t1.t()*diagmat(Expectile(Xn,Yn,Wn,tau,betaold,rhoold))*t1)/sigmanew -0.5*det((betaold-beta0).t()*B0.i()*(betaold-beta0));
    te2 =  -det(t2.t()*diagmat(Expectile(Xn,Yn,Wn,tau,betanew,rhoold))*t2)/sigmanew -0.5*det((betanew-beta0).t()*B0.i()*(betanew-beta0));
    if(te2>=te1){Te=2;}
    else{Te=exp(te2-te1);}
    Te = min(1.0,Te);
    U = randu();
    if(U<Te){
      BETA(i,span(2,d-1))=betanew.t();
      kbeta=kbeta+1;
    }
    else{BETA(i,span(2,d-1))=betaold.t();}

    rhonew = rhoold + D*rnorm;
    A1 = I-rhoold*Wn;
    A2 = I-rhonew*Wn;
    g1 = Epsilon(Xn,Yn,Wn,rhoold,betanew);
    g2 = Epsilon(Xn,Yn,Wn,rhonew,betanew);
    te1 =  log(det(A1))-det(g1.t()*diagmat(Expectile(Xn,Yn,Wn,tau,betanew,rhoold))*g1)/sigmanew -0.5*pow((rhoold-rho0),2)/10000;
    te2 =  log(det(A2))-det(g2.t()*diagmat(Expectile(Xn,Yn,Wn,tau,betanew,rhonew))*g2)/sigmanew -0.5*pow((rhonew-rho0),2)/10000;
    if(te2>=te1){Te=2;}
    else{Te=exp(te2-te1);}
    Te = min(1.0,Te);
    U = randu();
    if(U<Te){
      BETA(i,0)=rhonew;
      krho=krho+1;
    }
    else{BETA(i,0)=rhoold;}
  }
  // in armadillo ,0 stands for column and 1 stands for row cause the vector and restored by row
  Result(0,0)=mean(BETA.rows(m/2,m-1),0);//theta_hat
  Result(1,0)=ks;
  Result(2,0)=kbeta;
  Result(3,0)=krho;
  Result(4,0)=stddev(BETA.rows(m/2,m-1),0,0);//sdtheta_hat
  Result(5,0)=BETA;
  Result(5,0).save("mat_field");
  //stddev( M, norm_type, dim )
  //the default norm_type=0 performs normalisation using N-1(where N is the number of samples)
  //using norm_type=1 performs normalisation using N , which provides the second moment around the mean
  return Result;
}

// [[Rcpp::export]]
mat BETA(){
  mat temp;
  temp.load("mat_field");
  return temp;
}