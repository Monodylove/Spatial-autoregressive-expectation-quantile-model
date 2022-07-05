Type "sourceCpp("MH.cpp")" in R to import some functions in this program into R. 

Rexport：
①：Epsilon(independent variable matrix Xn,dependent variable vector Yn,parameter matrix Wn,double rho,parameter vector beta)
      Return the residual vector of a Spatial Autoregressive Quantile Model 

②：Expectile(independent variable matrix Xn,dependent variable vector Yn,parameter matrix Wn,double tau,parameter vector beta,double rho)
      Return the main diagonal elements of the loss function discriminant matrix,which is the $|\tau-I(u<0)|$ in $P_\tau(u)=u^2|\tau-I(u<0)|$.

③：MHS(independent variable matrix Xn,dependent variable vector Yn,parameter matrix Wn,double tau,int m,initial vector beta0,double S,double C,double D,double sigma0,double rho0,parameter matrix B0)
      This function returns a "field",seems like the "list" in R.The output of six lists in order is:the estimate value of beta,the iteration number of sigma,the iteration number of beta,the iteration number of rho,the standard error of beta,and the last is the estimation matrix of beta. 

④：BETA()
      This function needs no parameter and it will return the sparse matrix converted by "mat_field" which is loaded by function "MHS"'s sixth element.