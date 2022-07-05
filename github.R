library(Rcpp)
sourceCpp("MH.cpp")
betar=c(1,1)
R=100
m=5
n=m*R  
Bm=(rep(1,m)%*%t(rep(1,m))-diag(m))/(m-1)
Wn=kronecker(diag(R),Bm)
x1<-rnorm(n, 0, 2)
x2<-rnorm(n, 0, 2)
epsilon1 <- rnorm(n, 0, sqrt(2))
Xn<-cbind(x1,x2)
y=solve(diag(n)-0.5*Wn)%*%(Xn%*%betar+epsilon1)
Yn=y
epsilon125<-epsilon1-quantile(epsilon1,0.25)
y125=solve(diag(n)-0.5*Wn)%*%(Xn%*%betar+epsilon125)
S=0.1
C=0.1
D=0.1
sigmma0=1
B0=diag(c(1,1))*10000
beta0=c(0,0)
I<-diag(length(Yn))
rho_0=0
sa25 <- MHS(Xn,y125,Wn,0.25,1000,beta0,S,C,D,sigmma0,rho_0,B0)
sa25mat <- BETA()