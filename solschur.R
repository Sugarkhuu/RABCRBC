
# installed package QZ
# 
setwd("C:/forr/ABCRBC")
kbar<-12.6698
ybar<-1.2353
rbar<-0.0351
cbar<-0.9186
hbar<-1/3
delta<-0.025
theta<-.36
beta<-.99
gamma<-.95

B=matrix(c(kbar, 0, -ybar, 0, 0,
   0, 1, 0, 0, 0,
   0, -1, theta, 0, 0, 
   0, 0, 1, 0, 0,
   0, 0, 0, 1, -rbar*beta),nrow=5, ncol=5, byrow = T)

A=matrix(c((1-delta)*kbar, 0, 0, -cbar, 0,
   0, gamma, 0, 0, 0,
   theta, 0, 0, -(1-theta), 0,
   1, 0, 0, 0, 1,
   0, 0, 0, 1, 0
   ), nrow=5, ncol=5, byrow = T)
G<-matrix(c(0, 1, 0, 0, 0),nrow=5,ncol=1)

ret<-ordqz (A,B, keyword ="udi")
Qp<-t(ret$Q)
SS<-ret$S
Zp<-t(ret$Z)
TT<-ret$T

dim<-dim(Z)
a<-dim[1]
b<-dim[2]

# https://stackoverflow.com/questions/1826519/how-to-assign-from-a-function-which-returns-more-than-one-value
# assigning results to some symbols is problematic?

nx=2
N=solve(Zp[c((a-nx+1):a),c((a-nx+1):a)])%*% Zp[c((a-nx+1):a),c(1:(a-nx))]
L=solve(Zp[c((a-nx+1):a),c((a-nx+1):a)])%*%solve(SS[c((a-nx+1):a),c((a-nx+1):a)])
L=L%*%(Qp[c((a-nx+1):a),c(1:(a-nx))])%*%G[c(1:(a-nx)),1]+Qp[c((a-nx+1):a),c((a-nx+1):a)]%*%G[c((a-nx+1):a),1]
invBBN=solve(B[c(1:(a-nx)),c(1:(a-nx))]-B[c(1:(a-nx)),c((a-nx+1):a)]%*%N)
C=invBBN%*%(A[c(1:(a-nx)),c(1:(a-nx))]-A[c(1:(a-nx)),c((a-nx+1):a)]%*%N)
D=invBBN%*%(G[c(1:(a-nx)),1]-A[c(1:(a-nx)),c((a-nx+1):a)]%*%L)

print(N)
print(L)
print(C)
print(D)

