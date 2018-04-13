kbar<-12.6698
ybar<-1.2353
rbar<-0.0351
cbar<-0.9186
hbar<-1/3
delta<-0.025
theta<-.36
beta<-.99
gamma<-.95
A<-matrix(c(0, -kbar, 0, 0),nrow=4,ncol=1,byrow=T)
B<-matrix(c(0, (1-delta)*kbar, theta, -1),nrow=4,ncol=1,byrow=T)
C<-matrix(c(1, -1, -1, 0, ybar, -cbar, 0, 0, -1, 0, 1-theta, 0, 1, 0, 0, -1),nrow<-4,ncol <- 4,byrow <- T)
D<-matrix(c(0, 0, 1, 0),nrow=4,ncol=1,byrow=T)
F<-0
G<-F
H<-F
J<-matrix(c(0, -1, 0, beta*rbar),nrow=1,ncol=4,byrow=T)
K<-matrix(c(0, 1, 0, 0),nrow=1,ncol=4,byrow=T)
L<-F
M<-F
N<-.95
Cinv<-solve(C)
a<-F-J %*% Cinv %*% A
b=-(J%*%Cinv%*%B-G+K%*%Cinv%*%A)
c=-K%*%Cinv%*%B+H
P1=(-b+sqrt(b^2-4*a*c))/(2*a)
P2=(-b-sqrt(b^2-4*a*c))/(2*a)


if (P1<1) {
  P=P1
} else {
  P=P2
}

R=-Cinv%*%(A%*%P+B)
Q=(J%*%Cinv%*%D-L)%*%N+K%*%Cinv%*%D-M
QD=kronecker(t(N),(F-J%*%Cinv%*%A))+(J%*%R+F%*%P+G-K%*%Cinv%*%A)
Q=Q/QD
S=-Cinv%*%(A%*%Q+D)

print(P)
print(Q)
print(R)
print(S)

