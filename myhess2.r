myhess2 <- function (y,two.par = c(1,1))
{
k <- length(y)
sigma <- two.par[1]
xi <- two.par[2]
#print("Hello 1")
#print(k)
#print(sigma)
#print(xi)
#print(y)
#print("Goodbye 1")

g.00<- (sum((2*xi*y)/sigma^3)-sum((3*xi^2*y^2)/sigma^4)+sum((4*xi^3*y^3)/sigma^5)+sum((2*y)/sigma^3)-sum((3*xi*y^2)/sigma^4)+sum((4*xi^2*y^3)/sigma^5))-(k/sigma^2)

g.01<- -sum(y/sigma^2)+sum((2*xi*y^2)/sigma^3)-sum((3*xi^2*y^3)/sigma^4)+sum(y^2/sigma^3)-sum((2*xi*y^3)/sigma^4)

g.10<- -sum(y/sigma^2)+sum((2*xi*y^2)/sigma^3)-sum((3*xi^2*y^3)/sigma^4)+sum(y^2/sigma^3)-sum((2*xi*y^3)/sigma^4)

g.11<- -sum(y^2/sigma^2)+sum((2*xi*y^3)/sigma^3)+sum((2*y^3)/(3*sigma^3))


#print(g.00)
#print(g.01)
#print(g.10)
#print(g.11)
data<-c(g.00,g.01,g.10,g.11)
hessian<- matrix(data,nrow=2,byrow=T)
#cat("Hessian Matrix","\n")
#print(hessian)
#cat("","\n")
data2<-c((g.00),(g.01),(g.10),(g.11))
Information<-matrix(data2,nrow=2,byrow=T)

#cat("Information Matrix","\n")
#print(Information)
#cat("","\n")


det<-(g.00*g.11)-(g.01*g.10)
#print(det)
I.00<-(1/det)*g.11
I.01<-(1/det)*-g.01
I.10<-(1/det)*-g.10
I.11<-(1/det)*g.00
data3<-c(I.00,I.01,I.10,I.11)
Inv_Information<-matrix(data3,nrow=2,byrow=T)

#if(I.00 < 0 || I.11 < 0 )
#{
#cat("sigma=",sigma,"\n")
#cat("xi=",xi,"\n")
#cat("","\n")
#cat("Hessian Matrix","\n")
#print(hessian)
#cat("","\n")

cat("Inverse Information Matrix","\n")
print(Inv_Information)
#}
return(Inv_Information)

}
