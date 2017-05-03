myhess<-function (y,two.par = c(1,1))
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

g.00<- -(k-sum(y*(2*sigma+xi*y)/(sigma+xi*y)^2)*xi-sum(y*(2*sigma+xi*y)/(sigma+xi*y)^2))/sigma^2

#print((sigma+xi*y)^2)

g.01<- -(-sum(y/(sigma+xi*y))+sum(y/(sigma+xi*y)^2)*sigma*xi+sum(y/(sigma+xi*y)^2)*sigma)/xi/sigma

g.10<- -(-sum(y/(sigma+xi*y))+sum(y/(sigma+xi*y)^2)*sigma*xi+sum(y/(sigma+xi*y)^2)*sigma)/xi/sigma

g.11<- -(-2*sum(log((sigma+xi*y)/sigma))+2*sum(y/(sigma+xi*y))*xi+sum(y^2/(sigma+xi*y)^2)*xi^3+sum(y^2/(sigma+xi*y)^2)*xi^2)/xi^3



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


#det<-(g.00*g.11)-(g.01*g.10)
#print(det)
I.00<-(1/det(Information))*g.11
I.01<-(1/det(Information))*-g.01
I.10<-(1/det(Information))*-g.10
I.11<-(1/det(Information))*g.00
data3<-c(I.00,I.01,I.10,I.11)
Inv_Information<-matrix(data3,nrow=2,byrow=T)

#if(I.00 < 0 || I.11 < 0 )
#{
#cat("sigma=",sigma,"\n")
#cat("xi=",xi,"\n")
#cat("y=",sum(y),"\n")
#cat("","\n")
#cat("Hessian Matrix","\n")
#print(hessian)
#cat("","\n")

#cat("Inverse Information Matrix","\n")
#print(Inv_Information)
#}
return(Inv_Information)

}
