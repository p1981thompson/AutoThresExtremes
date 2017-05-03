paulgpd5<-function (data,threshold,maxit=10000,show=F,...)
{
z<-list()
  	n <- length(data)
##check for correct threshold input###

	if (is.function(threshold))
	stop("`threshold' cannot be a function")

##define exceedances over threshold##

    	exceed <- data[data > threshold]
   	index <- (1:n)[data > threshold]
	y <- exceed-threshold
	#y <- rgpd(10000,xi=0.25,beta=0.4,mu=1.5)

##set initial values for parameter estimates of xi and sigma##

	#initial_sig1 <- sqrt(6 * var(data))/pi
    	#initial_sig <- mean(data, na.rm = TRUE) - 0.57722 * initial_sig1

	initial_xi<- -0.2161221
	initial_sig<- 0.3425574
  
 	#initial_xi<- 0.1
	#initial_sig<- 0.1
      
	initials <- c(initial_sig,initial_xi)

	z$threshold <- threshold

	z$num_exceed <- length(exceed)
   	z$data <- exceed

## Define negative log likelihood for GPD for optimization##

gpd.lik<- function(parameters,y)
{
sigma <- parameters[1]
xi <- parameters[2]
# cat("","\n")
#print(paste("sigma=",sigma))
#print(paste("xi=",xi))

#xi: log derivative condition#
der_xi1 <- (-(-2*sum(log((sigma+xi*y)/sigma))+2*sum(y/(sigma+xi*y))*xi+sum(y^2/(sigma+xi*y)^2)*xi^3+sum(y^2/(sigma+xi*y)^2)*xi^2)/xi^3)

#sigma: log derivative condition#
der_sig1 <- (-(length(y)-sum(y*(2*sigma+xi*y)/(sigma+xi*y)^2)*xi-sum(y*(2*sigma+xi*y)/(sigma+xi*y)^2))/sigma^2)

#xi: taylor derivative condition#
der_xi2 <- (-sum(y^2/sigma^2)+sum((2*xi*y^3)/sigma^3)+sum((2*y^3)/(3*sigma^3)))

#sigma: taylor derivative condition#
der_sig2 <- (sum((2*xi*y)/sigma^3)-sum((3*xi^2*y^2)/sigma^4)+sum((4*xi^3*y^3)/sigma^5)+sum((2*y)/sigma^3)-sum((3*xi*y^2)/sigma^4)+sum((4*xi^2*y^3)/sigma^5))-(length(y)/sigma^2)

#print(paste("deriv xi 1=",der_xi1))
#print(paste("deriv xi 2=",der_xi2))
#print(paste("deriv sigma 1=",der_sig1))
#print(paste("deriv sigma 2=",der_sig2))


if( der_xi1 <= 0 || der_sig1 <= 0  || any(is.na(der_xi1)) || any(is.na(der_sig1)))
{

lik.total <- 10^6
}
else{
lik_with_log <- length(y)*logb(sigma) + sum(logb(1+y*xi/sigma)) * (1/xi + 1)
lik_with_taylor <- length(y)*logb(sigma) + (xi+1)*(sum(y)/sigma -xi*sum(y)^2/(2*sigma^2)+ xi^2*sum(y)^3/(3*sigma^3))
lik.total <- ifelse((xi < 10^(-10)) & (xi > -10^(-10)), lik_with_taylor, lik_with_log)
}
#cat("xi.parameter=",xi,"\n")
#cat("sigma.parameter=",sigma,"\n")
#cat("max(y)=",max(y),"\n")
#cat("min(y)=",min(y),"\n")

lik.total

}

## Define gradient function for negative log likelihood for GPD##

grad.lik <- function(parameters, y)
	{

	sigma <- parameters[1]
	xi <- parameters[2]

#	part1 <- 1+(y*xi)/sigma
#
#	if(any(part1 <= 0) || sigma <= 0 )
#		{
#		gradient <- c(NA,NA)
#		}
#	else{
##taylor derivatives##

	xi_derive1 <- sum((y/sigma)-((xi*y^2)/(2*sigma^2))+((xi^2*y^3)/(3*sigma^3)))+(xi+1)*(sum(-(y^2/(2*sigma^2))+((2*xi*y^3)/(3*sigma^3))))
	sig_derive1 <- (length(y)/sigma)+(xi+1)*sum((-y/sigma^2)+((xi*y^2)/sigma^3)-((xi^2*y^3)/sigma^4))

##log derivatives##

	xi_derive2 <- (1/xi^2)*(-sum(logb(1+(y*xi)/sigma))+sum(y/(sigma+y*xi))*xi^2+sum(y/(sigma+y*xi))*xi)
	sig_derive2 <- (length(y)/sigma) - sum((xi*y)/((sigma^2)*(1+(y*xi)/sigma))) * (1/xi + 1)

	if((xi < 10^(-10)) & (xi > -10^(-10)))
		{
		gradient<-c(sig_derive1,xi_derive1)
		}
	else
		{
		gradient<-c(sig_derive2,xi_derive2)
		}
#		}
	#cat("gradient=",gradient,"\n")
	return(gradient)

	}

##optimization##
	###two constaints###
	 
	x <- constrOptim2(initials, f = gpd.lik, grad = grad.lik, ui = rbind(c(1,0),c(1,max(y))),ci=c(0,0),y=y)

	###three constraints###
	#x <- constrOptim(initials, f = gpd.lik, grad = grad.lik, ui = matrix(c(1,0,1,max(y),0,1),3,2),ci=c(0,0,-0.5),y=y)



	if(any(1+(y*x$par[2])/x$par[1])<= 0 || x$par[1] <= 0 )
		{
		cat("query this sigma!",x$par[1],"\n")
		cat("query this xi!",x$par[2],"\n")
		}


	z$conv <- x$convergence
    	z$nllh <- x$value
	z$mle <-x$par
    	#z$rate <- length(data)/n

##Hessian calculations, varcov matrix & standard errors##
#if(z$mle[2] <= -0.5)
#{
	if(z$mle[2] < 10^(-10) & z$mle[2] > -10^(-10))
		{z$analytic_covar<- myhess2(y,c(z$mle[1],z$mle[2]))
		if(any(diag(z$analytic_covar) <= 0))
		{z$analytic_covar<-matrix(NA,2,2)}
 		 #z$analytic_se<-sqrt(diag(z$analytic_covar))
		}
	else{
		z$analytic_covar<- myhess(y,c(z$mle[1],z$mle[2]))
		if(any(diag(z$analytic_covar) <= 0))
		{z$analytic_covar<-matrix(NA,2,2)}
		#z$analytic_se<-sqrt(diag(z$analytic_covar))
		}
# }
##likelihood plot##

#	xi.seq <- seq(from = z$mle[2] - 3*z$analytic_se[2], to = z$mle[2] + 3*z$analytic_se[2],length = 1000)
#	sigma.seq <- seq(from = z$mle[1] - 3*z$analytic_se[1], to = z$mle[1] + 3*z$analytic_se[1],length = 1000)
#
#	#xi.seq <- seq(from = z$mle[2] - 0.1, to = z$mle[2] + 0.1,length = 1000)
#	#sigma.seq <- seq(from = 0 - 0.001, to = 0 + 0.001,length = 1000)
#
#	like.f1 <- matrix(NA, nrow = 1000, ncol = 1)
#	like.f2 <- matrix(NA, nrow = 1000, ncol = 1)
#
#	for(i in 1:1000)
#	{
#	like.f1[i] <- gpd.lik(cbind(sigma.seq[i],rep(z$mle[2],length.out=length(sigma.seq))),y=y)
#	like.f2[i] <- gpd.lik(cbind(z$mle[1],xi.seq[i]),y=y)
#	}
#	par(mfrow=c(2,2))
#	plot(sigma.seq,like.f1,type="l")
#	plot(xi.seq,like.f2,type="l")


##contour plot ##

#	xi.seq <- seq(from = z$mle[2] - 3*z$analytic_se[2], to = z$mle[2] + 3*z$analytic_se[2],length = 100)
#	sigma.seq <- seq(from = z$mle[1] - 3*z$analytic_se[1], to = z$mle[1] + 3*z$analytic_se[1],length = 100)
#
#	#xi.seq <- seq(from = z$mle[2] - 10, to = z$mle[2] + 100,length = 100)
##	#sigma.seq <- seq(from = z$mle[1] - 10, to = z$mle[1] + 100,length = 100)

#	like.f <- matrix(NA, nrow = 100, ncol = 100)

#	for(i in 1:100){
#	for(j in 1:100){
#	like.f[i,j] <- gpd.lik(cbind(sigma.seq[i],xi.seq[j]),y=y)
#	}
#	}
#	image(sigma.seq,xi.seq, like.f,xlab="sigma",ylab="xi",main="Negative log-likelihood")
#	contour(sigma.seq,xi.seq, like.f,add=TRUE)
#abline(h=-0.5)
	#persp(sigma.seq,xi.seq, like.f,theta=60,phi=30)
	#return(like.f)

################

#
# Plot the gradient
#

#sigma.xi.seq<-expand.grid(sigma=sigma.seq,xi=xi.seq)
#df.sigma.xi <- apply(sigma.xi.seq, 1, grad.lik, y = y)

#df.sigma.xi.mat <- matrix(df.sigma.xi[1,], nrow =  100)

#image(sigma.seq, xi.seq, df.sigma.xi.mat, xlab = "sigma", ylab = "xi", main = "Gradient in sigma direction")
#contour(sigma.seq, xi.seq, df.sigma.xi.mat, add = TRUE)
#points(z$mle[1], z$mle[2], cex = 2)
#df.sigma.xi.mat <- matrix(df.sigma.xi[2,], nrow =  100)

#image(sigma.seq, xi.seq, df.sigma.xi.mat, xlab = "sigma", ylab = "xi", main = "Gradient in xi direction")
#contour(sigma.seq, xi.seq, df.sigma.xi.mat, add = TRUE)
#points(z$mle[1], z$mle[2], cex = 2)
#par(mfrow=c(1,1))

#####################
	z$gradient <- grad.lik(z$mle,y)
	z$n <- n
  z$rate <- length(exceed)/n
    	z$xdata <- data
  z$npy<-365

#    if (show)
#	{
#       print(z[c(1,2,4,5,6,7,8)])
#    	}
   	invisible(z)

}
