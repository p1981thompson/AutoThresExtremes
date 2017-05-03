ppfitrange_idea5<-function (data, umin, umax, nint = 100,show=F)
{


#lowfindthresh<-function (data, ne)
#    {
#        data <- (sort(as.numeric(data)))
#        thresholds <- unique(data)
#        indices <- match(data[ne], thresholds)
#        indices <- pmin(indices + 1, length(thresholds))
#        thresholds[indices]
#    }

umax<-findthresh(data,50)

umin<-quantile(data,0.25)

m <- s <- paulvar <- up <- ul <- mnew <- matrix(NA, nrow = nint, ncol = 2)
paulcovar <- rep(NA, length = nint)
my.m <- m
u <- seq(umin, umax, length = nint)


  for (i in 1:nint)
  {
    #FITTING GPD AT RANGE OF THRESHOLDS U#
      z <- paulgpd5(data, u[i])
  #MAX LIKELIHOOD ESTIMATES & VARIANCES FOR EACH MODEL.
      m[i, ] <- z$mle
      mnew[i,] <- z$mle
      paulvar[i, ] <- diag(z$analytic_covar)
      paulcovar[i] <- z$analytic_covar[1,2]

    #REPARAMETERIZATION OF GPD SCALE PARAMETER.
        m[i, 1] <- m[i, 1] - m[i, 2] * u[i]
        d <- matrix(c(1, -u[i]), ncol = 1)

        names <- c("Differenced Modified Scale", "Differenced Shape")
        um <- max(up[, 2])
        ud <- min(ul[, 2])
        my.m[,2] <- m[,2]
        invisible()

}
lgth1 <- length(m[,1])-1
par(mfrow=c(1,1))
#par(mfrow=c(3,2))
# plot(u[1:lgth1], diff(mnew[,1]),xlab = expression(paste("Threshold ",  u[j-1])), ylab = expression(hat(sigma)[u[j]] - hat(sigma)[u[j-1]]))
 #abline(h=0.2*(umax-umin)/nint,col="blue")
 #abline(v=2.9,col="red")
# plot(u[1:lgth1], diff(m[,1]), xlab = expression(paste("Threshold ", u[j-1])), ylab = expression(hat(tau)[u[j]] - hat(tau)[u[j-1]]))
# abline(h=0)
 #abline(v=2.9,col="red")
# plot(u, mnew[,1],xlab = expression(paste("Threshold ", u[j-1])), ylab = expression(hat(sigma)[u[j-1]]))
 #abline(h=0.398634147485306+0.2*(2.678225-2.9),col="green")
 #abline(v=2.9,col="red")
# plot(u, ((m[,2])/sqrt(paulvar[,2])),  xlab = expression(paste("Threshold ", u[j-1])), ylab = expression((hat(xi)[u[j-1]]) / se(hat(xi)[u[j-1]])))
#abline(h=0,lty=2)
#abline(h=c(2,-2),lty=3)
#abline(v=2.9,col="red")
##plot(u, m[,2] ,xlab = expression(paste("Threshold ", u[j-1])), ylab = expression(hat(xi)[u[j-1]]))
#abline(v=2.9,col="red")
j.m <- m
j.paulvar <- paulvar

#return(j.m, j.paulvar, u)

#stop()


#SET UP BLANK GRIDS TO FILL#

diff<-matrix(NA,nrow = nint-1, ncol = 2)
diff2<-matrix(NA,nrow = nint-1, ncol = 2)
test<-vector(mode="numeric", length = (nint-1))


  for(g in 1:(length(m[,1])-1))
{
  diff[g,]<-m[(g+1),]-m[g,]
  diff2[g,]<-paulvar[(g+1),]-paulvar[g,]
}



p0old<-diff[,1]
uold<-u[2:length(u)]

diff.na<-is.na(diff[,1])
diff.na2<-is.na(diff2[,1])

p0 <- diff[!diff.na2,1]
u <- u[!diff.na2]
m <- m[-1,]

paulvar <- paulvar[-1,]
mmat <- m[!diff.na2,]
m <- m[!diff.na2,1]
paulvar <- paulvar[!diff.na2,1]

pointset<- length(p0)

   # quantilemin<-quantile(data, prob=0.98)
    newquant<-paulgpd5(data,umax)
    testex <- 0.25*newquant$num_exceed


 testcase <- function(pointset=pointset,u=u, m=m,paulvar=paulvar,umin=umin,umax=umax)
        { z<-list()
          test<-vector(mode="numeric", length = (pointset))
          my.test1<-vector(mode="numeric", length = (pointset))
          my.test2<-vector(mode="numeric", length = (pointset))
          my.test3<-vector(mode="numeric", length = (pointset))
          my.test4<-vector(mode="numeric", length = (pointset))
          my.test5<-vector(mode="numeric", length = (pointset))
          testend<-length(test)-7
                   for(ppit1 in 1:testend)
      {  z$thres<-NA

              # my.test1[ppit] <- pearson.test(p0[ppit:length(test)])$p.value
              # my.test2[ppit] <- sf.test(p0[ppit:length(test)])$p.value
              # my.test3[ppit] <- ad.test(p0[ppit:length(test)])$p.value
              #my.test4[ppit1] <- lillie.test(p0[ppit1:length(test)])$p.value

              #my.test5[ppit] <- cvm.test(p0[ppit:length(test)])$p.value

              #test[ppit]<-mypearson2(p0[ppit:length(test)],meanpear = m[ppit], sdpear = paulvar[ppit]*((u[ppit+1]-u[ppit])^2))
              test.mean <- 0
              test.var1 <- j.paulvar[ppit1,1] +   j.paulvar[ppit1,2]*u[ppit1]^2 - 2*u[ppit1]*paulcovar[ppit1]
              test.var2 <- j.paulvar[ppit1 + 1,1] +   j.paulvar[ppit1 + 1,2]*u[ppit1 + 1]^2 - 2*u[ppit1 + 1]*paulcovar[ppit1 + 1]
              test.var <- test.var1 + test.var2
              test.sd <- sqrt(test.var)
            #  print(test.sd)

              my.test4[ppit1]<-mypearson2(p0[ppit1:length(test)], meanpear = 0, sdpear = sd(p0[ppit1:length(test)]))
           #   print(my.test4[ppit1])


              #hist(p0[ppit1:length(test)])
              #rug(p0[ppit1:length(test)])
              #print( p0[ppit1:length(test)])
              #print(mean(p0[ppit1:length(test)]))
              #abline(v = test.mean)
              #abline(v = test.mean + 2*test.sd)
              #abline(v = test.mean - 2*test.sd)
              #scan("")
              #print(paste("my test value1=",my.test1[ppit]))
              # print(paste("my test value2=",my.test2[ppit]))
              # print(paste("my test value3=",my.test3[ppit]))
             # print(paste("my test value4=",my.test4[ppit1]))
              #print(paste("my test value5=",my.test5[ppit]))
              #cat("test value=",test[ppit],"\n")
              #cat("threshold",u[ppit],"\n")
             # cat(" ","\n")
              }
              #print(my.test4)
              for(ppit in 2:testend)
               {
               #if(my.test4[ppit] > 0.2 && my.test4[ppit-1] > 0.2 && my.test4[ppit-2] > 0.2)
               if(my.test4[ppit] > 0.05 )

# && my.test4[ppit-1] > 0.1
           {
                      #if(my.test4[ppit-1] > 0.2)
             # { print("testcase 5")
                        #  if(my.test4[ppit-2] > 0.2)
               #  {print("testcase 6")
                          z$thres<-u[ppit]
                            z$b<-ppit

                        #   }
                       # }
           }
     if(is.numeric(z$thres))
          {
                   break
                  }
               }
      z$thres<-z$thres
z$b<-z$b
invisible(z)
       }



loopy1<-function(pointset=pointset,u=u, m=m,paulvar=paulvar,umin=umin,umax=umax)
  {    a<-list()
       mytest<-testcase(pointset,u,m,paulvar,umin,umax)
        newtest<-paulgpd5(data,mytest$thres)

       testex<-newtest$num_exceed

      if(testex <= max(newquant$num_exceed))
        {pointset<-mytest$b}

      a$choice<-mytest$thres
      invisible(a)
  }



selecta<-loopy1(pointset,u,m,paulvar,umin,umax)
    ulgth<-length(u)
    u<-u[2:ulgth]

   threshold<-selecta$choice
   invisible(z)

#plot(uold, p0old,xlab = expression(paste("Threshold ", u[j-1])), ylab = expression(hat(tau)[u[j]] - hat(tau)[u[j-1]]), type = "b",ylim=c(-0.3,0.3),cex.lab=1.3, cex.axis=1.3)

#abline(v=selecta$choice,col=10,lwd=2)
#abline(h=0,lwd=2)
return(threshold)
}
