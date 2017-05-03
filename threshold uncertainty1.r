
        ##########################################################
        #          ~~UNCERTAINTY IN THRESHOLD CHOICE             #
        #    BOOTSTRAPPING PROCEDURE FOR CONFIDENCE INTERVALS    #
        #                                                        #
        ##########################################################
        
  ###LOAD DATA###

new.dat<-read.table("d://Local Data//p1thompson//University work//R work/Rana QR//my.samp.txt",header=T)
 #new.dat<-my.dist
  ###LOAD REQUIRED LIBRARY###
  library(evir)
  library(bootstrap)
     ######################
  ###FUNCTION FOR BOOTSTRAP###
     ######################

   ###THRESHOLD SELECTION FUNCTION###

 ###THRESHOLD SELECTION FUNCTION###

     ppfitrange_idea5<-function (data, umin=0, umax, nint = 100,show=F)
{
  library(nortest)
lowfindthresh<-function (data, ne)
{
    data <- (sort(as.numeric(data)))
    thresholds <- unique(data)
    indices <- match(data[ne], thresholds)
    indices <- pmin(indices + 1, length(thresholds))
    thresholds[indices]
}

umax<-findthresh(data,50)
#print(umax)
umin<-lowfindthresh(data,50)
#print(umin)

    	m <- s <- paulvar <- up <- ul <- matrix(NA, nrow = nint, ncol = 2)
    	my.m <- m

  u <- seq(umin, umax, length = nint)
  #print(u)
  for (i in 1:nint) {
        z <- paulgpd5(data, u[i])
	#MAX LIKELIHOOD EST FOR EACH MODEL.
        m[i, ] <- z$mle
	  paulvar[i, ] <- diag(z$analytic_covar)

	#REPARAMETERIZATION OF GPD SCALE PARAMETER.
     #   m[i, 1] <- m[i, 1] - m[i, 2] * u[i]
        d <- matrix(c(1, -u[i]), ncol = 1)
      #VARIANCE AND STD ERRORS FOR EACH MODEL.
	 # v <- t(d) %*% z$analytic_covar %*% d
#        s[i, ] <- z$se2
#	  s[i, 1] <- sqrt(v)
	#95% CONF INTERVAL FOR MLE.
   #     up[i, ] <- m[i, ] + 1.96 * s[i, ]
   #     ul[i, ] <- m[i, ] - 1.96 * s[i, ]


diff<-matrix(NA,nrow = nint-1, ncol = 2)
diff2<-matrix(NA,nrow = nint-1, ncol = 2)

test<-vector(mode="numeric", length = (nint-1))


   names <- c("Differenced Modified Scale", "Differenced Shape")
   # oldpar <- par(mfrow = c(2, 1))

        		um <- max(up[, 2])
     			ud <- min(ul[, 2])
			my.m[,2] <- m[,2]
			#par(oldpar)
			invisible()
			}

for(g in 1:(length(m[,1])-1))
		{
		diff[g,]<-m[(g+1),]-m[g,]
		diff2[g,]<-paulvar[(g+1),]-paulvar[g,]
		}
p0old<-diff[,1]
uold<-u[2:length(u)]
#m<-na.omit(m)
#paulvar<-na.omit(paulvar)
diff.na<-is.na(diff[,1])
diff.na2<-is.na(diff2[,1])

## diff ##
#p0 <- diff[!diff.na,2]
#if(any(is.na(p0)))
#{print(p0)}
#u <- u[!diff.na]

#if(is.na(m[1,]) || is.na(paulvar[1,]))
#{na.omit(m[1,])
#na.omit(paulvar[1,])}

#m <- m[-1,]
#paulvar <- paulvar[-1,]

#m <- m[!diff.na,2]
#paulvar <- paulvar[!diff.na,2]

####diff2####

p0 <- diff[!diff.na2,1]
#if(any(is.na(p0)))
#{print(p0)}
u <- u[!diff.na2]

#if(is.na(m[1,]) || is.na(paulvar[1,]))
#{na.omit(m[1,])
#na.omit(paulvar[1,])}

m <- m[-1,]
paulvar <- paulvar[-1,]

m <- m[!diff.na2,1]
paulvar <- paulvar[!diff.na2,1]

####

#cat("m dimension",dim(m),"\n")
#cat("paulvar dimension",dim(paulvar),"\n")
#print(paulvar[,2])

pointset<- length(p0)
#print(pointset)
 quantilemin<-quantile(data, prob=0.98)
 newquant<-paulgpd5(data,quantilemin)
  #print(newquant$num_exceed)
 testex <- 0.25*newquant$num_exceed

loopy1<-function(pointset=pointset,u=u, m=m,paulvar=paulvar,umin=umin,umax=umax)
{
a<-list()
#while(testex <= 100)
# {

    testcase <- function(pointset=pointset,u=u, m=m,paulvar=paulvar,umin=umin,umax=umax)
        {
          z<-list()
          test<-vector(mode="numeric", length = (pointset))
          my.test1<-vector(mode="numeric", length = (pointset))
          my.test2<-vector(mode="numeric", length = (pointset))
          my.test3<-vector(mode="numeric", length = (pointset))
          my.test4<-vector(mode="numeric", length = (pointset))
          my.test5<-vector(mode="numeric", length = (pointset))
          testend<-length(test)-7

        for(ppit in 1:testend)
				{ #cat("ppit, iteration=",ppit,"\n")
          z$thres<-NA
        # print(paste("mean part=",m[ppit]))
        #  print(paste("variance part=",paulvar[ppit]*((u[ppit+1]-u[ppit])^2)))


         # my.test1[ppit] <- pearson.test(p0[ppit:length(test)])$p.value
          my.test2[ppit] <- sf.test(p0[ppit:length(test)])$p.value
          my.test3[ppit] <- ad.test(p0[ppit:length(test)])$p.value
          my.test4[ppit] <- lillie.test(p0[ppit:length(test)])$p.value
          #my.test5[ppit] <- cvm.test(p0[ppit:length(test)])$p.value

          #test[ppit]<-mypearson2(p0[1:ppit],meanpear = (m[ppit]*(umax-umin))/length(m), sdpear = paulvar[ppit]/(length(m)^2))
          test[ppit]<-mypearson2(p0[ppit:length(test)],meanpear = m[ppit], sdpear = paulvar[ppit]*((u[ppit+1]-u[ppit])^2))
          #print(paste("my test value1=",my.test1[ppit]))
         ##print(paste("my test value2=",my.test2[ppit]))
          ##print(paste("my test value3=",my.test3[ppit]))
          print(paste("my test value4=",my.test4[ppit]))
          #print(paste("my test value5=",my.test5[ppit]))
          cat("test value=",test[ppit],"\n")
          ##cat("threshold",u[ppit],"\n")
          cat(" ","\n")

         if(my.test4[ppit] > 0.05)
					     {
                #if(my.test4[ppit-1] > 0.05)
					      #  {
                 #   if(my.test4[ppit-2] > 0.05)
	 				        #   {
	                     z$thres<-u[ppit]
                       z$b<-ppit
                   #  }
                  #}
					     }
					if(is.numeric(z$thres))
				  {
          #hist( p0[ppit:length(test)],10)
          #meanmine<-mean(p0[ppit:length(test)])
          #print(m[ppit])
         #abline(v=ppit,lwd=2,col="blue",lty=2)
          break
          }

        }
       # hist(my.test4,20)
        #plot(seq(c(1:length(my.test4))),my.test4,type="l",lty=1,xlab="normality test number",ylab="P values for normality tests")
        #lines(seq(c(1:length(my.test4))),my.test3,lty=2)
        #lines(seq(c(1:length(my.test4))),my.test2,lty=3)
        #abline(v=89,lwd=2,col="blue",lty=3)
         
  #      leg3.txt<-c("Shapiro-Francia","Lilliefors","Anderson-Darling")
  # legend(x=10,y=0.1,leg3.txt,cex=0.7,lty=c(3,2,1),lwd=1.5,bg="white")


         #print(z$thres)
	      z$thres<-z$thres
				z$b<-z$b
				invisible(z)
       }

       mytest<-testcase(pointset,u,m,paulvar,umin,umax)
    #   cat("b value=",mytest$b,"\n")
   #    cat("test thres=",mytest$thres,"\n")
      newtest<-paulgpd5(data,mytest$thres)
      testex<-newtest$num_exceed
  #     cat("test exceed=",testex,"\n")
      if(testex <= max(newquant$num_exceed))
        {pointset<-mytest$b}
 #     cat("pointset=",pointset,"\n")
    #rnd1<-round(mytest$thres,digits=3)
    a$choice<-mytest$thres
# }

 a$choice<-a$choice
 invisible(a)
}
selecta<-loopy1(pointset,u,m,paulvar,umin,umax)
    ulgth<-length(u)
    u<-u[2:ulgth]
plot(uold, p0old,xlab = expression(paste("Threshold ", u[j-1])), ylab = expression(hat(tau)[u[j]] - hat(tau)[u[j-1]]), type = "b",ylim=c(-0.2,0.2),cex.lab=1.3, cex.axis=1.3)
#abline(v=1.672323,col="red",lwd=2)
	abline(v=selecta$choice,col=10,lwd=2)
#  abline(v=umin,lty=3)
#  abline(v=umax,lty=3)
 # print(testcase$thres)
#cat("threshold choice=",selecta$choice)
#z<-list()
threshold<-selecta$choice
 #print(z$threshold)
#if (show) {
#           print(z[1])
#          }
invisible(z)
 #return(list(threshold=threshold))
 return(threshold)
 } 
 x<-my.dist2
  bboot<-vector(mode="numeric",length=100)
 for(i in 1:100)
 {
  bboot[i]<-ppfitrange_idea5(sample(x,10000,replace=T)) 
 # boot[i]<-boot1$threshold
 }
 print(bboot)
  sortboot3<-sort(bboot)
 
 
 
 
    #####################
 ### BOOTSTRAP PROCEDURE ###
    #####################
 #    x<-jitter(my.samp$Hs)
 #   bootfun3 <-function(x){ppfitrange_idea5(x)}
 #   
 #  my.boot3 <- bootstrap(x,nboot=10,theta=bootfun3)
 #  sortboot3<-sort(my.boot3$thetastar)


   ###EXTRACT CONFIDENCE INTERVALS AND MEDIAN###
   
   boot.thres <- median(sortboot3)
   
   boot.thres.up <- quantile(sortboot3,0.975)
   
   boot.thres.low <- quantile(sortboot3,0.025)
   
   hist(boot,main="Bootstrap Interval for Threshold Uncertainty",xlab="Threshold (m)")
   
   #abline(v=boot.thres)
   abline(v=boot.thres.up,col="blue", lty=2)
   abline(v=boot.thres.low,col="blue", lty=2)
   #abline(v=1.607,lty=3)