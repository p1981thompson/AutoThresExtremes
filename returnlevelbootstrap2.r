

rlgen<-function(){
no.of.years.of.data <- 27 # Or whatever

x <- sample(my.samp$Hs,10000,replace=T)

n <- length(x) # Same as length of original data

chosen.u <- ppfitrange_idea5(x)

gpdpaul<-gpd.fit(x,chosen.u,show=F)

gpd.mle <- gpdpaul$mle

no.of.points.above.threshold <- length(x[x > chosen.u])

my.rate <- no.of.points.above.threshold / n

number.of.observations.per.year <- n / no.of.years.of.data


return.period.in.years <- 1000 # Or whatever


#gpdq(gpd.mle , chosen.u , 1 / (my.rate * number.of.observations.per.year * return.period.in.years)  )


paul<-gpdq(gpd.mle , chosen.u , 1 / (return.period.in.years * no.of.points.above.threshold  / no.of.years.of.data)  )

return(paul)
}

   rlboot<-vector(mode="numeric",length=1000)

 for(i in 1:1000)
    {
      rlboot[i]<- rlgen()
    }

  print(rlboot)

  rlbootsort<-sort(rlboot2)

   rlboot.thres <- median(rlbootsort)

   rlboot.thres.up <- quantile(rlbootsort,0.975)

   rlboot.thres.low <- quantile(rlbootsort,0.025)

   hist(rlboot3,breaks=50,main="Bootstrap interval for 1000 year return level uncertainty",xlab="1000 year Return Level (m)")

   #abline(v=boot.thres)
   abline(v=rlboot.thres.up, lty=2,lwd=2)
   abline(v=rlboot.thres.low, lty=2,lwd=2)


   rlgenact<-function(){
no.of.years.of.data <- 27 # Or whatever

x <- my.samp$Hs

n <- length(x) # Same as length of original data

chosen.u <- 0.4878788

gpd.mle <-  c( 0.5756078, -0.2297250)

no.of.points.above.threshold <- length(x[x > chosen.u])

my.rate <- no.of.points.above.threshold / n

number.of.observations.per.year <- n / no.of.years.of.data


return.period.in.years <- 100 # Or whatever


#gpdq(gpd.mle , chosen.u , 1 / (my.rate * number.of.observations.per.year * return.period.in.years)  )


paul<-gpdq(gpd.mle , chosen.u , 1 / (return.period.in.years * no.of.points.above.threshold  / no.of.years.of.data)  )

return(paul)
}
actual_rtn <- rlgenact()
abline(v=actual_rtn,lwd=2)