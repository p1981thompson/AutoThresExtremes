vary_thres_paul_sampHs<-function (data=my.samp,n=40)
{
#old<-getTime()
pal<-list()
#### BASED ON ISMEV: GPD.FIT ##############


all.res <- numeric(0)
plot(data$coswav,data$Hs,xlab="cos(Wave direction)",ylab="Wave height (m)",cex=0.5,pch=20,col="skyblue", cex.lab=1.3, cex.axis=1.3)


#### CONSTANT THRESHOLD ####
#constant<-ppfitrange_idea2(data$Hs)
#abline(h=0.487878,col="blue",lwd=2,lty=2)
abline(v= -1,lwd=0.8)                                       

#### INITIAL BLOCKING OF DATA ####
borders <- seq(from = -1, to = 1, by = (2/n))
#abline(v=borders,lwd=0.7)
i_new<-1
start_i<-0
#### OPTIMISE BLOCKS ####
for(i in 1:n)
{	i<-i_new
	if(i == n+1)
	{break}
	last_i<-start_i
	start_i<-i
	new.block <- subset(data,coswav > borders[start_i] & coswav < borders[start_i+1])

	#model.thres<-ppfitrange_idea2(new.block[,1])
	#new.thres <- model.thres$threshold


## CONDITION 1: SUFFICIENT NUMBER OF DATA POINTS IN BLOCK##

  	if(length(new.block[,1]) < 500)
  	{
	if(borders[i+1] == borders[n+1])
		{
		new.block <- subset(data,coswav > borders[last_i] & coswav < borders[n+1])
		model.thres<-ppfitrange_idea2(new.block[,5])
		new.thres <- model.thres$threshold
		}
	   else
		{
		length.new.block <- 100
  		k <- 0
  		while(length.new.block < 500)
			{k <- k + 1
				if(borders[i+k] == borders[n+1])
				{
				new.block <- subset(data,coswav > borders[last_i] & coswav < borders[n+1])
				length.new.block <- length(new.block[,5])
				start_i<-last_i
				}
				else
				{
				new.block <- subset(data,coswav > borders[start_i] & coswav < borders[i+k])
				length.new.block <- length(new.block[,5])
				}
			last_k<-k
			}
		new.block<-subset(data,coswav > borders[start_i] & coswav < borders[i+k])
		i<-i+last_k
  		model.thres<-ppfitrange_idea2(new.block[,5])
		new.thres <- model.thres$threshold
		i_new<-i
    		}
	i_new<-i_new
	new_block<-new.block
	model.thres<-model.thres
	new.thres<-new.thres
	}

	model.thres<-ppfitrange_idea2(new.block[,5])
	new.thres <- model.thres$threshold


## CONDITION 2: SUFFICIENT NUMBER OF THRESHOLD EXCESSES FOR BLOCK ##

	test.block<-subset(new.block[,5],new.block[,5] > new.thres)
   quantilemin<-quantile(new.block[,5], prob=0.96)
 newquant<-paulgpd5(new.block[,5],quantilemin)

	if(length(test.block) < max(newquant$num_exceed,100))
 	{
	if(borders[i+1] == borders[n+1])
		{
		new.block <- subset(data,coswav > borders[last_i] & coswav < borders[n+1])
		model.thres<-ppfitrange_idea2(new.block[,5])
		new.thres <- model.thres$threshold
		start_i<-last_i
		}
	   else
		{
		length.test.block <- newquant$num_exceed*0.25
  		j <- 0
  		while(length.test.block < max(newquant$num_exceed,100))
			{j <- j + 1
				if(borders[i+j] == borders[n+1])
				{new.block <- subset(data,coswav > borders[last_i] & coswav < borders[n+1])
				model.thres<-ppfitrange_idea2(new.block[,5])
				new.thres <- model.thres$threshold
				test.block<-subset(new.block[,5],new.block[,5]>new.thres)
				length.test.block <- length(test.block)
				start_i<-last_i
				}
				else
				{
				new.block <- subset(data,coswav > borders[start_i] & coswav < borders[i+j])
				model.thres<-ppfitrange_idea2(new.block[,5])
				new.thres <- model.thres$threshold
				test.block<-subset(new.block[,5],new.block[,5]>new.thres)
				length.test.block <- length(test.block)
				}
			last_j<-j
			}
		new.block<-subset(data,coswav > borders[start_i] & coswav < borders[i+j])
		model.thres<-ppfitrange_idea2(new.block[,5])
		new.thres <- model.thres$threshold
		i<-i+last_j
		i_new<-i
		}
	i_new<-i_new
	new_block<-new.block
	model.thres<-model.thres
	new.thres<-new.thres
	}                                                              
	model.thres<-model.thres
	new.thres<-new.thres
 	if(length(new.block[,5]) >= 500 & length(test.block) >= max(newquant$num_exceed,100))
		{
		i_new<-i
		block.thres <- cbind(new.block,new.thres)
  		all.res <- rbind(all.res,block.thres)
  		}
   #indcomp<-gpd.fit(new.block[,5],new.thres,show=F)
	 #gpd.diag(indcomp)
  i_new<-i+1

	if(i_new > n+1)
		{
		lines(c(borders[start_i],borders[n+1]),c(new.thres,new.thres),lwd=2.5,col="black")
		#print(new.thres)
		#print(borders[n+1])
		abline(v=borders[n+1],lwd=0.5)
		}
	else
		{
		lines(c(borders[start_i],borders[i_new]),c(new.thres,new.thres),lwd=2.5,col="black")
		#print(new.thres)
		#print("border below")
		#print(borders[i_new])
		abline(v=borders[i_new],lwd=0.5)
		}

	#i_new<-i+1
if(borders[i_new-1] == borders[n+1])
{break}
		block.thres <- cbind(new.block,new.thres)
  		all.res <- rbind(all.res,block.thres)
}

############ ESTIMATED CONTINUOUS THRESHOLD FUNCTION FROM SMOOTHING DISCRETE THRESHOLD FUNCTION ############

#### LOESS MODEL ####
#interval <- (1/1000)
#x.seq<-seq(from = -1, to = 1, by = interval)
#fitlo1<-loess(all.res[,5]~all.res[,4],degree=2,span=0.5)
#pres<-predict(fitlo1,data.frame(x.seq))
#lines(x.seq,pres,col="blue",lwd=1.5)

#### SMOOTHING SPLINES ####
#lines(smooth.spline(all.res[,4],all.res[,5],df=5),col="yellow",lwd=1.5)
myfit<-smooth.spline(all.res[,8],all.res[,9],df=25)
#lines(myfit,col="red",lwd=3,lty=1)

#### LEGEND ####
leg.txt<-c("Constant threshold","Smoothed threshold","Piecewise constant threshold")
#legend(x=0,y=2.1,leg.txt,col=c("blue","red","black"),cex=1.2,lty=c(2,1,1),lwd=2,bg="white")

#timeElapsed(old,new=getTime())

}


pdf("newblock2.pdf",width=8,height=8)
vary_thres_paul_sampHs()
#contour(f2,nlevels=30,add=T,col="orangered")
dev.off()
