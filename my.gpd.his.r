
# pdf("histogramcompare2az.pdf",width=8,height=8)

#library(evir)
#library(ismev)

#my.samp<-read.table("d://Local Data//p1thompson//University work//R work/Rana QR//my.samp.txt",header=T)

my.gpd.his<-function (a1, a2, u1, u2, dat1,dat2)
{
    h1 <- hist(dat1, prob = TRUE, plot = FALSE)
    x1 <- seq(u1, max(h1$breaks), length = 100)
    y1 <- gpd.dens(a1, u1, x1)
    hist(dat1, prob = TRUE, ylim = c(0, max(y1)), xlab = "Daily Rainfall (mm)", ylab = "Density",
        main = " ")
    lines(x1, y1, col = 4,lwd=2)



    h2 <- hist(dat2, prob = TRUE, plot = FALSE)
    x2 <- seq(u2, max(h2$breaks), length = 100)

    one<-pgpd(u2,a2[2],u2,a2[1])
    two<-pgpd(u1,a2[2],u2,a2[1])
#print(one)
print(two)
#print(1-two-one)
	y2 <- gpd.dens(a2, u2, x2)
	y2<-y2/(1-two-one)
print(y2)
    lines(x2, y2, col = 4,lty=2,lwd=2)

 cap.txt<-c("Coles(2001) Threshold","Automated Threshold")
   legend("topright",cap.txt,lty=c(1,2),lwd=2,cex=1.1,bg="white")

}

 fit2<-gpd.fit(rain,20.31)
 fit1<-gpd.fit(rain, 30)


 my.gpd.his(fit1$mle,fit2$mle,fit1$threshold,fit2$threshold,fit1$data,fit2$data)
 
#  dev.off()