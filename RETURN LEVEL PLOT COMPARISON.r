 pdf("returnlevelcompare2az.pdf",width=8,height=8)

 ##RUN THESE TWO PROGRAMS IN SUCCESSION TO PRODUCE CONFIDENCE COMPARISON RETURN LEVEL PLOT###


gpd.rl<- function (a, u, la, n, npy, mat, dat, xdat)
{
    a <- c(la, a)
    eps <- 1e-06
    a1 <- a
    a2 <- a
    a3 <- a
    a1[1] <- a[1] + eps
    a2[2] <- a[2] + eps
    a3[3] <- a[3] + eps
    jj <- seq(-1, 3.75 + log10(npy), by = 0.1)
    m <- c(1/la, 10^jj)
    q <- gpdq2(a[2:3], u, la, m)
    d1 <- (gpdq2(a1[2:3], u, la, m) - q)/eps
    d2 <- (gpdq2(a2[2:3], u, la, m) - q)/eps
    d3 <- (gpdq2(a3[2:3], u, la, m) - q)/eps
    d <- cbind(d1, d2, d3)
    mat <- matrix(c((la * (1 - la))/n, 0, 0, 0, mat[1, 1], mat[1,
        2], 0, mat[2, 1], mat[2, 2]), nc = 3)
    v <- apply(d, 1, q.form, m = mat)
    plot(m/npy, q, log = "x", type = "n", xlim = c(0.1, max(m)/npy),
        ylim = c(u, max(xdat, 0.5+q[q > u - 1] + 1.96 * sqrt(v)[q >
            u - 1])), xlab = "Return period (years)", ylab = "Return level",
        main = " ")
    lines(m[q > u - 1]/npy, q[q > u - 1],lwd=2)
    lines(m[q > u - 1]/npy, q[q > u - 1] + 1.96 * sqrt(v)[q >
        u - 1], col = 4,lty=2,lwd=2)
    lines(m[q > u - 1]/npy, q[q > u - 1] - 1.96 * sqrt(v)[q >
        u - 1], col = 4,lty=2,lwd=2)
    nl <- n - length(dat) + 1
    sdat <- sort(xdat)
    points((1/(1 - (1:n)/(n + 1))/npy)[sdat > u], sdat[sdat > u])
}





gpd.rl2<-function (a, u, la, n, npy, mat, dat, xdat)
{
    a <- c(la, a)
    eps <- 1e-06
    a1 <- a
    a2 <- a
    a3 <- a
    a1[1] <- a[1] + eps
    a2[2] <- a[2] + eps
    a3[3] <- a[3] + eps
    jj <- seq(-1, 3.75 + log10(npy), by = 0.1)
    m <- c(1/la, 10^jj)
    q <- gpdq2(a[2:3], u, la, m)
    d1 <- (gpdq2(a1[2:3], u, la, m) - q)/eps
    d2 <- (gpdq2(a2[2:3], u, la, m) - q)/eps
    d3 <- (gpdq2(a3[2:3], u, la, m) - q)/eps
    d <- cbind(d1, d2, d3)
    mat <- matrix(c((la * (1 - la))/n, 0, 0, 0, mat[1, 1], mat[1,
        2], 0, mat[2, 1], mat[2, 2]), nc = 3)
    v <- apply(d, 1, q.form, m = mat)

    lines(m[q > u - 1]/npy, q[q > u - 1],lty=3,lwd=2)
    lines(m[q > u - 1]/npy, q[q > u - 1] + 1.96 * sqrt(v)[q >
        u - 1], col = 4,lty=4,lwd=2)
    lines(m[q > u - 1]/npy, q[q > u - 1] - 1.96 * sqrt(v)[q >
        u - 1], col = 4,lty=4,,lwd=2)
    nl <- n - length(dat) + 1
    sdat <- sort(xdat)
    #points((1/(1 - (1:n)/(n + 1))/npy)[sdat > u], sdat[sdat > u],pch=4)
}


fit1<-gpd.fit(new.dat$WHT,1.746)
fit2<-gpd.fit(new.dat$WHT, 1.970)

gpd.rl(fit1$mle, fit1$threshold, fit1$rate, fit1$n, fit1$npy, fit1$cov,  fit1$data, fit1$xdata)
gpd.rl2(fit2$mle, fit2$threshold, fit2$rate, fit2$n, fit2$npy, fit2$cov,  fit2$data, fit2$xdata)
rtn.leg<-c("Return Level curve (JOINSEA)","Return Level curve (Automated threshold)","Confidence Interval (JOINSEA Threshold)","Confidence Interval (Automated Threshold)")
legend("topleft",rtn.leg,cex=1.2,lty=c(3,1,4,2),lwd=2,col=c("black","black","blue","blue"),bg="white")

  dev.off()