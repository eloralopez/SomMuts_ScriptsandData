data<-read.delim("~/Documents/CurrentBiologySubmission/sizeandfreq.txt")

size<-data$sizecm2 
GoH<-data$GoH
LoH<-data$LoH
totalverified<-data$totalverified
totalbases<-data$totalbasesassayed
Nsamples<-data$Nsamples
freq<-totalverified/totalbases/Nsamples
lohfreq<-LoH/totalbases/Nsamples
denovofreq<-GoH/totalbases/Nsamples

lmdenovo<-lm(denovofreq~size)
lmtotal<-lm(freq~size)
lmLoH<-lm(lohfreq~size)

modelSummary <- summary(lmdenovo)  # capture model summary as an object
modelCoeffs <- modelSummary$coefficients  # model coefficients
beta.estimate <- modelCoeffs["size", "Estimate"]  # get beta estimate for speed
std.error <- modelCoeffs["size", "Std. Error"]  # get std.error for speed
t_value <- beta.estimate/std.error  # calc t statistic
p_value <- 2*pt(-abs(t_value), df=nrow(cars)-ncol(cars))  # calc p Value
f_statistic <- linearMod$fstatistic[1]  # fstatistic
f <- summary(linearMod)$fstatistic  # parameters for model p-value calc
model_p <- pf(f[1], f[2], f[3], lower=FALSE)

par(mfrow=c(1,1))
plot(size,denovofreq, col="blue",ylim=c(min(denovofreq),max(freq)), pch=16, cex=2, ylab="Mutations per nucleotide per sample
", xlab="Colony surface area (cm2)")
abline(lmdenovo, col="blue",lty="dashed")
points(size,freq, pch=16, cex=2)
abline(lmtotal,lty="dashed")
points(size,lohfreq,col="red",pch=16, cex=2)
abline(lmLoH, col="red",lty="dashed")

