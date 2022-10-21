library(tikzDevice)

clr1<-rgb(27/255,158/255,119/255)
clr2<-rgb(217/255,95/255,2/255)
clr3<-rgb(117/255,112/255,179/255)
clr4<-rgb(231/255,41/255,138/255)
clr5<-rgb(102/255,166/255,30/255)
clr6<-rgb(230/255,171/255,2/255)
clr7<-rgb(166/255,118/255,29/255)

clr1b<-rgb(27/255,158/255,119/255,0.2)
clr2b<-rgb(217/255,95/255,2/255,0.2)
clr3b<-rgb(117/255,112/255,179/255,0.2)
clr4b<-rgb(231/255,41/255,138/255,0.2)
clr5b<-rgb(102/255,166/255,30/255,0.2)
clr6b<-rgb(230/255,171/255,2/255,0.2)
clr7b<-rgb(166/255,118/255,29/255,0.2)

clr1c<-rgb(27/255,158/255,119/255,0.5)
clr2c<-rgb(217/255,95/255,2/255,0.5)
clr3c<-rgb(117/255,112/255,179/255,0.5)
clr4c<-rgb(231/255,41/255,138/255,0.5)
clr5c<-rgb(102/255,166/255,30/255,0.5)
clr6c<-rgb(230/255,171/255,2/255,0.5)
clr7c<-rgb(166/255,118/255,29/255,0.5)


###############################
#### Initialize parameters ####
###############################

args <-  commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)

n<-args[1]
sigma_0<-args[2]
sigma_1<-args[3]
delta<-args[4]
rho<-args[5]
ID<-args[6]


########################
#### Load libraries ####
########################

library(grf)

#######################
#### Simulate data ####
#######################

seed <- sample(1:2^15, 1)
set.seed(seed)

# Draw X1, X2, U0, U1, X3
x1 = rbinom(n,1,0.5) # Sex.
x2 = rnorm(n, mean=0, sd=1) # Normalized SBP.

u0 = rnorm(n, mean=0, sd=sigma_0) # define U0.
u1 = rnorm(n, mean=0, sd=sigma_1) # define U1.

x3<-sapply(u1,function(x){rnorm(1, mean=0+delta*rho*(x), sd=(delta*sigma_1*sqrt(1-rho^2)))} )

# Simulate Y^0:
y_0_intercept = 5.9
y_0_x1 = 0.8
y_0_x2 = 0.5
y_0 = y_0_intercept + y_0_x1*x1 + y_0_x2*x2 + u0

# Simulate Y^1:
y_1_intercept = 0.45
y_1_x1 = 0.1
y_1_x2 = 0.15
y_1 = y_0 + y_1_intercept + y_1_x1*x1 + y_1_x2*x2 + u1

# Treatment assignment
trt_prob_intercept = -1.7 
trt_prob_x1 = -0.1
trt_prob_x2 = 0.4

trt_prob = exp(trt_prob_intercept + trt_prob_x1*x1 + trt_prob_x2*x2)/
    (1 + exp(trt_prob_intercept + trt_prob_x1*x1 + trt_prob_x2*x2))

A = rbinom(n, 1, trt_prob)

data = data.frame(y_0,y_1,A,x1,x2,x3)
data$outcome = ifelse(data$A == 1, data$y_1, data$y_0)

#################
#### Fit grf ####
#################

X = data.frame(data$x1, data$x2,data$x3)
Y = data$outcome
W = data$A

#Run the causal forest on the data
if(n==200){
    tau.forest10 <- causal_forest(X,Y,W, num.trees = 2000, min.node.size = 1)
}

if(n>200){
    tau.forest10 <- causal_forest(X,Y,W, num.trees = 2000, min.node.size = 5)
}

tau.hat.oob10 <- predict(tau.forest10, estimate.variance=T)

alpha<-as.matrix(get_forest_weights(tau.forest10))

##Reproduce predict.forest()
#predict(tau.forest10)$predictions[1]
#lm(Y.centered ~ W.centered, weights = alpha[1,])$coefficients[[2]]
#solve(t(W.centered)%*%diag(alpha[1,])%*%t(t(W.centered)))%*%t(W.centered)%*%diag(alpha[1,])%*%t(t(Y.centered))
#sum(alpha[1,]*Y.centered*W.centered)/sum(alpha[1,]*W.centered^2)


#Predictions, see eq (19) grf paper, https://github.com/grf-labs/grf/issues/1095
Y.centered <- tau.forest10$Y.orig - tau.forest10$Y.hat
W.centered <- tau.forest10$W.orig - tau.forest10$W.hat

tau<-tau.hat.oob10[,1]
mu1<-tau.forest10$Y.hat+(1-tau.forest10$W.hat)*tau
mu0<-tau.forest10$Y.hat-tau.forest10$W.hat*tau


#########
#Fit regression forest
Y2.centered <- Y^2 - predict(regression_forest(X,Y^2))[,1]

#tau.tilde<-(alpha%*%(Y2.centered*W.centered)/ alpha%*%(W.centered^2))  
tau.tilde<-apply(alpha,1,function(x){lm(Y2.centered ~ W.centered, weights = x)$coefficients[[2]]})

tau2<-tau.tilde-2*tau*mu0
pred.var.neighbouring<-sapply(tau2-tau^2, function(x){max(0,x)})


#We assume var(Y(1))>var(Y(0)), otherwise exposure reduces variability and thus we can speak about heterogeneity of effect of absence of exposure
sd.custom<-sqrt(max(0,mean(tau2)-mean(tau.hat.oob10[,1])^2))
sd.custom.AIPW<-sqrt(max(0,mean(tau2)-(average_treatment_effect(tau.forest10)[1])^2))
P.custom<-mean(apply(cbind(tau.hat.oob10[,1],pred.var.neighbouring),1,function(x){pnorm(0,x[1],sqrt(x[2]))}))



########################
### Bootstrap for CI ###
########################

mean.B<-NULL
mean.dr.B<-NULL
sd.B<-NULL
P.B<-NULL

sd.custom.B<-NULL
sd.custom.AIPW.B<-NULL
P.custom.B<-NULL

for(i in 1:1000){
    boot<-sample(seq(1,n,1),n,replace=T)
    X.B<-X[boot,]
    Y.B<-Y[boot]
    W.B<-W[boot]
    
    if(n==200){
        tau.forest10.B <- causal_forest(X.B,Y.B,W.B, num.trees = 2000, min.node.size = 1)
    }
    
    if(n>200){
        tau.forest10.B <- causal_forest(X.B,Y.B,W.B, num.trees = 2000, min.node.size = 5)
    }
    tau.hat.oob10.B <- predict(tau.forest10.B, estimate.variance=T)
    
    #AIPW ATE estimate
    mean.dr.B<-c(mean.dr.B,average_treatment_effect(tau.forest10.B)[1])
    #Average of CATEs
    mean.B<-c(mean.B,mean(tau.hat.oob10.B[,1]))
    #SD
    sd.B<-c(sd.B,sd(tau.hat.oob10.B[,1]))
    #PEP
    P.B<-c(P.B,length(which(tau.hat.oob10.B[,1]>0))/n)
    
    
    alpha.B<-as.matrix(get_forest_weights(tau.forest10.B))
    
    Y.centered.B <- tau.forest10.B$Y.orig - tau.forest10.B$Y.hat
    W.centered.B <- tau.forest10.B$W.orig - tau.forest10.B$W.hat
    
    tau.B<-tau.hat.oob10.B[,1]
    mu1.B<-tau.forest10.B$Y.hat+(1-tau.forest10.B$W.hat)*tau.B
    mu0.B<-tau.forest10.B$Y.hat-tau.forest10.B$W.hat*tau.B
    
    #Fit regression forest
    Y2.centered.B <- Y.B^2 - predict(regression_forest(X.B,Y.B^2))[,1]
    tau.tilde.B<-apply(alpha.B,1,function(x){lm(Y2.centered.B ~ W.centered.B, weights = x)$coefficients[[2]]})
    
    tau2.B<-tau.tilde.B-2*tau.B*mu0.B
    pred.var.neighbouring.B<-sapply(tau2.B-tau.B^2, function(x){max(0,x)})
    
    sd.custom.B<-c(sd.custom.B,sqrt(max(0,mean(tau2.B)-mean(tau.hat.oob10.B[,1]))^2))
    sd.custom.AIPW.B<-c(sd.custom.AIPW.B,sqrt(max(0,mean(tau2.B)-(average_treatment_effect(tau.forest10.B)[1])^2)))
    P.custom.B<-c(P.custom.B,mean(apply(cbind(tau.hat.oob10.B[,1],pred.var.neighbouring.B),1,function(x){pnorm(0,x[1],sqrt(x[2]))})))
    
}

#####################
#### Save output ####
#####################

if(length(which(abs(mean.dr.B)%in%c(NaN,Inf)))>0){
    mean.dr.B<-mean.dr.B[-which(abs(mean.dr.B)%in%c(NaN,Inf))]
}

if(length(which(abs(mean.B)%in%c(NaN,Inf,NA)))>0){
    mean.B<-mean.B[-which(abs(mean.B)%in%c(NaN,Inf,NA))]
}

if(length(which(abs(sd.B)%in%c(NaN,Inf,NA)))>0){
    sd.B<-sd.B[-which(abs(sd.B)%in%c(NaN,Inf,NA))]
}

if(length(which(abs(P.B)%in%c(NaN,Inf,NA)))>0){
    P.B<-P.B[-which(abs(P.B)%in%c(NaN,Inf,NA))]
}

write.table(
    t(c(ID,seed,
        average_treatment_effect(tau.forest10)[1],average_treatment_effect(tau.forest10)[1]-1.96*average_treatment_effect(tau.forest10)[2],average_treatment_effect(tau.forest10)[1]+1.96*average_treatment_effect(tau.forest10)[2],
        average_treatment_effect(tau.forest10)[1],quantile(mean.dr.B,0.025),quantile(mean.dr.B,0.975),
        mean(tau.hat.oob10[,1]),quantile(mean.B,0.025),quantile(mean.B,0.975),
        sd(tau.hat.oob10[,1]),quantile(sd.B,0.025),quantile(sd.B,0.975),
        length(which(tau.hat.oob10[,1]>0))/n,quantile(P.B,0.025),quantile(P.B,0.975))),  
    file=paste("./GRFoutput/fit-",n,"-",sigma_0,"-",sigma_1,"-",delta,"-",rho,".txt", sep=""), sep=",", append = T,row.names = FALSE, col.names = FALSE, eol="\n")


#Save density at fixed points
c.trunc<-10
dat<-unlist(tau.hat.oob10[,1])
dat<-dat[which(abs(dat)<c.trunc)]
d<-density(dat)
keep.y0<-approx(d$x, d$y, xout = seq(-c.trunc,c.trunc,0.1))$y
keep.y0[which(is.na(keep.y0))]<-0
keep.y<-keep.y0/sum(keep.y0*1)

keep.y<-round(keep.y,digits=4)

write.table(t(c(ID,seed,keep.y)),  file=paste("./GRFoutput/DensityEstimate-",n,"-",sigma_0,"-",sigma_1,"-",delta,"-",rho,".txt", sep=""), sep=",", append = T,row.names = FALSE, col.names = FALSE, eol="\n")

write.table(t(c(ID,seed,coef(test_calibration(tau.forest10)),test_calibration(tau.forest10)[7],test_calibration(tau.forest10)[8])
),  file=paste("./GRFoutput/HeterogeneityTest-",n,"-",sigma_0,"-",sigma_1,"-",delta,"-",rho,".txt", sep=""), sep=",", append = T,row.names = FALSE, col.names = FALSE, eol="\n")



####Write results extended grf###

if(length(which(abs(sd.custom.B)%in%c(NaN,Inf,NA)))>0){
    sd.custom.B<-sd.custom.B[-which(abs(sd.custom.B)%in%c(NaN,Inf,NA))]
}

if(length(which(abs(sd.custom.AIPW.B)%in%c(NaN,Inf,NA)))>0){
    sd.custom.AIPW.B<-sd.custom.AIPW.B[-which(abs(sd.custom.AIPW.B)%in%c(NaN,Inf,NA))]
}

if(length(which(abs(P.custom.B)%in%c(NaN,Inf,NA)))>0){
    P.custom.B<-P.custom.B[-which(abs(P.custom.B)%in%c(NaN,Inf,NA))]
}


write.table(
    t(c(ID,seed,
        average_treatment_effect(tau.forest10)[1],average_treatment_effect(tau.forest10)[1]-1.96*average_treatment_effect(tau.forest10)[2],average_treatment_effect(tau.forest10)[1]+1.96*average_treatment_effect(tau.forest10)[2],
        average_treatment_effect(tau.forest10)[1],quantile(mean.B,0.025),quantile(mean.B,0.975),
        mean(tau.hat.oob10[,1]),quantile(mean.B,0.025),quantile(mean.B,0.975),
        sd.custom,quantile(sd.custom.B,0.025),quantile(sd.custom.B,0.975),
        sd.custom.AIPW,quantile(sd.custom.AIPW.B,0.025),quantile(sd.custom.AIPW.B,0.975),
        P.custom,quantile(P.custom.B,0.025),quantile(P.custom.B,0.975))),  
    file=paste("./GRFoutput/CUSTOMfit-",n,"-",sigma_0,"-",sigma_1,"-",delta,"-",rho,".txt", sep=""), sep=",", append = T,row.names = FALSE, col.names = FALSE, eol="\n")

c.trunc<-10

keep.y0<-sapply(seq(-c.trunc,c.trunc,0.05),function(x){mean(dnorm(x,mean(tau.hat.oob10[,1]), sqrt(pred.var.neighbouring)))})
keep.y0[which(is.na(keep.y0))]<-0

#When variance=0 for some individuals
if(sum(keep.y0*1)<0.1){
    dat<-unlist(apply(cbind(tau.hat.oob10[,1],pred.var.neighbouring),1,function(x){rnorm(10000,x[1],sqrt(x[2]))}))
    dat<-dat[which(abs(dat)<c.trunc)]
    d<-density(dat)
    keep.y0<-approx(d$x, d$y, xout = seq(-c.trunc,c.trunc,0.05))$y
    keep.y0[which(is.na(keep.y0))]<-0
    keep.y<-keep.y0/sum(keep.y0*1)    
} else{
    keep.y<-keep.y0/sum(keep.y0*1)}

keep.y<-round(keep.y,digits=4)


write.table(t(c(ID,seed,keep.y)),  file=paste("./GRFoutput/CUSTOMDensityEstimate-",n,"-",sigma_0,"-",sigma_1,"-",delta,"-",rho,".txt", sep=""), sep=",", append = T,row.names = FALSE, col.names = FALSE, eol="\n")

###################
### Load Output ###
###################
library(readr)




fit <- read_csv("~/HPC/GRFoutput/fit-2000-1.6-1.4-2-1.txt", 
                col_names = FALSE)

fit<-na.omit(fit)



################################
### Obtain bias MSE coverage ###
################################

ACE<-0.45+0.1*0.5+0.15*0+0
#SD.ICE<-sd(y_1-y_0)
SD.ICE<-sqrt(0.1^2*0.25+0.15^2*1+1.4^2)
#P.0<-length(which(y_1-y_0>0))/n
#P.0<-1-0.6386
P.0<-0.6386
P.0.c<-1-P.0


#Table 1

Table<-NULL

for (rho in c(0,0.25,0.5,0.75,1)){
    for (n in c(200,2000,20000)){
        
            fit  <- data.frame(read_csv(paste("~/HPC/GRFoutput/October2022/fit-",n,"-1.6-1.4-2-",rho,".txt",sep=""), 
                                         col_names = FALSE))
        
        fit<-na.omit(fit)

Table<-rbind(Table,round(c(rho,n,mean((fit$X3-ACE)),
        #mean((fit$X5-ACE))
        #mean((fit$X8-ACE))
        
        mean((fit$X12-SD.ICE)),
        
        mean((fit$X15-P.0)),
        
        
        mean((fit$X3-ACE)^2),
        #mean((fit$X5-ACE)^2)
        #mean((fit$X8-ACE)^2)
        
        mean((fit$X12-SD.ICE)^2),
        
        mean((fit$X15-P.0)^2),
        
        length(which(apply(fit,1,function(x){(ACE<=x[5] && ACE>=x[4])})))/length(unlist(fit[,1])),
        #length(which(apply(fit,1,function(x){(ACE<=x[7] && ACE>=x[6])})))/1000
        #length(which(apply(fit,1,function(x){(ACE<=x[10] && ACE>=x[9])})))/1000
        
        length(which(apply(fit,1,function(x){(SD.ICE<=x[14] && SD.ICE>=x[13])})))/length(unlist(fit[,1])),
        
        length(which(apply(fit,1,function(x){(P.0<=x[17] && P.0>=x[16])})))/length(unlist(fit[,1]))),digits=3)
)
    }
}

Table[,2]<-round(Table[,2],digits=0)
colnames(Table)<-c("rho","n","ATE.bias","SD.bias","PEP.bias","ATE.MSE","SD.MSE","PEP.MSE","ATE.cov","SD.cov","PEP.cov")
print(xtable(Table), include.rownames=FALSE)

#Table 2

Table.custom<-NULL

for (rho in c(0,0.25,0.5,0.75,1)){
    for (n in c(200,2000,20000)){
        
      
            cfit  <- data.frame(read_csv(paste("~/HPC/GRFoutput/October2022/CUSTOMfit-",n,"-1.6-1.4-2-",rho,".txt",sep=""), 
                                           col_names = FALSE))

       
        cfit<-na.omit(cfit)
        

Table.custom<-rbind(Table.custom,round(c(rho,n,mean((cfit$X3-ACE),na.rm=T),
        #mean((cfit$X5-ACE))
        #mean((cfit$X8-ACE))
        
        #mean((cfit$X11-SD.ICE))
        mean((cfit$X15-SD.ICE),na.rm=T),
        
        mean(((1-cfit$X18)-P.0),na.rm=T),
        
        mean((cfit$X3-ACE)^2,na.rm=T),
        #mean((cfit$X5-ACE)^2)
        #mean((cfit$X8-ACE)^2)
        
        #mean((cfit$X11-SD.ICE)^2)
        mean((cfit$X15-SD.ICE)^2,na.rm=T),
        
        mean(((1-cfit$X18)-P.0)^2,na.rm=T),
        
        
        #length(which(apply(cfit,1,function(x){(ACE<=unlist(x[4]) & ACE>=unlist(x[3]))})))/length(unlist(cfit[,1])),
        #length(which(apply(cfit,1,function(x){(ACE<=x[7] && ACE>=x[6])})))/1000
        #length(which(apply(cfit,1,function(x){(ACE<=x[10] && ACE>=x[9])})))/1000
        length(which(cfit$X5>=ACE & cfit$X4<=ACE))/length(unlist(cfit[,1])),
        
        #length(which(apply(cfit,1,function(x){(SD.ICE<=x[13] && SD.ICE>=x[12])})))/1000
        #length(which(apply(cfit,1,function(x){(SD.ICE<=x[16] && SD.ICE>=x[15])})))/length(unlist(cfit[,1])),
        length(which(cfit$X17>=SD.ICE & cfit$X16<=SD.ICE))/length(unlist(cfit[,1])),
        
        
        #length(which(apply(cfit,1,function(x){(P.0.c<=x[19] && P.0.c>=x[18])})))/length(unlist(cfit[,1]))
        length(which(cfit$X20>=(P.0.c) & cfit$X19<=(P.0.c)))/length(unlist(cfit[,1]))),digits=3)
   ) }}

Table.custom[,2]<-round(Table.custom[,2],digits=0)
colnames(Table.custom)<-c("rho","n","ATE.bias","SD.bias","PEP.bias","ATE.MSE","SD.MSE","PEP.MSE","ATE.cov","SD.cov","PEP.cov")
print(xtable(Table.custom), include.rownames=FALSE)

###############
### Figures ###
###############

#Fig1

# plot(density(y_1-y_0),ylim=c(0,2.5),xlim=c(-4,5),lwd=2,col=clr1,main="",xlab="Y",ylab="pdf")
# abline(v=0.5,lty=1,lwd=1)
# lines(density(y_1_intercept + y_1_x1*x1 + y_1_x2*x2),lwd=2,col=clr1c,lty=2)
# lines(density(y_1-y_0),ylim=c(0,2.5),xlim=c(-4,5),lwd=2,col=clr1,main="",xlab="Y",ylab="pdf")
# legend(1.5,2,c("Y^1-Y^0", "E[Y^1-Y^0|X]"),lty=c(1,2),col=c(clr1,clr1c),lwd=2)
# 
# 
# plot(density(y_1),ylim=c(0,0.25),xlim=c(0,14),lwd=2,col=clr2,main="",xlab="Y",ylab="pdf")
# lines(density(y_0),col=clr3,lwd=2)
# lines(density(Y[which(A==0)]),col=clr3c,lwd=2,lty=2)
# lines(density(Y[which(A==1)]),col=clr2c,lwd=2,lty=2)
# legend(9,0.225,c("Y^0","Y^1","Y|A=0","Y|A=1"),col=c(clr3,clr2,clr3c,clr2c),lwd=2)


#Fig2
tikz("Fig2b.tex",width=12,height=16)
par(mfrow = c(5,3),
    oma = c(6,6,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)

ymax<-c(2.5,2.0,1.5,1.0,1.0)

for (rho in c(0,0.25,0.5,0.75,1)){
    for (n in c(200,2000,20000)){
        #n<-20000
        if(n<20000){
            density <- data.frame(read_csv(paste("~/HPC/GRFoutput/DensityEstimate-",n,"-1.6-1.4-2-",rho,".txt",sep=""), 
                                           col_names = FALSE))
        }
        if(n==20000){
            density <- data.frame(read_csv(paste("~/HPC/GRFoutput/DensityEstimate-",n,"-1.6-1.4-2-",rho,".txt",sep=""), 
                                           col_names = FALSE)) 
        }
        density[] <- lapply(density, function(x) as.numeric(as.character(x)))

        
#         if(n<20000){
#         rgf.custom.mean<-10*colMeans(data.frame(density)[,-c(1,2)],na.rm=T)
#         
#         plot(seq(-19,19,0.1),rgf.custom.mean,col=clr4,type="l",xlim=c(-4,4),ylim=c(0,ymax[        which(c(0,0.25,0.5,0.75,1)==rho)
# ]),ylab="pdf",xlab="Y1-Y0",axes=F)
#         
#         ti<-seq(-19,19,0.1)
#         lb<-10*apply(density[,-c(1,2)],2,function(x){quantile(x,0.975,na.rm=T)})
#         ub<-10*apply(density[,-c(1,2)],2,function(x){quantile(x,0.025,na.rm=T)})
#         
#         polygon(x=c(rev(ti), ti), y=c(rev(lb), ub), col=clr4c, border=NA)
#         
#         
#         lines(seq(-19,19,0.1),rgf.custom.mean,col=clr4)
#         
#         lines(seq(-19,19,0.1),sapply(seq(-19,19,0.1),function(x){(0.5*dnorm(x,0.55,sqrt(0.15^2+1.4^2))+0.5*dnorm(x,0.45,sqrt(0.15^2+1.4^2)))}),col=clr1)
#         }
#         
#         if(n==20000){
            rgf.custom.mean<-10*colMeans(data.frame(density)[,-c(1,2)],na.rm=T)
            
            plot(seq(-6,6,0.1),rgf.custom.mean,col=clr4,type="l",xlim=c(-4,4),ylim=c(0,ymax[        which(c(0,0.25,0.5,0.75,1)==rho)
                                                                                                      ]),ylab="pdf",xlab="Y1-Y0",axes=F)
            
            ti<-seq(-6,6,0.1)
            lb<-10*apply(density[,-c(1,2)],2,function(x){quantile(x,0.975,na.rm=T)})
            ub<-10*apply(density[,-c(1,2)],2,function(x){quantile(x,0.025,na.rm=T)})
            
            polygon(x=c(rev(ti), ti), y=c(rev(lb), ub), col=clr4b, border=NA)
            
            
            lines(seq(-6,6,0.1),rgf.custom.mean,col=clr4)
            
            lines(seq(-6,6,0.1),sapply(seq(-6,6,0.1),function(x){(0.5*dnorm(x,0.55,sqrt(0.15^2+1.4^2))+0.5*dnorm(x,0.45,sqrt(0.15^2+1.4^2)))}),col=clr1)
     #   }
        
        
        if(n==200){
            axis(side = 2,
                 labels = TRUE, cex.axis=2)
        }
        if(n>200){
            axis(side = 2,
                 labels = FALSE)   
        }
        
        if(rho==1){
            axis(side = 1,
                 labels = TRUE, cex.axis=2)
        }
        if(rho<1){
            axis(side = 1,
                 labels = FALSE) 
        }
        box(which = "plot", bty = "l")
        
        
        
        
    }
}


title(xlab = "ICE",
      ylab = "pdf",
      outer = TRUE, line = 4, cex.lab=3)
dev.off()

#Fig3
tikz("Fig3b.tex",width=12,height=16)
par(mfrow = c(5,3),
    oma = c(6,6,0,0) + 0.1,
    mar = c(0,0,1,1) + 0.1)


    for (rho in c(0,0.25,0.5,0.75,1)){
        for (n in c(200,2000,20000)){
            
      
        
       
        if(n<20000){
            cdensity <- data.frame(read_csv(paste("~/HPC/GRFoutput/CUSTOMDensityEstimate-",n,"-1.6-1.4-2-",rho,".txt",sep=""), 
                                            col_names = FALSE))
        }
        if(n==20000){
            cdensity <- data.frame(read_csv(paste("~/HPC/GRFoutput/CUSTOMDensityEstimate-",n,"-1.6-1.4-2-",rho,".txt",sep=""), 
                                           col_names = FALSE)) 
        }
            
            cdensity[] <- lapply(cdensity, function(x) as.numeric(as.character(x)))        
        
      
        rgf.custom.mean<-10*colMeans(data.frame(cdensity)[,-c(1,2)],na.rm=T)
        
        plot(seq(-6,6,0.1),rgf.custom.mean,col=clr6,type="l",xlim=c(-4,4),ylim=c(0,0.75),ylab="pdf",xlab="Y1-Y0",axes=F)
        
        ti<-seq(-6,6,0.1)
        lb<-10*apply(cdensity[,-c(1,2)],2,function(x){quantile(x,0.975,na.rm=T)})
        ub<-10*apply(cdensity[,-c(1,2)],2,function(x){quantile(x,0.025,na.rm=T)})
        
        polygon(x=c(rev(ti), ti), y=c(rev(lb), ub), col=clr6b, border=NA)
        
        
        lines(seq(-6,6,0.1),rgf.custom.mean,col=clr6)
        
        lines(seq(-6,6,0.1),sapply(seq(-6,6,0.1),function(x){(0.5*dnorm(x,0.55,sqrt(0.15^2+1.4^2))+0.5*dnorm(x,0.45,sqrt(0.15^2+1.4^2)))}),col=clr1)
        #   }
        
        
        
        if(n==200){
        axis(side = 2,
             labels = TRUE, cex.axis=2)
        }
        if(n>200){
            axis(side = 2,
                 labels = FALSE)   
        }
        
        if(rho==1){
        axis(side = 1,
             labels = TRUE, cex.axis=2)
        }
        if(rho<1){
            axis(side = 1,
                 labels = FALSE) 
        }
        box(which = "plot", bty = "l")
        
       
        
        
    }
}


title(xlab = "ICE",
      ylab = "pdf",
      outer = TRUE, line = 4, cex.lab=3)
dev.off()
