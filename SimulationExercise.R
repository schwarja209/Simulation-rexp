library(plyr) #load libraries
library(dplyr)

lambda<-0.2 #assign default lambda
n<-40 #assign default number of samples
nsim1<-1000 #assign default number of simulations

set.seed(1987) #set seed for reproducibility

sim1<-matrix(data=rexp(n*nsim1,rate=lambda),nrow=nsim1,ncol=n) #create raw simulation matrix for rexp()
sim1_mean<-data.frame(Mean=rowMeans(sim1)) #take mean of each simulation, and store in data frame
sim1_mean<-mutate(sim1_mean,Index=row.names(sim1_mean))%>%select(Index,Mean) #tidy data frame

#graph simulation
hist(sim1_mean$Mean,breaks=100,freq=FALSE,col="lightblue",main="Simulation of rexp()",xlab="Sample Means")


mean_sample1<-mean(sim1_mean$Mean) #find sample mean
mean_theory1<-1/lambda #find population/theoretical mean

print(mean_sample1)
print(mean_theory1)

#simulation with means
hist(sim1_mean$Mean,breaks=100,freq=FALSE,col="lightblue",main="Simulation of rexp()",xlab="Sample Means")
abline(v=mean_sample1,col="orange",lwd=2)
abline(v=mean_theory1,col="green",lwd=2)


var_sample1<-var(sim1_mean$Mean) #find sample variance
var_theory1<-(1/lambda)^2/n #find population/theoretical variance

print(var_sample1)
print(var_theory1)

x1<-seq(min(sim1_mean$Mean),max(sim1_mean$Mean),length=100) #create x spread for population density curve
y1<-dnorm(x1,mean=mean_theory1,sd=sqrt(var_theory1)) #create y distribution for population density

#graph simulation with density curves
hist(sim1_mean$Mean,breaks=100,freq=FALSE,col="lightblue",main="Simulation of rexp()",xlab="Sample Means")
lines(density(sim1_mean$Mean),col="orange",lwd=3) #create sample density curve
lines(x1,y1,col="green",lwd=3)


conf95_sample1<-mean_sample1+c(-1,1)*qnorm(.975)*sqrt(var_sample1)/sqrt(n) #find 95% confidence intervals
conf95_theory1<-mean_theory1+c(-1,1)*qnorm(.975)*sqrt(var_theory1)/sqrt(n)

print(conf95_sample1)
print(conf95_theory1)


nsim2<-1000000 #increase number of simulations to demonstrate convergence toward normality
set.seed(1987) #reset seed for reproducibility

sim2<-matrix(data=rexp(n*nsim2,rate=lambda),nrow=nsim2,ncol=n) #recalculate simulation matrix
sim2_mean<-data.frame(Mean=rowMeans(sim2)) #recalculate mean of each simulation
sim2_mean<-mutate(sim2_mean,Index=row.names(sim2_mean))%>%select(Index,Mean)

x2<-seq(min(sim2_mean$Mean),max(sim2_mean$Mean),length=100) #recalculate population density curve
y2<-dnorm(x2,mean=mean_theory1,sd=sqrt(var_theory1))

#graph new simulation with density curves
hist(sim2_mean$Mean,breaks=100,freq=FALSE,col="lightblue",main="Simulation of rexp()",xlab="Sample Means")
lines(density(sim2_mean$Mean),col="orange",lwd=3) #create sample density curve
lines(x2,y2,col="green",lwd=3)