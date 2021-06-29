#install.packages("pls",lib="/home/x_garciaj/Rpackages/",repos = "http://cran.us.r-project.org")
#install.packages("MASS",lib="/home/x_garciaj/Rpackages/",repos = "http://cran.us.r-project.org")

library("pls", lib.loc="/home/x_garciaj/Rpackages/")
library("MASS", lib.loc="/home/x_garciaj/Rpackages/")

#library("pls"); library("MASS");
#read file with simulations
numComp<-35;
dir_in='/home/x_garciaj/Projects/SLiM_ABC/tempsubtrop.v3/';
directory<-'/home/x_garciaj/Projects/SLiM_ABC/tempsubtrop.v3/PLSfind/'; 
filename<-'tempsubtrop_sim_Obs0_sampling1.txt';
a<-read.table(paste(dir_in, filename, sep=''), header=T,nrows= 1000);
print(dim(a))
a <- a[ , colSums(is.na(a)) == 0]
a <- a[ , colSums(a) > 0]

print(dim(a));
stat<-a[,16:ncol(a)]; 
param<-a[,1:15]; 

print(a[1:5,1:20])
rm(a);
#standardize the params
for(i in 1:length(param)){param[,i]<-(param[,i]-mean(param[,i]))/sd(param[,i]);}
#force stat in [1,2]
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:length(stat)){
myMax<-c(myMax, max(stat[,i])); myMin<-c(myMin, min(stat[,i]));
stat[,i]<-1+(stat[,i] -myMin[i])/(myMax[i]-myMin[ i]);
}
#transform statistics via boxcox
for(i in 1:length(stat)){
d<-cbind(stat[,i], param);
mylm<-lm(as.formula(d), data=d);
myboxcox<-boxcox(mylm, lambda=seq(-20,100,1/10), interp=T, eps=1/50);
lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);
myGM<-c(myGM, mean(exp(log(stat[,i]))));
}
#standardize the BC-stat
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:length(stat)){
stat[,i]<-(stat[,i] ^lambda[i] - 1)/(lambda[i]*myGM[i] ^(lambda[i]-1));
myBCSDs<-c(myBCSDs, sd(stat[,i]));
myBCMeans<-c(myBCMeans, mean(stat[,i]));
stat[,i]<-(stat[,i] -myBCMeans[i])/myBCSDs[i];
}
#perform pls
#myPlsr<-plsr(as.matrix(stat) ~ as.matrix(param), scale=F, ncomp=numComp,
#validation='LOO');

print(dim(as.matrix(stat)))
myPCA<- prcomp(as.matrix(stat),scale= T, center= T);
print(dim(myPCA$rotation))
print(myPCA$sdev)
myloadings<-myPCA$rotation[,1:numComp] %*% diag(myPCA$sdev, numComp, numComp);
#print(names(myPCA))

print(dim(as.matrix(stat)))
print(dim(as.matrix(myPCA$x)))

myPlsr<-plsr(as.matrix(stat) ~ as.matrix(myPCA$x), scale=F, ncomp=numComp,
validation='LOO');
#write pls to a file
myPlsrDataFrame<-data.frame(comp1=myloadings[,1]);
table_name <- paste(directory, 'PLSfile_', filename, sep='');
print(table_name);
for(i in 2:numComp){myPlsrDataFrame<-cbind(myPlsrDataFrame, myloadings[,i]);}
write.table(cbind(names(stat), myMax, myMin, lambda, myGM, myBCMeans, myBCSDs,
myPlsrDataFrame), file=paste(directory, 'PLSfile_', filename, sep=''), col.names=F,
row.names=F, sep='\t', quote=F);
#make RMSEP plot
pdf(paste(directory, 'RMSEP_', filename, '.pdf', sep=''));
plot(RMSEP(myPlsr), legendpos = "topright");
dev.off();
