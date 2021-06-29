args = commandArgs(trailingOnly=TRUE)

library(eSMC)

##############
# Set Parameters for eSMC
##############
repeats= as.integer(args[1])
M= as.integer(args[2])
tag=args[3]

cat(paste("repeats: ", toString(repeats), ", haps: ", toString(M), sep= ""))

gen=1 # one year generation time for rice
mu=1.3*10^(-8) # mutation rate per position per generation 
r= 1*10^(-8)  # recombination rate per position per generation
rho=r/mu        # ratio recombination/mutation 

# set one to true (rest to false) 
ER=F # false to estimate recombination rate
SF=T # true to estimate selfing rate
SB=F # false to estimate germination rate

# set boundaries
BoxB=c(0.05,1) #  min and max value of germination rate 
Boxs=c(0,0.99) #  min and max value of selfing rate

#######################
path= "./"
timez=list()
sizes=list()
rep_list= list()

##
for (rep in 1:repeats){

filename= paste("rep", toString(rep),"_multihetsep.txt", sep= "")
data=Get_real_data(path,M,filename,delim="\t")

result=eSMC(n=30,rho=rho,data,BoxB=BoxB,Boxs=Boxs,SB=SB,SF=SF,Rho=ER,Check=F,NC=1)
rep_list[[rep]]= result

cat(names(result))
## for plotting and stats
time_rep=result$Tc * (result$mu/mu) * gen
size_rep=result$Xi * (result$mu/(2*mu))

timez<- c(timez, time_rep)
sizes<- c(sizes, size_rep)

}
##

Plot_esmc_results(rep_list,mu,WP=T,LIST=T)

for (res in rep_list) {
cat(paste("selfing: ", toString(res$Self), sep= ""))
}

save(rep_list,file=paste(tag,"_result_list.RData",sep=""))





