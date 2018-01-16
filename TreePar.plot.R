# codes for Zhang et al. 2017. Diversification of Rosaceae since the Late Cretaceous based on plastid phylogenomics
set.seed(12345)
library(TreePar)

"""update plot"""
bd.shifts.plot<-function(resall,shifts,outfile,timemax=100,ratemin=-1,ratemax=1,plotturnover=FALSE) {
	pick<-function(resall,shifts){
		resall[[shifts+1]]
	}
	estimatesall<-sapply(resall,pick,shifts=shifts)
	pdf(outfile, pointsize=12)
	plot(c(0,-timemax),c(ratemin,ratemax),col="white",xlab="time before the present",ylab="diversification rate")
	
	for (i in 1:length(estimatesall[1,])){
		estimates<-estimatesall[,i]
		rates<-length(estimates)/3
		estimates<-estimates[-1]
		if (rates>1){
			time<-estimates[(length(estimates)-rates+2):length(estimates)]
			time<-sort(c(time,time,0,timemax))
			turnover<-estimates[1]
			div<-estimates[rates+1]
			for (j in 1:(rates-1)){turnover<-c(turnover,estimates[j:(j+1)])
				div<-c(div,estimates[(rates+j):(rates+j+1)])
			}
			turnover<-c(turnover,estimates[rates])
			div<-c(div,estimates[2*rates])}
		else {time<-c(0,timemax)
			turnover<-c(estimates[1],estimates[1])
			div<-c(estimates[2],estimates[2])
		}
	if (plotturnover==TRUE){lines(-time,turnover,col="red")}
	lines(-time,div,col="grey")
	}
	for (i in 1:1){
		estimates<-estimatesall[,i]
		rates<-length(estimates)/3
		estimates<-estimates[-1]
		if (rates>1){
			time<-estimates[(length(estimates)-rates+2):length(estimates)]
			time<-sort(c(time,time,0,timemax))
			turnover<-estimates[1]
			div<-estimates[rates+1]
			for (j in 1:(rates-1)){turnover<-c(turnover,estimates[j:(j+1)])
				div<-c(div,estimates[(rates+j):(rates+j+1)])
			}
			turnover<-c(turnover,estimates[rates])
			div<-c(div,estimates[2*rates])}
		else {time<-c(0,timemax)
			turnover<-c(estimates[1],estimates[1])
			div<-c(estimates[2],estimates[2])
		}
	if (plotturnover==TRUE){lines(-time,turnover,col="red")}
	lines(-time,div,col="red")
	}
	dev.off()
}


#
con_tre <- read.nexus("ingroup.consensus.tre")
mul_tre <- read.nexus("ingroup.1000rep.trees")
results = bd.shifts.optim(sort(getx(con_tre), decreasing=TRUE), c(0.044, 1, 1, 1, 1, 1), 1, 0, 86)[[2]]
results_aic <- sapply(results, function(pars) 2*(length(pars) - 1)+2*pars[1])
pchisq(2*(results[[1]][1]-results[[2]][1]), 3)

results <- bd.shifts.optim(sort(getx(con_tre), decreasing=TRUE), c(0.044, 1), 1, 5, 80)[[2]]
rep_res <- list()
rep_res[[1]] <- results
for (i in 1:1000) {x <- sort(getx(mul_tre[[i]]), decreasing=TRUE); rep_res[[i+1]] <- bd.shifts.optim(x, c(0.044, 1), 1, 0, 90)[[2]]; print("===================="); print(i); print("====================")}

bd.shifts.plot(rep_res, 2, "~/Documents/20160930-Rosaceae-diversification/Rosaceae.diversification.pdf", ratemin=0, ratemax=0.5, timemax=90)  # plot of 1 shifts
