setwd("~/Documents/SomaticMutations/OfuAug/")
library("reshape2")
library(ggplot2)
library(stringr)
library(sciplot)
library(sinaplot)
library(ggforce)
library(gridExtra)
sessionInfo()


files<-list.files(path="~/Documents/SomaticMutations/OfuAug/WithVerifications/WithDepths", pattern="*FOR_R.txt", full.names=T, recursive=FALSE)

par(mfrow=c(2,4)) #if your output includes plots, this will put all of the plots in one image. change the (4,5) to different numbers depending on how large you want your grid of plots to be.

#then, start your loop:
lapply(files,function(x) {
	metadata<-read.delim(x) #load a file
	base<-basename(x)
	colony<-strsplit(base, "\\_")[[1]][2] #the colony name 
	
	genoanddepth<-(metadata$genotype)
	split<-str_split_fixed(genoanddepth, ",", 4)
	genotypes<-split[,1]
	totaldepth<-as.numeric(split[,2]) # total number of reads for the sample at this site
	refdepth<-as.numeric(split[,3]) # number of reads that match the REF allele
	altdepth<-as.numeric(split[,4]) # number of ALT allele reads
	
	
	metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, 	"sample"= metadata$sample, "ref" = metadata$ref, "alt" = 	metadata$alt, "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, 	"GoH_or_LoH"=metadata$DeNovo_LoH, "Ti/Tv"=metadata$TiTv, 	"WhattoWhat" = metadata$WhattoWhat, "PopulationPoly"= 			metadata$PopulationPoly, 	"VerifiedYesorNo"=metadata$VerifiedYesorNo) #turns the input file into a dataframe, adds the ref and alt depths as a column

	DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean) #outputs the mean totaldepth per locus


	metadatadf<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos") #adds the mean depth per site to the last column of the metadatadf table, under column name "totaldepth.y"
	
	####VERIFIED####
	verifiedlist<-metadatadf[ which(metadatadf$VerifiedYesorNo 			=="yes"),] #outputs all verified mutations
	
	uniqueverifiedlist<-												verifiedlist[match(unique(verifiedlist$chrom.pos), 					verifiedlist$chrom.pos),] #outputs just the first sample with a verified mutation at each site
	
	uniqueverifiedcount<-nrow(uniqueverifiedlist) # number of unique verified mutations

	### VERIFIED POP POLYS VS NONPOP (FIGURE 4) ###
	
	#verifiedpopulationpolys<-											verifiedlist[ which(verifiedlist$PopulationPoly=="TRUE"),]
		
	uniqueverifiedpopulationpolys<-									uniqueverifiedlist[ which(uniqueverifiedlist$PopulationPoly=="TRUE"),]

	uniqueverifiednonpop<-	uniqueverifiedlist[ which(uniqueverifiedlist$PopulationPoly=="FALSE"),]
	
	#uniqueverifiednonpopcount<-nrow(uniqueverifiednonpop)
	
	uniqueverifiedpopulationpolys<-uniqueverifiedlist[ which(uniqueverifiedlist$PopulationPoly=="TRUE"),]
	
	#uniqueverifiedpopulationpolyscount<-nrow(uniqueverifiedpopulationpolys)

	uniqueverifiedpopulationpolysGoH<-uniqueverifiedpopulationpolys[which(uniqueverifiedpopulationpolys$GoH_or_LoH=="DeNovo"),]

	uniqueverifiedpopulationpolysGoHcount<-nrow(uniqueverifiedpopulationpolysGoH)

	uniqueverifiedpopulationpolysLoH<-uniqueverifiedpopulationpolys[which(uniqueverifiedpopulationpolys$GoH_or_LoH=="LoH"),]

	uniqueverifiedpopulationpolysLoHcount<-nrow(uniqueverifiedpopulationpolysLoH)

	uniqueverifiednonpopGoH<-uniqueverifiednonpop[which(uniqueverifiednonpop$GoH_or_LoH=="DeNovo"),]

	uniqueverifiednonpopGoHcount<-nrow(uniqueverifiednonpopGoH)

	uniqueverifiednonpopLoH<-uniqueverifiednonpop[which(uniqueverifiednonpop$GoH_or_LoH=="LoH"),]

	uniqueverifiednonpopLoHcount<-nrow(uniqueverifiednonpopLoH)
	
	mat<-matrix(c(uniqueverifiedpopulationpolysGoHcount*100/uniqueverifiedcount, uniqueverifiednonpopGoHcount*100/uniqueverifiedcount, uniqueverifiedpopulationpolysLoHcount*100/uniqueverifiedcount, uniqueverifiednonpopLoHcount*100/uniqueverifiedcount),ncol=2,byrow=TRUE)

	colnames(mat)<-c("PopPoly", "NonPopPoly")
	rownames(mat)<-c("% GOH","% LOH")
	mat<-as.table(mat)
	#barplot(mat, main=colony, ylim=c(0,100)) # here is your stacked barplot with poppolys, nonpoppolys, and LoH/GoH - Figure 4


###GOH VS LOH for each colony - FIGURE 2 c,d,e,f ###
	verDeNovolist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="DeNovo"),] # outputs all GoH mutations from the uniquified list of mutations
	
	verDeNovocount<-nrow(verDeNovolist) # number of all unique GoH mutations 
	verDeNovoprop<-verDeNovocount/uniqueverifiedcount
	
	verLoHlist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="LoH"),] # outputs all LoH mutations from the uniquified list of mutations

	verLoHcount<-nrow(verLoHlist) # number of all unique LoH mutations
	
	verLoHprop<-verLoHcount/uniqueverifiedcount
	
	verifiedDeNovoLoHDF<-data.frame(Types=c("GoH","LoH"),Proportions=c(verDeNovoprop*100,verLoHprop*100))
	verifiedDeNovoLoHplot<-barplot(verifiedDeNovoLoHDF$Proportion,names.arg=verifiedDeNovoLoHDF$Types, ylim=c(0,100), main=colony) #Figure 3 b, c, d, e
	verifiedDeNovoLoHDF$Proportion

	### COMPARING DEPTHS (FIGURE S1) ###
	verifieddepth<-(uniqueverifiedlist$totaldepth.y)
	
	falsifiedlist<-metadatadf[ which(metadatadf$VerifiedYesorNo =="no"),]
	
	uniquefalsifiedlist<-falsifiedlist[match(unique(falsifiedlist$chrom.pos), falsifiedlist$chrom.pos),]
	
	falsifieddepth<-(uniquefalsifiedlist$totaldepth.y)

	allsites<-metadatadf[match(unique(metadatadf$chrom.pos), metadatadf$chrom.pos),]

	x<-c(verifieddepth, falsifieddepth, allsites$totaldepth.y)

	groups<-c(rep("Verified",nrow(uniqueverifiedlist)),rep("Falsified",nrow(uniquefalsifiedlist)),rep("AllSites",nrow(allsites)))
	df<-data.frame(groups, x)

	p<- ggplot(df, aes(groups,x))
	p + geom_violin(aes(fill = groups))
	depthsplots<- p +geom_boxplot() + geom_sina(aes(color=groups),size=1 ) + ylim(0,85) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Read Depth") # here is your depths plot, FIGURE S1

	
})

	verifieddepth<-(uniqueverifiedlist$totaldepth.y)
	averageverdepth<-mean(verifieddepth) #this gives the mean value of average depth at each site, not just the mutant depth at each site
	SEaverageverdepth<-se(verifieddepth) #finds standard deviation of verifieddepths

####FALSIFIED####
	falsifiedlist<-metadatadf[ which(metadatadf$VerifiedYesorNo =="no"),]
	uniquefalsifiedlist<-falsifiedlist[match(unique(falsifiedlist$chrom.pos), falsifiedlist$chrom.pos),]
	uniquefalsifiedcount<-nrow(uniquefalsifiedlist)

	uniquefalsifiednonpop<-uniquefalsifiedlist[ which(uniquefalsifiedlist$PopulationPoly=="FALSE"),]
	uniquefalsifiednonpopcount<-nrow(uniquefalsifiednonpop)
	uniquefalsifiedpopulationpolys<-uniquefalsifiedlist[ which(uniquefalsifiedlist$PopulationPoly=="TRUE"),]
	uniquefalsifiedpopulationpolyscount<-	nrow(uniquefalsifiedpopulationpolys)


	falsifieddepth<-(uniquefalsifiedlist$totaldepth.y)
	averagefalsdepth<-mean(falsifieddepth)
	SEaveragefalsdepth<-se(falsifieddepth)

	###PLOT DEPTHS###
	#names<-c("Verified Mean Depth","Falsified Mean Depth")
	#means<-c(averageverdepth, averagefalsdepth)
	#SEs<-c(SEaverageverdepth, SEaveragefalsdepth)
	#plotTop<-max(means+SEs*2)
	#barCenters<-barplot(means,names.arg=names, las=1, ylim=c(0,50))
	#title(colony)
	#arrows(barCenters, means-SEs*2, barCenters, means+SEs*2, lwd=1, angle=90, code=3)

	allsites<-metadatadf[match(unique(metadatadf$chrom.pos), metadatadf$chrom.pos),]

	x<-c(verifieddepth, falsifieddepth, allsites$totaldepth.y)

	groups<-c(rep("Verified",nrow(uniqueverifiedlist)),rep("Falsified",nrow(uniquefalsifiedlist)),rep("AllSites",nrow(allsites)))
	df<-data.frame(groups, x)

	#sinaplot(x,groups, col=2:4) #make a sinaplot with different colors for each class

	p<- ggplot(df, aes(groups,x))
	#p + geom_violin(aes(fill = groups))
	#depthsplots<- p +geom_boxplot() + geom_sina(aes(color=groups),size=1 ) + ylim(0,85) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Read Depth") # here is your depths plot
	uniqueverifiednonpopGoH<-uniqueverifiednonpop[which(uniqueverifiednonpop$GoH_or_LoH=="DeNovo"),]

	uniqueverifiednonpopGoHcount<-nrow(uniqueverifiednonpopGoH)

	uniqueverifiednonpopLoH<-uniqueverifiednonpop[which(uniqueverifiednonpop$GoH_or_LoH=="LoH"),]

	uniqueverifiednonpopLoHcount<-nrow(uniqueverifiednonpopLoH)


	mat<-matrix(c(uniqueverifiedpopulationpolysGoHcount*100/uniqueverifiedcount, uniqueverifiednonpopGoHcount*100/uniqueverifiedcount, uniqueverifiedpopulationpolysLoHcount*100/uniqueverifiedcount, uniqueverifiednonpopLoHcount*100/uniqueverifiedcount),ncol=2,byrow=TRUE)
	colnames(mat)<-c("PopPoly", "NonPopPoly")
	rownames(mat)<-c("% GOH","% LOH")
	mat<-as.table(mat)
	barplot(mat, main=colony, ylim=c(0,100)) # here is your stacked barplot with poppolys, nonpoppolys, and LoH/GoH - Figure S2
	
	

})	

##Figure 2 b ##
props<-read.delim("~/Documents/SomaticMutations/OfuAug/ColonyDeNovoLoHprops.txt")
frame<-data.frame(props)
se<-function(x) sd(x)/sqrt(length(x))
se1<-se(frame[,2])
se2<-se(frame[,3])
standarderrors<-c(se1,se2)
mean1<-(mean(frame[,2]))
mean2<-(mean(frame[,3]))
names<-c("DeNovo","LoH")
means<-c(mean1, mean2)
plotTop <- max(means+standarderrors*2)
barCenters <- barplot(means, names.arg=names, col="gray", las=1, ylim=c(0,100))
arrows(barCenters, means-standarderrors*2, barCenters, means+standarderrors*2, lwd=1, angle=90, code=3)

##Figure 1 ##
#first run this wrapper script:
#I got this from https://github.com/mrxiaohe/R_Functions/blob/master/functions/bar
bar <- function(dv, factors, dataframe, percentage=FALSE, errbar=!percentage, half.errbar=TRUE, conf.level=.95, 

		xlab=NULL, ylab=NULL, main=NULL, names.arg=NULL, bar.col="black", whisker=.015,args.errbar=NULL,

		legend=TRUE, legend.text=NULL, args.legend=NULL,legend.border=FALSE, box=TRUE, args.yaxis=NULL, 

		mar=c(5,4,3,2),...){

	axes=!percentage

	dv.name<-substitute(dv)

	if(length(dv.name)>1) stop("'dv' only takes one variable")

	dv.name<-as.character(dv.name)

	dv<-dataframe[[dv.name]]

	fnames<-substitute(factors)

	if(length(fnames)==1){

		factors<-as.character(fnames)

		nf<-1

	}else{

		factors<-as.character(fnames[-1L])

		nf<-length(factors)

	}

	if(nf>2) stop("This function accepts no more than 2 factors \n",

			"\t-i.e., it only plots one-way or two-way designs.")

	if(percentage & errbar){

		warning("percentage=TRUE; error bars were not plotted")

		errbar<-FALSE

	}

	if(!percentage) xbars<-tapply(dv, dataframe[,factors], mean, na.rm=TRUE)

	else {

		xbars<-tapply(dv, list(interaction(dataframe[,factors], lex.order=TRUE)), mean, na.rm=TRUE)

		if(sum(na.omit(dv)!=0&na.omit(dv)!=1)>0) 

			stop("Data points in 'dv' need to be 0 or 1 in order to set 'percentage' to TRUE")

		xbars<-rbind(xbars, 1-xbars)*100

	}

	if(errbar){

		se<-tapply(dv, dataframe[,factors], sd, na.rm=TRUE)/sqrt(tapply(dv, dataframe[,factors], length))

		conf.level=1-(1-conf.level)/2

		lo.bar<-xbars-se*qnorm(conf.level)

		hi.bar<-xbars+se*qnorm(conf.level)	

	}

	extras<-list(...)

	if(legend & !percentage){

		if(is.null(legend.text))

			legend.text<-sort(unique(dataframe[[factors[1]]]))

		args.legend.temp<-list(x="topright", bty=if(!legend.border)"n" else "o",

							   inset=c(0,0))

		if(is.list(args.legend))

			args.legend<-modifyList(args.legend.temp, args.legend)

		else 

			args.legend<-args.legend.temp

	} else if(legend & percentage){

		if(is.null(legend.text)) 

			legend.text<-c("1", "0")

		args.legend.temp<-list(x="topright", bty=if(!legend.border)"n" else "o",

							   inset=c(0,0))

		if(is.list(args.legend))

			args.legend<-modifyList(args.legend.temp, args.legend)

		else 

			args.legend<-args.legend.temp

	} else if(!legend){

		args.legend<-NULL

		legend.text<-NULL

	}

	if(errbar && legend && !percentage) ymax<-max(hi.bar)+max(hi.bar)/20

	else if(errbar && legend && percentage) ymax<-115

	else if(errbar && !legend) ymax <- max(xbars)

	else if(!errbar && legend && percentage) ymax<-110	

	else if(!errbar) ymax<-max(xbars) + max(xbars)/20

	if(!percentage){

		args.barplot<-list(beside=TRUE, height=xbars, ylim=c(0, ymax), main=main, names.arg=names.arg,

				col=hcl(h=seq(0,270, 270/(length(unique(dataframe[[factors[1]]]))))[-length(unique(dataframe[[factors[1]]]))]),

				legend.text=legend.text, args.legend=args.legend, xpd=TRUE,

				xlab=if(is.null(xlab)) factors[length(factors)] else xlab,

				ylab=if(is.null(ylab)) dv.name else ylab, axes=axes)

	}else{

		args.barplot<-list(beside=TRUE, height=xbars, ylim=c(0, ymax),  main=main, names.arg=names.arg,

				col=hcl(h=seq(0,270, 270/(length(unique(dataframe[[factors[1]]]))))[-length(unique(dataframe[[factors[1]]]))]),

				legend.text=legend.text, args.legend=args.legend, xpd=TRUE,

				xlab=if(is.null(xlab)) " "[length(factors)] else xlab,

				ylab=if(is.null(ylab)) "percentage" else ylab, axes=axes)		

	}

	args.barplot<-modifyList(args.barplot, extras)

	errbars = function(xvals, cilo, cihi, whisker, nc, args.errbar = NULL, half.errbar=TRUE) {

		if(half.errbar){

			cilo<-(cihi+cilo)/2

		}

		fixedArgs.bar = list(matlines, x=list(xvals), 

                       		 y=lapply(split(as.data.frame(t(do.call("rbind", 

                       		 list(cihi, cilo)))),1:nc),matrix, 

                       		 nrow=2, byrow=T))

	  	allArgs.bar = c(fixedArgs.bar, args.errbar)

 	 	whisker.len = whisker*(par("usr")[2] - par("usr")[1])/2

 	 	whiskers = rbind((xvals - whisker.len)[1,],

        	             (xvals + whisker.len)[1,])

  		fixedArgs.lo = list(matlines, x=list(whiskers), 	

  	         				y=lapply(split(as.data.frame(t(do.call("rbind", 

                      		list(cilo, cilo)))), 1:nc), matrix, nrow=2, byrow=T))

	  	allArgs.bar.lo = c(fixedArgs.lo, args.errbar)

		fixedArgs.hi = list(matlines, x=list(whiskers), 

		  					y=lapply(split(as.data.frame(t(do.call("rbind", 

	                      	list(cihi, cihi)))), 1:nc), matrix, nrow=2, byrow=T))

	  	allArgs.bar.hi = c(fixedArgs.hi, args.errbar)  

		invisible(do.call(mapply, allArgs.bar))

	  	if(!half.errbar) invisible(do.call(mapply, allArgs.bar.lo))

	  	invisible(do.call(mapply, allArgs.bar.hi))

	}

	par(mar=mar)

	errloc<-as.vector(do.call(barplot, args.barplot))

	if(errbar){

		errloc<-rbind(errloc, errloc)

		lo.bar<-matrix(as.vector(lo.bar))

		hi.bar<-matrix(as.vector(hi.bar))

		args.errbar.temp<-list(col=bar.col, lty=1)

		args.errbar<-if(is.null(args.errbar)|!is.list(args.errbar)) 

		                args.errbar.temp

                 	 else if(is.list(args.errbar)) 

                 	 	modifyList(args.errbar.temp, args.errbar)

		errbars(errloc, cilo=lo.bar, cihi=hi.bar, nc=1, whisker=whisker, 

				args.errbar=args.errbar, half.errbar=half.errbar)

	}

	if(box) box()

	if(percentage){

		args.yaxis.temp<-list(at=seq(0,100, 20), las=1)

		args.yaxis<-if(!is.list(args.yaxis)) args.yaxis.temp else modifyList(args.yaxis.temp, args.yaxis)

		do.call(axis, c(side=2, args.yaxis))

	}

}
#load file with colony information: 
stats<-read.delim("~/Documents/SomaticMutations/PrimersforVerification/ColonyMutationStats.txt")
#create factors for the four categories:
f1 <- factor(c("Putative","Inconclusive","Falsified","Verified"), levels = c("Putative","Inconclusive","Falsified","Verified"))
#turn the categories into a dataframe
my.data <- data.frame(Names = f1)
levels(my.data$Names)
stats$Type <-factor(c("Putative","Inconclusive","Falsified","Verified"), levels = c("Putative","Inconclusive","Falsified","Verified"))
#then apply the wrapper script to your data:
bar(dv = Number,
	factors = c(Type, Colony),
	dataframe = stats,
	errbar =FALSE,
	ylim = c(0, 240),
	col=c("firebrick","black","gray","blue"))

### to run an anlysis on all the combined data:###
files<-list.files(path="~/Documents/SomaticMutations/OfuAug/WithVerifications/WithDepths", pattern="*FOR_R.txt", full.names=T, recursive=FALSE)
metadata= NULL
for (i in 1:length(files)) {
	file =files[i]
	data<-read.delim(file)
	base<-basename(file)
	colony<-strsplit(base, "\\_")[[1]][2]
	metadata <- rbind(metadata, data.frame(data))
}

genoanddepth<-(metadata$genotype)
	split<-str_split_fixed(genoanddepth, ",", 4)
	genotypes<-split[,1]
	totaldepth<-as.numeric(split[,2])
	refdepth<-as.numeric(split[,3])
	altdepth<-as.numeric(split[,4])
	metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, 	"sample"= metadata$sample, "ref" = metadata$ref, "alt" = 	metadata$alt, "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, 	"GoH_or_LoH"=metadata$DeNovo_LoH, "Ti/Tv"=metadata$TiTv, 	"WhattoWhat" = metadata$WhattoWhat, "PopulationPoly"= 			metadata$PopulationPoly, 	"VerifiedYesorNo"=metadata$VerifiedYesorNo)

	DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean)


	metadatadf<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos")
	
	uniquemetadatadf<-metadatadf[match(unique(metadatadf$chrom.pos), metadatadf$chrom.pos),]
	poppolys<-uniquemetadatadf[ which(uniquemetadatadf$PopulationPoly 			=="TRUE"),]
	poppolyscount<-nrow(poppolys)
	nonpoppolys<-uniquemetadatadf[ which(uniquemetadatadf$PopulationPoly 			=="FALSE"),]
	nonpoppolyscount<-nrow(nonpoppolys)
	
	####VERIFIED####
	verifiedlist<-metadatadf[ which(metadatadf$VerifiedYesorNo 			=="yes"),]
	verifiedpopulationpolys<-											verifiedlist[ which(verifiedlist$PopulationPoly=="TRUE"),]
	uniqueverifiedlist<-												verifiedlist[match(unique(verifiedlist$chrom.pos), 					verifiedlist$chrom.pos),]
	uniqueverifiedpopulationpolys<-									uniqueverifiedlist[ which(uniqueverifiedlist$PopulationPoly=="TRUE"),]

	uniqueverifiedcount<-nrow(uniqueverifiedlist)

	uniqueverifiednonpop<-	uniqueverifiedlist[ which(uniqueverifiedlist$PopulationPoly=="FALSE"),]
	uniqueverifiednonpopcount<-nrow(uniqueverifiednonpop)
	uniqueverifiedpopulationpolys<-uniqueverifiedlist[ which(uniqueverifiedlist$PopulationPoly=="TRUE"),]
	uniqueverifiedpopulationpolyscount<-nrow(uniqueverifiedpopulationpolys)

	uniqueverifiedpopulationpolysGoH<-uniqueverifiedpopulationpolys[which(uniqueverifiedpopulationpolys$GoH_or_LoH=="DeNovo"),]

	uniqueverifiedpopulationpolysGoHcount<-nrow(uniqueverifiedpopulationpolysGoH)

	uniqueverifiedpopulationpolysLoH<-uniqueverifiedpopulationpolys[which(uniqueverifiedpopulationpolys$GoH_or_LoH=="LoH"),]

	uniqueverifiedpopulationpolysLoHcount<-nrow(uniqueverifiedpopulationpolysLoH)
	
	names<-c("GoH","LoH")
	proportions<-c(uniqueverifiedpopulationpolysGoHcount*100/uniqueverifiedpopulationpolyscount, uniqueverifiedpopulationpolysLoHcount*100/uniqueverifiedpopulationpolyscount)
	#barplot(proportions, names.arg=names, main=colony, ylim=c(0,100))



###GOH VS LOH###
	verDeNovolist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="DeNovo"),]
	verDeNovocount<-nrow(verDeNovolist)
	verDeNovoprop<-verDeNovocount/uniqueverifiedcount
	verLoHlist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="LoH"),]
	verLoHcount<-nrow(verLoHlist)
		verLoHprop<-verLoHcount/uniqueverifiedcount

	verifieddepth<-(uniqueverifiedlist$totaldepth.y)
	averageverdepth<-mean(verifieddepth) #this gives the mean value of average depth at each site, not just the mutant depth at each site
	SEaverageverdepth<-se(verifieddepth) #finds standard deviation of verifieddepths

###MUTATION SPECTRA###
verATGClist<-uniqueverifiedlist[ which(uniqueverifiedlist$WhattoWhat=="AtoG" | uniqueverifiedlist$WhattoWhat=="TtoC"),]
	verATGCcount<-nrow(verATGClist)
	verATGCprop<-verATGCcount/uniqueverifiedcount
	SEverATGCprop<-se(verATGCprop)
	verGCATlist<-uniqueverifiedlist[ which(uniqueverifiedlist$WhattoWhat=="GtoA" | uniqueverifiedlist$WhattoWhat=="CtoT"),]
	verGCATcount<-nrow(verGCATlist)
	verGCATprop<-verGCATcount/uniqueverifiedcount

	verATTAlist<-uniqueverifiedlist[ which(uniqueverifiedlist$WhattoWhat=="AtoT" | uniqueverifiedlist$WhattoWhat=="TtoA"),]
	verATTAcount<-nrow(verATTAlist)
	verATTAprop<-verATTAcount/uniqueverifiedcount

	verACTGlist<-uniqueverifiedlist[ which(uniqueverifiedlist$WhattoWhat=="AtoC" | uniqueverifiedlist$WhattoWhat=="TtoG"),]
	verACTGcount<-nrow(verACTGlist)
	verACTGprop<-verACTGcount/uniqueverifiedcount

	verGCCGlist<-uniqueverifiedlist[ which(uniqueverifiedlist$WhattoWhat=="CtoG" | uniqueverifiedlist$WhattoWhat=="GtoC"),]
	verGCCGcount<-nrow(verGCCGlist)
	verGCCGprop<-verGCCGcount/uniqueverifiedcount

	verGCTAlist<-uniqueverifiedlist[ which(uniqueverifiedlist$WhattoWhat=="GtoT" | uniqueverifiedlist$WhattoWhat=="CtoA"),]
	verGCTAcount<-nrow(verGCTAlist)
	verGCTAprop<-verGCTAcount/uniqueverifiedcount
	vertypesDF<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), Proportion=c(verATGCprop, 	verGCATprop, verATTAprop,verACTGprop, verGCCGprop, verGCTAprop))
	
#"A>G,T>C","G>A,C>T","A>T,T>A","A>C,T>G","G>C,C>G","G>T,C>A"
#FIGURE 4:
par(mfrow=c(1,4))
	vertypesDFplot<-barplot(vertypesDF$Proportion,names.arg=vertypesDF$Types, ylim=c(0,0.7), main="Mutation spectrum across 4 colonies", las=2) #Figure 4a
	
#humantypesDF<-data.frame(Types=c("ATtoGC","GCtoAT","ATtoTA","ATtoCG","GCtoCG","GCtoTA"), Proportion=c(.221, .408, .067, .078, .125, .110))
#humantypesDFplot<-barplot(humantypesDF$Proportion, names.arg=humantypesDF$Types, ylim=c(0,0.7), main="Human Germline Spectrum",las=2)# find the citation

#humantypesDF<-data.frame(Types=c("ATtoGC","GCtoAT","ATtoTA","ATtoCG","GCtoCG","GCtoTA"), Proportion=c(.02, .57, .02, .24, .08, .08))
#humantypesDFplot<-barplot(humantypesDF$Proportion, names.arg=humantypesDF$Types, ylim=c(0,0.7), main="Human Germline Spectrum",las=2)# Greenman et al 2007

humantypesDF<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), Proportion=c(.27, .24, .07, .06, .09, .09))
humantypesDFplot<-barplot(humantypesDF$Proportion, names.arg=humantypesDF$Types, ylim=c(0,0.7), main="Human Germline Spectrum- Rahbari",las=2)# Rahbari et al 2016 NATURE GENETICS	
#Figure 4b

humanesophagus<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), Proportion=c(.19, .42, .09, .06, .07, .16))
humanesophagusplot<-barplot(humanesophagus$Proportion, names.arg=humanesophagus$Types, ylim=c(0,0.7), main="Human Esophagus Spectrum",las=2)# Martincorena et al 2018 Science
#Figure 4c

#y<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"),Proportion=c(.324,.378, .081, .054, .081, .081))
#barplot(y$Proportion, names.arg=y$Types,col="gray", las=1,ylim=c(0,0.7), main="Distribution of SNP types for A. palmata", las=2)

UVtypesDF<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), Proportion=c(2500/35850, 24000/35850, 2400/35850, 2450/35850, 500/35850, 4000/35850))
UVtypesDFplot<-barplot(UVtypesDF$Proportion, names.arg=UVtypesDF$Types, ylim=c(0, 0.7), main="High UV Spectrum",las=2) #Pleasance et al. 2009
#Figure 4d
#y<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"),Proportion=c(.66,.21, .1, .009, .014, .009))
#barplot(y$Proportion, names.arg=y$Types,col="gray",ylim=c(0,0.7), main="Distribution of error types for Taq DNA Polymerase",las=2)	 #Potapov and Ong

##Transiton/Transversion ratio for all verified muts across 4 colonies###
Transtions<- verATGCprop + verGCATprop
Transversions<- verATTAprop+verACTGprop+verGCCGprop+verGCTAprop
TiTv<- Transtions/Transversions
