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

par(mfrow=c(1,4)) #if your output includes plots, this will put all of the plots in one image. change the (4,5) to different numbers depending on how large you want your grid of plots to be.

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
	
	DeNovos<-subset(metadatadf, GoH_or_LoH=="DeNovo") #pulls out just the DeNovo (GOH) mutations
	
	LoH<-subset(metadatadf, GoH_or_LoH =="LoH") #pulls out just the mutations labeled LOH
		
	trueLoH<-subset(LoH, refdepth =="0" | altdepth=="0") #pulls out just the "true" LOH mutations- that is, removes all sites labeled LOH where either the minor allele is in fact present but at less than 10% frequency
	
	metadatadf<-rbind( DeNovos, trueLoH)
	uniquemetadatadf<-											metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),] #outputs just the first sample with a verified mutation at each site
  trueLoHunique<-subset(uniquemetadatadf, GoH_or_LoH =="LoH")
  DeNovosunique<-subset(uniquemetadatadf, GoH_or_LoH=="DeNovo")

	####VERIFIED####
	verifiedlist<-metadatadf[ which(metadatadf$VerifiedYesorNo 			=="yes"),] #outputs all verified mutations
	
		DeNovos<-subset(verifiedlist, GoH_or_LoH=="DeNovo")
	LoH<-subset(verifiedlist, GoH_or_LoH =="LoH")# && refdepth =="0" | uniqueverifiedlist$altdepth=="0") )
	trueLoH<-subset(LoH, refdepth =="0" | altdepth=="0")
	verifiedlist<-rbind( DeNovos, trueLoH)
	
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
	#barplot(mat, main=colony, ylim=c(0,100), ylab =("Percent of Verified Mutations") ) # here is your stacked barplot with poppolys, nonpoppolys, and LoH/GoH - Figure 4


###GOH VS LOH for each colony - FIGURE 2 c,d,e,f ###
	
	verDeNovolist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="DeNovo"),] # outputs all GoH mutations from the uniquified list of mutations
	
	verDeNovocount<-nrow(verDeNovolist) # number of all unique GoH mutations 
	verDeNovoprop<-verDeNovocount/uniqueverifiedcount
	
	verLoHlist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="LoH"),] # outputs all LoH mutations from the uniquified list of mutations
	verLoHlist<-uniqueverifiedlist[ which(uniqueverifiedlist$refdepth =="0" | uniqueverifiedlist$altdepth=="0"),]
	verLoHcount<-nrow(verLoHlist) # number of all unique LoH mutations

	verLoHprop<-verLoHcount/uniqueverifiedcount
	
	verifiedDeNovoLoHDF<-data.frame(Types=c("GoH","LoH"),Proportions=c(verDeNovoprop*100,verLoHprop*100))
	#verifiedDeNovoLoHplot<-barplot(verifiedDeNovoLoHDF$Proportion,names.arg=verifiedDeNovoLoHDF$Types, ylim=c(0,100), main=colony, ylab="Percent of Verified Mutations") #Figure 3 b, c, d, e
	#verifiedDeNovoLoHDF$Proportion

	### COMPARING DEPTHS (FIGURE S1a-d) ###
	verifieddepth<-(uniqueverifiedlist$totaldepth.y)
	
	falsifiedlist<-metadatadf[ which(metadatadf$VerifiedYesorNo =="no"),]
	
	uniquefalsifiedlist<-falsifiedlist[match(unique(falsifiedlist$chrom.pos), falsifiedlist$chrom.pos),]
	
	falsifieddepth<-(uniquefalsifiedlist$totaldepth.y)

	allsites<-metadatadf[match(unique(metadatadf$chrom.pos), metadatadf$chrom.pos),]

	x<-c(verifieddepth, falsifieddepth, allsites$totaldepth.y)

	groups<-c(rep("Verified",nrow(uniqueverifiedlist)),rep("Falsified",nrow(uniquefalsifiedlist)),rep("AllSites",nrow(allsites)))
	df<-data.frame(groups, x)
	
	#test for signficance in difference between means:
	#allVSver<-t.test(allsites$totaldepth.y,verifieddepth)
	wilcox.test(allsites$totaldepth.y,verifieddepth) #use wilcoxon instead of t test
	#falseVSver<-t.test(falsifieddepth,verifieddepth)
	wilcox.test(falsifieddepth,verifieddepth) #use wilcoxon instead of t test
})
#	p<- ggplot(df, aes(groups,x))
#	p + geom_violin(aes(fill = groups))
#	depthsplots<- p +geom_boxplot() + ggtitle(colony) + geom_sina(aes(color=groups),size=1 ) + ylim(0,85) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Read Depth") # here is your depths plot, FIGURE S1

	


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

###pop vs nonpop poly Figure S2###
	mat<-matrix(c(uniqueverifiedpopulationpolysGoHcount*100/uniqueverifiedcount, uniqueverifiednonpopGoHcount*100/uniqueverifiedcount, uniqueverifiedpopulationpolysLoHcount*100/uniqueverifiedcount, uniqueverifiednonpopLoHcount*100/uniqueverifiedcount),ncol=2,byrow=TRUE)
	colnames(mat)<-c("PopPoly", "NonPopPoly")
	rownames(mat)<-c("% GOH","% LOH")
	mat<-as.table(mat)
	#barplot(mat, main=colony, ylim=c(0,100), legend.text = c("GoH","LoH")) # here is your stacked barplot with poppolys, nonpoppolys, and LoH/GoH - Figure S2
	
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
	vertypesDF<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), coralProportion=c(verATGCprop, 	verGCATprop, verATTAprop,verACTGprop, verGCCGprop, verGCTAprop))
  ###FIGURE S3 - MUTATION SPECTRUM FOR EACH COLONY
	vertypesDFplot<-barplot(vertypesDF$coralProportion,names.arg=vertypesDF$Types, ylim=c(0,0.7), main=colony, las=2) #Figure 4a
	

})	

##Figure 2b ##
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
barCenters <- barplot(means, names.arg=names, col="gray", las=1, ylim=c(0,100), ylab="Percent of Verified Mutations")
arrows(barCenters, means-standarderrors*2, barCenters, means+standarderrors*2, lwd=1, angle=90, code=3)
wilcox.test(frame[,2],frame[,3])

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
	col=c("firebrick","black","gray","blue"),
	ylab="# of SNPs")

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
	
	
	DeNovos<-subset(metadatadf, GoH_or_LoH=="DeNovo")
	LoH<-subset(metadatadf, GoH_or_LoH =="LoH")# && refdepth =="0" | uniqueverifiedlist$altdepth=="0") )
	trueLoH<-subset(LoH, refdepth =="0" | altdepth=="0")
	metadatadf<-rbind( DeNovos, trueLoH)
	uniquemetadatadf<-											metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),] #outputs just the first sample with a verified mutation at each site
	
	poppolys<-uniquemetadatadf[ which(uniquemetadatadf$PopulationPoly 			=="TRUE"),]
	poppolyscount<-nrow(poppolys)
	nonpoppolys<-uniquemetadatadf[ which(uniquemetadatadf$PopulationPoly 			=="FALSE"),]
	nonpoppolyscount<-nrow(nonpoppolys)
	
	####VERIFIED####
	verifiedlist<-metadatadf[ which(metadatadf$VerifiedYesorNo 			=="yes"),]
	
	DeNovos<-subset(verifiedlist, GoH_or_LoH=="DeNovo")
	LoH<-subset(verifiedlist, GoH_or_LoH =="LoH")# && refdepth =="0" | uniqueverifiedlist$altdepth=="0") )
	trueLoH<-subset(LoH, refdepth =="0" | altdepth=="0")
	verifiedlist<-rbind( DeNovos, trueLoH)
	
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
	verLoHlist<-uniqueverifiedlist[ which(uniqueverifiedlist$refdepth =="0" | uniqueverifiedlist$altdepth=="0"),]
	verLoHcount<-nrow(verLoHlist)
		verLoHprop<-verLoHcount/uniqueverifiedcount
		verifiedDeNovoLoHDF<-data.frame(Types=c("GoH","LoH"),Proportions=c(verDeNovoprop*100,verLoHprop*100))
		verifiedDeNovoLoHplot<-barplot(verifiedDeNovoLoHDF$Proportion,names.arg=verifiedDeNovoLoHDF$Types, ylim=c(0,100), main=colony) #Figure 3 b, c, d, e
		
		#verifiedDeNovoLoHDF$Proportion
		
	verifieddepth<-(uniqueverifiedlist$totaldepth.y)
	averageverdepth<-mean(verifieddepth) #this gives the mean value of average depth at each site, not just the mutant depth at each site
	SEaverageverdepth<-se(verifieddepth) #finds standard deviation of verifieddepths
	
	
##LOH: to REF or to ALT? ##

toREF<-verLoHlist[which(verLoHlist$altdepth=="0"),]
 toREFcount<-nrow(toREF)	
 

toALT<-verLoHlist[which(verLoHlist$refdepth=="0"),]
toALTcount<-nrow(toALT)


	### COMPARING DEPTHS FOR COMBINED COLONY DATA (FIGURE S1E ###
	verifieddepth<-(uniqueverifiedlist$totaldepth.y)
	
	falsifiedlist<-metadatadf[ which(metadatadf$VerifiedYesorNo =="no"),]
	
	uniquefalsifiedlist<-falsifiedlist[match(unique(falsifiedlist$chrom.pos), falsifiedlist$chrom.pos),]
	
	falsifieddepth<-(uniquefalsifiedlist$totaldepth.y)

	allsites<-metadatadf[match(unique(metadatadf$chrom.pos), metadatadf$chrom.pos),]

	x<-c(verifieddepth, falsifieddepth, allsites$totaldepth.y)

	groups<-c(rep("Verified",nrow(uniqueverifiedlist)),rep("Falsified",nrow(uniquefalsifiedlist)),rep("AllSites",nrow(allsites)))
	df<-data.frame(groups, x)
	
	#test for signficance in difference between means for COMBINED COLONY DATA:
#	allVSver<-t.test(allsites$totaldepth.y,verifieddepth)
	wilcox.test(allsites$totaldepth.y,verifieddepth)
#	falseVSver<-t.test(falsifieddepth,verifieddepth)
  wilcox.test(falsifieddepth,verifieddepth)
	p<- ggplot(df, aes(groups,x))
	#p + geom_violin(aes(fill = groups))
	depthsplotscombined<- p +geom_boxplot() + geom_sina(aes(color=groups),size=1 ) + ylim(0,85) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Read Depth") # here is your depths plot, FIGURE S1
  depthsplotscombined
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
	vertypesDF<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), coralProportion=c(verATGCprop, 	verGCATprop, verATTAprop,verACTGprop, verGCCGprop, verGCTAprop))

###MUTATION SPECTRA for LOH ###
	verLoHATGClist<-verLoHlist[ which(verLoHlist$WhattoWhat=="AtoG" | verLoHlist$WhattoWhat=="TtoC"),]
	verLoHATGCcount<-nrow(verLoHATGClist)
	verLoHATGCprop<-verLoHATGCcount/verLoHcount
	SEverLoHATGCprop<-se(verLoHATGCprop)
	
	verLoHGCATlist<-verLoHlist[ which(verLoHlist$WhattoWhat=="GtoA" | verLoHlist$WhattoWhat=="CtoT"),]
	verLoHGCATcount<-nrow(verLoHGCATlist)
	verLoHGCATprop<-verLoHGCATcount/verLoHcount

	verLoHATTAlist<-verLoHlist[ which(verLoHlist$WhattoWhat=="AtoT" | verLoHlist$WhattoWhat=="TtoA"),]
	verLoHATTAcount<-nrow(verLoHATTAlist)
	verLoHATTAprop<-verLoHATTAcount/verLoHcount

	verLoHACTGlist<-verLoHlist[ which(verLoHlist$WhattoWhat=="AtoC" | verLoHlist$WhattoWhat=="TtoG"),]
	verLoHACTGcount<-nrow(verLoHACTGlist)
	verLoHACTGprop<-verLoHACTGcount/verLoHcount

	verLoHGCCGlist<-verLoHlist[ which(verLoHlist$WhattoWhat=="CtoG" | verLoHlist$WhattoWhat=="GtoC"),]
	verLoHGCCGcount<-nrow(verLoHGCCGlist)
	verLoHGCCGprop<-verLoHGCCGcount/verLoHcount

	verLoHGCTAlist<-verLoHlist[ which(verLoHlist$WhattoWhat=="GtoT" | verLoHlist$WhattoWhat=="CtoA"),]
	verLoHGCTAcount<-nrow(verLoHGCTAlist)
	verLoHGCTAprop<-verLoHGCTAcount/verLoHcount
	vertypesLoHDF<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), LoHcoralProportion=c(verLoHATGCprop, 	verLoHGCATprop, verLoHATTAprop,verLoHACTGprop, verLoHGCCGprop, verLoHGCTAprop))

###MUTATION SPECTRA for DeNovo ###
	verDeNovoATGClist<-verDeNovolist[ which(verDeNovolist$WhattoWhat=="AtoG" | verDeNovolist$WhattoWhat=="TtoC"),]
	verDeNovoATGCcount<-nrow(verDeNovoATGClist)
	verDeNovoATGCprop<-verDeNovoATGCcount/verDeNovocount
	SEverDeNovoATGCprop<-se(verDeNovoATGCprop)
	
	verDeNovoGCATlist<-verDeNovolist[ which(verDeNovolist$WhattoWhat=="GtoA" | verDeNovolist$WhattoWhat=="CtoT"),]
	verDeNovoGCATcount<-nrow(verDeNovoGCATlist)
	verDeNovoGCATprop<-verDeNovoGCATcount/verDeNovocount

	verDeNovoATTAlist<-verDeNovolist[ which(verDeNovolist$WhattoWhat=="AtoT" | verDeNovolist$WhattoWhat=="TtoA"),]
	verDeNovoATTAcount<-nrow(verDeNovoATTAlist)
	verDeNovoATTAprop<-verDeNovoATTAcount/verDeNovocount

	verDeNovoACTGlist<-verDeNovolist[ which(verDeNovolist$WhattoWhat=="AtoC" | verDeNovolist$WhattoWhat=="TtoG"),]
	verDeNovoACTGcount<-nrow(verDeNovoACTGlist)
	verDeNovoACTGprop<-verDeNovoACTGcount/verDeNovocount

	verDeNovoGCCGlist<-verDeNovolist[ which(verDeNovolist$WhattoWhat=="CtoG" | verDeNovolist$WhattoWhat=="GtoC"),]
	verDeNovoGCCGcount<-nrow(verDeNovoGCCGlist)
	verDeNovoGCCGprop<-verDeNovoGCCGcount/verDeNovocount

	verDeNovoGCTAlist<-verDeNovolist[ which(verDeNovolist$WhattoWhat=="GtoT" | verDeNovolist$WhattoWhat=="CtoA"),]
	verDeNovoGCTAcount<-nrow(verDeNovoGCTAlist)
	verDeNovoGCTAprop<-verDeNovoGCTAcount/verDeNovocount
	vertypesDeNovoDF<-data.frame(Types=c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), DeNovocoralProportion=c(verDeNovoATGCprop, 	verDeNovoGCATprop, verDeNovoATTAprop,verDeNovoACTGprop, verDeNovoGCCGprop, verDeNovoGCTAprop))
	
#"A>G,T>C","G>A,C>T","A>T,T>A","A>C,T>G","G>C,C>G","G>T,C>A"

par(mfrow=c(1,3))	

	#LoH spectrum:

	vertypesLoHDFplot<-barplot(vertypesLoHDF$LoHcoralProportion,names.arg=vertypesLoHDF$Types, ylim=c(0,0.7), main="Verified LoH coral somatic mutations", las=2)
	
	#GOH spectrum:
	vertypesDeNovoDFplot<-barplot(vertypesDeNovoDF$DeNovocoralProportion,names.arg=vertypesDeNovoDF$Types, ylim=c(0,0.7), main="Verified GoH coral somatic mutations", las=2)

	
#"A>G,T>C","G>A,C>T","A>T,T>A","A>C,T>G","G>C,C>G","G>T,C>A"
#FIGURE 4:
par(mfrow=c(1,4))
	vertypesDFplot<-barplot(vertypesDF$coralProportion,names.arg=vertypesDF$Types, ylim=c(0,0.7), main="Verified coral somatic mutations", las=2) #Figure 4a
	
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
#for LoH
TranstionsLoH<- verLoHATGCprop + verLoHGCATprop
TransversionsLoH<- verLoHATTAprop+verLoHACTGprop+verLoHGCCGprop+verLoHGCTAprop
TiTvLoH<- TranstionsLoH/TransversionsLoH

#for DeNovo
TranstionsDeNovo<- verDeNovoATGCprop + verDeNovoGCATprop
TransversionsDeNovo<- verDeNovoATTAprop+verDeNovoACTGprop+verDeNovoGCCGprop+verDeNovoGCTAprop
TiTvDeNovo<- TranstionsDeNovo/TransversionsDeNovo

##make a contingency table for chi-square test##
types<-c("A>G/T>C","G>A/C>T","A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A")
coralcounts<-c(verATGCcount,verGCATcount, verATTAcount, verACTGcount, verGCCGcount, verGCTAcount)

coralcountsLoH<-c(verLoHATGCcount,verLoHGCATcount, verLoHATTAcount, verLoHACTGcount, verLoHGCCGcount, verLoHGCTAcount)

coralcountsDeNovo<-c(verDeNovoATGCcount,verDeNovoGCATcount, verDeNovoATTAcount, verDeNovoACTGcount, verDeNovoGCCGcount, verDeNovoGCTAcount)


humanesophaguscounts<-c(1200,700,3300,900,1500,450)
humangermcounts<-c(192, 321, 41, 58, 63, 72)

contingency<-data.frame("coralmuts"= coralcounts, "humangerm"= humangermcounts, "humanesophagus"=humanesophaguscounts,"humanUV"=(3500*UVtypesDF$Proportion) )

c.humangerm<-data.frame("coralmuts"=coralcounts, "humangerm"=humangermcounts)

c.humaneso<-data.frame("coralmuts"= coralcounts, "humaneso"=humanesophaguscounts)
	
c.humanUV<-data.frame("coralmuts"= coralcounts, "humanUV"=(3500*UVtypesDF$Proportion))

c.LoH<-data.frame("coralmuts"=coralcounts, "LoH"=coralcountsLoH)

c.DeNovo<-data.frame("coralmuts"=coralcounts, "DeNovo"=coralcountsDeNovo)

LoH.DeNovo<-data.frame("LoH"=coralcountsLoH, "DeNovo"=coralcountsDeNovo)

c.LoH.DeNovo<-data.frame("coralmuts"=coralcounts, "LoH"=coralcountsLoH, "DeNovo"=coralcountsDeNovo)

#info about contingency tables here: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r#data-format-contingency-tables

chi<- chisq.test(contingency)
chi2<- chisq.test(c.humangerm)
chi3<- chisq.test(c.humaneso)
chi4<-chisq.test(c.humanUV)
chi5<-chisq.test(c.LoH)
chi6<-chisq.test(c.DeNovo)
chi7<-chisq.test(LoH.DeNovo)
chi8<-chisq.test(c.LoH.DeNovo)

chi$observed
chi$expected
chi$residuals
chi$p.value
chi$estimate
chi$statistic
chi$parameter


#Let’s visualize Pearson residuals using the package corrplot:
library(corrplot)
corrplot(chi$residuals, is.cor = FALSE)

# Contibution in percentage (%)
contrib <- 100*chi$residuals^2/chi$statistic
round(contrib, 3)

# Visualize the contribution
corrplot(contrib, is.cor = FALSE)

#library(graphics)
#mosaicplot(contingency, shade = TRUE, las=2,
           main = "housetasks")
           
install.packages("vcd")
library("vcd")
# plot just a subset of the table
assoc(head(contingency), shade = TRUE, las=3)


file_path <- "http://www.sthda.com/sthda/RDoc/data/housetasks.txt"
housetasks <- read.delim(file_path, row.names = 1)
chisq.test(housetasks)


#figure 3#
data<-read.delim("~/Documents/CurrentBiologySubmission/sizeandfreq.txt")

size<-data$sizecm2 

freq<-data$totalverified.bp.sample
lohfreq<-data$LOH.bp.sample
denovofreq<-data$denovo.bp.samples

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
