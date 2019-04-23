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


###GOH VS LOH for each colony - FIGURE 3 b,c,d,e ###
	verDeNovolist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="DeNovo"),] # outputs all GoH mutations from the uniquified list of mutations
	
	verDeNovocount<-nrow(verDeNovolist) # number of all unique GoH mutations 
	verDeNovoprop<-verDeNovocount/uniqueverifiedcount
	
	verLoHlist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="LoH"),] # outputs all LoH mutations from the uniquified list of mutations

	verLoHcount<-nrow(verLoHlist) # number of all unique LoH mutations
	
	verLoHprop<-verLoHcount/uniqueverifiedcount
	
	verifiedDeNovoLoHDF<-data.frame(Types=c("GoH","LoH"),Proportions=c(verDeNovoprop*100,verLoHprop*100))
	verifiedDeNovoLoHplot<-barplot(verifiedDeNovoLoHDF$Proportion,names.arg=verifiedDeNovoLoHDF$Types, ylim=c(0,100), main=colony) #Figure 3 b, c, d, e
	verifiedDeNovoLoHDF$Proportion

	### COMPARING DEPTHS (FIGURE 2) ###
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
	depthsplots<- p +geom_boxplot() + geom_sina(aes(color=groups),size=1 ) + ylim(0,85) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Read Depth") # here is your depths plot, FIGURE 2

	
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
	barplot(mat, main=colony, ylim=c(0,100)) # here is your stacked barplot with poppolys, nonpoppolys, and LoH/GoH - Figure 4
	
	

})	

##Figure 3 A ##
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
