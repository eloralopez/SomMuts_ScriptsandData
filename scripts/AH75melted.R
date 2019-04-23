###AH75###
setwd("~/Documents/SomaticMutations/OfuAug/")
library("reshape2")
library(ggplot2)
library(stringr)
library(sciplot)
library(sinaplot)
library(ggforce)
library(gridExtra)
sessionInfo()

AH09<-read.delim("AH09datatable20181029.txt")
mdf<-melt(AH09, id.vars="chrom.pos", measure.vars=c("AH09_10A", "AH09_10B", "AH09_11A", "AH09_11B", "AH09_1A","AH09_1B","AH09_2B","AH09_2C","AH09_3A","AH09_3C","AH09_4A","AH09_4C","AH09_5A","AH09_5B","AH09_6A","AH09_6C","AH09_7A","AH09_7B","AH09_8B","AH09_8C","AH09_9A","AH09_9C"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedAH09datatable20181029.txt",sep="\t",quote=FALSE, row.name=FALSE)


AH75<-read.delim("AH75datatable20181029.txt")
mdf<-melt(AH75, id.vars="chrom.pos", measure.vars=c("AH75_1", "AH75_10", "AH75_11", "AH75_12", "AH75_13", "AH75_14", "AH75_15", "AH75_16", "AH75_17", "AH75_2", "AH75_3","AH75_4", "AH75_5", "AH75_6", "AH75_7", "AH75_8", "AH75_9"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedAH75datatable20181029.txt",sep="\t",quote=FALSE, row.name=FALSE)

AH06<-read.delim("AH06datatable20181029.txt")
mdf<-melt(AH06, id.vars="chrom.pos", measure.vars=c("AH06_1", "AH06_10", "AH06_11", "AH06_12", "AH06_13", "AH06_14", "AH06_15", "AH06_16", "AH06_17", "AH06_2", "AH06_3","AH06_4", "AH06_5", "AH06_6", "AH06_7", "AH06_8", "AH06_9"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedAH06datatable20181029.txt",sep="\t",quote=FALSE, row.name=FALSE)

AH88<-read.delim("AH88datatable20181029.txt")
mdf<-melt(AH88, id.vars="chrom.pos", measure.vars=c("AH88_1", "AH88_10", "AH88_11", "AH88_12", "AH88_13", "AH88_14", "AH88_15", "AH88_16", "AH88_17", "AH88_2", "AH88_3","AH88_4", "AH88_5", "AH88_6", "AH88_7", "AH88_8", "AH88_9"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedAH88datatable20181029.txt",sep="\t",quote=FALSE, row.name=FALSE)


####ready to analzyze###

AH75_meta<-read.delim("/Users/eloralopez/Documents/SomaticMutations/OfuAug/WithVerifications/WithDepths/melted_AH06_datatable20181029_withtransitionandpopinfo_withverifications.txt")
AH75_meta<-read.delim("/Users/eloralopez/Documents/SomaticMutations/OfuAug/WithVerifications/WithDepths/melted_AH88_datatable20181029_withtransitionandpopinfo_withverifications.txt")
genoanddepth<-(AH75_meta$genotype)
split<-str_split_fixed(genoanddepth, ",", 4)
genotypes<-split[,1]
totaldepth<-as.numeric(split[,2])
refdepth<-as.numeric(split[,3])
altdepth<-as.numeric(split[,4])
AH75df<-data.frame("chrom.pos" = AH75_meta$chrom.pos, "sample"= AH75_meta$sample, "ref" = AH75_meta$ref, "alt" = AH75_meta$alt, "genotype"= genotypes, "totaldepth"=totaldepth, "refdepth"=refdepth, "altdepth"=altdepth, "GoH_or_LoH"=AH75_meta$DeNovo_LoH, "Ti/Tv"=AH75_meta$TiTv, "WhattoWhat" = AH75_meta$WhattoWhat, "PopulationPoly"= AH75_meta$PopulationPoly, "VerifiedYesorNo"=AH75_meta$VerifiedYesorNo)

DepthMeansdf<-aggregate(totaldepth~chrom.pos, AH75df, FUN=mean)


AH75df<-merge(AH75df, DepthMeansdf[, c("chrom.pos", "totaldepth")], by="chrom.pos")
uniqueAH75df<-AH75df[match(unique(AH75df$chrom.pos), AH75df$chrom.pos),]
	poppolys<-uniqueAH75df[ which(uniqueAH75df$PopulationPoly 			=="TRUE"),]
	poppolyscount<-nrow(poppolys)
	nonpoppolys<-uniqueAH75df[ which(uniqueAH75df$PopulationPoly 			=="FALSE"),]
	nonpoppolyscount<-nrow(nonpoppolys)


####VERIFIED####
verifiedlist<-AH75df[ which(AH75df$VerifiedYesorNo =="yes"),]
verifiedpopulationpolys<-verifiedlist[ which(verifiedlist$PopulationPoly=="TRUE"),]
uniqueverifiedlist<-verifiedlist[match(unique(verifiedlist$chrom.pos), verifiedlist$chrom.pos),]

uniqueverifiedcount<-nrow(uniqueverifiedlist)

uniqueverifiednonpop<-uniqueverifiedlist[ which(uniqueverifiedlist$PopulationPoly=="FALSE"),]
uniqueverifiednonpopcount<-nrow(uniqueverifiednonpop)
uniqueverifiedpopulationpolys<-uniqueverifiedlist[ which(uniqueverifiedlist$PopulationPoly=="TRUE"),]
uniqueverifiedpopulationpolyscount<-nrow(uniqueverifiedpopulationpolys)

uniqueverifiedpopulationpolysGoH<-uniqueverifiedpopulationpolys[which(uniqueverifiedpopulationpolys$GoH_or_LoH=="DeNovo"),]

uniqueverifiedpopulationpolysGoHcount<-nrow(uniqueverifiedpopulationpolysGoH)

uniqueverifiedpopulationpolysLoH<-uniqueverifiedpopulationpolys[which(uniqueverifiedpopulationpolys$GoH_or_LoH=="LoH"),]

uniqueverifiedpopulationpolysLoHcount<-nrow(uniqueverifiedpopulationpolysLoH)

uniqueverifiednonpopGoH<-uniqueverifiednonpop[which(uniqueverifiednonpop$GoH_or_LoH=="DeNovo"),]

uniqueverifiednonpopGoHcount<-nrow(uniqueverifiednonpopGoH)

uniqueverifiednonpopLoH<-uniqueverifiednonpop[which(uniqueverifiednonpop$GoH_or_LoH=="LoH"),]

uniqueverifiednonpopLoHcount<-nrow(uniqueverifiednonpopLoH)

	names<-c("GoH","LoH")
	proportions<-c(uniqueverifiedpopulationpolysGoHcount*100/uniqueverifiedpopulationpolyscount, uniqueverifiedpopulationpolysLoHcount*100/uniqueverifiedpopulationpolyscount)
	barplot(proportions, names.arg=names, main="AH75", ylim=c(0,100))
	


mat<-matrix(c(uniqueverifiedpopulationpolysGoHcount*100/uniqueverifiedcount, uniqueverifiednonpopGoHcount*100/uniqueverifiedcount, uniqueverifiedpopulationpolysLoHcount*100/uniqueverifiedcount, uniqueverifiednonpopLoHcount*100/uniqueverifiedcount),ncol=2,byrow=TRUE)

colnames(mat)<-c("PopPoly", "NonPopPoly")
rownames(mat)<-c("% GOH","% LOH")
mat<-as.table(mat)
barplot(mat)

verifieddepth<-(uniqueverifiedlist$totaldepth.y)
averageverdepth<-mean(verifieddepth) #this gives the mean value of average depth at each site, not just the mutant depth at each site
SEaverageverdepth<-se(verifieddepth) #finds standard deviation of verifieddepths

###GOH VS LOH###
	verDeNovolist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="DeNovo"),]
	verDeNovocount<-nrow(verDeNovolist)
	verDeNovoprop<-verDeNovocount/verifiedcount
	verLoHlist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="LoH"),]
	verLoHcount<-nrow(verLoHlist)
		verLoHprop<-verLoHcount/verifiedcount

	names<-c("GoH","LoH")
	groups<-c(verDeNovocount*100/uniqueverifiedcount, verLoHcount*100/uniqueverifiedcount)
	barplot(groups, names.arg=names, ylim=c(0,70), main="AH09")
	
	

####FALSIFIED####
falsifiedlist<-AH75df[ which(AH75df$VerifiedYesorNo =="no"),]
uniquefalsifiedlist<-falsifiedlist[match(unique(falsifiedlist$chrom.pos), falsifiedlist$chrom.pos),]
uniquefalsifiedcount<-nrow(uniquefalsifiedlist)

uniquefalsifiednonpop<-uniquefalsifiedlist[ which(uniquefalsifiedlist$PopulationPoly=="FALSE"),]
uniquefalsifiednonpopcount<-nrow(uniquefalsifiednonpop)
uniquefalsifiedpopulationpolys<-uniquefalsifiedlist[ which(uniquefalsifiedlist$PopulationPoly=="TRUE"),]
uniquefalsifiedpopulationpolyscount<-nrow(uniquefalsifiedpopulationpolys)


falsifieddepth<-(uniquefalsifiedlist$totaldepth.y)
averagefalsdepth<-mean(falsifieddepth)
SEaveragefalsdepth<-se(falsifieddepth)

###PLOT DEPTHS###
names<-c("Average Depth for Verified","Average Depth for Falsified")
means<-c(averageverdepth, averagefalsdepth)
SEs<-c(SEaverageverdepth, SEaveragefalsdepth)
plotTop<-max(means+SEs*2)
barCenters<-barplot(means,names.arg=names, las=1, ylim=c(0,50))
title("AH75")
arrows(barCenters, means-SEs*2, barCenters, means+SEs*2, lwd=1, angle=90, code=3)

allsites<-AH75df[match(unique(AH75df$chrom.pos), AH75df$chrom.pos),]

x<-c(verifieddepth, falsifieddepth, allsites$totaldepth.y)
groups<-c(rep("Verified",uniqueverifiedcount),rep("Falsified",uniquefalsifiedcount),rep("AllSites",nrow(allsites)))
df<-data.frame(groups, x)

#sinaplot(x,groups, col=2:4) #make a sinaplot with different colors for each class

p<- ggplot(df, aes(groups,x))
#p + geom_violin(aes(fill = groups))
p +geom_boxplot() + geom_sina(aes(color=groups),size=1 ) + ylim(0,85) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Read Depth") + ggtitle("test")

	###PLOT POPPOLY vs NON-POPPOLY
	names<-c("Verified PopPoly","Falsified PopPoly","Verified New Variant","Falsified New Variant")
	values<-c(uniqueverifiedpopulationpolyscount, uniquefalsifiedpopulationpolyscount, uniqueverifiednonpopcount, uniquefalsifiednonpopcount)
	barplot(values,names.arg=names)



###here's a for loop that works (until you get to the sinaplot):###

files<-list.files(path="~/Documents/SomaticMutations/OfuAug/WithVerifications/WithDepths", pattern="*.txt", full.names=T, recursive=FALSE)

par(mfrow=c(1,4)) #if your output includes plots, this will put all of the plots in one image. change the (4,5) to different numbers depending on how large you want your grid of plots to be.

#then, start your loop:
lapply(files,function(x) {
	metadata<-read.delim(x) #load a file
	#then paste in whatever functions you would normally use.
	base<-basename(x)
	colony<-strsplit(base, "\\_")[[1]][2]
	
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
print(colony)
print("LOHporop is:")
print(verLoHprop)
print("GOHprop is:")
print(verDeNovoprop)
})
####FALSIFIED####
	falsifiedlist<-AH75df[ which(AH75df$VerifiedYesorNo =="no"),]
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
	barplot(mat, main=colony, ylim=c(0,100)) # here is your stacked barplot with poppolys, nonpoppolys, and LoH/GoH
	
	

})	
	###PLOT POPPOLY vs NON-POPPOLY
	names<-c("% of PopPolys verified","% of New Variants verified")
	values<-c(uniqueverifiedpopulationpolyscount*100/poppolyscount,  uniqueverifiednonpopcount*100/nonpoppolyscount)
	par(mar=c(15, 4.1, 4.1, 2.1))
	#barplot(values,names.arg=names, las=2, ylim=c(0,30), main=colony)
	names<-c("% of verified mutations that are poppoly","% of verified mutations that are new")
	groups<-c(uniqueverifiedpopulationpolyscount*100/uniqueverifiedcount, uniqueverifiednonpopcount*100/uniqueverifiedcount)
	par(mar=c(15, 4.1, 4.1, 2.1))
	#barplot(groups,names.arg=names, las=2, ylim=c(0,100), main=colony)

	verDeNovolist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="DeNovo"),]
	verDeNovocount<-nrow(verDeNovolist)

	verLoHlist<-uniqueverifiedlist[ which(uniqueverifiedlist$GoH_or_LoH =="LoH"),]
	verLoHcount<-nrow(verLoHlist)
	
	names<-c("GoH","LoH")
	groups<-c(verDeNovocount*100/uniqueverifiedcount, verLoHcount*100/uniqueverifiedcount)
	#barplot(groups, names.arg=names, ylim=c(0,100), main=colony)
	
})
grid.arrange(plots)

