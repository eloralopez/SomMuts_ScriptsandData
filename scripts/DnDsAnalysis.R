data  = NULL

setwd("~/Documents/GitHub/SomMuts/ScriptsandData/datafiles/dNdS/DnDsAnalysis.R-inputfiles/")
dnds<-read.delim("AH88_DnDs.txt")
#attach(dnds)
data <- rbind(data, data.frame(dnds))
##DN/DS CALCULATION WITH SNAP: https://www.hiv.lanl.gov/content/sequence/SNAP/SNAP.html


### it doesn't make sense to try to do syn/nonsyn for things not in the ORF. don't do this!

inORF<-dnds[ which(dnds$In_ORF == "yes"),]
	inORFcount<-nrow(inORF)
	
	ORFSyn<- inORF[ which(inORF$Syn_or_NonSyn=="Syn"),]
	countORFSyn<- nrow(ORFSyn)

	ORFNonSyn<- inORF[ which(inORF$Syn_or_NonSyn=="NonSyn"),]
	countORFNonSyn<- nrow(ORFNonSyn)

	DnDs_InORF<- countORFNonSyn/countORFSyn

	totalSyn<-dnds[ which(dnds$Syn_or_NonSyn == "Syn"),]
	counttotalSyn<- nrow(totalSyn)
	
	totalNonSyn<-dnds[ which(dnds$Syn_or_NonSyn == "NonSyn"),]
	counttotalNonSyn<- nrow(totalNonSyn)
	
	DnDs_Total<-counttotalNonSyn/counttotalSyn
	
	allrows<- nrow(dnds)
	
	100*countORFNonSyn/inORFcount
	100*counttotalNonSyn/allrows

notORF<-dnds[ which(dnds$In_ORF =="no"),]
notORFcount<-nrow(notORF)

###here's a for loop that works:###

files<-list.files(path="~/Documents/GitHub/SomMuts/ScriptsandData/datafiles/dNdS/DnDsAnalysis.R-inputfiles/", pattern="*DnDs.txt", full.names=T, recursive=FALSE)

par(mfrow=c(1,4)) #if your output includes plots, this will put all of the plots in one image. change the (4,5) to different numbers depending on how large you want your grid of plots to be.

data<-data.frame()

#then, start your loop:
lapply(files,function(x) {
	dnds<-read.delim(x) #load a file
	
	#then paste in whatever functions you would normally use:

	base<-basename(x)
	colony<-strsplit(base, "\\_")[[1]][1]

	#attach(dnds)

	inORF<-dnds[ which(dnds$In_ORF == "yes"),]
	inORFcount<-nrow(inORF)
	
	ORFSyn<- inORF[ which(inORF$Syn_or_NonSyn=="Syn"),]
	countORFSyn<- nrow(ORFSyn)

	ORFNonSyn<- inORF[ which(inORF$Syn_or_NonSyn=="NonSyn"),]
	countORFNonSyn<- nrow(ORFNonSyn)
	
	notORF<-dnds[ which(dnds$In_ORF == "no"),]
	notORFcount<-nrow(notORF)
	
	notORFSyn<- notORF[ which(notORF$Syn_or_NonSyn=="Syn"),]
	countnotORFSyn<- nrow(notORFSyn)

	notORFNonSyn<- notORF[ which(notORF$Syn_or_NonSyn=="NonSyn"),]
	countnotORFNonSyn<- nrow(notORFNonSyn)

	DnDs_InORF<- countORFNonSyn/countORFSyn

	totalSyn<-dnds[ which(dnds$Syn_or_NonSyn == "Syn"),]
	counttotalSyn<- nrow(totalSyn)
	
	totalNonSyn<-dnds[ which(dnds$Syn_or_NonSyn == "NonSyn"),]
	counttotalNonSyn<- nrow(totalNonSyn)
	
	DnDs_Total<-counttotalNonSyn/counttotalSyn
	
	notORF<-dnds[ which(dnds$In_ORF =="no"),]
	notORFcount<-nrow(notORF)

	allrows<- nrow(dnds)
	
	df_ORFornot<- data.frame(Types=c("In ORF","Not in ORF"), Proportion=c(inORFcount*100/allrows, notORFcount*100/allrows))
	#barplot(df_ORFornot$Proportion, names.arg=df_ORFornot$Types, main=colony, ylim=c(0,100), ylab="Percent of verified mutations")

	df<-data.frame(Types=c("In ORF","Total"), Proportion=c(DnDs_InORF, DnDs_Total))

	#barplot(df$Proportion, names.arg=df$Types, ylim=c(0,3),main=colony, ylab="Dn/Ds")
	
	df_percents<-data.frame(Types=c("in ORF","Total"), Proportion=c(100*countORFNonSyn/inORFcount, 100*counttotalNonSyn/allrows))
	
	#barplot(df_percents$Proportion, names.arg=df_percents$Types,main=colony,ylab="% of mutations that are Non-Synonymous",ylim=c(0,100))
	
	df_justORF<-data.frame(Types=c("Synonymous","Nonsynonymous"), Proportion=c(100*countORFSyn/inORFcount, 100*countORFNonSyn/inORFcount))
	
	barplot(df_justORF$Proportion, names.arg = df_justORF$Types, main=colony, ylab = "% of mutations in ORF",ylim=c(0,100))
})	
	print(colony)
	print("The total number of mutations in the ORF is:") 
	print(inORFcount)
	print("The total number of nonsynonymous mutations in the ORF is:")
	print(countORFNonSyn)
	print("The total number of synonymous mutations in the ORF is")
	print(countORFSyn)
	

})	

###here's another type of for loop:
data= NULL
files<-list.files(path="~/Documents/GitHub/SomMuts/ScriptsandData/datafiles/dNdS/DnDsAnalysis.R-inputfiles", pattern="*DnDs.txt", full.names=T, recursive=FALSE)
for (i in 1:length(files)) {
	file =files[i]
	base<-basename(file)
	colony<-strsplit(base, "\\_")[[1]][1]
	dnds<-read.delim(file)
	dnds<-data.frame(dnds)
	len<-nrow(dnds)
	colonyrep<-rep(colony, len)
	dndswcolony<-data.frame(dnds, colonyrep)
	data <- rbind(data, dndswcolony)
}

coding<-data[ which(data$In_ORF=="yes"),]
LOH<-coding[ which(coding$DeNovo_or_LoH=="LoH"),]
GOH<-coding[ which(coding$DeNovo_or_LoH=="DeNovo"),]

LOH.Length.dn<-LOH$Length.dn
LOH.Length.ds<-LOH$Length.ds
dNdSLOH<-mean(LOH.Length.dn)/mean(LOH.Length.ds)

GOH.Length.dn<-GOH$Length.dn
GOH.Length.ds<-GOH$Length.ds
dNdSGOH<-mean(GOH.Length.dn)/mean(GOH.Length.ds)

totalLength.dn<-coding$Length.dn
totalLength.ds<-coding$Length.ds
dNdStotal<-mean(totalLength.dn)/mean(totalLength.ds)

AH06<-subset(coding, colonyrep=="AH06")
AH06dNdS<-mean(AH06$Length.dn)/mean(AH06$Length.ds)

AH09<-subset(coding, colonyrep=="AH09")
AH09dNdS<-mean(AH09$Length.dn)/mean(AH09$Length.ds)

AH75<-subset(coding, colonyrep=="AH75")
AH75dNdS<-mean(AH75$Length.dn)/mean(AH75$Length.ds)

AH88<-subset(coding, colonyrep=="AH88")
AH88dNdS<-mean(AH88$Length.dn)/mean(AH88$Length.ds)

dNdS<- data.frame(Types=c("Total", "LOH", "GOH"), Proportion=c(dNdStotal, dNdSLOH, dNdSGOH))

barplot(dNdS$Proportion, names.arg = dNdS$Types, ylim = c(0,1.2))

