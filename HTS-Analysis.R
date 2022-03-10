# HTS-Analysis.R

# A script to generally analyze data from a high-throughput TREX1 screening assay

##########

rm(list=ls())

## Define these variables

path.to.file='~/Desktop/' # directory file is located in
file.name='HTS-Example-Data4' # name of file
pos.ctrls=c(1) # which plates have the positive controls (DMSO without enzyme)
neg.ctrls=c(2,3,4) # which plates have the negative controls (DMSO with enzyme)
candidate.threshold=0.5 # compounds with this level of relative activity deviation from controls will be identified as candidates; ex. 0.5 means compounds with <50% activity or >150% activity will be flagged
# script assumes that controls are in column 12 of plates, organized by plate
# script designates the 96-well reactions plates such that plates 1, 2, 3, and 4 are found in the topleft, topright, bottomleft, and bottomright staggers of the 384-well microplate, respectively
# in the graphs: purple dots are controls, inhibitors are red, activators are green, and all other compounds are grey
# Z-factor is a metric of assay precision: Z'≥0 is absolutely necessary when using a candidate threshold of 50% inhibition, Z'≥0.35 is decent, and Z'≥0.5 is excellent

#### BEGIN SCRIPT

## Import raw data

raw.data=as.matrix(read.csv(paste(path.to.file,file.name,'.csv',sep = ''),header = FALSE, sep = ','))

## Rearrange and tidy data

plates.fluorescence=array(0,dim = c(8,12,4))

plates.fluorescence[,,1]=raw.data[seq(1,nrow(raw.data),2),1:12]
plates.fluorescence[,,2]=raw.data[seq(1,nrow(raw.data),2),13:24]
plates.fluorescence[,,3]=raw.data[seq(2,nrow(raw.data),2),1:12]
plates.fluorescence[,,4]=raw.data[seq(2,nrow(raw.data),2),13:24]

## Normalize data

plates.activity=(plates.fluorescence-mean(plates.fluorescence[,12,pos.ctrls]))/(mean(plates.fluorescence[,12,neg.ctrls])-mean(plates.fluorescence[,12,pos.ctrls]))

## Calculate stats

Z.factor=1-3*(sd(plates.activity[,12,pos.ctrls])+sd(plates.activity[,12,neg.ctrls]))/abs(mean(plates.activity[,12,pos.ctrls])-mean(plates.activity[,12,neg.ctrls]))

screen.quality='N/A'
if (Z.factor<0){
  screen.quality='UNACCEPTABLE'
}
if (Z.factor>=0){
  screen.quality='POOR'
}
if (Z.factor>=0.1){
  screen.quality='ADEQUATE'
}
if (Z.factor>=0.35){
  screen.quality='GOOD'
}
if (Z.factor>=0.5){
  screen.quality='EXCELLENT'
}

num.inhibitors=sum(c(plates.activity)<=(1-candidate.threshold) & rep(rep(1:12,each=8),times=4)!=12)
num.activators=sum(c(plates.activity)>=(1+candidate.threshold) & rep(rep(1:12,each=8),times=4)!=12)

REPORT=paste('This assay was of',screen.quality,'quality, and identified',num.inhibitors,'candidate inhibitors and',num.activators,'candidate activators!',sep = ' ')

## Graph Results

y.axis=rep(c('A','B','C','D','E','F','G','H'),times=12)
y.axis.n=rep(8:1,times=12)
x.axis=rep(1:12,each=8)
colors=c(plates.activity); colors[]='grey'
colors[c(plates.activity)<=(1-candidate.threshold) & rep(x.axis,times=4)!=12]='red'
colors[c(plates.activity)>=(1+candidate.threshold) & rep(x.axis,times=4)!=12]='green'
colors[rep(x.axis,times=4)==12]='purple'
colors=array(colors,dim = c(8,12,4))

par(mfrow=c(2,2))

plot(x.axis,y.axis.n,pch=19,cex=1,col=colors[,,1],main='Plate 1',xlab = '',ylab = '',yaxt='n',xaxt='n')
axis(3,at=1:12,labels = 1:12,tick = FALSE,padj=1)
axis(2,at = 8:1,labels = c('A','B','C','D','E','F','G','H'),tick = FALSE,hadj=0,las=1)

plot(x.axis,y.axis.n,pch=19,cex=1,col=colors[,,2],main='Plate 2',xlab = '',ylab = '',yaxt='n',xaxt='n')
axis(3,at=1:12,labels = 1:12,tick = FALSE,padj=1)
axis(2,at = 8:1,labels = c('A','B','C','D','E','F','G','H'),tick = FALSE,hadj=0,las=1)

plot(x.axis,y.axis.n,pch=19,cex=1,col=colors[,,3],main='Plate 3',xlab = '',ylab = '',yaxt='n',xaxt='n')
axis(3,at=1:12,labels = 1:12,tick = FALSE,padj=1)
axis(2,at = 8:1,labels = c('A','B','C','D','E','F','G','H'),tick = FALSE,hadj=0,las=1)

plot(x.axis,y.axis.n,pch=19,cex=1,col=colors[,,4],main='Plate 4',xlab = '',ylab = '',yaxt='n',xaxt='n')
axis(3,at=1:12,labels = 1:12,tick = FALSE,padj=1)
axis(2,at = 8:1,labels = c('A','B','C','D','E','F','G','H'),tick = FALSE,hadj=0,las=1)

## Report plate activities

plate.1=as.data.frame(plates.activity[,,1])
plate.2=as.data.frame(plates.activity[,,2])
plate.3=as.data.frame(plates.activity[,,3])
plate.4=as.data.frame(plates.activity[,,4])
colnames(plate.1)=1:12
colnames(plate.2)=1:12
colnames(plate.3)=1:12
colnames(plate.4)=1:12
row.names(plate.1)=c('A','B','C','D','E','F','G','H')
row.names(plate.2)=c('A','B','C','D','E','F','G','H')
row.names(plate.3)=c('A','B','C','D','E','F','G','H')
row.names(plate.4)=c('A','B','C','D','E','F','G','H')

show('')
show('Plate 1:')
show(plate.1)
show('')
show('Plate 2:')
show(plate.2)
show('')
show('Plate 3:')
show(plate.3)
show('')
show('Plate 4:')
show(plate.4)

show('')
show('Assay Z-factor:')
show(Z.factor)

show('')
show('REPORT:')
show(REPORT)

## Save and tidy

RESULTS=list('Plate 1'=plate.1,'Plate 2'=plate.2,'Plate 3'=plate.3,'Plate 4'=plate.4,'Assay Z-factor'=Z.factor)
save(RESULTS,REPORT,file=paste(path.to.file,'RESULTS_',file.name,'.RData',sep = ''))

rm(list=setdiff(ls(),c('RESULTS','REPORT')))

#### END SCRIPT


