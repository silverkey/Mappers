all = read.delim(file='PR_RESULTS_FEATURE_REPEAT_ANNOTATED.xls')
sel.dir = '/Users/remo/ANALYSIS/RIKEN/FINAL_SELECTED_DE/MAPPING_ANALYSIS/SECOND/'
data = 'RESULTS_FEATURE_REPEAT_ANNOTATED.xls'
source.tab = read.delim(file=paste(sel.dir,data,sep=''))
pr = source.tab[grep('pr',source.tab$compname),]

pdf(file='PR_repeatclass_mapping.pdf',width=18,height=10)

pr.rep = table(pr$r.class)
all.rep = table(all$r.class)

pr.rep.sense = table(pr[pr$r.dir=='sense',]$r.class)
all.rep.sense = table(all[all$r.dir=='sense',]$r.class)
pr.rep.anti = table(pr[pr$r.dir=='antisense',]$r.class)
all.rep.anti = table(all[all$r.dir=='antisense',]$r.class)

par(mfrow=c(3,2),las=1,mar=c(5,12,5,5))

barplot(all.rep[names(all.rep) %in% names(pr.rep[pr.rep>14])]/sum(all.rep)*100,
        main=paste('Positions of All Clusters Overlapping Repeats (n = ',sum(all.rep),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,40),xlab='percentage')

barplot(pr.rep[pr.rep>14]/sum(pr.rep)*100,
        main=paste('Positions of PR DE Clusters Overlapping Repeats (n = ',sum(pr.rep),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,40),xlab='percentage')

barplot(all.rep.sense[names(all.rep.sense) %in% names(pr.rep.sense[pr.rep.sense>14])]/sum(all.rep.sense)*100,
        main=paste('Positions of All Clusters on Sense Strand Overlapping Repeats (n = ',sum(all.rep.sense),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,40),xlab='percentage')

barplot(pr.rep.sense[pr.rep.sense>14]/sum(pr.rep.sense)*100,
        main=paste('Positions of PR DE Clusters on Sense Strand Overlapping Repeats (n = ',sum(pr.rep.sense),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,40),xlab='percentage')

barplot(all.rep.anti[names(all.rep.anti) %in% names(pr.rep.anti[pr.rep.anti>14])]/sum(all.rep.anti)*100,
        main=paste('Positions of All Clusters on Antisense Strand Overlapping Repeats (n = ',sum(all.rep.anti),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,40),xlab='percentage')

barplot(pr.rep.anti[pr.rep.anti>14]/sum(pr.rep.anti)*100,
        main=paste('Positions of PR DE Clusters on Antisense Strand Overlapping Repeats (n = ',sum(pr.rep.anti),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,40),xlab='percentage')

dev.off()

c=cbind(pr.rep.sense,pr.rep.anti,sum(pr.rep.sense),sum(pr.rep.anti))
test = function(x) { prop.test(c(x[1],x[2]),c(x[3],x[4]),alternative='l')$p.value }
p = apply(c,1,test)
p.ad = p.adjust(p)
p.ad


pdf(file='PR_repeatfamily_mapping.pdf',width=18,height=10)

pr.rep = table(pr$r.family)
all.rep = table(all$r.family)

pr.rep.sense = table(pr[pr$r.dir=='sense',]$r.family)
all.rep.sense = table(all[all$r.dir=='sense',]$r.family)
pr.rep.anti = table(pr[pr$r.dir=='antisense',]$r.family)
all.rep.anti = table(all[all$r.dir=='antisense',]$r.family)

par(mfrow=c(3,2),las=1,mar=c(5,12,5,5))

barplot(all.rep[names(all.rep) %in% names(pr.rep[pr.rep>14])]/sum(all.rep)*100,
main=paste('Positions of All Clusters Overlapping Repeats (n = ',sum(all.rep),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(pr.rep[pr.rep>14]/sum(pr.rep)*100,
main=paste('Positions of PR DE Clusters Overlapping Repeats (n = ',sum(pr.rep),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(all.rep.sense[names(all.rep.sense) %in% names(pr.rep.sense[pr.rep.sense>14])]/sum(all.rep.sense)*100,
main=paste('Positions of All Clusters on Sense Strand Overlapping Repeats (n = ',sum(all.rep.sense),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(pr.rep.sense[pr.rep.sense>14]/sum(pr.rep.sense)*100,
main=paste('Positions of PR DE Clusters on Sense Strand Overlapping Repeats (n = ',sum(pr.rep.sense),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(all.rep.anti[names(all.rep.anti) %in% names(pr.rep.anti[pr.rep.anti>14])]/sum(all.rep.anti)*100,
main=paste('Positions of All Clusters on Antisense Strand Overlapping Repeats (n = ',sum(all.rep.anti),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(pr.rep.anti[pr.rep.anti>14]/sum(pr.rep.anti)*100,
main=paste('Positions of PR DE Clusters on Antisense Strand Overlapping Repeats (n = ',sum(pr.rep.anti),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

dev.off()

c=cbind(pr.rep.sense,pr.rep.anti,sum(pr.rep.sense),sum(pr.rep.anti))
test = function(x) { prop.test(c(x[1],x[2]),c(x[3],x[4]),alternative='l')$p.value }
p = apply(c,1,test)
p.ad = p.adjust(p)
p.ad[p.ad<=0.05 && !is.na(p.ad)]







pdf(file='PR_repeatname_mapping.pdf',width=18,height=10)

pr.rep = table(pr$r.name)
all.rep = table(all$r.name)

pr.rep.sense = table(pr[pr$r.dir=='sense',]$r.name)
all.rep.sense = table(all[all$r.dir=='sense',]$r.name)
pr.rep.anti = table(pr[pr$r.dir=='antisense',]$r.name)
all.rep.anti = table(all[all$r.dir=='antisense',]$r.name)

par(mfrow=c(3,2),las=1,mar=c(5,12,5,5))

barplot(all.rep[names(all.rep) %in% names(pr.rep[pr.rep>14])]/sum(all.rep)*100,
main=paste('Positions of All Clusters Overlapping Repeats (n = ',sum(all.rep),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(pr.rep[pr.rep>14]/sum(pr.rep)*100,
main=paste('Positions of PR DE Clusters Overlapping Repeats (n = ',sum(pr.rep),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(all.rep.sense[names(all.rep.sense) %in% names(pr.rep.sense[pr.rep.sense>14])]/sum(all.rep.sense)*100,
main=paste('Positions of All Clusters on Sense Strand Overlapping Repeats (n = ',sum(all.rep.sense),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(pr.rep.sense[pr.rep.sense>14]/sum(pr.rep.sense)*100,
main=paste('Positions of PR DE Clusters on Sense Strand Overlapping Repeats (n = ',sum(pr.rep.sense),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(all.rep.anti[names(all.rep.anti) %in% names(pr.rep.anti[pr.rep.anti>14])]/sum(all.rep.anti)*100,
main=paste('Positions of All Clusters on Antisense Strand Overlapping Repeats (n = ',sum(all.rep.anti),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

barplot(pr.rep.anti[pr.rep.anti>14]/sum(pr.rep.anti)*100,
main=paste('Positions of PR DE Clusters on Antisense Strand Overlapping Repeats (n = ',sum(pr.rep.anti),')',sep=''),
horiz=T,xlim=c(0,20),xlab='percentage')

dev.off()

c=cbind(pr.rep.sense,pr.rep.anti,sum(pr.rep.sense),sum(pr.rep.anti))
test = function(x) { prop.test(c(x[1],x[2]),c(x[3],x[4]),alternative='l')$p.value }
p = apply(c,1,test)
p.ad = p.adjust(p)
p.ad[p.ad<=0.05 && !is.na(p.ad)]
