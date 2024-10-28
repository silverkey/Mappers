all = read.delim(file='PR_RESULTS_FEATURE_REPEAT_ANNOTATED.xls')
sel.dir = '/Users/remo/ANALYSIS/RIKEN/FINAL_SELECTED_DE/MAPPING_ANALYSIS/SECOND/'
data = 'RESULTS_FEATURE_REPEAT_ANNOTATED.xls'

source.tab = read.delim(file=paste(sel.dir,data,sep=''))

pr = source.tab[grep('pr',source.tab$compname),]
write.table(unique(pr),file='PR_RES_ANNOTATED.xls',sep="\t",row.names=F,quote=F)

pr.gene = table(pr$f.name)[c('promoter','utr5','cds','exon','intron','utr3')]
pr.gene.sense = table(pr[pr$f.dir=='sense',]$f.name)[c('promoter','utr5','cds','exon','intron','utr3')]
pr.gene.anti = table(pr[pr$f.dir=='antisense',]$f.name)[c('promoter','utr5','cds','exon','intron','utr3')]

all.gene = table(all$f.name)[c('promoter','utr5','cds','exon','intron','utr3')]
all.gene.sense = table(all[all$f.dir=='sense',]$f.name)[c('promoter','utr5','cds','exon','intron','utr3')]
all.gene.anti = table(all[all$f.dir=='antisense',]$f.name)[c('promoter','utr5','cds','exon','intron','utr3')]

all.inter = table(all$f.name)['intergenic']
all.genic = sum(all.gene)
names(all.genic) = 'genic'

pr.inter = table(pr$f.name)['intergenic']
pr.genic = sum(pr.gene)
names(pr.genic) = 'genic'

pdf(file='PR_genemapping.pdf',width=18,height=10)

par(mfrow=c(4,2),las=1,mar=c(5,8,5,5))

barplot((c(all.inter,all.genic)/(all.inter+all.genic)*100),
main=paste('Positions of All Clusters (n = ',all.genic+all.inter,')',sep=''),
horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,90),xlab='percentage')

barplot((c(pr.inter,pr.genic)/(pr.inter+pr.genic)*100),
main=paste('Positions of PR DE Clusters (n = ',pr.genic+pr.inter,')',sep=''),
horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,90),xlab='percentage')

barplot(rev(all.gene/sum(all.gene)*100),
        main=paste('Positions of All Clusters (n = ',sum(all.gene),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,50),xlab='percentage')

barplot(rev(pr.gene/sum(pr.gene)*100),
        main=paste('Positions of PR DE Clusters (n = ',sum(pr.gene),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,50),xlab='percentage')

barplot(rev(all.gene.sense/sum(all.gene.sense)*100),
        main=paste('Positions of All Clusters on Sense Strand (n = ',sum(all.gene.sense),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,50),xlab='percentage')

barplot(rev(pr.gene.sense/sum(pr.gene.sense)*100),
        main=paste('Positions of PR DE Clusters on Sense Strand (n = ',sum(pr.gene.sense),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,50),xlab='percentage')

barplot(rev(all.gene.anti/sum(all.gene.anti)*100),
        main=paste('Positions of All Clusters on Antiense Strand (n = ',sum(all.gene.anti),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,50),xlab='percentage')

barplot(rev(pr.gene.anti/sum(pr.gene.anti)*100),
        main=paste('Positions of PR DE Clusters on Antisense Strand (n = ',sum(pr.gene.anti),')',sep=''),
        horiz=T,cex.names=1.5,cex.axis=1.5,cex.main=1.5,xlim=c(0,50),xlab='percentage')
 
dev.off()

c=cbind(pr.gene.sense,pr.gene.anti,sum(pr.gene.sense),sum(pr.gene.anti))
test = function(x) { prop.test(c(x[1],x[2]),c(x[3],x[4]),alternative='l')$p.value }
p = apply(c,1,test)
p.ad = p.adjust(p)
p.ad
