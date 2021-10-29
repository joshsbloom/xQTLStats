library(xQTLStats)
data(crosses.to.parents)

gmaps=readRDS(system.file('reference', 'yeast_gmaps.RDS', package='xQTLStats'))
#reference vcf
ref.vcf=system.file('reference', 'parents_w_svar_sorted.vcf.gz', package='xQTLStats')

#assuming input vcf with haploid parents and haploid specific variant calls 
# [1] "273614xa" "BYa"      "CBS2888a" "CLIB219x" "CLIB413a" "I14a"    
# [7] "M22"      "PW5a"     "RMx"      "Y10x"     "YJM145x"  "YJM454a" 
#[13] "YJM978x"  "YJM981x"  "YPS1009x" "YPS163a"

#yeast genome chromosome names sorted
#
#p1.name='YPS163a' #'BYa'
#p2.name='YJM145x' #'RMx'
sample.size=5e4 #2.5e5 #.5e5 #2.5e5 #2.5e5

cross='B'
#for(cross in names(crosses.to.parents)){
    p1.name=crosses.to.parents[[cross]][1]
    p2.name=crosses.to.parents[[cross]][2]
   
vcf.cross=getCrossVCF(ref.vcf,p1.name, p2.name)
    #  saveRDS(vcf.cross, file = paste0('/data/xQTL/', cross, '_vcf.RDS')) 

    #simulation bit ----------------------------------------------------------
    #    geno.matrix=simHaploidSegsFN(vcf.cross, gmaps[[cross]], sample.size, ngenerations=24)
geno.matrix=simHaploidSegsFN(vcf.cross, gmaps[[cross]], sample.size, ngenerations=2)

    #saveRDS(geno.matrix, file = paste0('/data/xQTL/', cross, '_250K_24.RDS')) 
#}
    #'/data/xQTL/B_48.RDS')
y=simPhenotypesForRefPop(geno.matrix, h2=.4, nQTL=40)

xqtl_reps=list()

for(i in 1:3){
    tf=tempfile()
    expected=simXQTLExperiment(y, geno.matrix, vcf.cross, sel.low=0.1, sel.high=0.1,depth.low=100, depth.high=100, vcf.out=tf)
    low.tail.expected=expected$low.tail.expected; high.tail.expected=expected$high.tail.expected

    #------------------------------------------------------------------------
    experiment.vcf=vcfR::read.vcfR(tf) #'/data/xQTL/sim.vcf.gz')

    #experiment names (should match col in vcf)
    low.tail.name='low.tail.sim'
    high.tail.name='high.tail.sim'

    low.tail=getBiallelicCounts(experiment.vcf, low.tail.name)
    high.tail=getBiallelicCounts(experiment.vcf, high.tail.name)

    low.tail$expected=low.tail.expected[low.tail$ID]
    high.tail$expected=high.tail.expected[high.tail$ID]


    #get simulated high tail
    sel.low  = 0.1
    sel.high = 0.1

    #phase everything 
    low.tail  = phaseCounts(vcf.cross, p1.name, low.tail)
    high.tail = phaseCounts(vcf.cross, p1.name, high.tail)


    #----------------------------------------------------------------------------
    low.tail  = calcAFD(low.tail, experiment.name=paste0('low',i),sample.size=sample.size, sel.strength=sel.low)
    high.tail = calcAFD(high.tail, experiment.name=paste0('high',i),sample.size=sample.size, sel.strength=sel.high)

    xqtl_reps[[as.character(i)]]$low.tail=low.tail
    xqtl_reps[[as.character(i)]]$high.tail=high.tail
}


results=list(xqtl_reps[[1]]$low.tail, xqtl_reps[[2]]$low.tail, xqtl_reps[[3]]$low.tail)
results.low=calcMetaAFD(results, meta_name='meta_low')
results=list(xqtl_reps[[1]]$high.tail, xqtl_reps[[2]]$high.tail, xqtl_reps[[3]]$high.tail)
results.high=calcMetaAFD(results, meta_name='meta_high')

results=calcContrastStats(list(results.high, results.low), L='meta_high', R='meta_low')


simulatedQTL=data.frame(chrom=results$chrom[abs(attr(y,"add.qtl"))], 
                    physical.position=results$physical.position[abs(attr(y,"add.qtl"))])





h1=plotIndividualExperiment(results, 'high1',simulatedQTL)
l1=plotIndividualExperiment(results, 'low1',simulatedQTL)
c1=plotContrast(results, 'high1', 'low1',simulatedQTL)

ggpubr::ggarrange(h1, l1, c1, nrow=3) 

x11()
plotSummary(results, simulatedQTL)




#library(ggplot2)
#high.plot= 
#    ggplot(results.data, aes(x=physical.position,y=p1.high/(p1.high+p2.high)))+
#    geom_point(size=0.3,alpha=0.6, color='gray21')+
#    geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')+
#    facet_grid(~chrom, scales='free_x', space='free_x')+
#    geom_ribbon(aes(ymin=afd.high-1.96*afd.high.se, ymax=afd.high+1.96*afd.high.se, fill='grey'), 
#                alpha=0.7,linetype='dashed', color='grey')+
#    geom_line(aes(x=physical.position, y=afd.high),color='red', size=2, alpha=1)+ 
#    geom_line(aes(x=physical.position, y=expected.af.high),color='black', size=.5)+
#    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none')+ggtitle('high tail')
#
#low.plot=
#    ggplot(results.data, aes(x=physical.position,y=p1.low/(p1.low+p2.low)))+
#    geom_point(size=0.3,alpha=0.6, color='gray21')+
#    geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')+
#    facet_grid(~chrom, scales='free_x', space='free_x')+
#    geom_ribbon(aes(ymin=afd.low-1.96*afd.low.se, ymax=afd.low+1.96*afd.low.se, fill='grey'), 
#                alpha=0.7,linetype='dashed', color='grey')+
#    geom_line(aes(x=physical.position, y=afd.low),color='red', size=2, alpha=1)+ 
#    geom_line(aes(x=physical.position, y=expected.af.low),color='black', size=.5)+
#    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none')+ggtitle('low tail')
#
##ggpubr::ggarrange(high.plot, low.plot, nrow=2) 
#
#contrast.plot=ggplot(results.data, aes(x=physical.position,y=(p1.high/(p1.high+p2.high))-(p1.low/(p1.low+p2.low))))+
#    geom_point(size=0.3,alpha=0.6, color='gray21')+
#    geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')+
#    facet_grid(~chrom, scales='free_x', space='free_x')+
#    geom_ribbon(aes(ymin=afd.contrast-1.96*afd.contrast.se, ymax=afd.contrast+1.96*afd.contrast.se, fill='grey'), 
#                     alpha=0.7,linetype='dashed', color='grey')+
#    geom_line(aes(x=physical.position, y=afd.contrast),color='red', size=2, alpha=1)+ 
#    geom_line(aes(x=physical.position, y=expected.af.high-expected.af.low),color='black', size=.5)+
#    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1), legend.position='none')+
#    ggtitle('allele frequence contrast')
#
##ggplot(results.data, aes(x=physical.position, y=afd.contrast.LOD))+geom_line(color='black', size=1.5)+ 
##        facet_grid(~chrom, scales='free_x', space='free_x')+
##        theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle('LOD')
#neglogp.plot=ggplot(results.data, aes(x=physical.position, y=-log10(afd.contrast.p)))+geom_line(color='black', size=1.5)+ 
#         geom_vline(data=simulatedQTL, aes(xintercept=physical.position), color='blue')+
#        facet_grid(~chrom, scales='free_x', space='free_x')+
#        geom_hline(aes(yintercept=-log10(.05/600)), color='red')+
#        theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))+ggtitle('-log10(p)')
#
# ggpubr::ggarrange(high.plot, low.plot, contrast.plot, neglogp.plot, nrow=4) 














