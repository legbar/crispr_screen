


library(biomaRt) ## in R   
library(rDGIdb) ##in Bioconductor ## https://rdrr.io/bioc/rDGIdb/f/inst/doc/vignette.pdf
library(mygene)
library(devtools)
library(rentrez)
library(calibrate)
library(ggplot2)
library(dplyr)


## set up the mart
mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",         
                    mart    = useMart("ENSEMBL_MART_ENSEMBL",       
                    host    = "www.ensembl.org"))    



## for final run using the controls as nointron list

#outf2 = paste0(inpath,"negcontrol_sgrna.list")
#write.table(paste(negcol_lib$NAME),file=outf2,quote=FALSE,row.names=FALSE,col.names=FALSE)

#outf22 = paste0(inpath,"negcontrol_sgrna_nointron.list")
#write.table(paste(negcol_lib2$NAME),file=outf22,quote=FALSE,row.names=FALSE,col.names=FALSE)


#########################################
## read in results
#########################################

###############################################################################################################################
## check the gRNA library to test whether all the gRNA are well and evenly distributed/represented
## Transfect the gRNA with lintovirus into the pool of cells
## Flow Cytometory to sort out the cell population, and prepare the top 5-10% of the cell population for gene expressed and 
## 	control population, in this case, 50% of the overall population
## sequence through PCR of particular primer for the gRNA (it will amplify wherever the guide will be)
## Analyze all the counts of gene and gRNA levels
## Significance is defined not only p-values but also, more guides within the gene shows the similar effect, both in p-values
## 	as well as directions (at least 50% of the gRNA within identified genes with sig p-values 
##	will have the similar significance and direction)
###############################################################################################################################

#DEFAULT COUNT, DEFAULT TEST
## screen I

gene_dat_scI = read.delim("test_data/defaultTest_defaultNormCount_screen1/defaultTest_defaultNormCount_screen1.gene_summary.txt",header=TRUE) #19292
sgrna_dat_scI = read.delim("test_data/defaultTest_defaultNormCount_screen1/defaultTest_defaultNormCount_screen1.sgrna_summary.txt",header=TRUE) #77736

## screen II

gene_dat_scII = read.delim("test_data/defaultTest_defaultNormCount_screen2/defaultTest_defaultNormCount_screen2.gene_summary.txt",header=TRUE) #19292
sgrna_dat_scII = read.delim("test_data/defaultTest_defaultNormCount_screen2/defaultTest_defaultNormCount_screen2.sgrna_summary.txt",header=TRUE) #77736


# raw_count_scII = read.delim(inf23,header=TRUE) ## sum(raw_count_scII$BL2) 142728229,sum(raw_count_scII$FXN2) 102431968
#norm_count = read.delim(inf24,header=TRUE) ##
# 
# raw_count_scII$BL2_LOG = ifelse(raw_count_scII$BL2>0,log10(raw_count_scII$BL2),0)
# raw_count_scII$FXN2_LOG = ifelse(raw_count_scII$FXN2>0,log10(raw_count_scII$FXN2),0)

 
# xx = intersect(raw_count_scI$sgRNA,raw_count_scII$sgRNA)
# 
# raw_count_scI2 = raw_count_scI[match(xx,raw_count_scI$sgRNA),]
# raw_count_scII2 = raw_count_scII[match(xx,raw_count_scII$sgRNA),]
# 
# raw_count_comb = data.frame(sgRNA=raw_count_scI2$sgRNA,Gene=raw_count_scI2$Gene,BL1=raw_count_scI$BL,BL2=raw_count_scII2$BL2,
#                             FXN1=raw_count_scI2$FXN,FXN2=raw_count_scII2$FXN2)
# 
# outf001 = paste0(outpath,"FXN_SCREEN_combined.count.txt")
# write.table(raw_count_comb,file=outf001,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

raw_count_comb = read.delim("count_data/default/count_20200409_.count.txt",header=TRUE) %>%
  mutate(baseline_1_log = log10(baseline_1 + 1), 
         fxn_1_log = log10(fxn_1 + 1), 
         baseline_2_log = log10(baseline_2 + 1), 
         fxn_2_log = log10(fxn_2 + 1))

## plot the raw counts for both screens

# rg_max = max(raw_count_comb[,c(7,8)])
# rg_min = min(raw_count_comb[,c(7,8)])
# raw_count = raw_count_comb
# 
# jpeg("screen_1_logCounts.jpeg",width = 600, height = 600, units = "px", pointsize = 12,quality = 100)
# plot(0,0,xlim=c(rg_min,rg_max),ylim=c(rg_min,rg_max),type="n",xlab="BL",ylab="FXN+",main="Screen I")
# points(raw_count$baseline_1_log,raw_count$fxn_1_log,col="red")
# abline(a=0,b=1,lty=2,col="gray")
# dev.off()
# 
# rg_max = max(raw_count_comb[,c(9,10)])
# rg_min = min(raw_count_comb[,c(9,10)])
# raw_count = raw_count_comb
# 
# outg12 = paste0(outpath,"sgRNA_scII_nointron_counts_FXN_BL.jpeg")
# 
# jpeg("screen_2_logCounts.jpeg",width = 600, height = 600, units = "px", pointsize = 12,quality = 100)
# plot(0,0,xlim=c(rg_min,rg_max),ylim=c(rg_min,rg_max),type="n",xlab="BL",ylab="FXN+",main="Screen II")
# points(raw_count$baseline_2_log,raw_count$fxn_2_log,col="red")
# abline(a=0,b=1,lty=2,col="gray")
# dev.off()


## process results
# 
# fdr_cf = 0.1 ## fdr cutoff 10%
# lfc_cf = log(1) ## log fold change (>0, positive change)
#  
# range(gene_dat_scI$pos.fdr)
# range(gene_dat_scII$pos.fdr)

 
#gene_dat_scI_ord_pos = gene_dat_scI[order(gene_dat_scI$pos.rank),]
#gene_dat_scII_ord_pos = gene_dat_scII[order(gene_dat_scII$pos.rank),]
#gene_dat_scI_ord_pos_top10 = head(gene_dat_scI_ord_pos,10)
#gene_dat_scII_ord_pos_top10 = head(gene_dat_scII_ord_pos,10)



## filtering the results
# gene_dat_scI_fdr_pos = subset(gene_dat_scI,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2) ## 474
# gene_dat_scII_fdr_pos = subset(gene_dat_scII,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2) ## 156
#  
# gene_dat_scI_fdr_pos_ord = gene_dat_scI_fdr_pos[order(gene_dat_scI_fdr_pos$pos.rank),]
# gene_dat_scII_fdr_pos_ord = gene_dat_scII_fdr_pos[order(gene_dat_scII_fdr_pos$pos.rank),]
#   
# gene_dat_scI_fdr_pos_ord_top10 = head(gene_dat_scI_fdr_pos_ord,10)
# gene_dat_scII_fdr_pos_ord_top10 = head(gene_dat_scII_fdr_pos_ord,10)
# 
# outt11 = paste0(outpath,"gene_dat_scI_fdr_pos_ord_top10.csv")
# outt21 = paste0(outpath,"gene_dat_scII_fdr_pos_ord_top10.csv")
# 
# write.csv(gene_dat_scI_fdr_pos_ord_top10,file=outt11)
# write.csv(gene_dat_scII_fdr_pos_ord_top10,file=outt21)
# 
# 
# 
# sgrna_dat_scI_fdr_pos = sgrna_dat_scI[sgrna_dat_scI$Gene%in%gene_dat_scI_fdr_pos$id,] ## 1896
# sgrna_dat_scII_fdr_pos = sgrna_dat_scII[sgrna_dat_scII$Gene%in%gene_dat_scII_fdr_pos$id,] ## 624
# 
# gene_dat_scI_fdr_pos$Gene = unlist(lapply(gene_dat_scI_fdr_pos$id,term_extr))
# gene_dat_scII_fdr_pos$Gene = unlist(lapply(gene_dat_scII_fdr_pos$id,term_extr))
# 
# 
# sh_gene = intersect(gene_dat_scI_fdr_pos$id,gene_dat_scII_fdr_pos$id) ## 18 genes
# 
# gene_dat_scI_fdr_pos_sh = gene_dat_scI_fdr_pos[match(sh_gene,gene_dat_scI_fdr_pos$id),]
# gene_dat_scII_fdr_pos_sh = gene_dat_scII_fdr_pos[match(sh_gene,gene_dat_scII_fdr_pos$id),]
# 
# 
# gene_col = c("red","green","hotpink","blue","orange","brown","plum","black","sienna","turquoise") ##sample(colors(),length(genein))
# topn = 10


## screen I
## only draw top 10

if(nrow(gene_dat_scI_fdr_pos)<topn)
{
  topin = nrow(gene_dat_scI_fdr_pos)
} else {
  topin = topn
}
 
genein_scI = paste((gene_dat_scI_fdr_pos[order(gene_dat_scI_fdr_pos$pos.rank),]$id))[1:topin] ## ordered by significance of genes

rg_max = max(raw_count_scI[,c(5,6)])
rg_min = min(raw_count_scI[,c(5,6)])
raw_count = raw_count_scI


outg21 = paste0(outpath,"sgRNA_nointron_raw_counts_scI_sig_top10.jpeg")

jpeg(outg21,width = 600, height = 600, units = "px", pointsize = 12,quality = 100)

plot(0,0,xlim=c(rg_min,rg_max),ylim=c(rg_min,rg_max),type="n",xlab="log(# of BL sequences)",ylab="log(# of FXN+ sequences)",main=paste("top ",topin," genes - Screen I",sep=""))
points(raw_count$BL_LOG,raw_count$FXN_LOG,col="grey")
abline(a=0,b=1,lty=2,col="black")
 
gen_col_tb = data.frame(GENE=genein_scI,COLOR=gene_col[1:length(genein_scI)])

for(i in 1:length(genein_scI))
{
  tmp = subset(sgrna_dat_scI_fdr_pos,Gene==genein_scI[i])
 
  tmp_raw = raw_count[match(tmp$sgrna,raw_count$sgRNA),]

  points(tmp_raw$BL_LOG,tmp_raw$FXN_LOG,col=gene_col[i],pch=18,cex=1.2)

}

legend(0,5,paste(gen_col_tb$GENE),col=paste(gen_col_tb$COLOR),pch=18,bty="n")

dev.off()






## screen II
## only draw top 10

if(nrow(gene_dat_scII_fdr_pos)<topn)
{
  topin = nrow(gene_dat_scII_fdr_pos)
} else {
  topin = topn
}
 
genein_scII = paste((gene_dat_scII_fdr_pos[order(gene_dat_scII_fdr_pos$pos.rank),]$id))[1:topin] ## ordered by significance of genes

rg_max = max(raw_count_scII[,c(5,6)])
rg_min = min(raw_count_scII[,c(5,6)])
raw_count = raw_count_scII


outg22 = paste0(outpath,"sgRNA_nointron_raw_counts_scII_sig_top10.jpeg")

jpeg(outg22,width = 600, height = 600, units = "px", pointsize = 12,quality = 100)

plot(0,0,xlim=c(rg_min,rg_max),ylim=c(rg_min,rg_max),type="n",xlab="log(# of BL sequences)",ylab="log(# of FXN+ sequences)",main=paste("top ",topin," genes - Screen II",sep=""))
points(raw_count$BL2_LOG,raw_count$FXN2_LOG,col="grey")
abline(a=0,b=1,lty=2,col="black")
 
gen_col_tb = data.frame(GENE=genein_scII,COLOR=gene_col[1:length(genein_scII)])

for(i in 1:length(genein_scII))
{
  tmp = subset(sgrna_dat_scII_fdr_pos,Gene==genein_scII[i])
 
  tmp_raw = raw_count[match(tmp$sgrna,raw_count$sgRNA),]

  points(tmp_raw$BL2_LOG,tmp_raw$FXN2_LOG,col=gene_col[i],pch=18,cex=1.2)

}

legend(0,5,paste(gen_col_tb$GENE),col=paste(gen_col_tb$COLOR),pch=18,bty="n")

dev.off()



######################################################
## subsampling approach
## randomly sampling 50 times, each time 
## 50% of the read from fastq of FXN+ and Baseline 
## will be sampled and analyzed (with the same seed)
######################################################



sim_num = 50 
gene_pos_scI_dat = data.frame()
sgrna_pos_scI_dat = data.frame()

gene_pos_scII_dat = data.frame()
sgrna_pos_scII_dat = data.frame()

## create a null distribution
gene_scI_null = gene_scII_null = data.frame()



for(i in 1:sim_num)
{
  
  inf110 = paste0(outpath,"subsm_ctrl_nointron_screenI_sgrna",i,".sgrna_summary.txt")
  inf120 = paste0(outpath,"subsm_ctrl_nointron_screenI_sgrna",i,".gene_summary.txt")

  inf210 = paste0(outpath,"subsm_ctrl_SCII_NEG_CTRL_nointron_sgrna",i,".sgrna_summary.txt")
  inf220 = paste0(outpath,"subsm_ctrl_SCII_NEG_CTRL_nointron_sgrna",i,".gene_summary.txt")  
 
  gene_sub_dat_scI = read.delim(inf120,header=TRUE)
  sgrna_sub_dat_scI = read.delim(inf110,header=TRUE)
 
  gene_sub_dat_scI_ord = gene_sub_dat_scI[order(gene_sub_dat_scI$pos.rank),]

  gene_sub_dat_scII = read.delim(inf220,header=TRUE)
  sgrna_sub_dat_scII = read.delim(inf210,header=TRUE)

  gene_sub_dat_scII_ord = gene_sub_dat_scII[order(gene_sub_dat_scII$pos.rank),]

  #subn2 = subset(gene_sub_dat,neg.fdr<=fdr_cf)
  #subp2 = subset(gene_sub_dat,pos.fdr<=fdr_cf) #7
 
  #gene_sub_dat_ord2 = subset(gene_sub_dat_ord,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2)
  #gene_sub_dat_null = subset(gene_sub_dat_ord,pos.fdr>fdr_cf)

  gene_sub_dat_scI_ord2 = subset(gene_sub_dat_scI_ord,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2)
  gene_sub_dat_scII_ord2 = subset(gene_sub_dat_scII_ord,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2)

  #gene_sub_dat_ord2$TRIAL = i
  #gene_sub_dat_null$TRIAL = i

  gene_sub_dat_scI_ord2$TRIAL = i
  gene_sub_dat_scII_ord2$TRIAL = i

  #sgrna_sub_dat_top = sgrna_sub_dat[sgrna_sub_dat$Gene%in%gene_sub_dat_ord2$id,] ## 18
  #sgrna_sub_dat_top2 = sgrna_sub_dat_top[order(sgrna_sub_dat_top$sgrna),]
  #sgrna_sub_dat_top2$TRIAL = i

  sgrna_sub_dat_scI_top = sgrna_sub_dat_scI[sgrna_sub_dat_scI$Gene%in%gene_sub_dat_scI_ord2$id,]
  sgrna_sub_dat_scII_top = sgrna_sub_dat_scII[sgrna_sub_dat_scII$Gene%in%gene_sub_dat_scII_ord2$id,]

  sgrna_sub_dat_scI_top2 = sgrna_sub_dat_scI_top[order(sgrna_sub_dat_scI_top$sgrna),]
  sgrna_sub_dat_scII_top2 = sgrna_sub_dat_scII_top[order(sgrna_sub_dat_scII_top$sgrna),]

  sgrna_sub_dat_scI_top2$TRIAL = i
  sgrna_sub_dat_scII_top2$TRIAL = i

  if(nrow(gene_pos_scI_dat)==0)
  {
    gene_pos_scI_dat = gene_sub_dat_scI_ord2    
  } else {
    gene_pos_scI_dat = rbind(gene_pos_scI_dat,gene_sub_dat_scI_ord2)
  }

  if(nrow(sgrna_pos_scI_dat)==0)
  {
    sgrna_pos_scI_dat = sgrna_sub_dat_scI_top2  
  } else {
    sgrna_pos_scI_dat = rbind(sgrna_pos_scI_dat,sgrna_sub_dat_scI_top2)
  }


  if(nrow(gene_pos_scII_dat)==0)
  {
    gene_pos_scII_dat = gene_sub_dat_scII_ord2    
  } else {
    gene_pos_scII_dat = rbind(gene_pos_scII_dat,gene_sub_dat_scII_ord2)
  }

  if(nrow(sgrna_pos_scII_dat)==0)
  {
    sgrna_pos_scII_dat = sgrna_sub_dat_scII_top2  
  } else {
    sgrna_pos_scII_dat = rbind(sgrna_pos_scII_dat,sgrna_sub_dat_scII_top2)
  }

  #if(nrow(gene_null)==0)
  #{
  #  gene_null = gene_sub_dat_null
  #} else {
  #  gene_null = rbind(gene_null,gene_sub_dat_null)
  #}

  print(i)

}

## 

outt12 = paste0(outpath,"sgrna_pos_scI_dat_subsampling_QCed.RData")
outt22 = paste0(outpath,"sgrna_pos_scII_dat_subsampling_QCed.RData")

save(sgrna_pos_scI_dat,file=outt12)
save(sgrna_pos_scII_dat,file=outt22)

outt13 = paste0(outpath,"gene_pos_scI_dat_subsampling_QCed.RData")
outt23 = paste0(outpath,"gene_pos_scII_dat_subsampling_QCed.RData")

save(gene_pos_scI_dat,file=outt13)
save(gene_pos_scII_dat,file=outt23)

load(outt13) 
load(outt23)





##

## get estimates of lfc from 50 subsampling


sim_gene_list_scI = unique(gene_pos_scI_dat$id)

lfc_sum_scI = data.frame(matrix(-1,length(sim_gene_list_scI),8))
sim_n_scI = c()

for(j in 1:length(sim_gene_list_scI))
{
  tmp_dat = subset(gene_pos_scI_dat,id==sim_gene_list_scI[j])$pos.lfc

  #tmp = log10(CI_gen(10^(tmp_dat),0.05))
  tmp = log(CI_gen(exp(tmp_dat),0.05))

  lfc_sum_scI[j,] = tmp

  sim_n_scI[j] = length(tmp_dat)
}


lfc_sum_scI2 = data.frame(Gene=sim_gene_list_scI,N=sim_n_scI,N_PERC=sim_n_scI/50,lfc_sum_scI)
colnames(lfc_sum_scI2)[-c(1:3)] = c("MIN","F_Qu","MEDIAN","MEAN","T_Qu","MAX","L95","U95")

gene_scI_sig_tb = table(gene_pos_scI_dat$TRIAL,paste(gene_pos_scI_dat$id))
gene_scI_sig_ratio = apply(gene_scI_sig_tb,2,sum)/nrow(gene_scI_sig_tb)
res_subsample_scI_list = list(lfc_sum_scI2,gene_scI_sig_tb,gene_scI_sig_ratio)


lfc_sum_scI3 = lfc_sum_scI2[match(names(gene_scI_sig_ratio),lfc_sum_scI2$Gene),]
lfc_sum_scI3$REPLI_RATIO = gene_scI_sig_ratio
lfc_sum_scI3_ranked = lfc_sum_scI3[with(lfc_sum_scI3,order(-N,MEDIAN)),]

dim(subset(lfc_sum_scI3_ranked,REPLI_RATIO>=0.8)) ## 228


sim_gene_list_scII = unique(gene_pos_scII_dat$id)

lfc_sum_scII = data.frame(matrix(-1,length(sim_gene_list_scII),8))
sim_n_scII = c()

for(j in 1:length(sim_gene_list_scII))
{
  tmp_dat = subset(gene_pos_scII_dat,id==sim_gene_list_scII[j])$pos.lfc

  #tmp = log10(CI_gen(10^(tmp_dat),0.05))
  tmp = log(CI_gen(exp(tmp_dat),0.05))

  lfc_sum_scII[j,] = tmp

  sim_n_scII[j] = length(tmp_dat)
}

lfc_sum_scII2 = data.frame(Gene=sim_gene_list_scII,N=sim_n_scII,N_PERC=sim_n_scII/50,lfc_sum_scII)
colnames(lfc_sum_scII2)[-c(1:3)] = c("MIN","F_Qu","MEDIAN","MEAN","T_Qu","MAX","L95","U95")

gene_scII_sig_tb = table(gene_pos_scII_dat$TRIAL,paste(gene_pos_scII_dat$id))
gene_scII_sig_ratio = apply(gene_scII_sig_tb,2,sum)/nrow(gene_scII_sig_tb)
res_subsample_scII_list = list(lfc_sum_scII2,gene_scII_sig_tb,gene_scII_sig_ratio)


lfc_sum_scII3 = lfc_sum_scII2[match(names(gene_scII_sig_ratio),lfc_sum_scII2$Gene),]
lfc_sum_scII3$REPLI_RATIO = gene_scII_sig_ratio
lfc_sum_scII3_ranked = lfc_sum_scII3[with(lfc_sum_scII3,order(-N,MEDIAN)),]

dim(subset(lfc_sum_scII3_ranked,REPLI_RATIO>=0.8)) ## 228




  
lfc_sum_scI3_ranked_sig = subset(lfc_sum_scI3_ranked,REPLI_RATIO>=0.8)
lfc_sum_scII3_ranked_sig = subset(lfc_sum_scII3_ranked,REPLI_RATIO>=0.8)

lfc_sum_scI3_ranked_sig$SCREEN = "I"
lfc_sum_scII3_ranked_sig$SCREEN = "II"

subsample_sig_sh = intersect(lfc_sum_scI3_ranked_sig$Gene,lfc_sum_scII3_ranked_sig$Gene)
subsample_sig_uniq = union(lfc_sum_scI3_ranked_sig$Gene,lfc_sum_scII3_ranked_sig$Gene)

lfc_sum_comb_ranked_sig = rbind(lfc_sum_scI3_ranked_sig,lfc_sum_scII3_ranked_sig)
lfc_sum_comb_ranked_sig2 = lfc_sum_comb_ranked_sig[order(lfc_sum_comb_ranked_sig$Gene),]

SH = rep("",dim(lfc_sum_comb_ranked_sig2)[1])
SH[lfc_sum_comb_ranked_sig2$Gene%in%subsample_sig_sh] = "Y"
lfc_sum_comb_ranked_sig2$SH = SH

  
gene_dat_scI_ol = gene_dat_scI[match(subsample_sig_sh,gene_dat_scI$id),]
gene_dat_scII_ol = gene_dat_scII[match(subsample_sig_sh,gene_dat_scII$id),]

lfc_sum_scI3_ranked_sig_sh = lfc_sum_scI3_ranked_sig[match(subsample_sig_sh,lfc_sum_scI3_ranked_sig$Gene),]
lfc_sum_scII3_ranked_sig_sh = lfc_sum_scII3_ranked_sig[match(subsample_sig_sh,lfc_sum_scII3_ranked_sig$Gene),]


ol_scI_comb = data.frame(gene_dat_scI_ol,lfc_sum_scI3_ranked_sig_sh)
x11 = round(ol_scI_comb$MEAN,2)
x21 = round(ol_scI_comb$L95,2)
x31 = round(ol_scI_comb$U95,2)
LFC_SUB3 = paste0(x11," (",x21,"-",x31,")")

ol_scI_comb2 = data.frame(Gene=ol_scI_comb$id,Rank=ol_scI_comb$pos.rank,sgrna_num=ol_scI_comb$pos.goodsgrna,FDR=ol_scI_comb$pos.fdr,Orig_LFC=ol_scI_comb$pos.lfc,
                          Repli_ratio=ol_scI_comb$N_PERC,SIM_LFC=LFC_SUB3)



ol_scII_comb = data.frame(gene_dat_scII_ol,lfc_sum_scII3_ranked_sig_sh)
x11 = round(ol_scII_comb$MEAN,2)
x21 = round(ol_scII_comb$L95,2)
x31 = round(ol_scII_comb$U95,2)
LFC_SUB4 = paste0(x11," (",x21,"-",x31,")")

ol_scII_comb2 = data.frame(Gene=ol_scII_comb$id,Rank=ol_scII_comb$pos.rank,sgrna_num=ol_scII_comb$pos.goodsgrna,FDR=ol_scI_comb$pos.fdr,Orig_LFC=ol_scII_comb$pos.lfc,
                          Repli_ratio=ol_scII_comb$N_PERC,SIM_LFC=LFC_SUB4)



outt15 = paste0(outpath,"gene_dat_scI_sim_ol_list.csv")
outt25 = paste0(outpath,"gene_dat_scII_sim_ol_list.csv")

write.csv(ol_scI_comb2,file=outt15,row.names=FALSE)
write.csv(ol_scII_comb2,file=outt25,row.names=FALSE)



gene_dat1 = read.csv(outt15,header=TRUE)
gene_dat2 = read.csv(outt25,header=TRUE)





## gene_dat_scI_fdr_pos_ord_top10
## generate the sim results for the top signals

lfc_sum_scI2_orig_top10 = lfc_sum_scI2[match(gene_dat_scI_fdr_pos_ord_top10$id,lfc_sum_scI2$Gene),]

lfc_sum_scII2_orig_top10 = lfc_sum_scII2[match(gene_dat_scII_fdr_pos_ord_top10$id,lfc_sum_scII2$Gene),]



x1 = round(lfc_sum_scI2_orig_top10$MEAN,2)
x2 = round(lfc_sum_scI2_orig_top10$L95,2)
x3 = round(lfc_sum_scI2_orig_top10$U95,2)

LFC_SUB1 = paste0(x1," (",x2,"-",x3,")")

gene_dat_scI_fdr_pos_ord_top10_sim = data.frame(Gene=gene_dat_scI_fdr_pos_ord_top10$id,Num=gene_dat_scI_fdr_pos_ord_top10$num,Pvalue=gene_dat_scI_fdr_pos_ord_top10$pos.p.value,
                                       FDR=gene_dat_scI_fdr_pos_ord_top10$pos.fdr,RANK=gene_dat_scI_fdr_pos_ord_top10$pos.rank,
						   Replicate_Ratio=lfc_sum_scI2_orig_top10$N_PERC,
                                       LFC_ORIG = gene_dat_scI_fdr_pos_ord_top10$pos.lfc,SUB_LFC_CI=LFC_SUB1)


x1 = round(lfc_sum_scII2_orig_top10$MEAN,2)
x2 = round(lfc_sum_scII2_orig_top10$L95,2)
x3 = round(lfc_sum_scII2_orig_top10$U95,2)

LFC_SUB2 = paste0(x1," (",x2,"-",x3,")")

gene_dat_scII_fdr_pos_ord_top10_sim = data.frame(Gene=gene_dat_scII_fdr_pos_ord_top10$id,Num=gene_dat_scII_fdr_pos_ord_top10$num,Pvalue=gene_dat_scII_fdr_pos_ord_top10$pos.p.value,
                                       FDR=gene_dat_scII_fdr_pos_ord_top10$pos.fdr,RANK=gene_dat_scII_fdr_pos_ord_top10$pos.rank,
						   Replicate_Ratio=lfc_sum_scII2_orig_top10$N_PERC,
                                       LFC_ORIG = gene_dat_scII_fdr_pos_ord_top10$pos.lfc,SUB_LFC_CI=LFC_SUB2)


outt14 = paste0(outpath,"gene_dat_scI_fdr_pos_ord_top10_sim.csv")
outt24 = paste0(outpath,"gene_dat_scII_fdr_pos_ord_top10_sim.csv")

write.csv(gene_dat_scI_fdr_pos_ord_top10_sim,file=outt14,row.names=FALSE)
write.csv(gene_dat_scII_fdr_pos_ord_top10_sim,file=outt24,row.names=FALSE)




#########################################################################################################
inf61 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron.gene_summary.txt") 
inf62 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron.sgrna_summary.txt")
 
fdr_cf = 0.1 ## fdr cutoff 10%
lfc_cf = log(1) ## log fold change (>0, positive change)

sum(raw_count_scI$BL) ## 88216257
sum(raw_count_scI$FXN) ## 73132444

sum(raw_count_scII$BL2) ## 142728229
sum(raw_count_scII$FXN2) ## 102431968

gene_dat_comb = read.table(inf61,header=TRUE)
sgrna_dat_comb = read.table(inf62,header=TRUE)

gene_dat_comb_ord = gene_dat_comb[order(gene_dat_comb$pos.rank),]
gene_dat_comb_fdr_pos = subset(gene_dat_comb,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2) ## 474
table(gene_dat_comb_ord$pos.goodsgrna)
 

gene_dat_scI_fdr_pos_ord = gene_dat_scI_fdr_pos[order(gene_dat_scI_fdr_pos$pos.rank),]
gene_dat_scII_fdr_pos_ord = gene_dat_scII_fdr_pos[order(gene_dat_scII_fdr_pos$pos.rank),]
  
gene_dat_scI_fdr_pos_ord_top10 = head(gene_dat_scI_fdr_pos_ord,10)
gene_dat_scII_fdr_pos_ord_top10 = head(gene_dat_scII_fdr_pos_ord,10)





inf71 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron_v2.gene_summary.txt") 
inf72 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron_v2.sgrna_summary.txt")
 
gene_dat_comb2 = read.table(inf71,header=TRUE)
sgrna_dat_comb2 = read.table(inf72,header=TRUE)

gene_dat_comb2_ord = gene_dat_comb2[order(gene_dat_comb2$pos.rank),]
gene_dat_comb2_fdr_pos = subset(gene_dat_comb2,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2) 
table(gene_dat_comb2_ord$pos.goodsgrna)
 











#############################################################################
## test 1
## sample both screen I and II data into the same number of sequences 
#############################################################################

## 85000000 # 88216257 for BL
## 73000000 # 73132444 for FXN
## for down sampling

## prepare the combined counts 

inf11 = paste0(outpath,"FXN_SCREEN_scII_downS.count.count.txt")
inf12 = paste0(outpath,"FXN_SCREEN_scI_downS.count.count.txt")
 
#raw_count_scII = read.delim(inf23,header=TRUE) ## sum(raw_count_scII$BL2) 142728229,sum(raw_count_scII$FXN2) 102431968
#norm_count = read.delim(inf24,header=TRUE) ##
#raw_count_scII$BL2_LOG = ifelse(raw_count_scII$BL2>0,log10(raw_count_scII$BL2),0)
#raw_count_scII$FXN2_LOG = ifelse(raw_count_scII$FXN2>0,log10(raw_count_scII$FXN2),0)

scI_ct_dat = read.delim(inf11,header=TRUE)
scII_ct_dat = read.delim(inf12,header=TRUE)

scI_ct_dat$BL_LOG = ifelse(scI_ct_dat$BL>0,log10(scI_ct_dat$BL),0)
scI_ct_dat$BL_LOG = ifelse(scI_ct_dat$FXN>0,log10(scI_ct_dat$FXN),0)

scII_ct_dat$BL_LOG = ifelse(scII_ct_dat$BL>0,log10(scII_ct_dat$BL),0)
scII_ct_dat$BL_LOG = ifelse(scII_ct_dat$FXN>0,log10(scII_ct_dat$FXN),0)
 
scI_ct_dat2 = scI_ct_dat[order(scI_ct_dat$sgRNA),]
scII_ct_dat2 = scII_ct_dat[order(scII_ct_dat$sgRNA),]

#sum(scII_ct_dat2$sgRNA==scI_ct_dat2$sgRNA)
#nrow(scI_ct_dat2)

scI_II_ct_dat = data.frame(sgRNA=scI_ct_dat2$sgRNA,Gene=scI_ct_dat2$Gene,
                           BL1=scI_ct_dat2$BL,FXN1=scI_ct_dat2$FXN,
                           BL2=scII_ct_dat2$BL,FXN2=scI_ct_dat2$FXN)

outf3 = paste0(outpath,"FXN_SCREEN_combined_downS.count.count.txt")
write.table(scI_II_ct_dat,file=outf3,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")



 

inf21 = paste0(outpath,"FXN_SCREEN_scII_downS_v2.count_normalized.txt")
inf22 = paste0(outpath,"FXN_SCREEN_scI_downS_v2.count_normalized.txt")
 
#raw_count_scII = read.delim(inf23,header=TRUE) ## sum(raw_count_scII$BL2) 142728229,sum(raw_count_scII$FXN2) 102431968
#norm_count = read.delim(inf24,header=TRUE) ##
#raw_count_scII$BL2_LOG = ifelse(raw_count_scII$BL2>0,log10(raw_count_scII$BL2),0)
#raw_count_scII$FXN2_LOG = ifelse(raw_count_scII$FXN2>0,log10(raw_count_scII$FXN2),0)

scI_ct2_dat = read.delim(inf21,header=TRUE)
scII_ct2_dat = read.delim(inf22,header=TRUE)

scI_ct2_dat$BL_LOG = ifelse(scI_ct2_dat$BL>0,log10(scI_ct2_dat$BL),0)
scI_ct2_dat$BL_LOG = ifelse(scI_ct2_dat$FXN>0,log10(scI_ct2_dat$FXN),0)

scII_ct2_dat$BL_LOG = ifelse(scII_ct2_dat$BL>0,log10(scII_ct2_dat$BL),0)
scII_ct2_dat$BL_LOG = ifelse(scII_ct2_dat$FXN>0,log10(scII_ct2_dat$FXN),0)
 
scI_ct2_dat2 = scI_ct2_dat[order(scI_ct2_dat$sgRNA),]
scII_ct2_dat2 = scII_ct2_dat[order(scII_ct2_dat$sgRNA),]

#sum(scII_ct2_dat2$sgRNA==scI_ct2_dat2$sgRNA)
#nrow(scI_ct_dat2)

scI_II_ct2_dat = data.frame(sgRNA=scI_ct2_dat2$sgRNA,Gene=scI_ct2_dat2$Gene,
                           BL1=scI_ct2_dat2$BL,FXN1=scI_ct2_dat2$FXN,
                           BL2=scII_ct2_dat2$BL,FXN2=scI_ct2_dat2$FXN)
 
outf32 = paste0(outpath,"FXN_SCREEN_combined_downS_v2.count_normalized.txt")
write.table(scI_II_ct2_dat,file=outf32,row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")




## handle results

inf71 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron_downS.gene_summary.txt") 
inf712 = paste0(outpath,"FXN_SCREENI_NEG_CTRL_nointron_downS.gene_summary.txt") 
inf713 = paste0(outpath,"FXN_SCREENII_NEG_CTRL_nointron_downS.gene_summary.txt") 
 
#inf81 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron_downS.sgrna_summary.txt") 
#inf812 = paste0(outpath,"FXN_SCREENI_NEG_CTRL_nointron_downS.sgrna_summary.txt") 
#inf813 = paste0(outpath,"FXN_SCREENII_NEG_CTRL_nointron_downS.sgrna_summary.txt") 


inf91 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron_downS_v2.gene_summary.txt") 
inf912 = paste0(outpath,"FXN_SCREENI_NEG_CTRL_nointron_downS_v2.gene_summary.txt") 
inf913 = paste0(outpath,"FXN_SCREENII_NEG_CTRL_nointron_downS_v2.gene_summary.txt") 
 
#inf101 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron_downS_v2.sgrna_summary.txt") 
#inf1012 = paste0(outpath,"FXN_SCREENI_NEG_CTRL_nointron_downS_v2.sgrna_summary.txt") 
#inf1013 = paste0(outpath,"FXN_SCREENII_NEG_CTRL_nointron_downS_v2.sgrna_summary.txt") 
 




sum(scI_ct_dat2$BL)
sum(scI_ct_dat2$FXN)

sum(scII_ct_dat2$BL)
sum(scII_ct_dat2$FXN)


gene_ds_dat_III = read.table(inf71,header=TRUE)
gene_ds_dat_I = read.table(inf712,header=TRUE)
gene_ds_dat_II = read.table(inf713,header=TRUE)

sgrna_ds_dat_III = read.table(inf81,header=TRUE)
sgrna_ds_dat_I = read.table(inf812,header=TRUE)
sgrna_ds_dat_II = read.table(inf813,header=TRUE)


table(gene_ds_dat_III$pos.goodsgrna)
table(gene_ds_dat_I$pos.goodsgrna)
table(gene_ds_dat_II$pos.goodsgrna)


gene_dsn_dat_III = read.table(inf91,header=TRUE)
gene_dsn_dat_I = read.table(inf912,header=TRUE)
gene_dsn_dat_II = read.table(inf913,header=TRUE)

sgrna_dsn_dat_III = read.table(inf101,header=TRUE)
sgrna_dsn_dat_I = read.table(inf1012,header=TRUE)
sgrna_dsn_dat_II = read.table(inf1013,header=TRUE)

table(gene_dsn_dat_III$pos.goodsgrna)
table(gene_dsn_dat_I$pos.goodsgrna)
table(gene_dsn_dat_II$pos.goodsgrna)




fdr_cf = 0.1 ## fdr cutoff 10%
lfc_cf = log(1) ## log fold change (>0, positive change)


gene_ds_dat_III_QC = subset(gene_ds_dat_III,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2)
gene_ds_dat_I_QC = subset(gene_ds_dat_I,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2)
gene_ds_dat_II_QC = subset(gene_ds_dat_II,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2)


gene_dsn_dat_III_QC = subset(gene_dsn_dat_III,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2)
gene_dsn_dat_I_QC = subset(gene_dsn_dat_I,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2)
gene_dsn_dat_II_QC = subset(gene_dsn_dat_II,pos.fdr<=fdr_cf&pos.lfc>=lfc_cf&pos.goodsgrna>=2)




## fdr_cf 0.1
## lfc_cf 0 (fc > 0)
## pos.goodsgrna 

outf11 = paste0(outpath,"screenI_orig_ranked_pos_QC.csv") ## 474
outf12 = paste0(outpath,"screenII_orig_ranked_pos_QC.csv") ## 156
outf13 = paste0(outpath,"screenI_DownSampling_ranked_pos_QC.csv") ## 137
outf14 = paste0(outpath,"screenII_DownSampling_ranked_pos_QC.csv") ## 161

write.csv(gene_dat_scI_fdr_pos_ord,file=outf11,row.names=FALSE)
write.csv(gene_dat_scII_fdr_pos_ord,file=outf12,row.names=FALSE)
write.csv(gene_dsn_dat_I_QC,file=outf13,row.names=FALSE)
write.csv(gene_dsn_dat_II_QC,file=outf14,row.names=FALSE)








subset(gene_ds_dat_III,id=="DDX6_a")
subset(gene_ds_dat_I,id=="DDX6_a")
subset(gene_ds_dat_II,id=="DDX6_a")

subset(sgrna_ds_dat_III,Gene=="DDX6_a")
subset(sgrna_ds_dat_I,Gene=="DDX6_a")
subset(sgrna_ds_dat_II,Gene=="DDX6_a")



subset(gene_ds_dat_III,id=="RPL14_a")
subset(gene_ds_dat_I,id=="RPL14_a")
subset(gene_ds_dat_II,id=="RPL14_a")

subset(sgrna_ds_dat_III,Gene=="RPL14_a")
subset(sgrna_ds_dat_I,Gene=="RPL14_a")
subset(sgrna_ds_dat_II,Gene=="RPL14_a")





#####################################################################################################
## how to combine two experiments
#####################################################################################################

## gene_dat_scII
## gene_dat_scI

library(metap)

gene_dat_scI2 = gene_dat_scI[order(gene_dat_scI$id),]
gene_dat_scII2 = gene_dat_scII[order(gene_dat_scII$id),]

## Assuming independence
## combine with Fisher's Method


pmat = data.frame(P1=gene_dat_scI2$pos.p.value,P2=gene_dat_scII2$pos.p.value)

plotp(pmat)

fisher_res = apply(pmat,1,sumlog)
 
p_vec = c()
for(i in fisher_res)
{
  p_vec = c(p_vec,i$p)
}

pmat$FisherP = p_vec

#fisher_res2 = do.call(rbind.data.frame,fisher_res)
 
scI_II_comb = data.frame(Gene=gene_dat_scI2$id,pmat,SCI_FDR=gene_dat_scI2$pos.fdr,SCI_RANK=gene_dat_scI2$pos.rank,SCI_GSGRNA=gene_dat_scI2$pos.goodsgrna,SCI_FC=gene_dat_scI2$pos.lfc,
                                                    SCII_FDR=gene_dat_scII2$pos.fdr,SCII_RANK=gene_dat_scII2$pos.rank,SCII_GSGRNA=gene_dat_scII2$pos.goodsgrna,SCII_FC=gene_dat_scII2$pos.lfc)

scI_II_comb_ord = scI_II_comb[order(scI_II_comb$FisherP),]
 
scI_II_comb_ord_sel = subset(scI_II_comb_ord,SCI_GSGRNA>=2&SCII_GSGRNA>=2)

outf001 = paste0(outpath,"scI_II_comb_ord_sel.csv")
write.csv(head(scI_II_comb_ord_sel,10),file=outf001,row.names=FALSE)




## Assuming less heterogeneity
## from combining the raw reads
inf16 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron_v2.gene_summary.txt")
inf17 = paste0(outpath,"FXN_SCREEN_combined_NEG_CTRL_nointron_v2.normalized.txt")

gene_dat_SCI_II = read.delim(inf16,header=TRUE) # 19292
raw_count_SCI_II = read.delim(inf17,header=TRUE) # 77736

gene_dat_SCI_II2 = gene_dat_SCI_II[order(gene_dat_SCI_II$pos.p.value),]

gene_dat_SCI_II2_sel = subset(gene_dat_SCI_II2,pos.goodsgrna>=2)
