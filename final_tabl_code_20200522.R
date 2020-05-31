library(metap)

#DEFAULT COUNT, DEFAULT TEST
## screen I

gene_dat_scI = read.delim("test_data/defaultTest_defaultNormCount_screen1/defaultTest_defaultNormCount_screen1.gene_summary.txt",header=TRUE) #19292
sgrna_dat_scI = read.delim("test_data/defaultTest_defaultNormCount_screen1/defaultTest_defaultNormCount_screen1.sgrna_summary.txt",header=TRUE) #77736

## screen II

gene_dat_scII = read.delim("test_data/defaultTest_defaultNormCount_screen2/defaultTest_defaultNormCount_screen2.gene_summary.txt",header=TRUE) #19292
sgrna_dat_scII = read.delim("test_data/defaultTest_defaultNormCount_screen2/defaultTest_defaultNormCount_screen2.sgrna_summary.txt",header=TRUE) #77736

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

library(tidyverse)

combined_table <- scI_II_comb_ord_sel %>%
  as_tibble() %>%
  group_by(Gene) %>%
  mutate(avg=mean(c(SCI_FC, SCII_FC), na.rm=T)) %>%
  ungroup() %>%
  arrange(desc(avg)) %>%
  mutate("Fold Enrichment" = 2^avg) %>%
  select(Gene,
         "Screen 1 P value" = P1,
         "Screen 2 P value" = P2,
         "Fisher P value" = FisherP,
         "Screen 1 log2 Fold Enrichment" = SCI_FC,
         "Screen 2 log2 Fold Enrichment" = SCII_FC,
         "Combined log2 Fold Enrichment" = avg,
         "Fold Enrichment",
         "Screen 1 'Good' sgRNAs" = SCI_GSGRNA,
         "Screen 2 'Good' sgRNAs" = SCII_GSGRNA) 

write_csv(combined_table, path = "combined_screens_entire_list_20200531.csv")

combined_table <- scI_II_comb_ord_sel %>%
  as_tibble() %>%
  group_by(Gene) %>%
  mutate(avg=mean(c(SCI_FC, SCII_FC), na.rm=T)) %>%
  ungroup() %>%
  filter(FisherP < 0.01) %>%
  arrange(desc(avg)) %>%
  mutate("Fold Enrichment" = 2^avg) %>%
  select(Gene,
         "Screen 1 P value" = P1,
         "Screen 2 P value" = P2,
         "Fisher P value" = FisherP,
         "Screen 1 log2 Fold Enrichment" = SCI_FC,
         "Screen 2 log2 Fold Enrichment" = SCII_FC,
         "Combined log2 Fold Enrichment" = avg,
         "Fold Enrichment",
         "Screen 1 'Good' sgRNAs" = SCI_GSGRNA,
         "Screen 2 'Good' sgRNAs" = SCII_GSGRNA) 

write_csv(combined_table, path = "combined_screens_signif_list_20200531.csv")

library(gt)

combined_table %>%
  head(50) %>%
  gt() %>%
  tab_header(
    title = md("Top Gene Targets from **Screens 1 and 2**")
  ) %>%
  gtsave("top50.png")


combined_table %>%
  select(Gene, 
         as.name("Enrichment"), 
         as.name("P value")) %>%
  # top_n(10) %>%
  head(50) %>%
  mutate(Rank = row_number()) %>%
  select(4, 1:3) %>%
  gt() %>%
  tab_header(
    title = md("Top Gene Targets from **Screens 1 and 2**")
  ) %>%
  fmt_number(columns = 3,
    decimals = 2,
    use_seps = TRUE
  ) %>%
  gtsave("top30.png")
  
combined_table %>%
  write_delim("final_combine.txt", delim = "\t")


library(formattable)
top_10_formattable <- combined_table %>%
  select(Gene, 
         as.name("Fold Enrichment"), 
         as.name("Fisher P value")) %>%
  # top_n(10) %>%
  head(10) %>%
  # mutate(Rank = row_number()) %>%
  # select(4, 1:3) %>%
  formattable(align = c("l", rep("r", 3)), 
              list(Gene = formatter("span", style = ~ style(color = "grey", font.weight = "bold", "font-family" = "arial")), 
                   area(col = 2) ~ color_tile("#DeF7E9", "#71CA97")))

#To get a png of the table, I viewed the html of the table, zoomed Chrome in to 250% and took a screenshot. Resized in Mac Preview



#TOTAL COUNT, TOTAL TEST
## screen I

gene_dat_scI_totalC_totalT = read.delim("test_data/totalNormTest_totalNormCount_screen1/totalNormTest_totalNormCount_screen1.gene_summary.txt",header=TRUE) #19292
sgrna_dat_scI_totalC_totalT = read.delim("test_data/totalNormTest_totalNormCount_screen1/totalNormTest_totalNormCount_screen1.sgrna_summary.txt",header=TRUE) #77736

## screen II

gene_dat_scII_totalC_totalT = read.delim("test_data/totalNormTest_totalNormCount_screen2/totalNormTest_totalNormCount_screen2.gene_summary.txt",header=TRUE) #19292
sgrna_dat_scII_totalC_totalT = read.delim("test_data/totalNormTest_totalNormCount_screen2/totalNormTest_totalNormCount_screen2.sgrna_summary.txt",header=TRUE) #77736



gene_dat_scI2_totalC_totalT = gene_dat_scI_totalC_totalT[order(gene_dat_scI_totalC_totalT$id),]
gene_dat_scII2_totalC_totalT = gene_dat_scII_totalC_totalT[order(gene_dat_scII_totalC_totalT$id),]

## Assuming independence
## combine with Fisher's Method


pmat = data.frame(P1=gene_dat_scI2_totalC_totalT$pos.p.value,P2=gene_dat_scII2_totalC_totalT$pos.p.value)

plotp(pmat)

fisher_res = apply(pmat,1,sumlog)

p_vec = c()
for(i in fisher_res)
{
  p_vec = c(p_vec,i$p)
}

pmat$FisherP = p_vec

#fisher_res2 = do.call(rbind.data.frame,fisher_res)

scI_II_comb_totalC_totalT = data.frame(Gene=gene_dat_scI2_totalC_totalT$id,
                                       pmat,
                                       SCI_FDR=gene_dat_scI2_totalC_totalT$pos.fdr,
                                       SCI_RANK=gene_dat_scI2_totalC_totalT$pos.rank,
                                       SCI_GSGRNA=gene_dat_scI2_totalC_totalT$pos.goodsgrna,
                                       SCI_FC=gene_dat_scI2_totalC_totalT$pos.lfc,
                                       SCII_FDR=gene_dat_scII2_totalC_totalT$pos.fdr,
                                       SCII_RANK=gene_dat_scII2_totalC_totalT$pos.rank,
                                       SCII_GSGRNA=gene_dat_scII2_totalC_totalT$pos.goodsgrna,
                                       SCII_FC=gene_dat_scII2_totalC_totalT$pos.lfc)

scI_II_comb_ord_totalC_totalT = scI_II_comb_totalC_totalT[order(scI_II_comb_totalC_totalT$FisherP),]

scI_II_comb_ord_sel_totalC_totalT = subset(scI_II_comb_ord_totalC_totalT,SCI_GSGRNA>=2&SCII_GSGRNA>=2)