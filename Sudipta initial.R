if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("universalmotif")

library(universalmotif)

BiocManager::install("MotifDb")
#install.packages("ape")
library("ape")
library(MotifDb)
#View(MotifDb)
motifs <- filter_motifs(MotifDb, organism = "Hsapiens")
#> motifs converted to class 'universalmotif'
# Compare the first motif with everything and return P-values
head(compare_motifs(motifs, 1))
#BiocManager::install("ggtree")
#BiocManager::install("rentrez")
library(rentrez)
library("ggtree")
invest <- read.FASTA("C:/Users/natal/Downloads/sequence.txt")
fetched <- entrez_fetch("")
?entrez_fetch

#invest$`lcl|NM_012156.2_cds_NP_036288.2_1 [gene=EPB41L1] [db_xref=CCDS:CCDS13271.1] [protein=band 4.1-like protein 1 isoform a] [protein_id=NP_036288.2] [location=172..2817] [gbkey=CDS]`
##Warning in compare_motifs(motifs, 1): Some comparisons failed due to low motif
#> IC
#> DataFrame with 6 rows and 8 columns
#> subject subject.i target target.i score logPval
#> <character> <numeric> <character> <integer> <numeric> <numeric>
#> 1 ORA59 1 ERF11 [duplicated #6.. 1371 0.991211 -13.5452
#> 2 ORA59 1 CRF4 [duplicated #566] 1195 0.990756 -13.5247
#> 3 ORA59 1 LOB 1297 0.987357 -13.3725
#> 4 ORA59 1 ERF15 618 0.977213 -12.9254
#> 5 ORA59 1 ERF2 [duplicated #294] 649 0.973871 -12.7804
#> 6 ORA59 1 ERF2 [duplicated #483] 1033 0.973871 -12.7804
#> Pval Eval
#> <numeric> <numeric>
#> 1 1.31042e-06 0.00359318
#> 2 1.33754e-06 0.00366754
#> 3 1.55744e-06 0.00427049
#> 4 2.43548e-06 0.00667809
#> 5 2.81553e-06 0.00772019
#> 6 2.81553e-06 0.00772019
motifs <- filter_motifs(MotifDb, family = c("AP2", "B3", "bHLH", "bZIP",
                                            "AT hook"))
MotifDb
#> motifs converted to class 'universalmotif'
motifs <- motifs[sample(seq_along(motifs), 100)]
tree <- motif_tree(motifs, linecol = "family")

## Make some changes to the tree in regular ggplot2 fashion:
# tree <- tree + ...
tree



library(scales)
library(ggplot2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(AnnotationDbi)



library(cowplot)
## Get our starting set of motifs:
motifs <- convert_motifs(MotifDb[1:10])
## Get the tree: make sure it's a horizontal type layout
tree <- motif_tree(motifs, layout = "rectangular", linecol = "none")

## Now, make sure we order our list of motifs to match the order of tips:
mot.names <- sapply(motifs, function(x) x["name"])
names(motifs) <- mot.names
new.order <- tree$data$label[tree$data$isTip]
new.order <- rev(new.order[order(tree$data$y[tree$data$isTip])])
motifs <- motifs[new.order]
## Plot the two together (finessing of margins and positions may be required):
plot_grid(nrow = 1, rel_widths = c(1, -0.15, 1),
          tree + xlab("try"), NULL,
          view_motifs(motifs, names.pos = "right") +
            ylab(element_blank()) +
            theme(
              axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.text = element_text(colour = "white")
            )
)

#install.packages("XLConnect")

library(XLConnect)

?readWorksheetFromFile

gene.list <- readWorksheetFromFile("C:/Users/natal/Downloads/SG_gene lists_for_5UTR.xlsx", sheet=1)

gene.list <- as.data.frame(gene.list)
typeof(gene.list)

gene.list <- mapIds(org.Hs.eg.db, keys = gene.list[,1], column = "ENTREZID", keytype = "SYMBOL")

View(gene.list)

#write.csv()

write.csv(gene.list, "C:/Users/natal/Downloads/gene_convert.csv")

library("Biostrings")


#BiocManager::install("biomaRt")

library("biomaRt")

s = readDNAStringSet("C:/Users/natal/Downloads/sequence_down.txt", )
View(s$`lcl|NM_012156.2_cds_NP_036288.2_1 [gene=EPB41L1] [db_xref=CCDS:CCDS13271.1] [protein=band 4.1-like protein 1 isoform a] [protein_id=NP_036288.2] [location=172..2817] [gbkey=CDS]`
)
subseq(s, start=10, end =100)

###NOTE FROM THE 8th, CHECK THIS FUNCTION OUT

getSequence_5UTR <- function(ensg, which.mart = mart, which.type = "ensembl_gene_id", which.seq.type = "5utr", fasta.output = "") {
  df.temp <- biomaRt::getSequence(id = ensg, type = which.type, seqType = which.seq.type, mart = which.mart)
  return(df.temp)
  print(paste("Retrieved", nrow(df.temp), "sequence(s) from getSequence. Input had", length(input), "entries."))
  
  if (fasta.output != "") {
    print("Moving any existing file to trash folder using trash command")
    system(paste("trash", fasta.output))
    exportFASTA(df.temp, fasta.output)
    fasta.temp <- gsub(".fasta", ".temp.fasta", fasta.output)
    print(fasta.out)
    print(fasta.output)
    command <- capture.output(cat("awk '", '/^>/{print s? s"\\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}', "' ", fasta.output, " > ", fasta.temp, sep = ""))
    system(command)
    system(paste("mv",  fasta.temp, fasta.output))
    system(paste("trash", fasta.temp))
  }
}




View(s)
s$`lcl|NM_001199849.1_cds_NP_001186778.1_1 [gene=DAP3] [db_xref=CCDS:CCDS1120.1] [protein=28S ribosomal protein S29, mitochondrial isoform 1] [protein_id=NP_001186778.1] [location=185..1381] [gbkey=CDS]`[0:50]

#can select top 100 nucleotides
values <- subseq(s, start=c(1), end=c(100))
#View(values)
library(R.utils)
library(pipeR)
library(rlist)
library(stringr)
library(dplyr)
g <- rtracklayer::import('C:/Users/natal/Downloads/gencode.v38.annotation.gtf')


gtf2 <- as.data.frame(g)

gtf2$gene_id


name <- read.delim("C:/Users/natal/Downloads/mart_export (2).txt")
name$Gene.stable.ID
gtf2$gene_id

?str_detect()

second_df <-  gtf2[gtf2$gene_id %in% "ENSG00000134882",]
View(gtf2)
df2 <- gtf2[stringr::str_detect(gtf2$gene_id, paste(name$Gene.stable.ID, collapse = "|")),]

dim(df2)

df_utr <- df2[str_detect(df2$type, "UTR"),]
df_utr  

df_utr_ex1 <- df_utr[df_utr$exon_number == "1",]

df_sud <- df_utr_ex1[df_utr_ex1$gene_type == "protein_coding" & df_utr_ex1$transcript_type =="protein_coding",]

df_quality <- df_sud[df_sud$transcript_support_level == "1",]


df_quality2 <- df_quality[df_quality$exon_number == "1",]
#View(df_quality22)
df_quality2 <- df_quality2[df_quality2$width > "3",]
df_quality2
df_quality2[10,]

typeof(df_quality2[1,])
View(df_quality2)
exon_number <- unique(df_quality2$exon_number[!is.na(df_quality2$exon_number)])
length(exon_number)
exon_id <- unique(df_quality2$exon_id[!is.na(df_quality2$exon_id)])

transcript_id <- unique(df_quality2$transcript_id[!is.na(df_quality2$transcript_id)])
length(transcript_id)
gene_id <- unique(df_quality2$gene_id[!is.na(df_quality2$gene_id)])
length(gene_id)
######ENDED HERE, NEED TO REMOVE NA SOMEHOW
###THEN TAKE THE WIDTH INFO AND FILTER THE SEQUENCES BASED ON THE id AND WIDTH 
###THEN SUBMIT TO XSTREAME 
###!!!!missing 91 entries 
length(transcript_id)-dim(example)

"""Alternative approach"""

#what attributes are available 

ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

#?getSequence()
selected_df <- getSequence(id = transcript_id, type = "ensembl_transcript_id_version" , seqType = "5utr", mart=ensembl )
View(selected_df)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

transc <- selected_df$ensembl_transcript_id_version
View(unique(listAttributes(mart)))
res <- getBM(attributes = c('external_gene_name', 
                            'refseq_mrna',
                            'ensembl_gene_id', 
                            'external_transcript_name',
                            'external_gene_name',
                            'ensembl_transcript_id_version'),
             filters = 'ensembl_transcript_id_version', 
             values = transc,
             mart = mart)
dim(selected_df)
dim(res)

merge <- merge(res, selected_df, by="ensembl_transcript_id_version")
View(merge)
n_occur <- data.frame(table(df_quality2$gene_id))
n_occur
n_occur[n_occur$Freq > 1,]
View(df_quality2[df_quality2$gene_id %in% n_occur$Var1[n_occur$Freq > 1],])
length(unique(selected_df$ensembl_transcript_id_version))
length((selected_df$ensembl_transcript_id_version))
selected_df$new <- selected_df$ensembl_transcript_id_version
selected_df$UTR <- selected_df$`5utr`
drop <- c("5utr", "ensembl_transcript_id_version")

selected_df <- selected_df[,!(names(selected_df) %in% drop)]
selected_df <-  as.data.frame(selected_df)
selected_df$new

org.Hs.eg.db


BiocManager::install("EnsDb.Hsapiens.v79")

library(seqinr)
library(EnsDb.Hsapiens.v79)
#?EnsDb.Hsapiens.v79
head(keys(EnsDb.Hsapiens.v79))

EnsDb.Hsapiens.v79@tables[["gene"]][3]

symbols <- ensembldb::select(EnsDb.Hsapiens.v79, keys = c(selected_df$new), column = "SYMBOL")

symbols

dim(selected_df)

selected_df$new
merge %>% distinct(merge$`5utr`, refseq_mrna, .keep_all = T)

length(unique(merge$`5utr`))
example <- merge[!duplicated(merge$`5utr`),]
dim(example)
View(example)
#?write.fasta()
##had to manually delete sequence unavail;able in notepad to make it work 
write.fasta(sequences = as.list(example$`5utr`), names=c(as.list(example$ensembl_transcript_id_version)),as.string = F, file.out = "C:/Users/natal/Downloads/UTR5test2.fasta", nbchar =30000)



length(transcript_id)

ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
View(ensembl@filters)

attributes = listAttributes(ensembl)
View(attributes)
unique(attributes)
?mapIds
sample_ids <- read.csv(file=#path to your csv within quote marks "myfile.csv"
)

library("org.Hs.eg.db")

library(tidyverse)

my_second_df <- filter(my_first_df,
                       SampleID %in% sample_ids$SampleID)

writeXStringSet(values, "C:/Users/natal/Downloads/100only.txt", append=FALSE,
                compress=FALSE, compression_level=NA, format="fasta")
