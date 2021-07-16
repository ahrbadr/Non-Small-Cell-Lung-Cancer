
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6
# # Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# library(BiocManager)
# 
# BiocManager::install("vsn")

library(readr)
library("vsn")
library(DESeq2)



##### load the mirna-Seq data #####
#data.path="E:/Bioinformatics/Diploma bioinformatics/Nile/integrative bioinformatics/data/microRNA/gdc_download_20181218_230142.517289"


data.path="E:/Integrative/miRNA/gdc_download_20210109_041328.916761"


files <- list.files(path=data.path,recursive=T, pattern = "mirnas.quantification.txt")

# read the first file for the first time
file=files[1]
file.id=strsplit(file,"/")[[1]][1]

#open a connection to your gz file and read the file
temp <- read.table(file.path(data.path,files[1]), header=T)



#create a storing object mirna.exp to save the whole read counts of each file read in an iteration
mirna.exp=cbind(temp[1],temp[2])

rownames(mirna.exp)=mirna.exp[,1]
mirna.exp=mirna.exp[-1]
colnames(mirna.exp)=c(file.id)

for(i in 2: length(files))
{
  
  ## refer to the next file (note that we start from index 2, bec we already read the first file)
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  
  # read the next file  
  temp <- read.table(file.path(data.path,files[i]), header=T)
  
  
  ## remove the first column, bec we had it already
  temp=temp[2]
  colnames(temp)=c(file.id)
  
  mirna.exp=cbind(mirna.exp,temp)
}


#pheno <-  read_delim("E:/Bioinformatics/Diploma bioinformatics/Nile/integrative bioinformatics/data/microRNA/gdc_sample_sheet.2018-12-18.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
pheno <-  read_delim("E:/Integrative/miRNA/gdc_sample_sheet.2021-01-09.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

table(pheno$`Sample Type`)

# rename column names : replace the spaces with dots 
pheno.names=names(pheno)
names(pheno)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))
table(pheno$Sample.Type)


#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to 

file.ids=colnames(mirna.exp)

file.ids.pheno=pheno$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mirna.exp)=pheno$Sample.ID[index.files]




# for simplifying the analysis (and for time considerations) we will consider only  20 sample from each type (normal and cancer)


all.normal.samples= pheno[ pheno$Sample.Type %in% c("Solid Tissue Normal"),]$Sample.ID


all.tumor.samples= pheno[ pheno$Sample.Type %in% c("Primary Tumor"),]$Sample.ID



# now we will retrieve  the exp data for the tumor and normal samples only and also the pheno data


normal.exp=mirna.exp[, names(mirna.exp)%in% all.normal.samples]
tumor.exp=mirna.exp[, names(mirna.exp)%in% all.tumor.samples]
pheno.sub=pheno[pheno$Sample.ID %in% c(all.normal.samples, all.tumor.samples), c("Sample.ID", "Sample.Type")]




exp.sub=cbind(normal.exp,tumor.exp)
exp.sub=apply (exp.sub, 2,as.integer)
rownames(exp.sub)=rownames(normal.exp)


save(mirna.exp,pheno, exp.sub,pheno.sub ,file="microRNA-seq.RDATA")

save.image( )

###### DO the differential EXP analysis using DeSeq2
cond1="Solid Tissue Normal" 
cond2="Primary Tumor"

dds = DESeqDataSetFromMatrix( countData = exp.sub , colData = pheno.sub , design = ~ Sample.Type)
dds.run = DESeq(dds)
### direct results or specifying the contrast (to make a res object based on two specific conditions/treatment)
res=results(dds.run)
res=results(dds.run, contrast = c("Sample.Type",cond1 ,cond2) )

# remove nulls
res=res[complete.cases(res), ]
summary(res)



res.df=as.data.frame(res)
write.table(res.df, file = "res.txt")

res.dems=res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange)>log2(2),]


#### get the normalized and loggedtransformed values of all exp data
#using the the variance stabilizing transformation. vsn package
ntd=normTransform(dds)
exp.sub.norm= assay(ntd)

plotMA(res, ylim=c(-1,1)) 
summary (res)
# Make a basic volcano plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-4,4)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))



library(xlsx)
write.xlsx(res.degs, file = "res_degs.xlsx")



