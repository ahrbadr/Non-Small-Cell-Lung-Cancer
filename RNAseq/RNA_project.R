###########################################
#  -Integrative bioinformatics            # 
#  - RNA-Seq                              #
#  - 2018- 12- 18                         #
#  - Copyright: Mohamed Hamed             #
###########################################
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6
# # Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
BiocManager::install("vsn")
library(readr)
library(org.Hs.eg.db)
library(vsn)
library(DESeq2)
library(pheatmap)

##### load the mRNA-Seq data #####
data.path="gdc_download_20210201_073404.771248"
files <- list.files(path=data.path,recursive=T, pattern = "gz")


# read the first file for the first time
file=files[1]
file.id=strsplit(file,"/")[[1]][1]
#open a connection to your gz file and read the file
gz.con=gzfile(file.path(data.path,files[1]))
temp <- read.table(gz.con, header=F)


#create a storing object mrna.exp to save the whole read counts of each file read in an iteration
mrna.exp=temp
rownames(mrna.exp)=mrna.exp[,1]
mrna.exp=mrna.exp[-1]
colnames(mrna.exp)=c(file.id)


# on other files
for(i in 3: length(files))
{
  ## refer to the next file (note that we start from index 2, bec we already read the first file)
  file=files[i]
  file.id=strsplit(file,"/")[[1]][1]
  
  # read the next file  
  gz.con=gzfile(file.path(data.path,files[i]))
  temp <- read.table(gz.con, header=F)
  
  ## remove the first column, bec we had it already
  temp=temp[-1]
  colnames(temp)=c(file.id)
  
  mrna.exp=cbind(mrna.exp,temp)
}
View(mrna.exp)

### do the mapping of ensembel.id to gene symbol###############

# prepare the ensembel id to be like the one in the database
ensemble.id=sapply(rownames(mrna.exp), function(x) strsplit(as.character(x),"\\.")[[1]][1])
View(ensemble.id)
mrna.exp=cbind(ensemble.id,mrna.exp)


mapper<- mapIds(org.Hs.eg.db, keys=ensemble.id, column="SYMBOL",keytype="ENSEMBL", multiVals="first") 
mapper.df=as.data.frame(mapper)
mapper.df=cbind(rownames(mapper.df), mapper.df)
names(mapper.df)=c("ensemble.id","symbol")


mrna.exp2=merge(mrna.exp,mapper.df,by="ensemble.id",all.x=T) 


# drop the first column (ensemble.id)
mrna.exp2=mrna.exp2[-1]
mrna.exp2=mrna.exp2[ ! is.na(mrna.exp2$symbol),]


# check duplciation of of gene symbols?  
x=duplicated(mrna.exp2$symbol)  
sum(x)

### yes .. why ? transcripts?  solutions : aggregation
mrna.exp.data=mrna.exp2[-dim(mrna.exp2)[2]]
mrna.exp.data=apply(mrna.exp.data,2, as.numeric)
####remove  duplication by aggregation
mrna.exp.data.agg= aggregate(mrna.exp.data, list(mrna.exp2$symbol),FUN=mean)


rownames(mrna.exp.data.agg)=mrna.exp.data.agg$Group.1
mrna.exp.data.agg=mrna.exp.data.agg[-1]
file.ids=colnames(mrna.exp.data.agg)


###### load the mrna sample sheets
pheno <- read_delim("gdc_sample_sheet.2021-02-01.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
View(pheno)
table(pheno$`Sample Type`)


# rename column names : replace the spaces with dots 
pheno.names=names(pheno)
names(pheno)= as.character( sapply ( pheno.names, function(x) gsub(" ",".",x)))
table(pheno$Sample.Type)


#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to 

file.ids.pheno=pheno$File.ID
index.files=match(file.ids,file.ids.pheno)
names(mrna.exp.data.agg)=pheno$Sample.ID[index.files]


 pheno=pheno[index.files,]
 names(mrna.exp.data.agg)= pheno$Sample.ID


# for simplifying the analysis (and for time considerations) we will consider only  20 sample from each type (normal and cancer)

all.normal.samples= pheno[ pheno$Sample.Type %in% c("Solid Tissue Normal"),]$Sample.ID 
normal.samples_ID=all.normal.samples[1:length(all.normal.samples)]

all.tumor.samples= pheno[ pheno$Sample.Type %in% c("Primary Tumor"),]$Sample.ID
tumor.samples_ID=all.tumor.samples[1:length(all.tumor.samples)]

# now we will retrieve  the exp data for these 20 samples only and also the pheno data
normal.exp=mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% normal.samples_ID]
tumor.exp=mrna.exp.data.agg[, names(mrna.exp.data.agg)%in% tumor.samples_ID]

pheno.sub=pheno[pheno$Sample.ID %in% c(normal.samples_ID,tumor.samples_ID), c("Sample.ID", "Sample.Type")]


exp.sub=cbind(normal.exp,tumor.exp)
exp.sub=apply (exp.sub, 2,as.integer)
rownames(exp.sub)=rownames(normal.exp)


save(mrna.exp.data.agg,pheno, exp.sub,pheno.sub ,file="RNA-seq.RDATA")


###### DO the differential EXP analysis using DeSeq2
cond1="Solid Tissue Normal" 
cond2="Primary Tumor"

dds = DESeqDataSetFromMatrix( countData = exp.sub , colData = pheno.sub , design = ~ Sample.Type)
dds.run = DESeq(dds)
### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)
res=results(dds.run)
res=results(dds.run, contrast = c("Sample.Type",cond1 ,cond2) )

# remove nulls
res=res[complete.cases(res), ]
summary(res)
res.df=as.data.frame(res)


plotMA(res, ylim=c(-1,1)) 
summary (res)

res.degs=res.df[res.df$padj< 0.01 & abs(res.df$log2FoldChange)>log2(2),]
#expression of these degs for plots
exp.degs= exp.sub[rownames(exp.sub) %in% rownames(res.degs), ]


#### get the normalized and loggedtransformed values of all exp data
#using the the variance stabilizing transformation. vsn package
ntd=normTransform(dds)
exp.sub.norm= assay(ntd)


# Make a basic volcano plot
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res.degs, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res.degs, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

#heatmap
### creating a heatmap
# 1- get the expression profiles of the degs only

colnames=colnames(exp.degs) #ht3rfo b3d kda 3mlt leh kda
case.vector=rep("Case", length(tumor.samples_ID))
ctrl.vector=rep("Ctrl", length(normal.samples_ID))
Sample=c( case.vector, ctrl.vector) 

annotation=as.data.frame(Sample)
rownames(annotation)=colnames
annotation

# Specify colors
Sample = c("lightgreen", "navy")
names(Sample) = c("Ctrl", "Case")
ann_colors = list(Sample = Sample) #????? sample 2li fo2 w 2li t7t

m2=scale(t(exp.degs),center=T,scale=T)

m2=t(m2)

pheatmap(m2, annotation = annotation, annotation_colors = ann_colors,fontsize_row = 5,fontsize_col = 9)##### ??????bytrsm m3aha alscale 


View(res.degs)
write.table(res.degs,file = "DEGs.csv")
