

BiocManager::install("affy")
 BiocManager::install("hgu133a.db")
 BiocManager::install("affyPLM")
 BiocManager::install("hgu133plus2.db")
 BiocManager::install("hgu133acdf")
 BiocManager::install("limma")
 BiocManager::install("hgu133plus2cdf")
 BiocManager::install("GSEABase")
 BiocManager::install("GOstats")
 BiocManager::install("genefilter")
 BiocManager::install("multtest")
 BiocManager::install("pheatmap")
 install.packages("ggplot2")
 install.packages("curl")
 install.packages("RCurl")
 install.packages("readr")

# load the reqiured libraries
 library(hgu133plus2.db)
 library(hgu133acdf)
 library(limma)
 library(hgu133plus2cdf)
 library(GSEABase)
 library(GOstats)
 library(ggplot2)
 library(curl)
 library(RCurl)
library(affy)
library(readr)
library(hgu133a.db)
library(genefilter)
library(multtest)

library(affyPLM)
library(pheatmap)
 
 
 '%ni%'=Negate('%in%')
 
 #### write  function for multiple test correction
 correctPvalueandReturnAll<-function(tt.pval,method)
 {
   mt=mt.rawp2adjp(tt.pval,proc=method)
   adjp=mt$adjp[order(mt$index),]
   return(adjp[,2])
 } 
 
 
 
 
 
 ### load the CEL files
 celFilesDirectory_GDS596="cell_files"
 cels_GDS596 = list.files(celFilesDirectory_GDS596, pattern = "CEL")
 cels_GDS596
 
 affyData_GDS596 <- ReadAffy(celfile.path=celFilesDirectory_GDS596)
 affyData_GDS596
 
 
 ### data exploration
 
 # explore data
 class(affyData_GDS596)
 sampleNames(affyData_GDS596)
 featureNames(affyData_GDS596)
 head(featureNames(affyData_GDS596))
 annotation(affyData_GDS596)
 dim(affyData_GDS596)
 
 # see how the RAW expression look like without processing : notice the large values
 head(exprs(affyData_GDS596))
 #another way to explore expressions in the first 3 genes/probesin the first 5 columns
 exprs(affyData_GDS596)[1:3, 1:5]
 
 # Exploratory analysis  1- histogram
 cols=seq(1:length(sampleNames(affyData_GDS596)))
 hist(affyData_GDS596,main = "Histogram affyData_GDS596",col=cols)
 legend(12,0.9, sampleNames(affyData_GDS596),col=cols,lty=1, lwd=2,cex=0.5)
 
 # Exploratory analysis  2- box plots
 boxplot(affyData_GDS596,main = "Box Plot GDS596",col=seq(1:23))
 
 # data pre-processing in one step : life is so easy !
 # threestep (background correction, normalization, summarization from probes to probesets)
 eset = threestep(affyData_GDS596,
                  background.method = "IdealMM",
                  normalize.method = "quantile",
                  summary.method = "median.polish")
 
 ### exercise for you :  can you do the same preprocessing  using another method
 ### hint:  search for "espresso"
 
 
 # export the expression data to a text file
 write.exprs(eset,file="expData.processed.txt")
 
 #### Now you can continue where we started last time. 
 
 
 ### load already processed data
 data <- read_delim("expData.processed.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
 
 ### or work directly from the eset object
 data.expression=exprs(eset)  #### notice that this object is a matrix (not a dataframe) 
 # and saves the probe ids already as rownames... TAKE CARE ..adjust below code accordingly if 
 # you wanna use this object
 
 
 names(data)
 names(data)[1]="probe_id"
 names(data)
 
 #### mapping the probe ids into gene symbole
 hgu133a()
 mapper = hgu133aSYMBOL
 mapper
 map.df = as.data.frame(mapper)
 head(map.df)
 
 # merge the two data frames to have the symbole annotation in the data object
 data2=merge(data,map.df,by="probe_id",all.x=T)
 head(data2)
 
 # do i need the probe id again?  no, then let's drop it
 data2=data2[,-1]
 
 
 # remove nulls : some probes were not mapped to any gene symbol
 data2=data2[ ! is.na(data2$symbol),]
 
 
 
 # check duplciation of of gene symbols?  
 x=duplicated(data2$symbol)  
 sum(x)
 
 ### yes .. why ? probesets?  solutions : aggregation
 exp.data=data2[-dim(data2)[2]]
 exp.data=apply(exp.data,2, as.numeric)
 
 ####remove duplication
 exp.data.agg= aggregate(exp.data, list(data2$symbol),FUN=mean)
 names(exp.data.agg)
 
 rownames(exp.data.agg)=exp.data.agg$Group.1
 exp.data.agg=exp.data.agg[- 1]
 
 
 ## log transformation of data
 exp.data.agg.log <- log2(exp.data.agg)
 
 #save the object in a RDATA file
 save(exp.data.agg,file="processed.RDATA")
 
 #### do the DIFF EXP Analysis
 exp=exp.data.agg
 case.indecies=c(1:8)
 ctrl.indecies=c(9:dim(exp)[2]) #### or  seq(from=9,to=dim(exp)[2])
 
 ## calculating LFC
 lfc.diff=apply(exp,1, function(x)  mean(x[case.indecies]) -mean(x[ctrl.indecies]))
 x=as.data.frame(lfc.diff)
 
 ## calcualting p values
 f=factor( c( rep(1, length(case.indecies)) , rep(2, length(ctrl.indecies)) ))
 t.pval=rowttests(as.matrix(exp),f)$p.value
 t.pval.adj=correctPvalueandReturnAll(t.pval,"BH")
 
 res=cbind(lfc.diff,t.pval,t.pval.adj)
 
 #####  selection criteria for identifying DEGS
 degs.res=res[t.pval<0.05,]  ##### identify DEGs based on teh significance level only
 degs.res=res[abs(lfc.diff) > log2(2),]  ##### identify DEGs based on the LFC only
 degs.res=res[abs(lfc.diff) > log2(2)   & t.pval <0.05,]  ##### identify DEGs based on both  LFC and the significane level 
 
 degs.genes= rownames(degs.res)
 
 #export them for further analysis in DAVID
 write.table(degs.genes,file = "DEGs.txt",row.names = F,col.names = F,quote = F)
 
 ### creating a heatmap
 # 1- get the expression profiles of the degs only
 
 exp.degs=exp[rownames(exp) %in% degs.genes, ]
 dsm=exp.degs
 
 colnames=colnames(dsm)
 case.vector=rep("Case", length(case.indecies))
 ctrl.vector=rep("Ctrl", length(ctrl.indecies))
 Sample=c( case.vector, ctrl.vector)
 
 annotation=as.data.frame(Sample)
 rownames(annotation)=colnames
 annotation
 
 # Specify colors
 Sample = c("lightgreen", "navy")
 names(Sample) = c("Ctrl", "Case")
 ann_colors = list(Sample = Sample)
 
 m2=scale(t(dsm),center=T,scale=T)
 
 m2=t(m2)
 
 pheatmap(m2, annotation = annotation, annotation_colors = ann_colors,fontsize_row = 5,fontsize_col = 9)
 
 
 
 
 
 
 
 
