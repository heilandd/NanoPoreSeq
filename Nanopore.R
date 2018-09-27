####Full Pipeline for Nanopor RNA-seq

#Define new Class
library(DESeq2)
library(rapport)
library(RColorBrewer)
library(AutoPipe)
#### Full Pipe 
n_por_seqPIPE=function(path,Samples_discription){
  
  #Einfügen der Samole Beschreibung
  Set=Poreseq(Samples_discription)
  
  #Angeben des Ortes
  Set@path=path
  
  # Zusammenfügen der fasq files
  GET_FASQ(Set)
  
  # Alignment der Sequences
  TRIM_Barcodes(Set)
  
  # Alignment der Sequences
  Set=NANOPORE_Aligner(Set)
  
  #Analyze Data
  Set=Analyzer(Set,Write_exp=T,filter_genes=5,MA_PLOT=T)
  
  
}

Poreseq <- setClass("Poreseq", slots = c(Samples_discription = "data.frame", path = "character", TRIM="character", Out_Counts = "list", counts="matrix", norm_Exp = "matrix", tsne = "data.frame", dds="DESeqDataSet", diff_exp="DESeqResults"))
#create input the files into the new class
setValidity("Poreseq", function(object) {
  msg <- NULL
  if ( ! is.data.frame(object@Samples_discription) ){
    msg <- c(msg, "input data must be data.frame")
  }else if ( nrow(object@Samples_discription) < 2 ){
    msg <- c(msg, "input data must have more than one row")
  }else if ( ncol(object@Samples_discription) < 2 ){
    msg <- c(msg, "input data must have more than one column")
  }
  if (is.null(msg)) TRUE
  else msg
})
setMethod("initialize",signature = "Poreseq", definition = function(.Object, Samples_discription ){
  .Object@Samples_discription <- Samples_discription
  validObject(.Object)
  return(.Object)
})
setGeneric("GET_FASQ", function(object) standardGeneric("GET_FASQ"))
setMethod("GET_FASQ",signature = "Poreseq", definition = function(object){
  
  path=object@path
  Samples_discription=object@Samples_discription
  
  print="################################## Start the Auto Nanopore Pipe ###################################"
  
  print(paste("#########################", nrow(Samples_discription), "Samples will be analyzed #############"))
  
  print(" ######################## Start to find fasq data ########################")
  
  setwd(path)
  ordner=dir()[file.info(dir())$isdir]
  
  cnames_fail=paste(ordner,"_combfail",sep="")
  cnames_pass=paste(ordner,"_combpass",sep="")
  cnames_comb=paste(ordner,"_combcomb",sep="")
  
  list_combfiles=c()
  list_combfiles_path=c()
  
  for (io in 1: length(ordner)){
    
    print(paste("############ Start with folder ->", ordner[io], "   ######################"))
    setwd(paste(path,"/",ordner[io], sep=""))
    subfile=dir()    
    test=subfile %in% c("fail", "pass")
    test=test[!duplicated(test)]
    if(test==F){print(" Error Structure of files is not corret, pleas correct the orginal data from Nanopore")}
    else{
      
      ##Start with pass and fail data 
      path_fail=paste(path,"/",ordner[io], "/fail", sep="")
      setwd(path_fail)
      system(paste("cd ", path_fail, sep="" ))
      system(paste("mkdir ", "Combined_file", sep="" ))
      system(paste("cat ", "*fastq",">", "combine.fastq", sep="" ))
      file.rename("combine.fastq", "Combined_file/combine_fail.fastq")
      
      #Now for pass reads
      path_fail=paste(path,"/",ordner[io], "/pass", sep="")
      setwd(path_fail)
      system(paste("cd ", path_fail, sep="" ))
      system(paste("mkdir ", "Combined_file", sep="" ))
      system(paste("cat ", "*fastq",">", "combine.fastq", sep="" ))
      file.rename("combine.fastq", "Combined_file/combine_pass.fastq")
      
      #create the comb file
      
      setwd(paste(path,"/",ordner[io], sep=""))
      system(paste("mkdir ", "comb", sep="" ))
      path_fail=paste(path,"/",ordner[io], "/fail/Combined_file", sep="")
      system(paste("cd ", path_fail))
      file.rename(paste(path_fail,"/combine_fail.fastq", sep=""), paste(path,"/",ordner[io],"/comb/",cnames_fail[io],".fastq", sep=""))
      
      path_fail=paste(path,"/",ordner[io], "/pass/Combined_file", sep="")
      system(paste("cd ", path_fail))
      file.rename(paste(path_fail,"/combine_pass.fastq", sep=""), paste(path,"/",ordner[io],"/comb/",cnames_pass[io],".fastq", sep=""))
      
      #combine fail and pass
      
      path_n=paste(path,"/",ordner[io], "/comb", sep="")
      
      setwd(path_n)
      system(paste("cd ", path_n, sep="" ))
      system(paste("mkdir ", "Combined_file", sep="" ))
      system(paste("cat ", "*fastq",">", "combine.fastq", sep="" ))
      file.rename("combine.fastq", paste(path_n,"/Combined_file/",cnames_comb[io],".fastq", sep="") )
      
      list_combfiles=c(list_combfiles, paste(path_n,"/Combined_file/",cnames_comb[io],".fastq", sep=""))
      list_combfiles_path=c(list_combfiles_path, paste(path_n,"/Combined_file/", sep=""))
      print("################## FASTQ files were combined #########################################")
      
      
      
    }
    
    
    
  }
  
  
  #make new files
  setwd(path)
  system(paste("cd ", path, sep="" ))
  system(paste("mkdir ", "Analysis", sep="" ))
  setwd(paste(path,"/Analysis", sep=""))
  system(paste("cd ", paste(path,"/Analysis", sep=""), sep="" ))
  system(paste("mkdir ", "FASTQ", sep="" ))
  
  #Copy combined files into fasq folder
  
  file.copy(list_combfiles, to=paste(path,"/Analysis/FASTQ", sep=""))
  system(paste("cd ", paste(path,"/Analysis/FASTQ", sep=""), sep="" ))
  setwd(paste(path,"/Analysis/FASTQ", sep=""))
  system(paste("cat ", "*fastq",">", "combine.fastq", sep="" ))
  
  
})


setGeneric("NANOPORE_Aligner", function(object) standardGeneric("NANOPORE_Aligner"))
setMethod("NANOPORE_Aligner",signature = "Poreseq", definition = function(object){
  
  print(paste("########## Go with Alignment ########################"))
  
  Samples_discription=object@Samples_discription
  path=paste(object@path,"/Analysis", sep="")
  #remove non trimmer reads
  setwd(paste(path,"/TRIM_FASQ", sep=""))
  
  rm_files=dir()[dir() %in% Samples_discription$File==F]
  
  file.remove(rm_files)
  files=dir()
  
  
  
  
  ######Alignment
  #Alignment with minimap2
  setwd(path)
  system(paste("mkdir ", "Minimap2", sep="" ))
  setwd(paste(path,"/Minimap2", sep=""))
  ref_fa="/media/daka/Data/RNA_Seq/ref.fa"
  
  for (i in 1:length(files)){
    
    print(paste("Sample: ",files[i], " will be aligned", sep=""))   
    outname=paste(path, "/Minimap2/",paste(gsub(".fastq","",files[i]),".sam", sep=""), sep="" ) 
    inputfile=paste(path, "/TRIM_FASQ/",files[i], sep="" ) 
    
    system(paste("cd ", "
                 ","cd minimap2","
                 ", "./minimap2 -ax splice", ref_fa, " ",inputfile, " > ", outname), show.output.on.console = F)
  }
  
  
  #samtools
  #samtools view -S -b sample.sam > sample.bam
  
  
  
  for (i in 1:length(files)){
    
    
    inputfile=paste(path, "/Minimap2/",paste(gsub(".fastq","",files[i]),".sam", sep=""), sep="" ) 
    outputfile=paste(path, "/Minimap2/",paste(gsub(".fastq","",files[i]),".bam", sep=""), sep="" )
    Fusion=paste(path, "/Minimap2/",paste(gsub(".fastq","",files[i]),"_fusion.sam", sep=""), sep="" )
    
    system(paste("samtools view -S -b  ",inputfile, " > ", outputfile))
    
    system(paste("samtools view -f 0x800 ",outputfile, " > " , Fusion))
    
    
  }
  
  path_analyze=paste(path, "/Minimap2/", sep="" ) 
  setwd(paste(path,"/Minimap2", sep="")) 
  
  Out_counts=lapply(1:length(files), function(i){
    
    #Read counts
    sam=paste(path_analyze,paste(gsub(".fastq","",files[i]),".sam", sep=""), sep="" ) 
    bam=paste(path_analyze,paste(gsub(".fastq","",files[i]),".bam", sep=""), sep="" ) 
    bam.sort=paste(path_analyze,paste(gsub(".fastq","",files[i]),"SORT", sep=""), sep="" ) 
    bam.sort2=paste(path_analyze,paste(gsub(".fastq","",files[i]),"SORT.bam", sep=""), sep="" ) 
    
    readcount.txt=paste(path_analyze,paste(gsub(".fastq","",files[i]),"readcount.txt", sep=""), sep="" ) 
    
    #system(paste("samtools view ", bam.index, sep=""))
    system(paste("samtools flagstat ", bam, sep=""))
    #system(paste("htseq-count -f sam ", sam, gtf, " > counts.file", sep=""))
    print("---------- Sort and Index files---------------------------------")
    system(paste("samtools sort ", bam," ", bam.sort,  sep=""))
    system(paste("samtools index ", bam.sort2,  sep=""))
    #system(paste("bam-readcount -f ", ref_fa, " ", bam.sort, " -w 1 > readcount.txt", sep=" "))
    print("---------- Generate Counts ---------------------------------")
    system(paste("samtools idxstats ", bam.sort2, " > ",readcount.txt, sep=" "))
    print("---------- Generate Counts -> Align Gene Regions ---------------------------------")
    
    #Load data to back to R
    reads=read.table(readcount.txt)
    names(reads)=c("ENST","Raw-Reads", "Mapped", "Unmapped")
    reads$genes=apply(reads, 1, function(x){return(unlist(strsplit(as.character(x[1]), "[.]"))[1])})
    reads=reads[reads$Mapped>=1, ]
    library(annotables)
    libraryx=grch38_tx2gene
    library2=grch38
    HUGO=lapply(1:nrow(reads), function(i2){
      ent=as.character(reads[i2, ]$genes )
      #print(ent)
      getEN=libraryx[libraryx$enstxp==ent, ]$ensgene
      
      #print(i2)
      return(getEN[1])
    })
    HUGOx=as.data.frame(do.call(rbind, HUGO))
    reads$genes=HUGOx$V1
    reads$HUGO=apply(reads,1,function(x){
      yy=library2[library2$ensgene==x[5], ]$symbol 
      return(yy[1])
    })
    #Remove duplicates
    
    reads=reads[!is.na(reads$HUGO), ]
    
    genes=unlist(reads[!duplicated(reads$HUGO), ]$HUGO)
    
    reads_new=data.frame(do.call(rbind, lapply(1:length(genes), function(ix){
      sym=genes[ix]
      n_reads=sum(reads[reads$HUGO==sym, ]$Mapped)
      out=data.frame(sym,n_reads)
      names(out)=c("HUGO", "Counts")
      return(out)
    })))
    
    
    
  })
  saveRDS(Out_counts, file = "Out_counts.R")
  
  
  object@Out_Counts=Out_counts
  
  
  
  
  })


setGeneric("TRIM_Barcodes", function(object) standardGeneric("TRIM_Barcodes"))
setMethod("TRIM_Barcodes",signature = "Poreseq", definition = function(object){
  
  print(paste("########## Go with analyis ########################"))
  
  path=paste(object@path,"/Analysis", sep="")
  
  
  
  # Start with porchop
  
  #Seperate barcodes into new fastq files
  
  setwd(path)
  system(paste("mkdir ", "TRIM_FASQ", sep="" ))
  
  
  system("cd 
         cd Porechop ")
  input=paste(path, "/FASTQ/combine.fastq", sep="" )
  output=paste(path, "/TRIM_FASQ", sep="" )
  
  print(paste("########## Porechop Analysis  ########################"))
  
  system(paste("porechop ","-i ",input," -b ", output, " --threads 10", sep=""),show.output.on.console = F)
  
  Set@TRIM="DONE"
  return(print("Done"))
  
})

setGeneric("Analyzer", function(object,filter_genes=5,MA_PLOT=T) standardGeneric("Analyzer"))
setMethod("Analyzer",signature = "Poreseq", definition = function(object,Write_exp=T, filter_genes=5,MA_PLOT=T){
  
  print(paste("########## Analyze Counts ########################"))
  
  Out_counts=object@Out_Counts
  Samples_discription=object@Samples_discription
  samples=nrow(Samples_discription)
  All_genes=unlist(lapply(1:samples,function(i){Out_counts[[i]]$HUGO}))
  All_genes=as.data.frame(as.character(All_genes[!duplicated(All_genes)]))
  names(All_genes)="genes"
  for(i in 1:length(files)){
    l=Out_counts[[i]]
    l[,1]=as.character(l[,1])
    All_genes[,i+1]=0
    for(ix in 1:nrow(All_genes)){
      count=sum(l[l[,1]==as.character(All_genes[ix,1]),2])
      All_genes[ix,i+1]=count
      
    }
    
  }
  input=All_genes
  
  
  input_m=as.matrix((input[2:(samples+1)]))
  rownames(input_m)=input[,1]
  class(input_m)="integer"
  colnames(input_m)=Samples_discription$Sample
  
  object@counts=input_m
  
  condition=factor(Samples_discription$treatment)
  dds <- DESeqDataSetFromMatrix(countData = input_m ,DataFrame(condition), ~ condition)
  
  saveRDS(dds, file = "DDS.R")
  object@dds=dds
  
  ####Filter####
  dds <- dds[ rowSums(counts(dds)) >= filter_genes, ]
  ##### Normalize #######
  rld <- rlog(dds, blind=FALSE)
  head(assay(rld), 3)
  Expression_data=(assay(rld))
  colnames(Expression_data)=Samples_discription$Sample
  
  object@norm_Exp=Expression_data
  
  if(Write_exp==T){write.table(Expression_data, file="Expression_data.csv", sep=";")}
  
  dds=DESeq(dds)
  res=results(dds)
  
  object@diff_exp=res
  
  
  if(MA_PLOT==T){
    plotMA(res, bty="n")
    res=data.frame(res)
    plot(res$log2FoldChange, -log(res$padj), bty="n", pch=19, cex=0.5)
    top_genes=rownames(res[order(res$pvalue, decreasing = F), ])[1:20]
    text(res[top_genes, ]$log2FoldChange, -log(res[top_genes, ]$padj), labels = top_genes)
    
    
  }
  
  
  
  top_genes2=c(rownames(res[order(res$log2FoldChange, decreasing = F), ])[1:10],rownames(res[order(res$log2FoldChange, decreasing = T), ])[1:10])
  heatmap(Expression_data[top_genes2, ], col=brewer.pal(11, "RdBu"))
  
  
  heatmap(cor(t(Expression_data[top_genes2, ])), col=brewer.pal(11, "Spectral"))
  
  
  
  data_out=Expression_data
  
  
  
  
  
  
  
  new = clusterProfiler::bitr(rownames(data_out), fromType = "SYMBOL", 
                              toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  new = new[!duplicated(new$ENTREZID), ]
  rownames(new) = new$ENTREZID
  data_out = data_out[new$SYMBOL, ]
  rownames(data_out) = new$ENTREZID
  me=data_out
  
  
  
  ## calculate best number of clusters 
  res<-TopPAM(me, max_clusters = 2, TOP=1000)
  
  me_TOP=res[[1]]
  number_of_k=res[[3]]
  
  ## Compute top genes of each cluster, with "TRw" samples with a negative Silhouette widths could be cut-off
  
  File_genes=Groups_Sup(me_TOP, me=me, number_of_k,TRw=-2)
  
  groups_men=File_genes[[2]]
  me_x=File_genes[[1]]
  
  groups_men[as.character(Samples_discription$Sample), ]$cluster=Samples_discription$treatment
  groups_men=groups_men[order(groups_men$cluster, decreasing = F), ]
  
  Supervised_Cluster_Heatmap(groups_men = groups_men, gene_matrix=me_x, db="c2", method="PAMR",TOP_Cluster=100,genes_to_print=10, show_sil=F,print_genes=TRUE, TOP=1000,GSE=TRUE,plot_mean_sil=F)
  
  
  
})






# Easy Workflow for Nanopore Alignment

Samples_discription=read.csv("samp_desc.csv", row.names=1)

#Einfügen der Samole Beschreibung
Set=Poreseq(Samples_discription)

#Angeben des Ortes
Set@path=c("path_to_file")

# Zusammenfügen der fasq files
GET_FASQ(Set)

# Alignment der Sequences
TRIM_Barcodes(Set)

# Alignment der Sequences
Set=NANOPORE_Aligner(Set)

#Analyze Data
Set=Analyzer(Set,Write_exp=T,filter_genes=5,MA_PLOT=T)


