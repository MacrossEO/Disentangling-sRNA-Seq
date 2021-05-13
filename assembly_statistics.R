#02/may/18
rm(list = ls())
suppressMessages(require(protr))

computer="production"
if(computer=="test"){
  md = "/Users/macross/Documents/data/sRNAs_assembly/assembly_eval/hp_mm/mapped/hpmm"
  dir = paste(md,"assemblies/trim_chimera",sep="/") #contigs folder
  all_fq_file = paste(md,"all.fq",sep="/");all_fq_file
  output = paste(dir,"assembly_statistics.tab",sep="/")
  bowtie_bin = "~/Documents/bin/bowtie-1.2.2/bowtie"
  ref_genome = "/Users/macross/Documents/data/sRNAs_assembly/hpmm_genome/hpmm"
  ref_genome_len = 3480323312#hpmm
  threads=7
}
if(computer=="production"){
  args = commandArgs(trailingOnly=TRUE)
  md =  args[1]
  dir = paste(md,args[2],sep="/")#contigs folder
  contigs_db = paste(md,args[3],sep="/")
  output = paste(md,args[4],sep="/")
  all_fq_file = paste(md,args[5],sep="/");#all_fq_file
  ref_genome = args[6]#"/LUSTRE/usuario/oramirez/assembly_eval/at_bc/at_genome/at_genome"
  ref_genome_len = as.numeric(args[7])
  threads = as.numeric(args[8])
  bowtie_bin = args[9]
  
  #ref_genome_len = 121182535#at
  #ref_genome_len = 43272801#bc
  #ref_genome_len = 164455336#atbc
  #ref_genome_len = 121182535#at
  #ref_genome_len = 471217518#cc
  #ref_genome_len = 592400053#atcc
  #ref_genome_len = 121182535#at
  #ref_genome_len = 79367360#ha
  #ref_genome_len = 200549895#atha
  #ref_genome_len = 703936975#hp
  #ref_genome_len = 2776391127#mm
  #ref_genome_len = 3480328102#hpmm
  #ref_genome_len = 2554703818#mu
  #ref_genome_len = 65897716#ls
  #ref_genome_len = 2620601534#muls
}

#1-9)contigs lenght
ext <- ".fa$"
files <- list.files(path =dir, pattern = ext, full.names = T, recursive = F);files
len = length(files) 
zeros = rep(0,len)
statistics = matrix(nrow=len,ncol=13);statistics
names = array()
i=1
for (i in 1:len){
  file = files[i]
  bsfile = basename(file); print(bsfile)
  names[i] = bsfile
  file_fa = readFASTA(file)
  statistics[i,1] = length(names(file_fa))
  file_len = unlist(lapply(file_fa,nchar))
  statistics[i,c(2:7)] = summary(file_len);statistics
  
  #contigsCov
  cov_contig_sum = sapply(names(file_fa),function(x)unlist(strsplit(x,","))[2])
  cov_contig_sum
  statistics[i,9] = sum(as.numeric(cov_contig_sum)) 
}

rownames(statistics) = names

#8)contigsCoverage
for(i in 1:len){
  file = files[i]
  bsfile = basename(file); print(bsfile)
  
  sam_file = paste(contigs_db,".sam",sep="")
  bowtie_call = paste(bowtie_bin,"-S -f -v 0 -p",threads,contigs_db,all_fq_file,sam_file,"2>&1");print(bowtie_call)
  bowtie_res = system(bowtie_call,intern=T);print(bowtie_res)
  
  a = unlist(strsplit(bowtie_res[2],"[(]"))
  b = unlist(strsplit(a[2],"[%]"))
  statistics[i,8]=as.double(b[1]);statistics
  
  rm_call=paste("rm ",dir,"/*.ebwt",sep = "");rm_call
  system(rm_call)
  system(paste("rm ",sam_file,sep = ""))
}

cp.statistics = statistics
statistics = cp.statistics

#10) GenomeMap (mapping contigs to reference genome)
i=1
for(i in 1:len){
  file = files[i]
  bs = basename(file);print(bs)
  sam_file = sub(pattern = ".fa",".sam",file)
  
  bowtie_call=paste(bowtie_bin,"-S -v 0 -f -p",threads,ref_genome,file,sam_file," 2>&1");print(bowtie_call)
  bowtie_res = system(bowtie_call,intern = T);print(bowtie_res)
  
  a = unlist(strsplit(bowtie_res[2],"[ ]"));a
  statistics[i,10] = as.double(a[9]);statistics
  
  system(paste("rm",sam_file))
}
cp.statistics = statistics

#11) Genome_Mapping% = ((number of aligned contigs * median length)/genome_len)*100
genome_Mapping_f = function(x){
  return( ((x[10]*x[4])/ref_genome_len) *100 )
}
genomeMapping = apply(statistics,1,genome_Mapping_f)
statistics[,11] = t(genomeMapping)

#12)Chimeras (number of NOT alignment contigs to genome/number of contigs)
statistics[,12] = statistics[,1]-statistics[,10] 

#13)Chimeras% 
statistics[,13] = (statistics[,11]*100)/statistics[,1]
statistics

colnames(statistics) = c("contigsNum","Min","1stQu","Median","Mean","3rdQu","Max","libMapping","contigsCov","GenomeMap",
                         "GenomeMap%","Chimeras","Chimeras%")#15
write.table(statistics, file = output, quote = FALSE, sep = "\t", col.names = TRUE,row.names = T)

