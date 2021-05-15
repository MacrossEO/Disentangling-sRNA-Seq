#21/dic/2020
#find and add origin (host/par/amb) to fasta sequences
rm(list=ls())
suppressMessages(library(seqinr))

computer="production"
if(computer=="production"){
  setwd("/home/macross/compartida/data/miRNA/mm_hp")
  
  input = "dea_contigs/contigs.fa" 
  host_genome = "/home/macross/compartida/data/miRNA/mmusculus_genome/Mus_musculus"
  par_genome = "/home/macross/compartida/data/miRNA/hp_genome/nHp.2.0"
  host_genome_bt2 = "/home/macross/compartida/data/miRNA/mmusculus_genome/bt2/mm_genome"
  par_genome_bt2 = "/home/macross/compartida/data/miRNA/hp_genome/bt2/hp_genome"
  
  output = "dea_contigs/contigs_tag.fa"
  threads=7
  
  bowtie_bin = "~/Documentos/bin/bowtie-1.2.2/bowtie"
  bowtie2_bin = "~/Documentos/bin/bowtie2-2.4.2-linux-x86_64/bowtie2"
}
if(computer=="production"){
  args = commandArgs(trailingOnly=TRUE)
  input = args[1]
  host_genome = args[2]
  par_genome = args[3]
  host_genome_bt2 = args[4]
  par_genome_bt2 = args[5]
  
  output = args[6]
  threads = args[7]
  
  bowtie_bin = args[8]
  bowtie2_bin = args[9]
}

chimera_assignment = function(chimeras_fa,chimeras_tab){
  #align chimeras to host
  host_sam_file = paste(tmpdir,"host.sam",sep="/")
  bowtie2_call = paste(bowtie2_bin,"-p",threads,"-f -x",host_genome_bt2,"-U",chimeras_fa,"-S",host_sam_file);print(bowtie2_call)
  system(bowtie2_call)
  
  #align chimeras to parasite
  par_sam_file = paste(tmpdir,"par.sam",sep="/")
  bowtie2_call = paste(bowtie2_bin,"-p",threads,"-f -x",par_genome_bt2,"-U",chimeras_fa,"-S",par_sam_file);print(bowtie2_call)
  system(bowtie2_call)
  
  #get chimeras
  chimeras = read.fasta(chimeras_fa,as.string=T,forceDNAtolower=F)
  chimeras = data.frame(contig_id=names(chimeras),tag="a")
  
  #assign chimeras
  i=1
  host_par = array()
  for(i in 1:nrow(chimeras)){
    id = chimeras$contig_id[i]
    
    #to host
    grep_call = paste("grep ",id,host_sam_file);grep_call
    sam_row = system(grep_call,intern = T);sam_row
    if(length(grep("XM:i",sam_row))==0){
      host_par[1]= -1
    }else{
      grep_call = paste("grep",id,host_sam_file, "| sed 's/^.*NM:i://' | awk '{print $1}'")
      host_par[1] = as.numeric(system(grep_call,intern = T))
    }
    
    #to par
    grep_call = paste("grep ",id,par_sam_file);grep_call
    sam_row = system(grep_call,intern = T);sam_row
    if(length(grep("XM:i",sam_row)) == 0){
      host_par[2] = -1
    }else{
      grep_call = paste("grep",id,par_sam_file, "| sed 's/^.*NM:i://' | awk '{print $1}'")
      host_par[2] = as.numeric(system(grep_call,intern = T))
    }
    #print(host_par)
    
    if(host_par[1]==host_par[2]){
      chimeras$tag[i] = "a"
    }else{
      idx = which.max(host_par)
      if(idx==1){
        if(host_par[1]>10){
          chimeras$tag[i] = "a"
        }else{
          chimeras$tag[i] = "h"#host
        }
      }else{
        if(host_par[2]>10){
          chimeras$tag[i] = "a"
        }else{
          chimeras$tag[i] = "p"#par
        }
      }
    }
  }
  write.table(chimeras, file = chimeras_tab, quote = FALSE, sep = "\t", col.names = F,row.names = F)
  return(chimeras)
}

tmpdir = paste(dirname(output),"tmp_dir",sep="/")#switch when necessary
dir.create(tmpdir)

#intermediate_files
fa_file = input
base_name = unlist(strsplit(basename(fa_file),"[.]"))[1]
map_to_host_file = paste(tmpdir,paste(base_name,"_mHost.fa",sep=""),sep="/")
map_to_par_file = paste(tmpdir,paste(base_name,"_mPar.fa",sep=""),sep="/")
sam_file = paste(tmpdir,"lib.sam",sep="/")

#map to host 
bowtie_call = paste(bowtie_bin,"-S -v 0 -f -p",threads,"-f",host_genome,fa_file,sam_file,"--al ",map_to_host_file);print(bowtie_call)
system(bowtie_call)
awk_call = paste("grep \">\" ", map_to_host_file ," | awk -F',' '{print $1}' | sed 's/[>]//' 2>&1 ");print(awk_call)
map_to_host = system(awk_call,intern=T)

#map to parasite
bowtie_call = paste(bowtie_bin,"-S -v 0 -f -p",threads,"-f",par_genome,fa_file,sam_file,"--al ",map_to_par_file);print(bowtie_call);
system(bowtie_call)
awk_call = paste("grep \">\" ", map_to_par_file ," | awk -F',' '{print $1}' | sed 's/[>]//' 2>&1 ");print(awk_call)
map_to_par = system(awk_call,intern=T)

#intersection 
ambiguous_reads = intersect(map_to_host,map_to_par)

#get ids
fasta_file =  read.fasta(input,as.string=T,forceDNAtolower=F)
fasta_ids = data.frame(matrix(unlist(strsplit(names(fasta_file),",")),ncol=3,byrow=T))
colnames(fasta_ids) = c("id","counts","length")
fasta_ids$id = as.character(fasta_ids$id)

fasta_ids$category = "chimera"
fasta_ids$category[fasta_ids$id %in% map_to_host] = "h"
fasta_ids$category[fasta_ids$id %in% map_to_par] = "p"
fasta_ids$category[fasta_ids$id %in% ambiguous_reads] = "a"
table(fasta_ids$category)

#bt2
if(sum(fasta_ids$category=="chimera")>0){
  chimeras_fa = paste(dirname(output),"chimeras.fa",sep="/")
  chimeras_idx = which(fasta_ids$category=="chimera")
  chimeras_tab = paste(dirname(output),"chimeras.tab",sep="/")
  
  write.fasta(fasta_file[chimeras_idx],names(fasta_file[chimeras_idx]),chimeras_fa)
  chimeras = chimera_assignment(chimeras_fa,chimeras_tab)
  chimeras = cbind(chimeras,data.frame(matrix(unlist(strsplit(chimeras$contig_id,",")),ncol=3,byrow=T)))
  fasta_ids$category[fasta_ids$id %in% chimeras$X1] = chimeras$tag 
}
table(fasta_ids$category)

fasta_ids$newID = apply(fasta_ids,1,function(x)paste(x,collapse=","))
write.fasta(fasta_file,fasta_ids$newID,output)
system(paste("rm -r",tmpdir))
