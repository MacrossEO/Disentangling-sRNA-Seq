#02/oct/18
#chimera assignment 
rm(list=ls())
library(protr)

computer="mazorka"
if(computer=="mac"){
  md = "/Users/macross/Documents/data/sRNAs_assembly/assembly_eval/hp_mm/mapped/hpmm/test"
  fa_file = "chimeric_in_DEA.fasta" #
  host_par_genome = "~/Documents/data/sRNAs_assembly/hpmm_genome/hpmm_genome"
  host_genome_bt2 = "/Users/macross/Documents/data/sRNAs_assembly/Mus_musculus_genome/bt2/mm_genome"
  par_genome_bt2 = "/Users/macross/Documents/data/sRNAs_assembly/hp_genome/bt2/hp_genome"
   
  bowtie_bin = "~/Documents/bin/bowtie-1.2.2/bowtie"
  bowtie2_bin = "~/Documents/bin/bowtie2-2.3.3/bowtie2"
  tally_to_fasta_bin = "~/Dropbox/proyectos/miRNA/bin/tally_to_fasta.pl"
}
if(computer=="mazorka"){
  args = commandArgs(trailingOnly=TRUE)
  md = args[1]
  fa_file = args[2]
  host_par_genome = args[3]
  host_genome_bt2 = args[4]
  par_genome_bt2 = args[5]
  
  bowtie_bin = args[6]
  bowtie2_bin = args[7]
  tally_to_fasta_bin = args[8]
}

dir = paste(md,"tmp_dir",sep="/")#switch when necessary
dir.create(dir)

#files
query =  paste(md,fa_file,sep="/")
sam_file = paste(dir,"tmp.sam",sep="/")
base_name = paste(unlist(strsplit(fa_file,"[.]"))[1],"_un.fa",sep="")
chimera_file = paste(md,base_name,sep="/");chimera_file

#get chimeras
bowtie_call = paste(bowtie_bin," -S -v 0 -p 8 -f ",host_par_genome,query,sam_file," --un", chimera_file);bowtie_call
system(bowtie_call)

#align chimeras to host
host_sam_file = paste(dir,paste(base_name,"_toHost_bt2.sam",sep=""),sep="/")
bowtie2_call = paste(bowtie2_bin,"-p 8 -f -x",host_genome_bt2,"-U",chimera_file,"-S",host_sam_file);bowtie2_call
system(bowtie2_call)

#align chimeras to parasite
par_sam_file = paste(dir,paste(base_name,"_toPar_bt2.sam",sep=""),sep="/")
bowtie2_call = paste(bowtie2_bin,"-p 8 -f -x",par_genome_bt2,"-U",chimera_file,"-S",par_sam_file);bowtie2_call
system(bowtie2_call)

#get chimeras
chimeras = readFASTA(chimera_file)
chimeras = data.frame(unlist(chimeras),stringsAsFactors = F)
chimeras$contig_id = rownames(chimeras)
names(chimeras) = c("seq","contig_id")
head(chimeras)

#asign chimeras
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
  print(host_par)
  
  if(host_par[1]==host_par[2]){
    chimeras$contig_id[i] = gsub("c","a",id)
  }else{
    idx = which.max(host_par)
    if(idx==1){
      if(host_par[1]>2){
        chimeras$contig_id[i] = gsub("c","a",id)
      }else{
        chimeras$contig_id[i] = gsub("c","h",id)#host
      }
    }else{
      if(host_par[2]>2){
        chimeras$contig_id[i] = gsub("c","a",id)
      }else{
        chimeras$contig_id[i] = gsub("c","p",id)#par
      }
    }
  }
}

##write to a dataframe
tally_out_f = sub(pattern = ".fa$",replacement = "_f.tally",chimera_file);tally_out_f
write.table(chimeras, file = tally_out_f, quote = FALSE, sep = "\t", col.names = F,row.names = F)

#tally to fasta
reads_fasta = sub(pattern = ".tally",replacement = ".fa",tally_out_f)
tally_to_fasta_call = paste(tally_to_fasta_bin,"-f",tally_out_f,">",reads_fasta);tally_to_fasta_call
system(tally_to_fasta_call)

grep_call = paste("grep -c ',h'",reads_fasta);grep_call
print("un_al to host:"); print(system(grep_call,intern = T))
awk_call = paste("grep ',h'",reads_fasta, "| awk -F',' '{print $3}'  | awk '{s+=$1}END{print s}' 2>&1");awk_call
print("un_al counts to host :");print(system(awk_call,intern = T))

grep_call = paste("grep -c ',p'",reads_fasta);grep_call
print("un_al to parasite:"); print(system(grep_call,intern = T))
awk_call = paste("grep ',p'",reads_fasta, "| awk -F',' '{print $3}'  | awk '{s+=$1}END{print s}' 2>&1");awk_call
print("un_al counts to parasite :");print(system(awk_call,intern = T))

grep_call = paste("grep -c ',a'",reads_fasta);grep_call
print("un_al ambiguous:"); print(system(grep_call,intern = T))
awk_call = paste("grep ',a'",reads_fasta, "| awk -F',' '{print $3}'  | awk '{s+=$1}END{print s}' 2>&1");awk_call
print("un_al counts ambiguous :");print(system(awk_call,intern = T))


#rm temp files
rm_files = paste(dir,"*",sep="/")
#system(paste("rm",rm_files))


