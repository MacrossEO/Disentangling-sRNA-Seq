#21/nov/18
#get amiguous reads
rm(list=ls())
suppressMessages(library(protr))
suppressMessages(library(Rsamtools))

computer="mazorka"
if(computer=="mac"){
  md = "/Users/macross/Documents/data/sRNAs_assembly/assembly_eval/hp_mm/mapped/hpmm/cov_cor2"#contigs dir
  fa_file = "contigs_unitigs.fa" # fastq_fasta_file (libray)
  host_genome = "/Users/macross/Documents/data/sRNAs_assembly/Mus_musculus_genome/mm_genome"
  par_genome = "/Users/macross/Documents/data/sRNAs_assembly/hp_genome/hp_genome"
  
  bowtie_bin = "~/Documents/bin/bowtie-1.2.2/bowtie"
  samtools_bin = "~/Documents/bin/samtools-1.5/samtools"
}
if(computer=="mazorka"){
  args = commandArgs(trailingOnly=TRUE)
  md = args[1]#contigs dir
  fa_file = args[2]#fastq_fasta_file (library)
  host_genome = args[3]
  par_genome = args[4]
  
  bowtie_bin = args[5]
  samtools_bin = args[6]
}

tmpdir = paste(md,"tmp_dir",sep="/")#switch when necessary
dir.create(tmpdir)

#out_files
base_name = unlist(strsplit(fa_file,"[.]"))[1];base_name
map_to_host_file = paste(md,paste(base_name,"_mHost.fa",sep=""),sep="/");map_to_host_file
map_to_par_file = paste(md,paste(base_name,"_mPar.fa",sep=""),sep="/");map_to_par_file
sam_file = paste(tmpdir,"lib.sam",sep="/");sam_file
fa_file = paste(md,fa_file,sep="/")
ambiguous_reads_file = paste(md,paste(base_name,"_amb.txt",sep=""),sep="/")

#map to host 
bowtie_call = paste(bowtie_bin,"-S -v 0 -p 8 -f",host_genome,fa_file,sam_file,"--al ",map_to_host_file);bowtie_call
system(bowtie_call)
map_to_host = readFASTA(map_to_host_file)

#map to parasit
bowtie_call = paste(bowtie_bin,"-S -v 0 -p 8 -f",par_genome,fa_file,sam_file,"--al ",map_to_par_file);bowtie_call
system(bowtie_call)
map_to_par = readFASTA(map_to_par_file)

#intersection 
ambiguous_reads = intersect(names(map_to_host),names(map_to_par))
write.table(data.frame(ambiguous_reads),ambiguous_reads_file,quote=F,sep="\t",row.names = F,col.names = F)

awk_call = paste("awk -F',' '{print $2}' ", ambiguous_reads_file," | awk '{s+=$1}END{print s}' 2>&1");awk_call
print(paste("ambiguous contigs_unitigs:",length(ambiguous_reads)))
print("ambiguous contigs_unitigs_counts:");system(awk_call)

