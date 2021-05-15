#22/feb/19
#===1===#
#Tabulate individual contigs_counts
#Concatenate contigs_counts and reads_counts
#22g quantification

rm(list=ls())
computer="production"
if(computer=="test"){
  option = "contigs"
  dir1 = "/home/macross/compartida/data/miRNA/assembly_eval/mm_hp/mapped_18-50/diSrna/DEA_contigs_CEM/mmhp"
  dir2 = NA 
  dir3 = NA
  
}
if(computer=="production"){
  args = commandArgs(trailingOnly=TRUE)
  option = args[1]
  dir1 = args[2]
  dir2 = args[3]
  dir3 = args[4]
}

#===individual contig_counts files===#
files_from_dirs = function(dirs,pattern){
  for(i in 1:length(dirs)){
    fi = list.files(dirs[i],pattern = pattern,full.names = T,recursive = F)
    if(i==1){ft = fi}else{ft=c(ft,fi)}
  }
  return(ft)
}


get_big_matrix_counts = function(count_files){
  for(i in 1:length(count_files)){
    ccounts_f = count_files[i];print(ccounts_f)
    ccounts = read.table(ccounts_f,header=F)
    if(i==1){ccountsF = data.frame(ccounts$V2,row.names=ccounts$V1)}else{ccountsF=cbind(ccountsF,ccounts$V2)}
  }
  colnames(ccountsF) = gsub(".counts","",basename(count_files))
  return(ccountsF)
}

paste_contigs_clusters_reads = function(c_files,r_files,dir1){
  for(i in 1:length(c_files)){
    c = read.table(c_files[i],header = T,stringsAsFactors = F)
    r = read.table(r_files[i],header=T,stringsAsFactors = F)
    cr = rbind(c,r)
    
    lib = sub(".contigs","",basename(c_files[i]))
    if(option=="contigs"){cr_file = paste(dir1,paste(lib,"contigsreads",sep="."),sep="/")}
    if(option=="clusters"){cr_file = paste(dir1,paste(lib,"clustersreads",sep="."),sep="/")}
    write.table(cr,cr_file,row.names = F,col.names = T,quote = F,sep = "\t")
  }
}

write_22g = function(files_22g,dir1){
  for(i in 1:length(files_22g)){
    cu = read.table(files_22g[i],header=T,sep="\t",stringsAsFactors = F)
    
    ids = unique(cu$contig_id)
    len = length(ids)
    g22 = rep(0,len)
    u22 = rep(0,len)
    otherT = rep(0,len)
    other = rep(0,len)
    u_cu = data.frame(g22,u22,otherT,other,row.names = ids)
    
    cu$length = nchar(cu$sequence)
    cu$nucletide = substr(cu$sequence,0,1)
    
    cu_22G = cu[cu$length>=21 & cu$length<=24 & cu$nucletide=="G",]
    g22 = tapply(cu_22G$final_counts,factor(cu_22G$contig_id),sum)
    u_cu[names(g22),"g22"] = g22
    
    cu_22u = cu[cu$length>=21 & cu$length<=24 & cu$nucletide=="T",]
    u22 = tapply(cu_22u$final_counts,factor(cu_22u$contig_id),sum)
    u_cu[names(u22),"u22"] = u22
    
    otherT = tapply(cu$final_counts,factor(cu$contig_id),sum)
    u_cu[names(otherT),"otherT"] = otherT
    
    u_cu$other = u_cu$otherT - (u_cu$g22 + u_cu$u22)
    
    lib = basename(files_22g[i])
    lib = gsub(".cr","",lib)
    lib = gsub(".reads","",lib)
    u_cu_file = paste(dir1,paste(lib,"22g",sep="."),sep="/")
    t = apply(u_cu[,c(1,2,4)],1,sum)
    print(paste(lib,"counts:",sum(t)))
    write.table(u_cu[,c(1,2,4)],u_cu_file,row.names = T,col.names = T,quote = F,sep="\t")
  }
}

##################################3MAIN PROGRAM#########################
if(option=="reads"){
  r_files = list.files(dir1,pattern = ".reads$",full.names = T,recursive = F)
  write_22g(r_files,dir1)
}else{
  #counts
  contig_counts_f = files_from_dirs(c(dir1,dir2,dir3),".counts$");contig_counts_f
  ccounts = get_big_matrix_counts(contig_counts_f)
  
  cc_file = paste(dir1,"c_alg_unitigs.tbl",sep="/")
  write.table(ccounts,cc_file,row.names = T,col.names = T,quote = F,sep="\t")
  
  ##22g
  cr_files = files_from_dirs(c(dir1,dir2,dir3),".cr$")
  write_22g(cr_files,dir1)
}

