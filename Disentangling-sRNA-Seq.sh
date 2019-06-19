#binaries and scripts
tally_bin="/Users/macross/Documents/bin/reaper-15-065/src/tally"
tally_to_fasta_bin="/Users/macross/Dropbox/proyectos/miRNA/bin/tally_to_fasta.pl"
bowtie_build_bin="/Users/macross/Documents/bin/bowtie-1.2.2/bowtie-build"
bowtie_bin="/Users/macross/Documents/bin/bowtie-1.2.2/bowtie"
bowtie2_bin="/Users/macross/Documents/bin/bowtie2-2.3.3"
samtools_bin="/Users/macross/Documents/bin/samtools-1.5/samtools"
Rscripts_dir="/Users/macross/Dropbox/proyectos/miRNA/bin/r"
cci_bin="$Rscripts_dir/cov_correction_isoforms_0.7.R"
dea_counts_22g_bin="Rscripts_dir/dea_counts_22g_v2.R"
get_ambiguous_reads_bin="Rscripts_dir/get_ambiguous_reads.t.R"
chimera_assignment_bin="Rscripts_dir/chimera_assignment.r"
tag_ids_DEA_tab_bin="Rscripts_dir/tag_ids_DEA_tab.R"

##parameters
interaction="mm_hp" #interaction_model alias libraries
host="mm" # host alias libraries
par="hp" # parasite alias libraries
base_dir="/Users/macross/Documents/data/sRNAs_assembly/assembly_eval"
output_dir="DEA_contigs_informed_0.1" ###########
working_dir="${base_dir}/${interaction}/mapped/diSrna/${output_dir}"
threshold=0.1 ###########


#genome folders
host_genome="${base_dir}/${interaction}/${host}_genome/${host}_genome"
par_genome="${base_dir}/${interaction}/${par}_genome/${par}_genome"
host_par_genome="${base_dir}/${interaction}/${host}${par}_genome/${host}${par}_genome"
host_genome_bt2="${base_dir}/${interaction}/${host}_genome/bt2/${host}_genome"
par_genome_bt2="${base_dir}/${interaction}/${par}_genome/bt2/${par}_genome"

#library folders
#fastq files must be converted to fasta using tally + tally_to_fasta.pl
infected_libs="${base_dir}/${interaction}/mapped/${host}${par}"
host_libs="${base_dir}/${interaction}/mapped/${host}_combined"
contigs="${infected_libs}/trinity/inchworm.K19.L19_f.fa"

#obtain unique sequences from host_infected, host and parasite libraries 
mkdir $output_dir
cd $working_dir
mkdir tmp_dir
#cat ${infected_libs}/all.fq ${host_libs}/all.fq > tmp_dir/${host}${par}_${host}.fq
#$tally_bin -i tmp_dir/${host}${par}_${host}.fq -o tmp_dir/${host}${par}_${host}.tmp --nozip -format '%R%t%X%n'
#awk '{print $1 "\t"$2  "\t"length($1)} ' tmp_dir/${host}${par}_${host}.tmp > tmp_dir/${host}${par}_${host}.tally
#perl $tally_to_fasta_bin -f tmp_dir/${host}${par}_${host}.tally > ${host}${par}_${host}.fa

#get unitgs and alg_seqs: a) Unitigs are those reads that do not map to contigs b)alg_seqs are reads that map to contigs
#$bowtie_build_bin -f $contigs tmp_dir/contig_db
#$bowtie_bin -S -v 0 -p 8 -f tmp_dir/contig_db ${host}${par}_${host}.fa tmp_dir/tmp.sam --un tmp_dir/unitigs.fa --al alg.fa
#sed 's/>s/>u/' tmp_dir/unitigs.fa > unitigs.fa
#sed 's/>s/>c/' $contigs > contigs.fa
#cat contigs.fa alg.fa unitigs.fa  > contigs_alg_unitigs.fa

#Get global support counts (uniquely-mapping reads)
#Rscript $cci_bin ${working_dir} ${host}${par}_${host}.fa alg.fa unitigs.fa contigs.fa $threshold tmp_dir write_counts_hostpar empty $bowtie_build_bin $bowtie_bin $tally_to_fasta_bin $samtools_bin   

#Distribute multi-mapping reads among contigs based on global support counts
mkdir ${host}${par}
#for lib in "$infected_libs"/*.fa;
#do Rscript $cci_bin ${working_dir}/${host}${par} $lib alg.fa unitigs.fa contigs.fa $threshold tmp_dir informed_SpliT base_counts_hostpar.txt $bowtie_build_bin $bowtie_bin $tally_to_fasta_bin $samtools_bin
#done;

mkdir $host
#for lib in "$host_libs"/*.fa;
#do Rscript $cci_bin ${working_dir}/${host} $lib alg.fa unitigs.fa contigs.fa $threshold tmp_dir informed_SpliT base_counts_hostpar.txt $bowtie_build_bin $bowtie_bin $tally_to_fasta_bin $samtools_bin
#done;

#Merge library counts and calculate 22g frecuencies
Rscript $dea_counts_22g_bin contigs ${working_dir}/${host}${par} ${working_dir}/${host}${par} ${working_dir}

#Classify contigs_alg_reads sequences (bowtie1): host, parasite or ambiguous
Rscript $get_ambiguous_reads_bin $working_dir $contigs.fa2 $host_genome $par_genome

#Some chimeric contigs_alg_reads do not map perfectly, To classify those chimeric sequences we used Bowtie2 
Rscript $chimera_assignment_bin $working_dir ${contigs}.fa2 $host_par_genome $host_genome_bt2 $par_genome_bt2

final="${g22_opt}_alg_unitigs"
#add sequence classification () column to DEA matrix
Rscript $tag_ids_DEA_tab_bin $working_dir $working_dir c_alg_unitigs.tbl ${final}_mPar.fa ${final}_mHost.fa ${final}_amb.txt ${final}.fa ${final}_un_f.fa F

