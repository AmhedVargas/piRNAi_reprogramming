##Bowtie mapping of short RNA libraries
#SG0919_lib29,SG0919_lib30,SG0919_lib31,SG0919_lib32,SG1020_lib10,SG1020_lib11,SG1020_lib12,SG1020_lib13,SG1020_lib14,SG1020_lib8,SG1020_lib9

##Set working directory where all compressed files are stored
#cd #
##FOr each library do
for file in `ls *.gz`; do 
#Mention library
echo ${file%*_*_*.gz}; 
#Make a directory
mkdir ${file%*_*_*.gz}; 
#Place file into that directory
mv ${file} ${file%*_*_*.gz};
#Go to that directory
cd ${file%*_*_*.gz};
#Perform basic quality assessment
fastqc ${file}
#To be consisstent with previous analysis, do not remove sequencing adapter 
#cutadapt3 -u 4 -m 15 -M 45 ${file} -o ${file%*_*_*.gz}_TrimmednFiltered.fastq
#But do keep only reads within 15 to 45 bp
cutadapt3 -m 15 -M 45 ${file} -o ${file%*_*_*.gz}_TrimmednFiltered.fastq
#Get Fasta reference containing ce11/WS235 and piRNA sequences
cp ../Fragments_reference.fasta .
#Index the reference
bowtie-build Fragments_reference.fasta Fragments_bowtie
#Map files to reference
bowtie -S -5 4 -v 0 -a Fragments_bowtie ${file%*_*_*.gz}_TrimmednFiltered.fastq ${file%*_*_*.gz}_TrimmednFiltered_bowtie.sam
#Use awk to obtain only mapped reads with a 21bp length and starting with T (21U reads)
awk -F"\t" '{if(length($10) == 21){print $0}}' ${file%*_*_*.gz}_TrimmednFiltered_bowtie.sam | awk -F"\t" '{if($2 != 4){print $0}}' | awk -F"\t" '{if($2 == 16){if(substr($10,21,1) == "A"){print $0}}else{if(substr($10,1,1) == "T"){print $0}}}' > ${file%*_*_*.gz}_TrimmednFiltered.21-U.bowtie.sam
#Convert sam into bam
grep "@" ${file%*_*_*.gz}_TrimmednFiltered_bowtie.sam | cat - ${file%*_*_*.gz}_TrimmednFiltered.21-U.bowtie.sam | samtools view -bS - | samtools sort - -o ${file%*_*_*.gz}_TrimmednFiltered.21-U.sort.bowtie.bam
#Get coverage by position by using mpileup
samtools mpileup -f Fragments_reference.fasta ${file%*_*_*.gz}_TrimmednFiltered.21-U.sort.bowtie.bam | awk -F"\t" '{print $1"\t"$3"\t"$2"\t"$4}' > ${file%*_*_*.gz}.coverage.tab
#PLot a preliminary figure
Rscript --vanilla ../simpleggplot.R ${file%*_*_*.gz}.coverage.tab ${file%*_*_*.gz}
#DOne analysis, return to main directory
cd ..;
done


