##Main script to get the coverage of piRNAi fragments on DNA sequencing libraries
#Libraries are grouped depending on the injection mix and hence the fasta reference used, i.e.:
#Set1:CFJ91_G0,SCFJ91_G6
#Set2:S72_2_G10,S72_5_G10
#Set3:S247_1,S247_5,S247_8,S248_4,S248_5,S248_8
#See suplemental table for more info

##Main script

##Add reference genome to each set
cd Set1
cat ../c_elegans.PRJNA13758.WS275.genomic_softmasked.fa Background.fasta Plasmids.fasta > Reference.fasta
cd ..

cd Set2
cat ../c_elegans.PRJNA13758.WS275.genomic_softmasked.fa Background.fasta Plasmids.fasta > Reference.fasta
cd ..

cd Set3
cat ../c_elegans.PRJNA13758.WS275.genomic_softmasked.fa Plasmids.fasta > Reference.fasta
cd ..

##Also note that each directory contains the DNA libraries uncompressed within folders named as them

##Main parameters for data processing
#Minimim mapping per base in bam file
minMaqQforBam=1		
#Number of threads to use
Ncpu=28
#Ram memory for Java
ramG=40

#Additional samtools instalation for samtools coverage
sami=../samtools-1.12/bin/bin/samtools

##Start with simplest set one
cd Set3
#For each library
for lib in `ls -d S*`; do
#Specify library
echo ${lib}
#Go to directory
cd ${lib}
#Copy reference file to map
cp ../Reference.fasta .
#Index the file by samtools, bwa, and picard
samtools faidx Reference.fasta
bwa index Reference.fasta
java -jar /home/velazqam/Downloads/Sonia_WGS/Software/picard2.23.6/picard.jar CreateSequenceDictionary R=Reference.fasta O=Reference.dict

#Perform fastqc
fastqc -t $Ncpu ${lib}*.gz;

##Perform mapping
bwa mem -t $Ncpu -aM Reference.fasta ${lib}*.gz | samtools view -buS - |samtools sort - -o ${lib}_rawmap.bam

#Samtools sort
samtools sort -@ $Ncpu ${lib}_rawmap.bam > ${lib}.sort.bam
#rm ${lib}_rawmap.bam

#Samtools index
samtools index ${lib}.sort.bam

##Remove minimal quality for Bam
samtools view -@ $Ncpu -q $minMaqQforBam -bh ${lib}.sort.bam > ${lib}.UM.bam
#rm ${lib}.sort.bam ${lib}.sort.bam.bai

#Picard add replace groups ## Explain its importance for combining multiple libraries
java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/picard2.23.6/picard.jar AddOrReplaceReadGroups I=${lib}.UM.bam O=${lib}.RG.bam RGID=${lib} RGLB=LB RGPL=illumina RGPU=PU RGSM=${lib}
#rm ${lib}.UM.bam

#Picard remove duplicates
java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/picard2.23.6/picard.jar MarkDuplicates I=${lib}.RG.bam O=${lib}.deDup.bam M=${lib}.dedupMetrics REMOVE_DUPLICATES=true
#rm ${lib}.RG.bam ${lib}.dedupMetrics

##Samtools
samtools index ${lib}.deDup.bam

${sami} coverage ${lib}.deDup.bam > ${lib}.coverage.tsv

cd ..;

done

cd ..

##Set2
cd Set2
for lib in `ls -d S*`; do
echo ${lib}
cd ${lib}
cp ../Reference.fasta .
samtools faidx Reference.fasta
bwa index Reference.fasta
java -jar /home/velazqam/Downloads/Sonia_WGS/Software/picard2.23.6/picard.jar CreateSequenceDictionary R=Reference.fasta O=Reference.dict

#Perform fastqc
fastqc -t $Ncpu ${lib}*.gz;

##Perform mapping
bwa mem -t $Ncpu -aM Reference.fasta ${lib}*.gz | samtools view -buS - |samtools sort - -o ${lib}_rawmap.bam

#Samtools sort
samtools sort -@ $Ncpu ${lib}_rawmap.bam > ${lib}.sort.bam
#rm ${lib}_rawmap.bam

#Samtools index
samtools index ${lib}.sort.bam

##Remove minimal quality for Bam
samtools view -@ $Ncpu -q $minMaqQforBam -bh ${lib}.sort.bam > ${lib}.UM.bam
#rm ${lib}.sort.bam ${lib}.sort.bam.bai

#Picard add replace groups ## Explain its importance for combining multiple libraries
java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/picard2.23.6/picard.jar AddOrReplaceReadGroups I=${lib}.UM.bam O=${lib}.RG.bam RGID=${lib} RGLB=LB RGPL=illumina RGPU=PU RGSM=${lib}
#rm ${lib}.UM.bam

#Picard remove duplicates
java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/picard2.23.6/picard.jar MarkDuplicates I=${lib}.RG.bam O=${lib}.deDup.bam M=${lib}.dedupMetrics REMOVE_DUPLICATES=true
#rm ${lib}.RG.bam ${lib}.dedupMetrics

##Samtools
samtools index ${lib}.deDup.bam

${sami} coverage ${lib}.deDup.bam > ${lib}.coverage.tsv

##GATK Indel realigner
java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/gatk-3.8.1.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt $Ncpu -R Reference.fasta -I ${lib}.deDup.bam -o ${lib}.forIndelRealigner.intervals

java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/gatk-3.8.1.0/GenomeAnalysisTK.jar -T IndelRealigner -R Reference.fasta -I ${lib}.deDup.bam -targetIntervals ${lib}.forIndelRealigner.intervals -o ${lib}.realigned.bam

#rm ${lib}.deDup.bam ${lib}.deDup.bam.bai ${lib}.forIndelRealigner.intervals

#Unified genotyper
##DATA NOT SHOWN NOR USED IN PAPER
cd ..;

done

cd ..


##Set1
cd Set1
for lib in `ls -d *_*`; do
echo ${lib}
cd ${lib}
cp ../Reference.fasta .
samtools faidx Reference.fasta
bwa index Reference.fasta
java -jar /home/velazqam/Downloads/Sonia_WGS/Software/picard2.23.6/picard.jar CreateSequenceDictionary R=Reference.fasta O=Reference.dict

#Perform fastqc
fastqc -t $Ncpu ${lib}*.gz;

##Perform mapping
bwa mem -t $Ncpu -aM Reference.fasta ${lib}*.gz | samtools view -buS - |samtools sort - -o ${lib}_rawmap.bam

#Samtools sort
samtools sort -@ $Ncpu ${lib}_rawmap.bam > ${lib}.sort.bam
#rm ${lib}_rawmap.bam

#Samtools index
samtools index ${lib}.sort.bam

##Remove minimal quality for Bam
samtools view -@ $Ncpu -q $minMaqQforBam -bh ${lib}.sort.bam > ${lib}.UM.bam
#rm ${lib}.sort.bam ${lib}.sort.bam.bai

#Picard add replace groups ## Explain its importance for combining multiple libraries
java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/picard2.23.6/picard.jar AddOrReplaceReadGroups I=${lib}.UM.bam O=${lib}.RG.bam RGID=${lib} RGLB=LB RGPL=illumina RGPU=PU RGSM=${lib}
#rm ${lib}.UM.bam

#Picard remove duplicates
java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/picard2.23.6/picard.jar MarkDuplicates I=${lib}.RG.bam O=${lib}.deDup.bam M=${lib}.dedupMetrics REMOVE_DUPLICATES=true
#rm ${lib}.RG.bam ${lib}.dedupMetrics

##Samtools
samtools index ${lib}.deDup.bam

${sami} coverage ${lib}.deDup.bam > ${lib}.coverage.tsv

##GATK Indel realigner
java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/gatk-3.8.1.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt $Ncpu -R Reference.fasta -I ${lib}.deDup.bam -o ${lib}.forIndelRealigner.intervals

java -Xmx${ramG}g -jar /home/velazqam/Downloads/Sonia_WGS/Software/gatk-3.8.1.0/GenomeAnalysisTK.jar -T IndelRealigner -R Reference.fasta -I ${lib}.deDup.bam -targetIntervals ${lib}.forIndelRealigner.intervals -o ${lib}.realigned.bam

#rm ${lib}.deDup.bam ${lib}.deDup.bam.bai ${lib}.forIndelRealigner.intervals

#Unified genotyper
##DATA NOT SHOWN NOR USED IN PAPER
cd ..;

done

cd ..
