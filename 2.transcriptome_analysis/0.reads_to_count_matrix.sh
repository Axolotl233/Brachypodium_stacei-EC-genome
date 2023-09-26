# align
hisat2 -x genome_index -p 20 -X 500 --fr --min-intronlen 20 --max-intronlen 500000 --dta -1 sample_1.clean.fq.gz -2 sample_2.clean.fq.gz 2>> sample.txt |samtools view -bS - | samtools sort -@ 20 -o sample.sort.bam

# assembly and quantification 
stringtie -e -G genome.gff -p 20 -o sample.gtf sample.sort.bam

# gtf to martix
Using "prepDE.py3" provided by stringtie website
