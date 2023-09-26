# population structure
Using https://github.com/Axolotl233/Simple_Script/blob/master/Vcf.structure.pl for population structure analysis (ML tree, PCA , admixture)

# demography
Using https://github.com/dportik/dadi_pipeline for demography simulation

# Transposable Element Polymorphism
TEP called by TEPID (https://github.com/ListerLab/TEPID)
TEP summarized using Perl script （z.util/TEP_stat*.pl）
PCA using R script (z.util/TEP_plot.R)
Gene_expression using R script (z.util/TEP_plot.R)

# structural variants analysis 
configManta.py --bam sample.realn.bam --referenceFasta genome.fasta --runDir manta_out
Using https://github.com/Axolotl233/Simple_Script/blob/master/Vcf.structure.pl for population structure analysis (NJ tree)