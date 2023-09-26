# Pi
vcftools --gzvcf Pop.final.SNP.filter.vcf.gz --keep ../AS.pop --out Pi.AS.lst.50K.12.5K --window-pi 50000 --window-pi-step 12500
vcftools --gzvcf Pop.final.SNP.filter.vcf.gz --keep ../ES.pop --out Pi.ES.lst.50K.12.5K --window-pi 50000 --window-pi-step 12500

#LD
PopLDdecay -InVCF ./Pop.final.SNP.filter.vcf.gz -OutStat Pop.ld.ES -SubPop ./ES.pop -MaxDist 5000
PopLDdecay -InVCF ./Pop.final.SNP.filter.vcf.gz -OutStat Pop.ld.AS -SubPop ./AS.pop -MaxDist 5000

#Fis and S
Using https://github.com/Axolotl233/Simple_Script/blob/master/Command.Plink.F.pl

#DXY
pixy --stat pi dxy --vcf vcf.file --populations pop.lst  --n_cores 10 --bypass_invariant_check yes --output_folder pixy.out

#FST
vcftools --gzvcf vcf.gz --fst-window-size 10000 --fst-window-step 0 --weir-fst-pop pop1.lst --weir-fst-pop pop2.lst  --out Fst.10K.0K
vcftools --gzvcf vcf.gz --weir-fst-pop pop1.lst --weir-fst-pop pop2.lst  --out Fst.single
python3 z.Util/island_select.py  Fst.10K.0K.fst  Fst.single.fst  10000 0 --top=0.05 --per=5000000 --T=60

#HKA
python3 z.Util/HKA.from_fst_bed_get_fixed.py Fst.10K.0K.fst gene.bed vcf.file pop1.lst 20
Rscript HKA.test.R gene_fixed_and_polymorphic
python3 z.Util/HKA.from_fst_bed_get_fixed.py Fst.10K.0K.fst gene.bed vcf.file pop2.lst 20
Rscript HKA.test.R gene_fixed_and_polymorphic

#Recombination
python3 z.Util/FastEPRR_v1.2.py O=pop2 V=pop2.vcf.gz R=path/to/genome.fasta.fai W=10 S=10 T=30 Phase=T
python3 z.Util/FastEPRR_v1.2.py O=pop1 V=pop1.vcf.gz R=path/to/genome.fasta.fai W=10 S=10 T=30 Phase=T

#Gene_anno
Using https://github.com/Axolotl233/Simple_Script/blob/master/Vcf.anno.bed.window.pl

