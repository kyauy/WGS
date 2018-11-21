# Comparative study of 3rd generation sequencing technologies

## Introduction

## Materials and methods

### Dataset
We select a presume de novo ID patient without diagnostics after Illumina HiSeq WGS and CGH array analysis.
We sequence this trio with Illumina Novaseq, 10x Genomics Chromium, PacBio, Oxford Nanopore and bionano.

### Protocol

#### NGS data processing

PacBio data have been aligned with minimap2 in hg38, SV calling have been made with pbsv2 and SNV calling with samtools mpileup.
Illumina data have been aligned with BWA-mem in hg38, SV calling have been made with Manta, Lumpy and Delly and SNV calling have been made with xAtlas.

Get clean SNV calling from Illumina Novaseq WGS :
```
{ ~/WGS/jobs }-> sbatch xatlas_child.job
{ ~/WGS/jobs }-> sbatch xatlas_father.job
{ ~/WGS/jobs }-> sbatch xatlas_mother.job
```

Split by sample all SV vcf :
```
##### PacBio
## pbsv2
{ ~/WGS/data/pbsv2 }-> bcftools view -s DNA17-06166 -O v -c 1  hg38.PBRT01+03+04+05+07.pbsv2.20180823.vcf > hg38.DNA17-06166.pbsv2.vcf
{ ~/WGS/data/pbsv2 }-> bcftools view -s DNA17-06167 -O v -c 1  hg38.PBRT01+03+04+05+07.pbsv2.20180823.vcf > hg38.DNA17-06167.pbsv2.vcf
{ ~/WGS/data/pbsv2 }-> bcftools view -s DNA17-06168 -O v -c 1  hg38.PBRT01+03+04+05+07.pbsv2.20180823.vcf > hg38.DNA17-06168.pbsv2.vcf

##### Illumina
## manta
15:16:18 kevin::login02 { ~/WGS/data/manta }-> bcftools view -O v -c 1 -s BvB41_child diploidSV.filtered.vcf.gz > hg38.DNA17-06166.manta.vcf
15:19:36 kevin::login02 { ~/WGS/data/manta }-> bcftools view -O v -c 1 -s BvB41_father diploidSV.filtered.vcf.gz > hg38.DNA17-06167.manta.vcf
15:19:55 kevin::login02 { ~/WGS/data/manta }-> bcftools view -O v -c 1 -s BvB41_mother diploidSV.filtered.vcf.gz > hg38.DNA17-06168.manta.vcf

## lumpy
15:21:35 kevin::login02 { ~/WGS/data/lumpy }-> bcftools view -O v -c 1 -s BvB41_child Trio5.Lumpy.svtyper2.GTalt+GQ20.vcf > hg38.DNA17-06166.lumpy.vcf
15:22:00 kevin::login02 { ~/WGS/data/lumpy }-> bcftools view -O v -c 1 -s BvB41_father Trio5.Lumpy.svtyper2.GTalt+GQ20.vcf > hg38.DNA17-06167.lumpy.vcf
15:22:11 kevin::login02 { ~/WGS/data/lumpy }-> bcftools view -O v -c 1 -s BvB41_mother Trio5.Lumpy.svtyper2.GTalt+GQ20.vcf > hg38.DNA17-06168.lumpy.vcf

## delly
15:24:30 kevin::login02 { ~/WGS/data/delly }-> bcftools view -O v -c 1 -s BvB41_child Trio5.Delly.GTalt+GQ20.vcf > hg38.DNA17-06166.delly.vcf
15:24:44 kevin::login02 { ~/WGS/data/delly }-> bcftools view -O v -c 1 -s BvB41_father Trio5.Delly.GTalt+GQ20.vcf > hg38.DNA17-06167.delly.vcf
15:24:55 kevin::login02 { ~/WGS/data/delly }-> bcftools view -O v -c 1 -s BvB41_mother Trio5.Delly.GTalt+GQ20.vcf > hg38.DNA17-06168.delly.vcf

## merge these 3 SV caller with maximum allowed distance of 500pb

# Patient 1
{ echo "/ifs/home/kevin/WGS/data/VCF/novaseq/manta/hg38.DNA17-06166.manta.vcf";
echo "/ifs/home/kevin/WGS/data/VCF/novaseq/lumpy/hg38.DNA17-06166.lumpy.vcf";
echo "/ifs/home/kevin/WGS/data/VCF/novaseq/delly/hg38.DNA17-06166.delly.vcf"; } > manta+lumpy+delly_DNA17-06166.fof

SURVIVOR merge manta+lumpy+delly_DNA17-06166.fof 500 1 1 0 0 20 manta+lumpy+delly_DNA17-06166.merged500pb.vcf

# Patient 2
{ echo "/ifs/home/kevin/WGS/data/VCF/novaseq/manta/hg38.DNA17-06167.manta.vcf";
echo "/ifs/home/kevin/WGS/data/VCF/novaseq/lumpy/hg38.DNA17-06167.lumpy.vcf";
echo "/ifs/home/kevin/WGS/data/VCF/novaseq/delly/hg38.DNA17-06167.delly.vcf"; } > manta+lumpy+delly_DNA17-06167.fof

SURVIVOR merge manta+lumpy+delly_DNA17-06167.fof 500 1 1 0 0 20 manta+lumpy+delly_DNA17-06167.merged500pb.vcf

# Patient 3
{ echo "/ifs/home/kevin/WGS/data/VCF/novaseq/manta/hg38.DNA17-06168.manta.vcf";
echo "/ifs/home/kevin/WGS/data/VCF/novaseq/lumpy/hg38.DNA17-06168.lumpy.vcf";
echo "/ifs/home/kevin/WGS/data/VCF/novaseq/delly/hg38.DNA17-06168.delly.vcf"; } > manta+lumpy+delly_DNA17-06168.fof

SURVIVOR merge manta+lumpy+delly_DNA17-06168.fof 500 1 1 0 0 20 manta+lumpy+delly_DNA17-06168.merged500pb.vcf
```

#### CNVs comparative processing

According to CLAMMS, deletions must satisfy no heterozygous SNPs and at least one homozygous SNP are called in the CNV region.
Duplications at least one heterozygous SNP is called in the CNV region and the average allele balance across all heterozygous SNPs in the region is in the range [0.611,0.723], corresponding to the 15th and 85th percentiles of inlier duplication calls.

Select only deletion and duplication CNV > 20pb and DP > 5 (assume that all DUP are INS) :
(beware that there could be a chrMT and chrM problem that need to be fixed)
```
################ Pacbio pbsv2 (according to https://github.com/PacificBiosciences/pbsv)
###### DELETION from 20pb to 100kb
14:01:52 kevin::login02 { ~/WGS/data/pbsv2 }-> bcftools view --output-type v --include 'INFO/SVLEN <= -20 && INFO/SVTYPE == "DEL" && DP>5 && GT="alt"' hg38.DNA17-06166.pbsv2.vcf > hg38.DNA17-06166.pbsv2.DEL.ONLY.vcf
14:02:45 kevin::login02 { ~/WGS/data/pbsv2 }-> bcftools view --output-type v --include 'INFO/SVLEN <= -20 && INFO/SVTYPE == "DEL" && DP>5 && GT="alt"' hg38.DNA17-06167.pbsv2.vcf > hg38.DNA17-06167.pbsv2.DEL.ONLY.vcf
14:04:04 kevin::login02 { ~/WGS/data/pbsv2 }-> bcftools view --output-type v --include 'INFO/SVLEN <= -20 && INFO/SVTYPE == "DEL" && DP>5 && GT="alt"' hg38.DNA17-06168.pbsv2.vcf > hg38.DNA17-06168.pbsv2.DEL.ONLY.vcf

###### DUPLICATION = PBSV call them INS ONLY = from 20bp to 5 kb / with DP filter
14:25:27 kevin::login02 { ~/WGS/data/pbsv2 }-> bcftools view --output-type v --include 'INFO/SVLEN >= 20 && INFO/SVTYPE == "INS" && DP>5 && GT="alt"' hg38.DNA17-06168.pbsv2.vcf > hg38.DNA17-06168.pbsv2.DUP.ONLY.vcf
14:25:50 kevin::login02 { ~/WGS/data/pbsv2 }-> bcftools view --output-type v --include 'INFO/SVLEN >= 20 && INFO/SVTYPE == "INS" && DP>5 && GT="alt"' hg38.DNA17-06167.pbsv2.vcf > hg38.DNA17-06167.pbsv2.DUP.ONLY.vcf
14:25:58 kevin::login02 { ~/WGS/data/pbsv2 }-> bcftools view --output-type v --include 'INFO/SVLEN >= 20 && INFO/SVTYPE == "INS" && DP>5 && GT="alt"' hg38.DNA17-06166.pbsv2.vcf > hg38.DNA17-06166.pbsv2.DUP.ONLY.vcf

################# Novaseq manta, lumpy and delly
###### DELETION of already filtered SV file
11:06:13 kevin::login02 { ~/WGS/data/VCF/novaseq }->  bcftools view --output-type v --include 'INFO/AVGLEN <= -20 && INFO/SVTYPE == "DEL" && GT="alt"' manta+lumpy+delly_DNA17-06166.merged500pb.vcf > manta+lumpy+delly_DNA17-06166.merged500pb.DEL.ONLY.vcf
11:09:25 kevin::login02 { ~/WGS/data/VCF/novaseq }->  bcftools view --output-type v --include 'INFO/AVGLEN <= -20 && INFO/SVTYPE == "DEL" && GT="alt"' manta+lumpy+delly_DNA17-06167.merged500pb.vcf > manta+lumpy+delly_DNA17-06167.merged500pb.DEL.ONLY.vcf
11:09:40 kevin::login02 { ~/WGS/data/VCF/novaseq }->  bcftools view --output-type v --include 'INFO/AVGLEN <= -20 && INFO/SVTYPE == "DEL" && GT="alt"' manta+lumpy+delly_DNA17-06168.merged500pb.vcf > manta+lumpy+delly_DNA17-06168.merged500pb.DEL.ONLY.vcf

###### DUPLICATION, with filtering of DEL& DUP call at same position
11:45:40 kevin::login02 { ~/WGS/data/VCF/novaseq }->  bcftools view --output-type v --include ' INFO/SVTYPE == "DUP" || INFO/SVTYPE == "INS" && GT="alt"' manta+lumpy+delly_DNA17-06168.merged500pb.vcf | grep -v DEL > manta+lumpy+delly_DNA17-06168.merged500pb.DUP.INS.vcf
11:46:26 kevin::login02 { ~/WGS/data/VCF/novaseq }->  bcftools view --output-type v --include ' INFO/SVTYPE == "DUP" || INFO/SVTYPE == "INS" && GT="alt"' manta+lumpy+delly_DNA17-06167.merged500pb.vcf | grep -v DEL > manta+lumpy+delly_DNA17-06167.merged500pb.DUP.INS.vcf
11:46:34 kevin::login02 { ~/WGS/data/VCF/novaseq }->  bcftools view --output-type v --include ' INFO/SVTYPE == "DUP" || INFO/SVTYPE == "INS" && GT="alt"' manta+lumpy+delly_DNA17-06166.merged500pb.vcf | grep -v DEL > manta+lumpy+delly_DNA17-06166.merged500pb.DUP.INS.vcf
## Change all DUP to INS
11:50:09 kevin::login02 { ~/WGS/data/VCF/novaseq }-> sed 's/DUP/INS/g' manta+lumpy+delly_DNA17-06168.merged500pb.DUP.INS.vcf > manta+lumpy+delly_DNA17-06168.merged500pb.ALL.INS.vcf
11:50:36 kevin::login02 { ~/WGS/data/VCF/novaseq }-> sed 's/DUP/INS/g' manta+lumpy+delly_DNA17-06167.merged500pb.DUP.INS.vcf > manta+lumpy+delly_DNA17-06167.merged500pb.ALL.INS.vcf
11:50:45 kevin::login02 { ~/WGS/data/VCF/novaseq }-> sed 's/DUP/INS/g' manta+lumpy+delly_DNA17-06166.merged500pb.DUP.INS.vcf > manta+lumpy+delly_DNA17-06166.merged500pb.ALL.INS.vcf
```

Sort and index VCF to be ready for filtering and annotation :
```
bcftools sort -O z manta+lumpy+delly.merged500bp.DEL.ONLY.vcf.gz -o manta+lumpy+delly.merged500bp.DEL.ONLY.sorted.vcf.gz
bcftools index manta+lumpy+delly.merged500bp.DEL.ONLY.sorted.vcf.gz
bcftools sort -O z manta+lumpy+delly.merged500bp.DUP.only.vcf.gz -o manta+lumpy+delly.merged500bp.DUP.ONLY.sorted.vcf.gz
bcftools index manta+lumpy+delly.merged500bp.DUP.ONLY.sorted.vcf.gz

bcftools index IlluminaSNV/calling-110918/BvB41_child.haplotypecaller.vcf.gz

bcftools sort -O z PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.DEL.ONLY.vcf.gz -o PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.DEL.ONLY.sorted.vcf.gz
bcftools index PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.DEL.ONLY.sorted.vcf.gz
bcftools sort -O z PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.DUP.ONLY.vcf.gz -o PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.DUP.ONLY.sorted.vcf.gz
bcftools index PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.DUP.ONLY.sorted.vcf.gz
```

Filter and annotate VCF by AnnotSV :

```
./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/AnnotSV/manta+lumpy+delly.merged500bp.DEL.ONLY.sorted.vcf -bedtools /cm/shared/apps/bioinf/bedtools/2.25.0/bin/bedtools -outputDir ~/WGS/data/AnnotSV/ -vcfFiles ~/WGS/data/AnnotSV/BvB41_child.haplotypecaller.vcf

./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/AnnotSV/PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.DEL.ONLY.sorted.vcf -bedtools /cm/shared/apps/bioinf/bedtools/2.25.0/bin/bedtools -outputDir ~/WGS/data/AnnotSV/ -vcfFiles ~/WGS/data/AnnotSV/BvB41_child.haplotypecaller.vcf
```

#### Haplotype by read and pedigree based phasing and genotype

We apply read-based phasing WhatsHap (Martin et al. 2018) with pedigree-phasing algorithm (Garg et al. 2018)
```
usage: whatshap phase
output-read-list FILE
                        Write reads that have been used for phasing to FILE.
```

### Results

#### SV calling performance among technologies


#### Phasing to resolve SNV calling

We reanalyse phased SNV data looking for unrecognized variants.


### References

```
Marcel Martin, Murray Patterson, Shilpa Garg, Sarah O. Fischer, Nadia Pisanti, Gunnar W. Klau, Alexander Schoenhuth, Tobias Marschall.
WhatsHap: fast and accurate read-based phasing
bioRxiv 085050
doi: 10.1101/085050

Shilpa Garg, Marcel Martin, Tobias Marschall.
Read-based phasing of related individuals
Bioinformatics 2016; 32 (12): i234-i242.
doi: 10.1093/bioinformatics/btw276
```
