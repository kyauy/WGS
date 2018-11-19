# Comparative study of 3rd generation sequencing technologies

## Introduction

## Materials and methods

### Dataset
We select a presume de novo ID patient without diagnostics after Illumina HiSeq WGS and CGH array analysis.
We sequence this trio with Illumina Novaseq, 10x Genomics Chromium, PacBio, Oxford Nanopore and bionano.

### Protocol

#### NGS data processing

PacBio data have been aligned with minimap2 in hg38 and SV calling have been made with pbsv2.
Illumina data have been aligned with BWA-mem in hg38 and SV calling have been made with Manta, Lumpy and Delly.

#### Filter CNVs

According to CLAMMS, deletions must satisfy no heterozygous SNPs and at least one homozygous SNP are called in the CNV region.
Duplications at least one heterozygous SNP is called in the CNV region and the average allele balance across all heterozygous SNPs in the region is in the range [0.611,0.723], corresponding to the 15th and 85th percentiles of inlier duplication calls.

Get clean SNV calling from Illumina Novaseq WGS :
```
{ ~/WGS/jobs }-> sbatch xatlas_child.job
{ ~/WGS/jobs }-> sbatch xatlas_father.job
{ ~/WGS/jobs }-> sbatch xatlas_mother.job
```

Select only deletion and duplication CNV :
```
{ ~/WGS/data }-> bcftools view -O v -i 'SVTYPE="DEL"' manta+lumpy+delly.merged500bp.vcf  > manta+lumpy+delly.merged500bp.DEL.only.vcf
{ ~/WGS/data }-> bcftools view -O v -i 'SVTYPE="DUP"' manta+lumpy+delly.merged500bp.vcf  > manta+lumpy+delly.merged500bp.DUP.only.vcf

bcftools view -O v -i 'SVTYPE="DEL"' PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.vcf > PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.DEL.ONLY.vcf
bcftools view -O v -i 'SVTYPE="DUP"' PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.vcf > PacBio+Illumina.merged500bp.SURVIVOR-1.0.4.DUP.ONLY.vcf
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
