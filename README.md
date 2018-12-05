# Comparative study of 3rd generation sequencing technologies

## Introduction

## Materials and methods

### Dataset
We select a presume de novo ID patient without diagnostics after Illumina HiSeq WGS and CGH array analysis.
We sequence this trio with Illumina Novaseq, 10x Genomics Chromium, PacBio, Oxford Nanopore and bionano.

### Protocol

### Requirements

bedtools 2.27.1 (or higher) is needed for AnnotSV, since they fixed an unfortunate bug.
AnnotSV 1.2
SURVIVOR 1.0.5

#### NGS data processing

PacBio data have been aligned with minimap2 in hg38, SV calling have been made with pbsv2 and SNV calling with samtools mpileup.
Illumina data have been aligned with BWA-mem in hg38, SV calling have been made with Manta, Lumpy and Delly and SNV calling have been made with xAtlas.

Get clean SNV calling from Illumina Novaseq WGS and create special selection for DUP filtering :
```
#### BATCH
{ ~/WGS/jobs }-> sbatch xatlas_child.job
{ ~/WGS/jobs }-> sbatch xatlas_father.job
{ ~/WGS/jobs }-> sbatch xatlas_mother.job

### DUPLICATION process
# VCF with heterozygous SNV & AB >= 0.611
17:01:02 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view --output-type v --include 'FMT/GT="0/1" && ((FMT/VR)/(FMT/DP)>=0.611)' BvB41_child_xAtlas.recode.sorted.vcf.gz >  BvB41_child_xAtlas.recode.sorted.DUP.vcf
17:03:04 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view --output-type v --include 'FMT/GT="0/1" && ((FMT/VR)/(FMT/DP)>=0.611)' BvB41_father_xAtlas.recode.sorted.vcf.gz >  BvB41_father_xAtlas.recode.sorted.DUP.vcf
17:04:50 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view --output-type v --include 'FMT/GT="0/1" && ((FMT/VR)/(FMT/DP)>=0.611)' BvB41_mother_xAtlas.recode.sorted.vcf.gz >  BvB41_mother_xAtlas.recode.sorted.DUP.vcf

# VCF with heterozygous SNV & AB [0.45-0.55]
bcftools view --output-type v --include 'FMT/GT="0/1" && ((FMT/VR)/(FMT/DP)>=0.45) && ((FMT/VR)/(FMT/DP)<=0.55)' BvB41_child_xAtlas.recode.sorted.vcf.gz > BvB41_child_xAtlas.recode.sorted.notDUP.vcf
sed -i 's/\t\/ifs/\tnotdup/g' BvB41_child_xAtlas.recode.sorted.notDUP.vcf
bcftools view --output-type v --include 'FMT/GT="0/1" && ((FMT/VR)/(FMT/DP)>=0.45) && ((FMT/VR)/(FMT/DP)<=0.55)' BvB41_father_xAtlas.recode.sorted.vcf.gz > BvB41_father_xAtlas.recode.sorted.notDUP.vcf
sed -i 's/\t\/ifs/\tnotdup/g' BvB41_father_xAtlas.recode.sorted.notDUP.vcf
bcftools view --output-type v --include 'FMT/GT="0/1" && ((FMT/VR)/(FMT/DP)>=0.45) && ((FMT/VR)/(FMT/DP)<=0.55)' BvB41_mother_xAtlas.recode.sorted.vcf.gz > BvB41_mother_xAtlas.recode.sorted.notDUP.vcf
sed -i 's/\t\/ifs/\tnotdup/g' BvB41_mother_xAtlas.recode.sorted.notDUP.vcf
### BGZIP all
for i in *DUP.vcf ; do bgzip $i ; done
for i in *notDUP.vcf ; do bgzip $i ; done
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

Merging sequencing technologies :

```
###### DELETION
{ echo "/ifs/home/kevin/WGS/data/VCF/pbsv2/hg38.DNA17-06166.pbsv2.DEL.ONLY.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/novaseq/manta+lumpy+delly_DNA17-06166.merged500pb.DEL.ONLY.vcf"; } > pacbio+novaseq_DNA17-06166.DEL.fof

SURVIVOR merge pacbio+novaseq_DNA17-06166.DEL.fof 500 1 1 0 0 20 pacbio+novaseq_DNA17-06166.DEL.vcf

{ echo "/ifs/home/kevin/WGS/data/VCF/pbsv2/hg38.DNA17-06167.pbsv2.DEL.ONLY.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/novaseq/manta+lumpy+delly_DNA17-06167.merged500pb.DEL.ONLY.vcf"; } > pacbio+novaseq_DNA17-06167.DEL.fof

SURVIVOR merge pacbio+novaseq_DNA17-06167.DEL.fof 500 1 1 0 0 20 pacbio+novaseq_DNA17-06167.DEL.vcf

{ echo "/ifs/home/kevin/WGS/data/VCF/pbsv2/hg38.DNA17-06168.pbsv2.DEL.ONLY.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/novaseq/manta+lumpy+delly_DNA17-06168.merged500pb.DEL.ONLY.vcf"; } > pacbio+novaseq_DNA17-06168.DEL.fof

SURVIVOR merge pacbio+novaseq_DNA17-06168.DEL.fof 500 1 1 0 0 20 pacbio+novaseq_DNA17-06168.DEL.vcf

###### DUPLICATION

{ echo "/ifs/home/kevin/WGS/data/VCF/pbsv2/hg38.DNA17-06166.pbsv2.DUP.ONLY.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/novaseq/manta+lumpy+delly_DNA17-06166.merged500pb.ALL.INS.vcf"; } > pacbio+novaseq_DNA17-06166.DUP.fof

SURVIVOR merge pacbio+novaseq_DNA17-06166.DUP.fof 500 1 1 0 0 20 pacbio+novaseq_DNA17-06166.DUP.vcf

{ echo "/ifs/home/kevin/WGS/data/VCF/pbsv2/hg38.DNA17-06167.pbsv2.DUP.ONLY.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/novaseq/manta+lumpy+delly_DNA17-06167.merged500pb.ALL.INS.vcf"; } > pacbio+novaseq_DNA17-06167.DUP.fof

SURVIVOR merge pacbio+novaseq_DNA17-06167.DUP.fof 500 1 1 0 0 20 pacbio+novaseq_DNA17-06167.DUP.vcf

{ echo "/ifs/home/kevin/WGS/data/VCF/pbsv2/hg38.DNA17-06168.pbsv2.DUP.ONLY.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/novaseq/manta+lumpy+delly_DNA17-06168.merged500pb.ALL.INS.vcf"; } > pacbio+novaseq_DNA17-06168.DUP.fof

SURVIVOR merge pacbio+novaseq_DNA17-06168.DUP.fof 500 1 1 0 0 20 pacbio+novaseq_DNA17-06168.DUP.vcf
```

Sort and index VCF to be ready for filtering and annotation :
```
#### Files have to be bgzipped, sorted and indexed
bcftools sort -O z BvB41_child_xAtlas.recode.vcf.gz -o BvB41_child_xAtlas.recode.sorted.vcf.gz
bcftools sort -O z BvB41_father_xAtlas.recode.vcf.gz -o BvB41_father_xAtlas.recode.sorted.vcf.gz
bcftools sort -O z BvB41_mother_xAtlas.recode.vcf.gz -o BvB41_mother_xAtlas.recode.sorted.vcf.gz

12:42:26 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools index BvB41_child_xAtlas.recode.sorted.vcf.gz
12:43:26 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools index BvB41_father_xAtlas.recode.sorted.vcf.gz
12:43:36 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools index BvB41_mother_xAtlas.recode.sorted.vcf.gz

17:05:15 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools index BvB41_child_xAtlas.recode.sorted.DUP.vcf.gz
17:06:16 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools index BvB41_father_xAtlas.recode.sorted.DUP.vcf.gz
17:06:29 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools index BvB41_mother_xAtlas.recode.sorted.DUP.vcf.gz

14:58:52 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools index BvB41_child_xAtlas.recode.sorted.notDUP.vcf.gz
15:01:16 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools index BvB41_father_xAtlas.recode.sorted.notDUP.vcf.gz
15:01:23 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools index BvB41_mother_xAtlas.recode.sorted.notDUP.vcf.gz

12:46:33 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq }-> bcftools sort -O z pacbio+novaseq_DNA17-06166.DEL.vcf.gz -o pacbio+novaseq_DNA17-06166.DEL.sorted.vcf.gz
12:46:33 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq }-> bcftools sort -O z pacbio+novaseq_DNA17-06167.DEL.vcf.gz -o pacbio+novaseq_DNA17-06167.DEL.sorted.vcf.gz
12:46:33 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq }-> bcftools sort -O z pacbio+novaseq_DNA17-06168.DEL.vcf.gz -o pacbio+novaseq_DNA17-06168.DEL.sorted.vcf.gz

12:47:21 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq }-> bcftools sort -O z pacbio+novaseq_DNA17-06166.DUP.vcf.gz -o pacbio+novaseq_DNA17-06166.DUP.sorted.vcf.gz
12:47:21 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq }-> bcftools sort -O z pacbio+novaseq_DNA17-06167.DUP.vcf.gz -o pacbio+novaseq_DNA17-06167.DUP.sorted.vcf.gz
12:47:21 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq }-> bcftools sort -O z pacbio+novaseq_DNA17-06168.DUP.vcf.gz -o pacbio+novaseq_DNA17-06168.DUP.sorted.vcf.gz

for i in *.sorted.vcf.gz ; do bcftools index $i ; done
```

Filter and annotate VCF by AnnotSV :

```
##### DELETION

./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06166.DEL.sorted.vcf.gz -genomeBuild GRCh38 -typeOfAnnotation full -bedtools /cm/shared/apps/bioinf/bedtools/2.25.0/bin/bedtools -outputDir  ~/WGS/data/VCF/pacbio-novaseq/DEL -vcfFiles ~/WGS/data/VCF/xatlasSNV/BvB41_child_xAtlas.recode.sorted.vcf.gz

./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06167.DEL.sorted.vcf.gz -genomeBuild GRCh38 -typeOfAnnotation full -bedtools /cm/shared/apps/bioinf/bedtools/2.25.0/bin/bedtools -outputDir  ~/WGS/data/VCF/pacbio-novaseq/DEL -vcfFiles ~/WGS/data/VCF/xatlasSNV/BvB41_father_xAtlas.recode.sorted.vcf.gz

./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06168.DEL.sorted.vcf.gz -genomeBuild GRCh38 -typeOfAnnotation full  -bedtools /cm/shared/apps/bioinf/bedtools/2.25.0/bin/bedtools -outputDir  ~/WGS/data/VCF/pacbio-novaseq/DEL -vcfFiles ~/WGS/data/VCF/xatlasSNV/BvB41_mother_xAtlas.recode.sorted.vcf.gz

##### DUPLICATION

./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06166.DUP.sorted.vcf.gz -bedtools /ifs/home/kevin/bedtools2/bin/bedtools -SVinputInfo 1 -outputDir  ~/WGS/data/VCF/pacbio-novaseq/DUP2 -vcfFiles "~/WGS/data/VCF/xatlasSNV/BvB41_child_xAtlas.recode.sorted.DUP.vcf.gz ~/WGS/data/VCF/xatlasSNV/BvB41_child_xAtlas.recode.sorted.notDUP.vcf.gz"

./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06167.DUP.sorted.vcf.gz -genomeBuild GRCh38 -typeOfAnnotation full -bedtools /ifs/home/kevin/bedtools2/bin/bedtools -SVinputInfo 1 -outputDir  ~/WGS/data/VCF/pacbio-novaseq/DUP2 -vcfFiles "~/WGS/data/VCF/xatlasSNV/BvB41_father_xAtlas.recode.sorted.DUP.vcf.gz ~/WGS/data/VCF/xatlasSNV/BvB41_father_xAtlas.recode.sorted.notDUP.vcf.gz"

./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06168.DUP.sorted.vcf.gz -genomeBuild GRCh38 -typeOfAnnotation full -bedtools /ifs/home/kevin/bedtools2/bin/bedtools -SVinputInfo 1 -outputDir  ~/WGS/data/VCF/pacbio-novaseq/DUP2 -vcfFiles "~/WGS/data/VCF/xatlasSNV/BvB41_mother_xAtlas.recode.sorted.DUP.vcf.gz ~/WGS/data/VCF/xatlasSNV/BvB41_mother_xAtlas.recode.sorted.DUP.vcf.gz"
```

Add SV size to annotated tsv from AnnotSV, yet respect excel format for analysis used for DELETION

```
### DUP
16:40:54 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> echo "SVSize" > child_INS_SVsize.tsv
16:43:02 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> echo "SVSize" > father_INS_SVsize.tsv
16:43:09 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> echo "SVSize" > child_INS_SVsize.tsv
16:43:13 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> grep -o -P '(?<=AVGLEN\=).*(?=\;SVTYPE)' pacbio+novaseq_DNA17-06166.DUP.sorted.annotated.tsv >> child_INS_SVsize.tsv
16:43:24 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> grep -o -P '(?<=AVGLEN\=).*(?=\;SVTYPE)' pacbio+novaseq_DNA17-06167.DUP.sorted.annotated.tsv >> father_INS_SVsize.tsv
16:43:55 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> grep -o -P '(?<=AVGLEN\=).*(?=\;SVTYPE)' pacbio+novaseq_DNA17-06168.DUP.sorted.annotated.tsv >> mother_INS_SVsize.tsv

17:10:43 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> paste -d '\t' pacbio+novaseq_DNA17-06166.DUP.sorted.annotated.tsv child_INS_SVsize.tsv | cut -f-3,5-6,9- | grep -v "split" > pacbio+novaseq_DNA17-06166.DUP.sorted.annotated.fullonly.SVsize.tsv
17:24:46 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> paste -d '\t' pacbio+novaseq_DNA17-06167.DUP.sorted.annotated.tsv father_INS_SVsize.tsv | cut -f-3,5-6,9- | grep -v "split" > pacbio+novaseq_DNA17-06167.DUP.sorted.annotated.fullonly.SVsize.tsv
17:25:11 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> paste -d '\t' pacbio+novaseq_DNA17-06168.DUP.sorted.annotated.tsv mother_INS_SVsize.tsv | cut -f-3,5-6,9- | grep -v "split" > pacbio+novaseq_DNA17-06168.DUP.sorted.annotated.fullonly.SVsize.tsv

### DUP2
16:19:21 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP2 }-> echo "SVSize" > child_INS_SVsize.tsv
16:19:23 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP2 }-> echo "SVSize" > father_INS_SVsize.tsv
16:19:30 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP2 }-> echo "SVSize" > child_INS_SVsize.tsv
16:19:36 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP2 }-> grep -o -P '(?<=AVGLEN\=).*(?=\;SVTYPE)' pacbio+novaseq_DNA17-06166.DUP.sorted.annotated.tsv >> child_INS_SVsize.tsv
16:19:42 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP2 }-> grep -o -P '(?<=AVGLEN\=).*(?=\;SVTYPE)' pacbio+novaseq_DNA17-06167.DUP.sorted.annotated.tsv >> father_INS_SVsize.tsv
16:19:47 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP2 }-> grep -o -P '(?<=AVGLEN\=).*(?=\;SVTYPE)' pacbio+novaseq_DNA17-06168.DUP.sorted.annotated.tsv >> mother_INS_SVsize.tsv

16:36:43 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP2 }-> paste -d '\t' pacbio+novaseq_DNA17-06166.DUP.sorted.annotated.tsv child_INS_SVsize.tsv | cut -f-3,5-6,9-15,18-  | grep -v "split"  > pacbio+novaseq_DNA17-06166.DUP.sorted.annotated.fullonly.SVsize.tsv
16:39:08 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP2 }-> paste -d '\t' pacbio+novaseq_DNA17-06167.DUP.sorted.annotated.tsv father_INS_SVsize.tsv | cut -f-3,5-6,9-15,18-  | grep -v "split"   > pacbio+novaseq_DNA17-06167.DUP.sorted.annotated.fullonly.SVsize.tsv
16:39:36 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP2 }-> paste -d '\t' pacbio+novaseq_DNA17-06168.DUP.sorted.annotated.tsv mother_INS_SVsize.tsv | cut -f-3,5-6,9-15,18-  | grep -v "split" > pacbio+novaseq_DNA17-06168.DUP.sorted.annotated.fullonly.SVsize.tsv

```


Number of heterozygous vs homozygous variant SNV :

```
####
17:25:54 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view -g hom -c1  BvB41_child_xAtlas.recode.sorted.vcf.gz | grep -v "#" - | wc -l
[W::bcf_hdr_check_sanity] PL should be declared as Number=G
1811424
17:26:11 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view -g het -c1  BvB41_child_xAtlas.recode.sorted.vcf.gz | grep -v "#" - | wc -l
[W::bcf_hdr_check_sanity] PL should be declared as Number=G
3088763
17:26:33 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view -c1  BvB41_child_xAtlas.recode.sorted.vcf.gz | grep -v "#" - | wc -l
[W::bcf_hdr_check_sanity] PL should be declared as Number=G
4900187
####
17:26:55 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view -g hom -c1  BvB41_father_xAtlas.recode.sorted.vcf.gz | grep -v "#" - | wc -l
[W::bcf_hdr_check_sanity] PL should be declared as Number=G
1798137
17:37:42 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view -g het -c1  BvB41_father_xAtlas.recode.sorted.vcf.gz | grep -v "#" - | wc -l
[W::bcf_hdr_check_sanity] PL should be declared as Number=G
3104193
17:38:01 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view -c1  BvB41_father_xAtlas.recode.sorted.vcf.gz | grep -v "#" - | wc -l
[W::bcf_hdr_check_sanity] PL should be declared as Number=G
4902330
####
17:38:23 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view -g hom -c1  BvB41_mother_xAtlas.recode.sorted.vcf.gz | grep -v "#" - | wc -l
[W::bcf_hdr_check_sanity] PL should be declared as Number=G
1789717
17:39:24 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view -g het -c1  BvB41_mother_xAtlas.recode.sorted.vcf.gz | grep -v "#" - | wc -l
[W::bcf_hdr_check_sanity] PL should be declared as Number=G
3119161
17:39:41 kevin::login02 { ~/WGS/data/VCF/xatlasSNV }-> bcftools view -c1  BvB41_mother_xAtlas.recode.sorted.vcf.gz | grep -v "#" - | wc -l
[W::bcf_hdr_check_sanity] PL should be declared as Number=G
4908878
```

Merge sample for WGS Trio mode analysis, sort and index :

```
###### DELETION Trio
{ echo "/ifs/home/kevin/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06166.DEL.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06167.DEL.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06168.DEL.vcf";} > pacbio+novaseq_trio.DEL.fof

SURVIVOR merge pacbio+novaseq_trio.DEL.fof 500 1 1 0 0 20 pacbio+novaseq_trio.DEL.vcf
bgzip pacbio+novaseq_trio.DEL.vcf
bcftools sort -O z pacbio+novaseq_trio.DEL.vcf.gz -o pacbio+novaseq_trio.DEL.sorted.vcf.gz
bcftools index pacbio+novaseq_trio.DEL.sorted.vcf.gz

###### DUPLICATION Trio
{ echo "/ifs/home/kevin/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06166.DUP.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06167.DUP.vcf";
  echo "/ifs/home/kevin/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_DNA17-06168.DUP.vcf";} > pacbio+novaseq_trio.DUP.fof

SURVIVOR merge pacbio+novaseq_trio.DUP.fof 500 1 1 0 0 20 pacbio+novaseq_trio.DUP.vcf
bgzip pacbio+novaseq_trio.DUP.vcf
bcftools sort -O z pacbio+novaseq_trio.DUP.vcf.gz -o pacbio+novaseq_trio.DUP.sorted.vcf.gz
bcftools index pacbio+novaseq_trio.DUP.sorted.vcf.gz
```

Filter Trio with only good CNV :

```
###### DELETION
### ADD CHR to chromosome number as AnnotSV removed them
10:49:58 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DEL }-> sed 's/^/chr/g' good_del_cnv.txt > good_del_cnv_chr.txt

10:55:50 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq }-> bcftools view ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DEL.sorted.vcf.gz -T ~/WGS/data/VCF/pacbio-novaseq/DEL/good_del_cnv_chr.txt | grep -v "SUPP_VEC=0" > ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DEL.sorted.filtered.vcf
bgzip ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DEL.sorted.filtered.vcf
bcftools index ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DEL.sorted.filtered.vcf.gz

./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DEL.sorted.filtered.vcf.gz -bedtools /ifs/home/kevin/bedtools2/bin/bedtools -outputDir  ~/WGS/data/VCF/pacbio-novaseq/DEL -vcfFiles ~/WGS/data/VCF/xatlasSNV/BvB41_child_xAtlas.recode.sorted.vcf.gz

##### DUPLICATION
12:13:11 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> bcftools view ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DUP.sorted.vcf.gz -T ~/WGS/data/VCF/pacbio-novaseq/DUP/good_dup_cnv_chr.txt | grep -v "SUPP_VEC=0" > ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DUP.sorted.filtered.vcf
12:13:29 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> bgzip ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DUP.sorted.filtered.vcf
12:14:06 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> bcftools index ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DUP.sorted.filtered.vcf.gz

12:17:00 kevin::login02 { ~ }-> ./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DUP.sorted.filtered.vcf.gz -bedtools /ifs/home/kevin/bedtools2/bin/bedtools -SVinputInfo 1 -outputDir  ~/WGS/data/VCF/pacbio-novaseq/DUP -vcfFiles ~/WGS/data/VCF/xatlasSNV/BvB41_child_xAtlas.recode.sorted.DUP.vcf.gz
```

Bcftools filter for de novo is not working (don't know why), so I did it manually again
```
10:10:49 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq }-> bcftools view --output-type v  --include 'INFO/SUPP_VEC="100"' pacbio+novaseq_trio.DUP.sorted.vcf.gz > pacbio+novaseq_trio.DUP.sorted.denovo.vcf
10:11:13 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq }-> bgzip pacbio+novaseq_trio.DUP.sorted.denovo.vcf
./AnnotSV_1.2/bin/AnnotSV -SVinputFile ~/WGS/data/VCF/pacbio-novaseq/pacbio+novaseq_trio.DUP.sorted.vcf.gz -bedtools /ifs/home/kevin/bedtools2/bin/bedtools -SVinputInfo 1 -outputDir  ~/WGS/data/VCF/pacbio-novaseq/DUP -vcfFiles ~/WGS/data/VCF/xatlasSNV/BvB41_child_xAtlas.recode.sorted.DUP.vcf.gz

15:51:31 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> echo "SVSize" > denovo_INS_SVsize.tsv
15:51:38 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> grep -o -P '(?<=AVGLEN\=).*(?=\;SVTYPE)' pacbio+novaseq_trio.DUP.sorted.denovo.annotated.tsv >> denovo_INS_SVsize.tsv
15:52:06 kevin::login02 { ~/WGS/data/VCF/pacbio-novaseq/DUP }-> paste -d '\t' pacbio+novaseq_trio.DUP.sorted.denovo.annotated.tsv  denovo_INS_SVsize.tsv | cut -f-3,5-6,9- > pacbio+novaseq_trio.DUP.sorted.denovo.annotated.SVSize.tsv
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
