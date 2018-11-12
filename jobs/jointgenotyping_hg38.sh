java -Dspring.config.location=/ifs/software/inhouse/Buzzzz/PRODUCTION/buzzzRecipeSender/PRODUCTION/config/application-TURBO.properties -Dspring.config.name=TURBO -cp "/ifs/software/inhouse/Buzzzz/PRODUCTION/buzzzRecipeSender/PRODUCTION/lib/*" org.umcn.ngs.recipesender.SampleFinder --type=joint_genotyping --partition=RESEARCH_NORMAL --outFolder=jointGenotype --include=/ifs/data/research/novaseq/genomes/old_samples/BvB41_child/remapping-040918/jointGenotyping/samples.txt --fastaFile=/ifs/data/lib/genomes/human/GRCh38/hs_ref_GRCh38.p2_all_contigs.fa --runType=research --callers=haplotypecaller --callerMerging=none 