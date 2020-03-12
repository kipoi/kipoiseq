"""
skript which creates vcf file from csv file with the mutations from gnomAD
the vcf file is only for testing!
after that the file should be preprocessed in order to be usable with the data loader:
bgzip file.vcf
tabix file.vcf.gz
ready!
"""





vcf_head = """##fileformat=VCFv4.0
##phasing=partial
##contig=<ID=1,length=16>
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA00001	NA00002	NA00003
"""

csv_file_path = '../../file.vcf'

with open("vcf_file_for_testing_missense_mutations.vcf", 'w+') as vcf_file:
    vcf_file.write(vcf_head)
    with open(csv_file_path, 'r+') as csv_file:
        for line in csv_file:
            if 'Position' in line:
                continue
            line = line.split(',')
            vcf_file.write(f"{line[0]}	{line[1]}	.	{line[3]}	{line[4]}	9.6	.	.	GT:HQ	./.:51,51	0|0:51,51	1/1:51,51\n")