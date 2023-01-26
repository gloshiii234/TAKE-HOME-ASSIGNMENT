##KIRABO GLORIA
##2022/HD07/2043U
##BASH ASSIGNMENT


##QUESTION ONE-1. Describe the format of the file and the data stored
    The file is a VCF file and it is a text file that contains infoemation about genetic variations in a genome.
    The type of data stored includes information about the position of the variant on the genome, the reference,
    the alternate alleles and qualitynmetrics such as read depths and quality scores.

##QUESTION TWO-2. What does the header section of the file contain
    -meta information about the file(format version, sample information,reference genome)


##QUESTION THREE-3. How many samples are in the file
    bcftools query -l sample.vcf.gz 
    bcftools query -l sample.vcf.gz  | wc -l
    *there are 6 samples in the file.


##QUESTION FOUR-4. How many variants are in the file

    bcftools query -f '%ALT\n' sample.vcf.gz | wc -l
    The file contains 398246 variants.


##QUESTION FIVE-5. How would you extract the chromosome, position, QualByDepth and#RMSMappingQuality fields? Save the output to a tab-delimited file

    bcftools query -f '%CHROM\t%POS[\t%QD;%MQ]\n' sample.vcf > sample_1.vcf


##QUESTION SIX-6. Extract data that belongs to chromosomes 2,4 and MT
    awk '$1=="2" || $1=="4" || $1=="MT"' sample.vcf


##QUESTION SEVEN-7.Print out variants that do not belong to chr20:1-30000000

    awk '$1 != "20" || ($1 == "chr20" && ($2 < 1 || $2 > 30000000)) {print $1, $2, $4, $5}' sample.vcf > sample_2_variants.vcf


##QUESTION EIGHT-8. Extract variants that belong to SRR13107019

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' -s SRR13107019 sample.vcf > SRR13107019_output.vcf 


##QUESTION NINE-9. Filter out variants with a QualByDepth above 7
    awk -F '\t' '{if ($6>=7) print $0}' sample.vcf > QBDoutput.vcf


##QUESTION TEN-10. How many contigs are referred to in the file. Check the header section

    grep -c "^##contig" sample.vcf
    There are 2211 contigs in the file.


##QUESTION ELEVEN-11. Comment on the eighth and ninth columns of the file
        The eighth column is made up of genotype format.
        the ninth column contains the genotype data for each sample in the file. 



##QUESTION TWELVE-12. Extract data on the read depth of called variants for sample SRR13107018

bcftools querry -f '%DP\n' -s SRR13107018 sample.vcf > SRR13107018_output.vcf


##QUESTION THIRTEEN-13. Extract data on the allele frequency of alternate alleles. Combine this data with the #-chromosome and position of the alternate allele

bcftools query -f '%CHROM\t%POS\t%AF\n' sample.vcf > alternate_alleles.vcf



######## THE SAM FILE ########

##QUESTION ONE-1. Describe the format of the file and the data stored
    The file is a SAM file and it is a text based alignment format which supports single- and paired-end reads produced by different sequencing
platforms which saves alignment information of short reads mapped against reference sequences

##QUESTION TWO-2. What does the header section of the file contain

    The header contains metadata about the file, such as the reference genome used, and information about the samples and sequencing runs in the file



##QUESTION THREE-3. How many samples are in the file

    grep -E '^@RG' sample.sam | awk '{print $2}' | awk -F ':' '{print $2}' | sort | uniq | wc -l
    There are 249


##QUESTION FOUR-4. How many alignments are in the file

    samtools view -F 4 sample.sam | wc -l
    There are 36142


##QUESTION FIVE-5. Get summary statistics for the alignments in the file

    samtools flagstat sample.sam > sam5.sam



##QUESTION SIX_6. Count the number of fields in the file
    awk '{print NF}' sample.sam | sort -nu | wc -l  
    There are 4


##QUESTION SEVEN-7. Print all lines in the file that have @SQ and sequence name tag beginning with NT_

    grep "@SQ.*NT_" sample.sam 


##QUESTION EIGHT-8. Print all lines in the file that have @RG and LB tag beginning with Solexa

    grep "@RG.*LB:Solexa" sample.sam

##QUESTION NINE-9. Extract primarily aligned sequences and save them in another file

    awk '$1 !~ /^@/ && $2 == "99" || $2 == "83"' sample.sam > primary_sequences.sam


##QUESTION TEN-10. Extract alignments that map to chromosomes 1 and 3. Save the output in #BAM format

    awk '$1 !~ /^@/ && ($3 == "1" || $3 == "3")' sample.sam | samtools view -Sb - > alignments_chr1_3.bam


##QUESTION ELEVN-11. How would you obtain unmapped reads from the file
    samtools view -f 4 sample.sam > unmapped_reads.sam


##QUESTION TWELVE-12. How many reads are aligned to chromosome 4

    grep -c "^4\t" sample.sam


##QUESTION THIRTEEN-13. Comment of the second and sixth column of the file
    The second column of the file refers to the read name.
    The sixth column refers to the read flag which is an integer value that encodes information about the alignment of the read(mapped or unmapped).



##14. Extract all optional fields of the file and save them in “optional_fields.txt”

    awk '{for(i=11;i<=NF;i++) print $i}' sample.sam > optional_fields.txt

