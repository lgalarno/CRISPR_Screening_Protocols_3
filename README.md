###### CRISPR Screening_Protocols

Adapted from :
Scripts from the Joung et al Nature Protocols 2016 manuscript on knockout and transcriptional activation screening
https://github.com/fengzhanglab/Screening_Protocols_manuscript

Python packages required are:
biopython
numpy
scipy
twobitreader

##### design_library
Download the genome 2bit file that the target gene coordinates corresponds to from the UCSC Genome Browser (http://hgdownload.cse.ucsc.edu/downloads.html). 

For seqmap, install the version 1.0.13 source code for all platforms (http://www-personal.umich.edu/~jianghui/seqmap/ or available here) and compile with g++ -O3 -m64 -o seqmap match.cpp. Place the seqmap folder in the same folder as the python script design_library.py.

About seqmap:
Author

SeqMap was started in the Wong Lab and is developed and maintained by Hui Jiang.

Reference

Jiang, H., Wong, W.H. (2008) SeqMap: Mapping Massive Amount of Oligonucleotides to the Genome, Bioinformatics, 24(20). [online]




example of genes.csv

"

name	chrom	start	end

EGFR	chr7	55086525	55086725

LPAR5	chr12	6745297	6745497

GPR35	chr2	241544625	241544825

"


parameters:

-o	Output csv file with names for target genes, spacer sequences, spacer orientations, chromosome locations, cleavage site locations, off-target scores, and oligo library sequences in columns from left to right,	default:final_guides.csv

-i	Prefix of input genome 2bit file,	default:hg38

-g	Target genes csv file with the gene name, chromosome, start of the targeted region, and end of the targeted region in columns from left to right,	default:genes.csv

-gc	Minimum GC content required for an sgRNA spacer sequence,	default:25

-s	Minimum spacing required between cleavage sites of sgRNAs targeting the same genomic region,	default:20

-n	Maximum number of guides selected targeting each gene in the target genes csv file,	default:3

-db	Use an existing off-target database constructed from a previous custom library design for a new library, default:False

-gecko or -sam	Specify the type of library and add the respective flanking sequences to the spacers for the oligo library synthesis, default: none


##### design_targeted_library

parameters:

-o	Output csv file with names for target genes, corresponding spacer sequences, and oligo library sequences in columns from left to right,	default:oligos.csv

-l	Annotated library csv file with names in the first column and corresponding spacer sequences in the second column,	default:annotated_library.csv

-g	Target genes csv file with names of target genes,	target_genes.csv

-gecko or -sam	Specify the type of library and add the respective flanking sequences to the spacers for the oligo library synthesis, default: none


##### count_spacers

parameters:

Flag	Description	Default

-f	Fastq file containing NGS data for analysis,	default:NGS.fastq

-o	Output csv file with guide spacer sequences in the first column and respective read counts in the second column,	default:library_count.csv

-i	Input csv file with guide spacer sequences, default:library_sequences.csv

-no-g	Indicate absence of a guanine before the guide spacer sequence, 	default:guanine is present


##### calculate_indel

example of a sample sheet:

"

example1	example1.fastq	GCCCGATCGCTATATCCACG	TGTATATACCTCGCGCCTAACTGCCAGCTGACCACGCCGTACAGGTTCTGTTGTCTACTGATGCATTACATCTCCTTAGGGTACTTCCGCTGAGATCATTGCCCGATCGCTATATCCACGCGTTTGGCTCGTTCAACCAGTCCAACCGACTCGTTGGTCTGGTAATGTGCCCAAGTTAAGGTGAGTATGGACATGGCGGG	Experimental

example2	example2.fastq	CGAGATAAGTCAGCAGGGGC	CTCTTCTGCTCAAGCGAGTTCCCAGAGGTCCTTGCCGAGGGGGTTATATCGATCCACGAAACATAGTATGTAATACGAAAGTCATCGGCGCCTATGCCCTCGAGATAAGTCAGCAGGGGCTTTTCGTACATTTTCCAAGATTCGGGATTGACGTTGCATCGCAAGCTAATTGGTTACCATTAGACCCCAGTCCTCAGCCC	Experimental

example3	example3.fastq	CACCCACACCAACCGCAGAA	CTGGGTTTAACCGAGCTAGTCCTGAAGATCTTGAGTAACTGCACATGTAAATA

"

parameters:


-f	Indicates input file is fasta format,	default:Fastq file format

-a	Uses alternative hashing algorithm for calculation74,	default:Ratcliff-Obershelp based algorithm73

-o	Output file with calculated indel rates and statistics,	default:calc_indel_out.csv

-i	Input file with sample name, fasta or fastq file name, guide sequence, PCR target amplicon, and experimental or control in columns from left to right,	default:sample_sheet.csv

-v or -q	Increase or decrease reporting as script runs,	default:Standard reporting

-no-m	Does not perform MLE correction,	default:MLE correction is performed

