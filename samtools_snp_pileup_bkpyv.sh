#/bin/bash


#  Pipeline for the samtools workflow for calling variants.
#
#
#  Run within root level of samples analysis dir.
#  No multi mode yet.
#



##################################################################
#  SOFTWARES
#  Provide locations of the softwares to be used. 

PICARD="java -jar /usr/local/src/picard-tools-1.119"
SAMTOOLS="samtools"
GATK="java -jar /usr/local/src/gatk/GenomeAnalysisTK.jar"


#################################################################
#   FILES 
# Aligned bam file and reference fasta file

mode=$1
reference=$2
bam=$3


if [ $# -ne 3 ] ; then
    echo 'Input: [run mode: single/multi] [reference.fasta] [sample.bam]'
    exit 0
fi

# Clip any filepath or suffix naming from sample
readGroupSMTag=`echo $bam | sed 's/_L001_R1.*.bam//' | sed 's/.*\///' `

# Set folder output
sampleDir=`echo $bam | sed 's/\/.*//'`


######################################################################
## If mode flag used, run as single or multimode
if [ $mode == "single" ]
then

	#####################################################################
	#   Check if the bam file satisfies the requirements of the GATK

	$SAMTOOLS view -H $bam > $sampleDir/$readGroupSMTag"_alpha.txt"
	head_file=$readGroupSMTag"_alpha.txt"
	grep '^@SQ' $head_file > $sampleDir/$readGroupSMTag"_beta.txt"
	file=$readGroupSMTag"_beta.txt"

	echo "Created $file"

	#less $file
	# check if the file is sorted.
	T="sort -c -t':' -nk2 $file"
	if [ "$T" ]; then
    	echo "The file is sorted!"

	else
    	echo "The file is not sorted."
    	echo "Sorting........."
    	$PICARD/SortSam.jar INPUT=$bam OUTPUT="${bam%.bam}.sorted.bam" SORT_ORDER=coordinate

	fi


	# Check if the file contains RG information.
	if grep -q '^@RG' $head_file; then
		echo "The file contains RG information."
	else
		echo "The file does not contain RG information and the GATK will not work!"
		echo "Attempting fix..."
		echo "Copying original bam input as $bam.original"
	
		cp $bam $bam".original"	

		$PICARD/AddOrReplaceReadGroups.jar \
		INPUT= $bam \
		OUTPUT=$bam".temp" RGID=1 RGLB=Library1 RGPL=ILLUMINA RGPU=1 RGSM=$readGroupSMTag RGCN=AA RGDS=AD

		# Renaming bam variable to new RG bam
		mv $bam".temp" $bam
	fi


	###################################################################
	#   Prepping reference geonome fasta file for GATK
	echo "Prepping reference geonome fasta file for GATK....."
	
	# Create sequence dictionary using Picard Tools.
	# the following command produces a SAM-style header file describing the contents of our fasta file.
	$PICARD/CreateSequenceDictionary.jar \
	reference=$reference \
	OUTPUT=$reference".dict"

	# Looking for reference without fasta suffix
	newLib=`echo $reference".dict" | sed 's/.fasta//'`
	mv $reference".dict" $newLib

	echo "Created sequence dictionary for the reference genome."
	echo "Indexing the reference genome...."

	# Create the fasta index file.
	# The index file describes byte offset in the fasta file for each contig. It is a text file with one record
	# per line for each of the fasta contigs. Each record is of the type -
	# contig, size, location, basePerLine, bytesPerLine
	$SAMTOOLS faidx $reference

	echo "Reference genome is now ready for GATK..."

	###############################################################
	## Summary Statistics
	## Used in debug

	$PICARD/MeanQualityByCycle.jar \
	INPUT=$bam \
	CHART_OUTPUT="$sampleDir/$readGroupSMTag_mean_quality_by_cycle.pdf" \
	OUTPUT="$sampleDir/$readGroupSMTag_read_quality_by_cycle.txt" \
	reference_SEQUENCE=$reference


	$PICARD/QualityScoreDistribution.jar \
	INPUT=$bam \
	CHART_OUTPUT="$sampleDir/$readGroupSMTag_mean_quality_overall.pdf" \
	OUTPUT="$sampleDir/$readGroupSMTag_read_quality_overall.txt" \
	reference_SEQUENCE=$reference


	$PICARD/CollectWgsMetrics.jar \
	INPUT=$bam OUTPUT="$sampleDir/$readGroupSMTag_stats_picard.txt" \
	reference_SEQUENCE=$reference \
	MINIMUM_MAPPING_QUALITY=20 \
	MINIMUM_BASE_QUALITY=20

	
	#############################################################
	# Mark duplicate reads.
	echo "Mark the duplicates in the bam file."

	$PICARD/MarkDuplicates.jar INPUT=$bam OUTPUT="${bam%.bam}_dups_marked.bam" \
	METRICS_FILE="${bam%.bam}_dups_metrics.txt" REMOVE_DUPLICATES=false

	echo "Index the dup-marked bam file,${bam%.bam}_dups_marked.bam"

	$SAMTOOLS index "${bam%.bam}_dups_marked.bam"

	newBam="${bam%.bam}_dups_marked.bam"


	#############################################################
	# GATK Data Pre-Processing

	# Step 1 - Local realignment around indels.
	# Create a target list of intervals to be realigned.
	echo "Creating a target list of intervals to be realigned...."

	$GATK \
	-T RealignerTargetCreator \
	-R $reference \
	-I $newBam \
	-o "${bam%.bam}_target_intervals.list"

	# Run the local realignment
	echo "Local realignment..."

	$GATK \
	-T IndelRealigner \
	-R $reference \
	-I $newBam \
	-targetIntervals "${bam%.bam}_target_intervals.list" \
	-o "${bam%.bam}_realigned_reads.bam"

	echo "Indexing the realigned bam file..."

	# Tidy up files
	rm "${bam%.bam}_dups_marked.bam"

	# Create a new index file.
	$SAMTOOLS index "${bam%.bam}_realigned_reads.bam"

	# Step 2 - Base recalibration (fixes them so they better reflect the probability of mismatching the genome).
	# Analyze patterns of covariation in the sequence.
	echo "Base recalibration...skipped as not needed"

	#$GATK \
	#-T BaseRecalibrator \
	#-R $reference \
	#-I "${bam%.bam}_realigned_reads.bam" \
	#-o "${bam%.bam}_recal_data.table"

	# Apply recalibration to the sequence data.
	echo "Recalibrating the sequence data...skipped as not needed"

	#$GATK \
	#-T PrintReads \
	#-R $reference \
	#-I "${bam%.bam}_realigned_reads.bam" \
	#-BQSR "${bam%.bam}_recal_data.table" \
	#-o "${bam%.bam}_recal_reads.bam"

	
	###########################################################################
	# Samtools Variant Calling
	
	echo "Running mpileup..."

	$SAMTOOLS mpileup -ugf $reference "${bam%.bam}_realigned_reads.bam" | bcftools call -vmO z -o "${bam%.bam}_realigned_reads.vcf.gz"

	echo "Running tabix..."

	/usr/local/bin/tabix -p vcf "${bam%.bam}_realigned_reads.vcf.gz"
	
	echo "Plotting snps..."

	bcftools stats -F $reference -s - "${bam%.bam}_realigned_reads.vcf.gz" > "${bam%.bam}_realigned_reads.vcf.gz.stats"

	# Filtering VCF
	echo "Generating soft filter VCF..."

	bcftools filter -O z -o "${bam%.bam}_realigned_reads.filtered.vcf.gz" -s LOWQUAL -i'%QUAL>30 && DP>5' "${bam%.bam}_realigned_reads.vcf.gz"

	# Create mask fasta file based on bed coverage cutoff X reads.
   	echo "Masking reference fasta file with low cov across sequenced genome (<=5 read depths).."
	
	awk '($3<=5) {print $1"\t"$2"\t"$2}' "${bam}.bed_coverage" > "${bam}.lowcov.bed"
	bedtools maskfasta -fi $reference -bed "${bam}.lowcov.bed" -fo $sampleDir/"${reference%.fasta}.masked.fasta"
	
	# Convert to fasta file with gentypes inserted, needing bcf files of the filtered vcr
	echo "Generate fasta with filtered genotype positions inserted..."

	bcftools view -Ov -f .,PASS "${bam%.bam}_realigned_reads.filtered.vcf.gz" > "${bam%.bam}_realigned_reads.filtered_passed.vcf"
	cat "${bam%.bam}_realigned_reads.filtered_passed.vcf" | bcftools convert -O b -o "${bam%.bam}_realigned_reads.filtered_passed.bcf"

	bcftools index "${bam%.bam}_realigned_reads.filtered_passed.bcf"

	cat "$sampleDir/${reference%.fasta}.masked.fasta" | bcftools consensus "${bam%.bam}_realigned_reads.filtered_passed.bcf" > "${bam%.bam}_realigned_reads.filtered_passed.fasta"

	echo "Generating VP1 region only consensus..."
	samtools faidx "$sampleDir/${reference%.fasta}.masked.fasta" gi\|9627180\|ref\|NC_001538.1\|:1564-2652 | bcftools consensus "${bam%.bam}_realigned_reads.filtered_passed.bcf" > "${bam%.bam}_realigned_reads.filtered_passed.VP1.fasta"
	
	# Rename fasta header to sample names
	echo "Renaming genotype consensus fasta headers using $readGroupSMTag prefix of the bam file..."

	sed "s/>.*/>$readGroupSMTag complete/" "${bam%.bam}_realigned_reads.filtered_passed.fasta" > "${bam%.bam}_realigned_reads.filtered_passed.fasta.tmp"
 	sed "s/>.*/>$readGroupSMTag VP1/" "${bam%.bam}_realigned_reads.filtered_passed.VP1.fasta" > "${bam%.bam}_realigned_reads.filtered_passed.VP1.fasta.tmp"

	mv "${bam%.bam}_realigned_reads.filtered_passed.fasta.tmp" "${bam%.bam}_realigned_reads.filtered_passed.fasta"
	mv "${bam%.bam}_realigned_reads.filtered_passed.VP1.fasta.tmp" "${bam%.bam}_realigned_reads.filtered_passed.VP1.fasta"

	echo "File complete..."
fi

printf "Finished..."

