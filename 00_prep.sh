# convert FASTA+QUAL to FASTQ.GZ with BBMAP
for a in *.fasta
do
	~/sharedscratch/apps/bbmap/reformat.sh in=$a qfin=$a.qual out=${a/fasta/fastq.gz} qin=33
done

# MID FASTA editing
sed -i -e 's/MID/>MID/g' -e 's/\s/\n/g' MIDs.fa

# primer checks
# 450X950

F: ATCAGACACG#GATAGAATCGTACGTGCATA#AATTACCCAATCCTGACACAGG
R: GTCCCTATTAATCATTACCCTGG (rev compl:  CCAGGGTAATGATTAATAGGGAC)
