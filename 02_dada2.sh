#
# dada2 workflow including taxonomy annotation
# needs input directory, taxonomy database and output directory
#

module load bioconductor	# dada2_1.22.0 

# 450X950
sbatch -c 16 --mem 64G --wrap "Rscript 02_dada2.R trimmed_450X950 dada2_450X950 silva_132.18s.99_rep_set.dada2.fa.gz 1"
sbatch -c 16 --mem 64G --wrap "Rscript 02_dada2.R trimmed_450X950 dada2_450X950 silva_132.18s.99_rep_set.dada2.fa.gz 2"


# Phylogenetic tree of ASVs
conda activate trimal-1.4.1
for fasta in dada2_*/seqtab.*.fasta
do
	sbatch -c 8 --mem 16G --wrap "
		trimal -in $fasta.mafft -out $fasta.mafft.trimal.phy -automated1 -phylip
		~/sharedscratch/apps/iqtree-2.1.1-Linux/bin/iqtree2 -s $fasta.mafft.trimal.phy -nt 8 -bb 10000"
done

# ~/sharedscratch/apps/mafft/bin/ginsi --thread 8 $fasta > $fasta.mafft

#
# additional taxonomy with BLASTN
#

# NCBI search for 18S or SSU sequenecs
# 18s[All Fields] OR SSU[All Fields] AND (fungi[filter] OR protists[filter]) AND ("1"[SLEN] : "3000"[SLEN]) AND (is_nuccore[filter] AND ("1"[SLEN] : "3000"[SLEN]))
# download as FASTA and GenBank format (for taxonomy info)

# get GI numbers for all sequences
#
wget -O ncbi/tax.esearch.out 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=18s[All%20Fields]%20OR%20SSU[All%20Fields]%20AND%20(fungi[filter]%20OR%20protists[filter])%20AND%20(%221%22[SLEN]%20:%20%223000%22[SLEN])%20AND%20(is_nuccore[filter]%20AND%20(%221%22[SLEN]%20:%20%223000%22[SLEN]))&usehistory=y'

# parse query key and web accession
key=$(grep -m 1 -o "<QueryKey>[^<]*" ncbi/tax.esearch.out | sed 's/<QueryKey>//g')
web=$(grep -m 1 -o "<WebEnv>[^<]*" ncbi/tax.esearch.out | sed 's/<WebEnv>//g')
total=$(grep -m 1 -o "<Count>[^<]*" ncbi/tax.esearch.out | sed 's/<Count>//g')

# fetch FASTA records
for ((retstart = 0; retstart < $total; retstart += 9999)) {
	url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&WebEnv=$web&query_key=$key&retstart=$retstart&retmax=9999&rettype=fasta&retmode=text"
	wget -O ncbi/efetch.fasta.$retstart $url
}
cat ncbi/efetch.fasta.* > ncbi/18S_SSU_NCBI.fasta
rm ncbi/efetch.fasta.*

# fetch GB records for taxonomy
for ((retstart = 459954; retstart < $total; retstart += 9999)) {
	url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&WebEnv=$web&query_key=$key&retstart=$retstart&retmax=9999&rettype=gb&retmode=xml"
	wget -O ncbi/efetch.xml.$retstart $url
}
cat ncbi/efetch.xml.* | grep -e "<GBSeq_locus>" -e "<GBSeq_organism>" -e "<GBSeq_taxonomy>" | sed -e 's@</.*>@@g' -e 's/<.*[my]>//g' -e 's/^ *//g' | paste -s | sed 's/<GBSeq_locus>/\n/g' | sort -k1,1 > ncbi/18S_SSU_NCBI.tax.txt
cat ncbi/efetch.xml.* | gzip > ncbi/efetch.xml.txt.gz
rm ncbi/efetch.xml.[0-9]*

# we probably want to remove uncultured environemntal samples because they have poor taxonomy
grep -v -e "environmental samples" -e "uncultured" ncbi/18S_SSU_NCBI.tax.txt | awk '$0!=""{print ">"$1}' > ncbi/18S_SSU_NCBI.tax.txt.keep
seqtk seq -l 0 ncbi/18S_SSU_NCBI.fasta | grep -A1 -Fw -f ncbi/18S_SSU_NCBI.tax.txt.keep > ncbi/18S_SSU_NCBI.fasta.keep

# need to follow dada2 format conventions: https://benjjneb.github.io/dada2/training.html
#
# >ID Genus species
# >AB001448.1.1538 Pseudomonas savastanoi pv. phaseolicola
#
# >Level1;Level2;Level3;Level4;Level5;Level6;

# make BLAST database
module load blast-plus
srun makeblastdb -in ncbi/18S_SSU_NCBI.fasta -dbtype nucl
srun makeblastdb -in ncbi/18S_SSU_NCBI.fasta.keep -dbtype nucl

# BLASTN
for fasta in dada2_*/seqtab.*.fasta
do
	sbatch -c 8 --mem 16G --wrap "blastn -db ncbi/18S_SSU_NCBI.fasta.keep -query $fasta -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs' -evalue 1e-6 -num_threads 8 > $fasta.ncbi.blastn"
done

# filter BLASTN results
for bl in dada2_*/seqtab.*.fasta.ncbi.blastn
do
	# species: 99% identity and 99% query coverage
	awk '$3>=99 && $13>=99' $bl | sort -u -k1,1 > $bl.species

	# genus: 90% identity and 90% query coverage
	awk '$3>=90 && $13>=90' $bl | sort -u -k1,1 > $bl.genus
done

# make taxonomy tables
# species first
for bl in dada2_*/seqtab.*.species
do
	join -t $'\t' -1 2 -2 1 <(cut -f 1,2 $bl | sed 's/\.[0-9]$//g' | sort -k 2,2) ncbi/18S_SSU_NCBI.tax.txt | awk -F'\t' '{print $2, $3, $4}' OFS="\t" | sed 's/; /\t/g' > $bl.taxa
done
# genus 
for bl in dada2_*/seqtab.*.genus
do
	join -t $'\t' -1 2 -2 1 <(cut -f 1,2 $bl | sed 's/\.[0-9]$//g' | sort -k 2,2) ncbi/18S_SSU_NCBI.tax.txt | awk -F'\t' '{print $2, $3, $4}' OFS="\t" | sed 's/; /\t/g' > $bl.taxa
done
# keep genus only for those that don't have species annotation
for bl in dada2_*/seqtab.*.ncbi.blastn
do
	cut -f 1,2 $bl.species.taxa > $bl.taxatable
	cut -f 1 $bl.species.taxa | grep -v -Fw -f - $bl.genus.taxa | awk -F'\t' '{gsub(" .*", "", $2); print $1, $2}' OFS="\t" >> $bl.taxatable
	# add empty rows for ASVs without any annotation
	cut -f 1 $bl.taxatable | grep -v -Fw -f - <(grep ">" ${bl/.ncbi.blastn/}) | tr -d ">" | awk '{print $1,"NA"}' OFS="\t" >> $bl.taxatable
done


