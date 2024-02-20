conda activate trimal-1.4.1

for fasta in dada2_*/seqtab.*Gachon.fasta
do
	sbatch -c 8 --mem 16G --wrap "
		~/sharedscratch/apps/mafft/bin/ginsi --thread 8 $fasta > $fasta.mafft
		trimal -in $fasta.mafft -out $fasta.mafft.trimal.phy -gt 0.75 -phylip
		~/sharedscratch/apps/iqtree-2.1.1-Linux/bin/iqtree2 -s $fasta.mafft.trimal.phy -nt 8 -bb 10000 -o Cafeteria_roenbergensis,Cafeteria_sp"
done


