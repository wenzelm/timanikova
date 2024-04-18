module load seqtk
# sed -i 's/Elachista [^ ]*/|&|/g' raw/Elachista.fas
# sed -i 's/Ahnfeltia [^ ]*/|&|/g' raw/Ahnfeltia.fas
# sed -i 's/Blidingia [^ ]*/|&|/g' raw/Blidingia.fas
# sed -i 's/Palmaria [^ ]*/|&|/g' raw/Palmaria\ palmata.fas
# sed -i 's/Ulva [^ ]*/|&|/g' raw/Ulva\ flexuosa.fas
# sed -i 's/Pylaiella [^ ]*/|&|/g' raw/Pylaiella\ littoralis.fas

cat raw/*.fas | seqtk seq -l 0 | awk '!/>/{gsub("-", "", $0)}{print}' > allseqs.fa
grep -A 1 "rbcL" allseqs.fa > allseqs.rbcL.fa

# keep longest sequence per species
awk -F'|' '/^>/{s=$2; l=$0}!/>/{print s, length($0), l}' OFS=":" allseqs.rbcL.fa | awk -F':' '$2>500' | sort -k1,1 -k2,2gr -t ':' | tee allseqs.rbcL.lengths.txt | sort -u -k 1,1 -t ':' > rbcL.best.species.txt

# we need to remove some accessions that represent geographically inappropriate isolates:
# Elachista antarctica (AJ439841) is Antarctic (we need Baffin)
# Palmaria decipiens (MF543838) is Antarctic (we need Baffin)
# Porphyra linearis is European (we need P. capensis or mumfordi instead)

# longest per genus
grep -v -e "MF543838" -e "AJ439841" -e "AF055398.1" rbcL.best.species.txt | awk -F':' '{gsub(" .*", "", $1); print $1, $2, $3}' OFS=":" | sort -k1,1 -k2,2gr -t ':' | sort -u -k 1,1 -t ':' | grep -v -e "aculeata" -e "Fucus" > rbcL.best.genus.txt

# we need additional barcodes for some genera/species:
# particularly Pylaiella, for which we need 3 additional sequences
# add to best genera file:
# -e "Fucus gardneri" -e "Fucus serratus" # Fucus is no longer needed
# -e "Porphyra umbilicalis" # No longer needed
# -e "Scytothamnus australis" # No longer needed
# -e "Desmarestia dudresnayi" # No longer needed

grep -Fw \
-e "Chorda asiatica" \
-e "Desmarestia chordalis" \
-e "Ectocarpus siliculosus" \
-e "Elachista scutulata" \
-e "Porphyra linearis" \
-e "Palmaria hecatensis" \
-e "Ulva flexuosa" \
rbcL.best.species.txt | cat - rbcL.best.genus.txt | grep -v -e "Chorda" -e "P. littoralis" | cut -f 3 -d ":" > rbcL.best.txt
grep "Pylaiella" allseqs.rbcL.lengths.txt | head -n 4 | tail -n 3 | cut -f 3 -d ":" >> rbcL.best.txt


# simplify names
# give Pylaiella sequences a unique name
grep --no-group-separator -A1 -Fw -f rbcL.best.txt allseqs.rbcL.fa | sed 's/>[^|]*|/>/g' | sed 's/|.*//g' | tr ' ' '_' \
	| awk '/Pylaiella/{$0=$0""NR}{print}' > rbcL.best.fa

# align
module load mafft
srun -c 8 --mem 16G ginsi --thread 8 --preservecase --adjustdirectionaccurately rbcL.best.fa > rbcL.best.mafft.fa

conda activate trimal-1.4.1
trimal -in rbcL.best.mafft.fa -out rbcL.best.mafft.trimal.phy -automated1 -phylip

# infer tree
srun -c 8 --mem 8G ~/sharedscratch/apps/iqtree-2.1.1-Linux/bin/iqtree2 -s rbcL.best.mafft.trimal.phy -bb 10000 -nt 8



