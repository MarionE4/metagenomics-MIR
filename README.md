# MIR zone metagenomic analysis

Iron rich microbial mat at MIR, a low-temperature active deep-sea hydrothermal zone, was collected during 2 oceanographic campaign: HERMINE 2 (2022; DOI: [10.17600/18001851](https://doi.org/10.17600/18001851)) and BICOSE 3 (2023; [10.17600/18002399](https://doi.org/10.17600/18002399)). These samples are composed of sediments (HER2_PL2064-22\_**PBT2** and BIC3-PL2090-08-**PBT3**) or microbial mats (BIC3-PL2090-08-**PLUME4** and BIC3-PL2090-08-**PLUME6**). A co-assembly of these samples allowed the reconstruction of 625 Metagenome-assembled genomes (MAGs) before my master 2 internship. My goal during my internship was to start metagenomic analysis of microbial communities in MIR zone with DATARMOR server (Ifremer) and RSutdio.

## Rename and reformat MAGs

First, we had to rename the MAGs so that they had a name with no punctuation other than “\_”. Each MAG is named “bin” followed by a number. In the repertory where all MAGs are present, these lines of code have been launched.

``` bash
# Loop for renaming bins
for i in *.fa
do name=$(basename $i .fa | tr '[:punct:]' '_') # replace every punctuations by _
mv $i "$name".fa
done
```

The MAGs were then reformatted to avoid any future problems (notably with Anvi'o). In a new repertory these lines are run.

```{bash}
# Reformat fasta files : in fasta file, in contig name, replace old bin name by the new bin name
for i in FASTA/*.fa
do prefix=$(basename $i .fa)
anvi-script-reformat-fasta $i -o REFORMATED/"$prefix"_reformated.fa --prefix $prefix --simplify-names -r reformated.txt
done 
```

## MAG quality and taxonomy assignment

The quality of the MAGs was checked using the CheckM2 tool (Chklovski et al., 2023). This tool is a singularity on the DATARMOR server, so the following lines of code run the singularity after loading it.

``` bash
# Load singularity (code not shown)

# Checkm2 singularity container
CHECKM2_SIF=path/to/CHECKM2/SIF

# Create variable for work directory pathway
inputdir=path/to/FASTA/REFORMATED
outputdir=path/to/output/CheckM2

# Run Checkm2
singularity run -B /appli -B path/to/tool/checkm2/1.0.1/ -B $DATAWORK  ${CHECKM2_SIF} checkm2 predict --thread 30 --input $inputdir -x fa --output-directory $outputdir 
```

Taxonomic assignment was performed on all 625 MAGs using the GTDB-Tk v.2.4.0 tool (Chaumeil et al., 2020). The *-x* argument adds the file extension (here fasta), *--cpus* adds the number of threads and\
*--skip_ani_screen* allow to skip comparison genome ANI with Mash database, a time-consuming stage and unecessary for us.

```{bash}
gtdbtk classify_wf --genome_dir 01_FASTA/REFORMATED --out_dir 03_TAXONOMY/GTDB-Tk -x fa --cpus 56 --skip_ani_screen 
```

After that, a quality filter was set to \>70% completion and \<5% redundancy and a symbolic link was created in order to work with good quality MAGs. Two files with *Zetaproteobacteria* (a class of neutrophilic iron-oxidizing bacteria) was created too.

```{bash}
cd 04_MAGs

# Create a file with all Zetaproteobacteria MAGs
grep "Zetaproteobacteria" MAGs_all.tsv > ../zeta_all.tsv

# Create a file with medium quality Zetaproteobacteria MAGs
grep "Zetaproteobacteria" MAGs_MQ.tsv > ../zeta_MQ.tsv
# MAGs MQ c'est complétion >70% et redondance <5%

# Loop to make symbolic link of MAGs
cd 04_MAGs/MAGs_SELECTED
cut -f1 ../MAGs_MQ.tsv | while read i ; do ln -s ../../01_FASTA/REFORMATED/"$i".fa . ; done

cd 04_MAGs/ZETA_MQ
cut -f1 ../zeta_MQ.tsv | while read i ; do ln -s ../../01_FASTA/REFORMATED/"$i".fa . ; done
```

## MAGs annotations

Selected MAGs were annotated using Prokka v.1.14.6 (Seemann, 2014) and KEGG Kofams (Kanehisa and Goto, 2000) using Anvi'o v.8_dev (Eren et al., 2021). The FeGenie v.1.2 tool (Garber et al., 2020) was used to annotate genes involved in the iron cycle.

### Prokka v.1.14.6

These argument was used: \
*-o* displays only the part of the text corresponding to the selection that follows; without *-o*, the entire line containing the selection would be displayed.\
*-P* activates Perl-compatible regular expressions mode, enabling you to use the *\\d+* (digits) expression, among others.\
*\^bin* searches for a name starting with bin and *\\d+* retrieves all digits following bin.

```{bash}
# Workflow Prokka
for i in ../04_MAGs/MAGs_SELECTED/*.fa
do 
BASENAME=$(basename "$i") # to get only the file name and not its path with the extension
BIN=$(echo "$BASENAME" | grep -oP '^bin\d+') 
prokka $i --outdir prokka/$BIN --prefix $BIN
done
```

### KEGG Kofams with Anvi'o v.8_dev

We need to create the *contigs.db* (database) for each MAG before run loop annotation.

```{bash}
# Creation of contigs.db
for i in MAGs_SELECTED/*.fa
do 
BIN=$(basename $i | cut -f1 -d '_') # Get the bin with its number.
anvi-gen-contigs-database -f $i -o MAGs_SELECTED_db/$BIN.db
done 

# Workflow KEGG Kofams
for j in MAGs_SELECTED_db/*.db
do
anvi-run-kegg-kofams -c $j -T 20
done 
```

### FeGenie v.1.2

These argument was used:\
*-bin_ext* for extension argument (here fasta)*\
-T* for number of threads

```{bash}
# Workflow FeGenie
FeGenie.py -bin_dir ../04_MAGs/MAGs_SELECTED/ -bin_ext fa -out fegenie -T 24
```
