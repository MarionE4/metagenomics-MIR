# MIR zone metagenomic analysis

Iron rich microbial mat at MIR, a low-temperature active deep-sea hydrothermal zone, was collected during 2 oceanographic campaign: HERMINE 2 (2022; DOI: [10.17600/18001851](https://doi.org/10.17600/18001851)) and BICOSE 3 (2023; [10.17600/18002399](https://doi.org/10.17600/18002399)). These samples are composed of sediments (HER2_PL2064-22\_**PBT2** and BIC3-PL2090-08-**PBT3**) or microbial mats (BIC3-PL2090-08-**PLUME4** and BIC3-PL2090-08-**PLUME6**). A co-assembly of these samples allowed the reconstruction of 625 Metagenome-assembled genomes (MAGs) before my master 2 internship. My goal during my internship was to start metagenomic analysis of microbial communities and *Zetaproteobacteria* (a class of neutrophilic iron-oxidizing bacteria) in MIR zone with DATARMOR server (Ifremer) and RStudio.

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

``` bash
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

``` bash
gtdbtk classify_wf --genome_dir 01_FASTA/REFORMATED --out_dir 03_TAXONOMY/GTDB-Tk -x fa --cpus 56 --skip_ani_screen 
```

After that, to define the quality filter, a step in RStudio was done to test a bunch of quality threshold (completion and redundancy).

```{r}
# Packages
library(dplyr)
library(tidyr)
library(formatR)
library(tibble)
```

We need to load the quality file and taxonomy files.

```{r}
# Load quality file, retrieve columns of interest and sort lines by MAGs name
checkm2 <- read.table("quality_report.tsv", header=TRUE, sep="\t")
checkm2_simplify <- select(checkm2, Name, Completeness, Contamination)
checkm2_simplify <- arrange(checkm2_simplify, Name)

# Load the taxonomy files and merge the 2 tables (archeae and bacteria) before arranging the lines according to the MAGs names.
taxo_b <- read.table(file = "gtdbtk.bac120.summary.tsv", header = TRUE, sep = "\t")
taxo_a <- read.table(file = "gtdbtk.ar53.summary.tsv", header = TRUE, sep = "\t")
taxo <- add_row(taxo_b, taxo_a)
taxo <- arrange(taxo, user_genome)
```

And then merged quality file, with MAGs name, completion and redundancy, and taxonomy file.

```{r}
Classification <- taxo$classification
MAGs_all <- add_column(checkm2_simplify, Classification)
```

Now we can test different quality filter and see how much *Zetaproteobacteria* we can work with.

```{r}
ninety_5 <- filter(MAGs_all, Completeness > 90, Contamination < 5) 
eighty_5 <- filter(checkm2_simplify, Completeness > 80, Contamination < 5)
seventy_5 <- filter(MAGs_all, Completeness > 70, Contamination < 5)
fivety_5 <- filter(checkm2_simplify, Completeness > 50, Contamination < 5)
sixty_5 <- filter(checkm2_simplify, Completeness > 60, Contamination < 5)

ninety_10 <- filter(checkm2_simplify, Completeness > 90, Contamination < 10)
eighty_10 <- filter(checkm2_simplify, Completeness > 80, Contamination < 10)
seventy_10 <- filter(checkm2_simplify, Completeness > 70, Contamination < 10)
fivety_10 <- filter(checkm2_simplify, Completeness > 50, Contamination < 10)
sixty_10 <- filter(checkm2_simplify, Completeness > 60, Contamination < 10)
```

The quality filter was set to \>70% completion and \<5% redundancy for "good quality" MAGs (MQ) and a file with high quality MAGs (HQ) was created too.

```{r}
MAGs_HQ <- filter(MAGs_all, Completeness > 90, Contamination < 5) #ninety_5
MAGs_MQ <- filter(MAGs_all, Completeness > 70, Contamination < 5) #seventy_5
```

And files were downloaded.

```{r}
write.table(MAGs_all, file='MAGs_all.tsv', quote = FALSE, sep = "\t", col.names = NA)
write.table(MAGs_MQ, file='MAGs_MQ.tsv', quote = FALSE, sep = "\t", col.names = NA)
write.table(MAGs_HQ, file='MAGs_HQ.tsv', quote = FALSE, sep = "\t", col.names = NA)
```

A symbolic link in DATARMOR was created in order to work with "good quality" MAGs (MQ). Two files with *Zetaproteobacteria* was created too.

``` bash
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

These argument was used:\
*-o* displays only the part of the text corresponding to the selection that follows; without *-o*, the entire line containing the selection would be displayed.\
*-P* activates Perl-compatible regular expressions mode, enabling you to use the *\\d+* (digits) expression, among others.\
*\^bin* searches for a name starting with bin and *\\d+* retrieves all digits following bin.

``` bash
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

``` bash
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

``` bash
# Workflow FeGenie
FeGenie.py -bin_dir ../04_MAGs/MAGs_SELECTED/ -bin_ext fa -out fegenie -T 24
```

### *grcJ* research (HMM profile)

*grcJ* is a potential periplasmic electron transporter involved in iron-oxidation (energetic metabolism) discovered in *Ghiorsea bivora* (*Zetaproteobacteria*) by Barco *et al.* in 2024. This gene was not in FeGenie so i decided to create a HMM profile for *grcJ* with HMMER v.3.3.2 (Eddy, 2023) in order to research it in "good quality" MAGs (MQ).
