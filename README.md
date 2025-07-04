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

The **quality filter** was set to **\>70% completion** and **\<5% redundancy** for "good quality" MAGs (MQ) and a file with high quality MAGs (HQ) was created too.

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

First, the gene sequence in Barco *et al.* (2024) is recovered and blasted into NCBI to obtain other *grcJ* sequences. Fasta sequences with an e-value less to 1e-50 were all paste in a unique file text. These sequences was align with MAFFT ([https://mafft.cbrc.jp/alignment/server/](#0)). Then *grcJ* HMM profile was created.

``` bash
hmmbuild grcJ.hmm grcJ_alignement.fasta 
```

The research of *grcJ* can now be run on *.faa* created by Prokka annotation.

``` bash
cd 05_ANNOTATION
# Do the HMM search from .faa (prokka output)
for ID in prokka/*/*.faa
    do
    name=$(basename $ID .faa)
    hmmsearch ../09_HMM/grcJ.hmm $ID
done
```

In the output file, an e-value less than 1e-30 means that *grcJ* is present in the MAG.

## DiTing

DiTing v.0.95 (Xue *et al*., 2021) was used to compare the metabolic pathways potentially present between sample. This tool use the .*fastq* file (reads), it's a gene-centered analysis. First, we need to create a symbolic link and rename the *.fastq* files.

``` bash
# Path to working directory
cd path/to08_DiTing/01_clean_reads &&

# Create a symbolic link 
ln -s path/to/LUDIVINE/01_metagenomic_mat/02_anvio/01_QC/*F*.fastq.gz . &&

# Loop to rename fastq
for i in *.fastq.gz
  do name=$(echo $i | sed 's/F-QUALITY_PASSED_R//g')
    mv $i $name
done  
```

Then we can load DiTing from Github and run it.

``` bash
# We need to load DiTing from GitHub because the database don't work on DATARMOR server
cd 08_DITING/ &&
git clone https://github.com/xuechunxu/DiTing.git
cp -r /appli/bioinfo/diting/0.95/DiTing/kofam_database DiTing/
# -r for repository

python DiTing/diting.py -r 01_clean_reads -o 02_Output_DiTing -n 38
# diting.py is a script, we need to tell DATARMOR server to open and read this script in python language
```

## MIR phylogenomic tree

### Phylogenomic bacterial tree construction

The MAG alignment file “gtdbtk.bac120.**user_msa**.fasta.gz” created by GTDB-Tk's *align* function is used to build the tree. This file was cleaned using TrimAl v.1.4.1 (Capella-Gutiérrez *et al.*, 2009). Nucleotide positions with more than 50% gaps were removed using the *-gt 0.50* argument.

``` bash
gunzip gtdbtk.bac120.user_msa.fasta.gz
trimal -in 03_TAXONOMY/GTDB-Tk/align/gtdbtk.bac120.user_msa.fasta -out 07_TREE/01_general/bacteria_trim.fasta -gt 0.50 
```

IQ-TREE v.2.0.3 (Minh *et al.*, 2020) was used to reconstruct the phylogenomic tree using the maximum likelihood method, the WAG model and a bootstrap with 1000 repetitions. \
The *-nt AUTO* argument is used to build the tree with the most optimal number of threads (don't forget to add the maximum number of threads available with *-ntmax*).

``` bash
cd 07_TREE/01_general
iqtree -s bacteria_trim.fasta -nt AUTO -ntmax 40 -m WAG -bb 1000 
#-redo argument allow to rerun the command if it was already run a first time
```

The tree was then visualized using Figtree v.1.4.4. Various information, such as completion, redundancy, abundance in samples and the number of genes involved in the iron cycle for each MAG, was added to the tree.

### Mapping (read recruitment)

A mapping or read recruitment is necessary to know abundance of each MAGs in each sample (MAGs were co-assembled so we didn't know from which sample MAGs came from). On DATARMOR there is a workflow for mapping with Anvi'o v.8_dev, we just need to modify the "config.json"file and path to fit with our needs.

``` bash
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG="S08_mapping"_"$TIMESTAMP"".log"

anvi-run-workflow -w metagenomics \
                  -c path/to/config.json \
                  --additional-params \
                      --directory path/to/06_MAPPING \
                      --jobs 20 \
                      --latency-wait 420 \
                      --cluster-config path/to/cluster.yml \
		      --rerun-triggers mtime \
                      --cluster \
                          "qsub -N {rule} \
                                -m n \
                                -q {cluster.queue} \
                                -l walltime={cluster.walltime} \
                                -l mem={cluster.memory} \
                                -l ncpus={threads} \
                                -V" >> path/to/00_LOGS/$LOG 2>&1
```

Then, we need to create a collection. A collection is an Anvi'o object that we'll need for next step. Here we need to create a "MAG-collection.txt" file with the names of contigs in each MAG, with the name of the MAG to which each contig belongs.

``` bash
cd 06_MAPPING
for split_name in `sqlite3 03_CONTIGS/MAG-contigs.db 'select split from splits_basic_info'`
do 
	MAG=`echo $split_name | awk 'BEGIN{FS="\_000000"}{print $1}'` # Get contig name (e.g. bin163_re_assembled_contigs_000000) and MAG name (e.g. bin163_re_assembly)
	echo -e "$split_name\t$MAG" >> MAG-collection.txt # Write tab-separated names to a collection.txt (1 line per contig)
done 
```

Now let's create our collection.

``` bash
anvi-import-collection MAG-collection.txt -p 06_MERGED/MAG/PROFILE.db -c 03_CONTIGS/MAG-contigs.db -C MAG_collection
```

Finally we can do the mapping summarize. This step will generate a lot of interesting file including our "abundance.txt" file (in bins_across_samples repertory).

``` bash
anvi-summarize -p 06_MERGED/MAG/PROFILE.db -c 03_CONTIGS/MAG-contigs.db -C MAG_collection
```

#### Estimate metabolism

A command to estimate metabolic pathway completion (KEGG) of each MAG was run after the mapping.

``` bash
anvi-estimate-metabolism -p 06_MAPPING/06_MERGED/MAG/PROFILE.db -c 06_MAPPING/03_CONTIGS/MAG-contigs.db -C MAG_collection 
```

### Layers

Now that we have all the file needed, we can create a layers.txt file with informations we want to add to our phylogenomic tree in RStudio. First, load packages and files.

```{r}
library(tidyr)
library(dplyr)
library(data.table)
library(stringr)
library(formatR)
library(imputeTS)
```

```{r}
checkm2 <- read.table("quality_report.tsv", sep = "\t", header = T)
gtdb_bac <- read.table("gtdbtk.bac120.summary.tsv", sep = "\t", header = T)
fegenie <- read.table("FeGenie-heatmap-data.csv", sep = ",", header = T)
abundance <- read.table("abundance.txt", sep = "\t", header = T )
```

Then select columns.

```{r}
checkm2 <- checkm2 %>% select(Name, Completeness, Contamination)

gtdb_bac_1 <- gtdb_bac %>%
  separate(classification, into = c("domain", "phylum", "classe"), sep = ";", extra = "drop") %>%
  mutate(Phylum = sub("p__", "", phylum), 
         Classe = sub("c__", "", classe)) %>%
  select(user_genome, Phylum, Classe)

abundance$bins <- paste0(abundance$bins, "_reformated")
```

For genes number in iron energetic metabolism the heatmap generate by FeGenie was use. The heatmap must be set in the right way to retrieve the information.

```{r}
fegenie <- pivot_longer(fegenie, cols=2:259, names_to="Name", values_to="values")
fegenie <- pivot_wider(fegenie, names_from = X, values_from = values)

fegenie$Name <- sub(".fa", "", fegenie$Name) # Name must be the same in all table we need to delete .fa at the end here
fegenie <- select(fegenie, Name, iron_oxidation, iron_reduction)
```

Now, we can create our layers table by merging and adding table one by one with *all.x* argument. It allows to merge tables on table *x* lines.

```{r}
super_tableau <- merge(abundance, checkm2, by.x = "bins", by.y = "Name")
super_tableau <- merge(super_tableau, gtdb_bac_1, by.x = "bins", by.y = "user_genome", all.x = TRUE)
super_tableau <- merge(super_tableau, fegenie, by.x = "bins", by.y = "Name", all.x = TRUE)
super_tableau <- na_replace(super_tableau, 0)
super_tableau$iron_oxidation <- as.numeric(super_tableau$iron_oxidation)
super_tableau$iron_reduction <- as.numeric(super_tableau$iron_reduction)
```

Finally we can download our layers.txt and add layers to our phylogenomic tree with Anvi'o v.8.

```{r}
write.table(super_tableau, "couches.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
```

``` bash
# Create tree profile.db 
anvi-interactive -p phylogenomic-profile.db -t general.treefile --manual-mode 

# Import layers into this profile.db
anvi-import-misc-data couches.tsv -p phylogenomic-profile.db --target-data-table items

# Visualize and modify the tree
anvi-interactive -p phylogenomic-profile.db -t general.treefile --manual-mode
```

## *Zetaproteobacteria*

### Phylogenomic tree

A *Zetaproteobacteria* phylogenomic tree was constructed with reference genomes and references MAGs from GTDB. First, the identification and alignment with GTDB-Tk using *--taxa_filter* to include only *Zetaproteobacteria* and the outgroup (*Magnetococcia* class).

``` bash
cd 07_TREE/02_zeta
gtdbtk identify --cpus 33 -x fa --genome_dir MAGs_zeta --out_dir identify
gtdbtk align --cpus 33 --identify_dir identify --out_dir align --taxa_filter c__Zetaproteobacteria,c__Magnetococcia
```

Then trim the alignment file. Here we use "gtdbtk.bac120.**msa**.fasta" because we want our MAGs [AND]{.underline} *Zetaproteobacteria* reference genomes and MAGs.

``` bash
cd align/align &&
gunzip gtdbtk.bac120.msa.fasta.gz &&
cd path/to/02_zeta &&

trimal -in align/align/gtdbtk.bac120.msa.fasta -out zeta_trim.fasta -gt 0.50
```

And now construct the tree and root it.

``` bash
iqtree -s zeta_trim.fasta -nt AUTO -ntmax $NCPUS -m WAG -bb 1000 -redo &&

gtdbtk root --input_tree zeta_trim.fasta.treefile --outgroup_taxon c__Magnetococcia --output_tree zeta.treefile
```

Now we are happy but not to much because the name of MAGs and reference in the tree are ugly... No problem we'll fix that ! First, download number accession and taxonomy of reference genome and MAGs from GTDB.

```{r}
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
```

```{r}

```

```{r}
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
```

```{r}
gtdb <- read.csv("Annotations_zeta_gtdbtk.csv", sep = "\t")
```
