# Commands used in the study

Generic commands for each analysis step. Replace the placeholders in `<angle brackets>` with your own files, directories and values.

| Placeholder | Meaning |
|---|---|
| `<threads>` | Number of CPU threads |
| `<genome.fasta>` | Genome assembly FASTA file |
| `<db.dmnd>` | DIAMOND database |

---

## 1. Genome annotation

**Funannotate v1.8.17**

```bash
funannotate clean -i <genome.fasta> -o <genome_clean.fasta>

funannotate sort -i <genome_clean.fasta> -o <genome_sort.fasta>

funannotate mask -i <genome_sort.fasta> -o <genome_masked.fasta>

funannotate predict --cpus <threads> -i <genome_masked.fasta> --species "<species_name>" --isolate <isolate_id> --transcript_evidence <transcript_evidence.fasta> --protein_evidence <protein_evidence.fasta> -o <output_dir>
```

---

## 2. Enterotoxin identification

**InterProScan v5.76-107.0**

```bash
./interproscan.sh -appl pfam -t p -i <proteins.fa> -f TSV --disable-precalc
```

---

## 3. GC content

**SeqKit v2.10.0**

```bash
seqkit fx2tab -n -g <archive.fasta> 
```

---

## 4. Statistical analysis

**R (base) and FSA package**

`tabgc` is a data frame with the columns `GC` (numeric) and `Group` (categorical).

```r
library(FSA)

# Shapiro-Wilk normality test
tapply(tabgc$GC, tabgc$Group, shapiro.test)

# Kruskal-Wallis test
kruskal.test(GC ~ Group, data = tabgc)

# Dunn's post hoc test (Holm correction)
dunnTest(GC ~ as.factor(Group), data = tabgc, method = "holm")
```

---

## 5. Redundancy removal

**CD-HIT v4.8.1**

```bash
cd-hit -i <input.fasta> -o <output.fasta> -c 1.00 -n 5 -M <memory> -d 0 -T <threads>
```

---

## 6. Phylogenetic inference

**IQ-TREE v3.1.3**

```bash
iqtree3 -s <alignment.fasta> -B 1000 -m TEST -nt <threads>
```

---

## 7. Enterotoxin conservation

**BUSCO v6.1.0**

```bash
busco -i <genomes_dir/> -l <lineage>_odb12 -o <output_name> -m genome -c <threads>
```

**OrthoFinder v3.1.5**

```bash
orthofinder -f <proteomes_dir/> -t <threads> -a <threads> -M msa -S diamond
```

**DIAMOND v2.1.13**

```bash
for file in *.fna; do
    genome=$(basename "$file" .fna)
    diamond blastx --threads <threads> -q "$file" -d <db.dmnd> \
        --outfmt 6 qseqid sseqid pident qcovhsp scovhsp \
        | awk -v genome="$genome" '{print genome "\t" $0}' >> <results.tsv>
done
```

---

## 8. RNA-seq quality control

**FastQC v0.12.1**

```bash
fastqc <reads.fastq> -t <threads>
```

**MultiQC v0.4**

```bash
multiqc .
```

**BBTools v39.15**

```bash
bbduk.sh in=<input.fastq> out=<output.fastq> minlen=30
```

---

## 9. Marker identification

**DIAMOND v2.1.13**

```bash
diamond blastx --query <query.fa> --db <db.dmnd> --outfmt 6 --max-target-seqs 5 --threads <threads> --out <results.tsv>
```

---

## 10. Expression quantification

**Salmon v0.14.1**

```bash
salmon index -t <transcripts.fa> -i <index_dir> -k 31

salmon quant -i <index_dir> -l A -1 <reads_1.fastq> -2 <reads_2.fastq> --validateMappings -o <output_dir>
```

---

## 11. Transposable element analysis

**RepeatModeler2 v2.0.9**

```bash
BuildDatabase -name <db_name> <genome.fasta>

RepeatModeler -database <db_name> -threads <threads> -LTRStruct
```

**RepeatMasker v4.2.4**

```bash
RepeatMasker -pa <threads> -lib <db_name>-families.fa -gff -dir <output_dir> <genome.fasta>
```

**BEDTools v2.31.1**

```bash
bedtools window -a <mRNA_loci.bed> -b <filtered_TEs.bed> -w 2500 > <output.tsv>
```

---

## 12. Detection of horizontal gene transfer

**AvP v1.0.10**

Taxon groups file (ingroup and outgroup, by NCBI taxid):

```bash
cat << 'EOF' > groups.yaml
Ingroup:
    <ingroup_taxid>: <ingroup_name>
EGP:
    <outgroup_taxid>: <outgroup_name>
EOF
```

Alien index calculation:

```bash
python <path/to/AvP>/aux_scripts/calculate_ai.py -i <blast_results.out> -x groups.yaml
```

Configuration file:

```bash
cat << 'EOF' > config.yaml
---
max_threads: <threads>

# DB path
blast_db_path: <path/to/blast_db>
fasta_path: <path/to/fasta_dir>
mode: blast
data_type: AA

## Algorithm options
# prepare
ai_cutoff: -200
ahs_cutoff: -1
outg_pct_cutoff: 0
selection: ai or ahs
percent_identity: 100
cutoffextend: 20
number_hits_noingroup: 50
trimal: false
min_num_hits: 4
percentage_similar_hits: 0.7

# detect, clasify, evaluate
fastml: true
node_support: 0
complex_per_ingroup: 20
complex_per_donor: 80
complex_per_node: 90

# Program specific options
mafft_options: '--anysymbol --auto'
trimal_options: '-automated1'

# IQ-Tree
iqmodel: '-mset WAG,LG,JTT -AICc -mrate E,I,G,R'
ufbootstrap: 1000
iq_threads: <threads>
EOF
```

Prepare step:

```bash
python3 <path/to/AvP>/avp prepare -a <ai_results.out> -o <output_dir> -f <query.fasta> -b <blast_results.out> -x groups.yaml -c config.yaml
```

Detect step:

```bash
python3 <path/to/AvP>/avp detect -i <output_dir>/mafftgroups/ -o <detect_results_dir> -g <output_dir>/groups.tsv -t <output_dir>/tmp/taxonomy_nexus.txt -c config.yaml
```
