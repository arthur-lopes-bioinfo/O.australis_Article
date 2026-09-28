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

## 13. Molecular Dynamics Simulation Commands for OaHLEα16
**Gromacs v.2023**

### `minim.sh`

```bash
#!/bin/bash

#SBATCH --job-name=Oa16_min
#SBATCH --partition=main
#SBATCH --constraint=gromacs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=min-%j.out
#SBATCH --error=min-%j.err

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"

if [[ "$(hostname)" == bellatrix* ]]; then
    source /home/node03/anaconda3/etc/profile.d/conda.sh
    conda activate gmx2023
else
    source /usr/local/gromacs/bin/GMXRC
fi

echo "======================================"
echo "Etapa: minimização"
echo "Nó: $(hostname)"
echo "Início: $(date)"
echo "GROMACS: $(command -v gmx)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-N/A}"
echo "SLURM_JOB_GPUS: ${SLURM_JOB_GPUS:-N/A}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-N/A}"
gmx --version
if command -v nvidia-smi >/dev/null 2>&1; then
    echo "GPUs visíveis:"
    nvidia-smi -L || true
fi
echo "======================================"

gmx grompp \
    -f step4.0_minimization.mdp \
    -o step4.0_minimization.tpr \
    -c step3_input.gro \
    -r step3_input.gro \
    -p topol.top \
    -n index.ndx

gmx mdrun \
    -v \
    -ntomp "$SLURM_CPUS_PER_TASK" \
    -deffnm step4.0_minimization \
    -nb gpu

test -s step4.0_minimization.gro

echo "SUCESSO: minimização concluída"
echo "Fim: $(date)"

```

### `step4.0_minimization.mdp`

```ini
define                  = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator              = steep
emtol                   = 1000.0
nsteps                  = 5000
nstlist                 = 10
cutoff-scheme           = Verlet
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
constraints             = h-bonds
constraint_algorithm    = LINCS

```

### `equil.sh`

```bash
#!/bin/bash

#SBATCH --job-name=Oa16_eq
#SBATCH --partition=main
#SBATCH --constraint=gromacs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=eq-%j.out
#SBATCH --error=eq-%j.err

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"

if [[ "$(hostname)" == bellatrix* ]]; then
    source /home/node03/anaconda3/etc/profile.d/conda.sh
    conda activate gmx2023
else
    source /usr/local/gromacs/bin/GMXRC
fi

echo "======================================"
echo "Etapa: equilibração"
echo "Nó: $(hostname)"
echo "Início: $(date)"
echo "GROMACS: $(command -v gmx)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-N/A}"
echo "SLURM_JOB_GPUS: ${SLURM_JOB_GPUS:-N/A}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-N/A}"
gmx --version
if command -v nvidia-smi >/dev/null 2>&1; then
    echo "GPUs visíveis:"
    nvidia-smi -L || true
fi
echo "======================================"

if [[ ! -s step4.0_minimization.gro ]]; then
    echo "ERRO: step4.0_minimization.gro não encontrado."
    exit 1
fi

gmx grompp \
    -f step4.1_equilibration.mdp \
    -o step4.1_equilibration.tpr \
    -c step4.0_minimization.gro \
    -r step3_input.gro \
    -p topol.top \
    -n index.ndx

gmx mdrun \
    -v \
    -ntomp "$SLURM_CPUS_PER_TASK" \
    -deffnm step4.1_equilibration \
    -nb gpu

test -s step4.1_equilibration.gro

echo "SUCESSO: equilibração concluída"
echo "Fim: $(date)"

```

### `step4.1_equilibration.mdp`

```ini
define                  = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator              = md
dt                      = 0.001
nsteps                  = 125000
nstxout-compressed      = 5000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = v-rescale
tc_grps                 = SOLU SOLV
tau_t                   = 1.0 1.0
ref_t                   = 309.65 309.65
;
constraints             = h-bonds
constraint_algorithm    = LINCS
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SOLU SOLV
;
gen-vel                 = yes
gen-temp                = 309.65
gen-seed                = -1

```

### `prod.sh`

```bash
#!/bin/bash

#SBATCH --job-name=Oa19_prod
#SBATCH --partition=main
#SBATCH --constraint=gromacs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=14-00:00:00
#SBATCH --output=prod-%j.out
#SBATCH --error=prod-%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

if [[ "$(hostname)" == bellatrix* ]]; then
    source /home/node03/anaconda3/etc/profile.d/conda.sh
    conda activate gmx2023
else
    source /usr/local/gromacs/bin/GMXRC
fi

echo "======================================"
echo "Nó: $(hostname)"
echo "Início: $(date)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-não definido}"
echo "======================================"

if [[ ! -s step4.1_equilibration.gro ]]; then
    echo "ERRO: step4.1_equilibration.gro não encontrado."
    exit 1
fi

if [[ ! -s step5_production.tpr ]]; then

    gmx grompp \
        -f step5_production.mdp \
        -o step5_production.tpr \
        -c step4.1_equilibration.gro \
        -p topol.top \
        -n index.ndx

fi

if [[ -s step5_production.cpt ]]; then

    echo "Checkpoint encontrado."
    echo "Retomando produção..."

    gmx mdrun \
        -v \
        -ntomp "$SLURM_CPUS_PER_TASK" \
        -deffnm step5_production \
        -cpi step5_production.cpt \
        -nb gpu

else

    echo "Iniciando nova produção..."

    gmx mdrun \
        -v \
        -ntomp "$SLURM_CPUS_PER_TASK" \
        -deffnm step5_production \
        -nb gpu

fi

echo "Fim: $(date)"

```

### `step5_production.mdp`

```ini
integrator              = md
dt                      = 0.002
nsteps                  = 100000000
nstxout-compressed      = 50000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
rlist                   = 1.2
rcoulomb                = 1.2
coulombtype             = PME
;
tcoupl                  = v-rescale
tc_grps                 = SOLU SOLV
tau_t                   = 1.0 1.0
ref_t                   = 309.65 309.65
;
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 5.0
compressibility         = 4.5e-5
ref_p                   = 1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SOLU SOLV
;


```

### `submit_all.sh`

```bash
#!/bin/bash

J1=$(sbatch --parsable minim.sh)

J2=$(sbatch --parsable \
    --dependency=afterok:$J1 \
    equil.sh)

J3=$(sbatch --parsable \
    --dependency=afterok:$J2 \
    prod.sh)

echo "Minimização: $J1"
echo "Equilibração: $J2"
echo "Produção:     $J3"

squeue -u $USER

```

## Post-Molecular Dynamics Analysis: plot_postMD_OaHLEa16.py

```python
#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

DEFAULT_DATA_DIR = Path('postMD_OaHLEa16/data')
DEFAULT_OUTPUT_DIR = Path('postMD_OaHLEa16/figures_replotted_v3')
DEFAULT_RECEPTOR_RESIDUES = 545
DPI = 300

# Okabe-Ito colorblind-friendly palette
COLORS = {
    'blue': '#0072B2',
    'orange': '#E69F00',
    'green': '#009E73',
    'vermillion': '#D55E00',
    'purple': '#CC79A7',
    'skyblue': '#56B4E9',
    'yellow': '#F0E442',
    'black': '#000000',
}

SERIES_STYLE = {
    'Complex': {'color': COLORS['blue'], 'linestyle': '-', 'linewidth': 1.35},
    'Octopamine receptor': {'color': COLORS['vermillion'], 'linestyle': '-', 'linewidth': 1.35},
    'OaHLEa16': {'color': COLORS['green'], 'linestyle': '-', 'linewidth': 1.35},
}


def read_xvg(path: Path) -> np.ndarray:
    rows = []
    with path.open('r', encoding='utf-8', errors='ignore') as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(('#', '@', '&')):
                continue
            try:
                values = [float(x) for x in line.split()]
            except ValueError:
                continue
            if len(values) >= 2:
                rows.append(values)
    if not rows:
        raise RuntimeError(f'No numerical data found in: {path}')
    return np.asarray(rows, dtype=float)


def require_file(data_dir: Path, filename: str) -> Path:
    path = data_dir / filename
    if not path.is_file():
        raise FileNotFoundError(f'Required file not found: {path}')
    return path


def load_all_data(data_dir: Path) -> dict[str, np.ndarray]:
    names = {
        'hbonds': 'hbonds_receptor_toxin.xvg',
        'rmsd_complex': 'rmsd_complex.xvg',
        'rmsd_receptor': 'rmsd_receptor.xvg',
        'rmsd_toxin': 'rmsd_toxin.xvg',
        'rg_complex': 'rg_complex.xvg',
        'rg_receptor': 'rg_receptor.xvg',
        'rg_toxin': 'rg_toxin.xvg',
        'rmsf_complex': 'rmsf_complex_CA.xvg',
        'rmsf_receptor': 'rmsf_receptor_CA.xvg',
        'rmsf_toxin': 'rmsf_toxin_CA.xvg',
    }
    return {key: read_xvg(require_file(data_dir, filename)) for key, filename in names.items()}


def time_ns(data: np.ndarray) -> np.ndarray:
    return data[:, 0] / 1000.0


def angstrom(values: np.ndarray) -> np.ndarray:
    return values * 10.0


def configure_matplotlib() -> None:
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 15,
        'legend.fontsize': 11,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'axes.linewidth': 1.0,
        'figure.dpi': 150,
        'savefig.dpi': DPI,
    })


def save_figure(fig: plt.Figure, output_base: Path) -> None:
    output_base.parent.mkdir(parents=True, exist_ok=True)

    # Raster formats
    fig.savefig(
        output_base.with_suffix('.png'),
        dpi=DPI,
        bbox_inches='tight'
    )

    fig.savefig(
        output_base.with_suffix('.tiff'),
        dpi=DPI,
        bbox_inches='tight'
    )

    # Vector format
    fig.savefig(
        output_base.with_suffix('.svg'),
        format='svg',
        bbox_inches='tight'
    )

    plt.close(fig)

def legend_inside_upper_right(ax: plt.Axes) -> None:
    ax.legend(
        loc='upper right',
        bbox_to_anchor=(0.985, 0.985),
        frameon=True,
        facecolor='white',
        edgecolor='0.75',
        framealpha=0.94,
        fancybox=False,
        borderpad=0.6,
    )


def add_legend_headroom(ax: plt.Axes, y_arrays, fraction: float = 0.28, force_zero: bool = False) -> None:
    arrays = [np.asarray(y, dtype=float) for y in y_arrays]
    finite = [a[np.isfinite(a)] for a in arrays if a.size > 0]
    if not finite:
        return
    values = np.concatenate(finite)
    ymin = float(np.min(values))
    ymax = float(np.max(values))

    if force_zero:
        ymin_plot = 0.0
    else:
        data_range = ymax - ymin
        if data_range <= 0:
            data_range = max(abs(ymax), 1.0)
        ymin_plot = ymin - 0.04 * data_range

    data_range = ymax - ymin_plot
    if data_range <= 0:
        data_range = max(abs(ymax), 1.0)

    ax.set_ylim(ymin_plot, ymax + fraction * data_range)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.10, 1.055, label, transform=ax.transAxes,
            fontsize=15, fontweight='bold', va='top', ha='left')


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(direction='out', length=4, width=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def plot_hbonds(data, output_dir):
    d = data['hbonds']
    y = d[:, 1]
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    ax.plot(time_ns(d), y, color=COLORS['purple'], linestyle='-', linewidth=1.25,
            label='Receptor–OaHLEa16')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Number of intermolecular hydrogen bonds')
    ax.set_title('Intermolecular Hydrogen Bonds')
    add_legend_headroom(ax, [y], fraction=0.22, force_zero=True)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'intermolecular_hydrogen_bonds')


def plot_rmsd(data, output_dir):
    fig, ax = plt.subplots(figsize=(12.5, 7.5))
    plotted_y = []
    for label, d in [
        ('Complex', data['rmsd_complex']),
        ('Octopamine receptor', data['rmsd_receptor']),
        ('OaHLEa16', data['rmsd_toxin']),
    ]:
        y = angstrom(d[:, 1])
        plotted_y.append(y)
        ax.plot(time_ns(d), y, label=label, **SERIES_STYLE[label])
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('RMSD')
    add_legend_headroom(ax, plotted_y, fraction=0.34)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'RMSD_complex_octopamine_receptor_OaHLEa16')


def plot_rg(data, output_dir):
    fig, ax = plt.subplots(figsize=(12.5, 7.5))
    plotted_y = []
    for label, d in [
        ('Complex', data['rg_complex']),
        ('Octopamine receptor', data['rg_receptor']),
        ('OaHLEa16', data['rg_toxin']),
    ]:
        y = angstrom(d[:, 1])
        plotted_y.append(y)
        ax.plot(time_ns(d), y, label=label, **SERIES_STYLE[label])
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Radius of gyration (Å)')
    ax.set_title('Radius of Gyration')
    add_legend_headroom(ax, plotted_y, fraction=0.34)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'radius_of_gyration_complex_octopamine_receptor_OaHLEa16')


def plot_rmsf_receptor(data, output_dir):
    d = data['rmsf_receptor']
    x = np.arange(1, len(d) + 1)
    y = angstrom(d[:, 1])
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    ax.plot(x, y, label='Octopamine receptor', **SERIES_STYLE['Octopamine receptor'])
    ax.set_xlabel('Residue index')
    ax.set_ylabel('Cα RMSF (Å)')
    ax.set_title('Cα RMSF – Octopamine Receptor')
    add_legend_headroom(ax, [y], fraction=0.26, force_zero=True)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'RMSF_octopamine_receptor')


def plot_rmsf_toxin(data, output_dir):
    d = data['rmsf_toxin']
    x = np.arange(1, len(d) + 1)
    y = angstrom(d[:, 1])
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    ax.plot(x, y, label='OaHLEa16', **SERIES_STYLE['OaHLEa16'])
    ax.set_xlabel('Residue index')
    ax.set_ylabel('Cα RMSF (Å)')
    ax.set_title('Cα RMSF – OaHLEa16')
    add_legend_headroom(ax, [y], fraction=0.26, force_zero=True)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'RMSF_OaHLEa16')


def plot_rmsf_complex(data, output_dir, receptor_residues):
    d = data['rmsf_complex']
    x = np.arange(1, len(d) + 1)
    y = angstrom(d[:, 1])
    if receptor_residues >= len(d):
        raise ValueError(f'Boundary {receptor_residues} is outside the complex RMSF range ({len(d)} residues).')
    fig, ax = plt.subplots(figsize=(12.0, 7.0))
    ax.plot(x, y, label='Complex', **SERIES_STYLE['Complex'])
    ax.axvline(receptor_residues + 0.5, color=COLORS['yellow'], linestyle='--',
               linewidth=2.0, label='Receptor–toxin boundary')
    ax.set_xlabel('Residue index')
    ax.set_ylabel('Cα RMSF (Å)')
    ax.set_title('Cα RMSF – Receptor–OaHLEa16 Complex')
    add_legend_headroom(ax, [y], fraction=0.30, force_zero=True)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'RMSF_complex')


def plot_summary(data, output_dir, receptor_residues):
    fig, axes = plt.subplots(2, 2, figsize=(20.0, 14.5))
    ax_a, ax_b, ax_c, ax_d = axes.flatten()

    d = data['hbonds']
    y_hb = d[:, 1]
    ax_a.plot(time_ns(d), y_hb, color=COLORS['purple'], linestyle='-', linewidth=1.20,
              label='Receptor–OaHLEa16')
    ax_a.set_xlabel('Time (ns)')
    ax_a.set_ylabel('Number of intermolecular hydrogen bonds')
    ax_a.set_title('Intermolecular Hydrogen Bonds')
    add_legend_headroom(ax_a, [y_hb], fraction=0.24, force_zero=True)
    legend_inside_upper_right(ax_a)
    style_axes(ax_a)
    add_panel_label(ax_a, 'A')

    rmsd_y = []
    for label, d in [
        ('Complex', data['rmsd_complex']),
        ('Octopamine receptor', data['rmsd_receptor']),
        ('OaHLEa16', data['rmsd_toxin']),
    ]:
        y = angstrom(d[:, 1])
        rmsd_y.append(y)
        ax_b.plot(time_ns(d), y, label=label, **SERIES_STYLE[label])
    ax_b.set_xlabel('Time (ns)')
    ax_b.set_ylabel('RMSD (Å)')
    ax_b.set_title('RMSD')
    add_legend_headroom(ax_b, rmsd_y, fraction=0.38)
    legend_inside_upper_right(ax_b)
    style_axes(ax_b)
    add_panel_label(ax_b, 'B')

    rg_y = []
    for label, d in [
        ('Complex', data['rg_complex']),
        ('Octopamine receptor', data['rg_receptor']),
        ('OaHLEa16', data['rg_toxin']),
    ]:
        y = angstrom(d[:, 1])
        rg_y.append(y)
        ax_c.plot(time_ns(d), y, label=label, **SERIES_STYLE[label])
    ax_c.set_xlabel('Time (ns)')
    ax_c.set_ylabel('Radius of gyration (Å)')
    ax_c.set_title('Radius of Gyration')
    add_legend_headroom(ax_c, rg_y, fraction=0.38)
    legend_inside_upper_right(ax_c)
    style_axes(ax_c)
    add_panel_label(ax_c, 'C')

    d = data['rmsf_complex']
    x = np.arange(1, len(d) + 1)
    y_rmsf = angstrom(d[:, 1])
    ax_d.plot(x, y_rmsf, label='Complex', **SERIES_STYLE['Complex'])
    ax_d.axvline(receptor_residues + 0.5, color=COLORS['yellow'], linestyle='--',
                 linewidth=2.0, label='Receptor–toxin boundary')
    ax_d.set_xlabel('Residue index')
    ax_d.set_ylabel('Cα RMSF (Å)')
    ax_d.set_title('Cα RMSF – Receptor–OaHLEa16 Complex')
    add_legend_headroom(ax_d, [y_rmsf], fraction=0.34, force_zero=True)
    legend_inside_upper_right(ax_d)
    style_axes(ax_d)
    add_panel_label(ax_d, 'D')

    fig.suptitle('Post-MD Structural Analysis of the Octopamine Receptor–OaHLEa16 Complex',
                 fontsize=18, y=0.985)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.07, top=0.92,
                        wspace=0.24, hspace=0.30)
    save_figure(fig, output_dir / 'postMD_summary_4panel')


def main():
    parser = argparse.ArgumentParser(description='Replot OaHLEa16 post-MD XVG data.')
    parser.add_argument('--data-dir', type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument('--receptor-residues', type=int, default=DEFAULT_RECEPTOR_RESIDUES)
    args = parser.parse_args()

    configure_matplotlib()
    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 76)
    print('OaHLEa16 post-MD figure generator — v3')
    print('=' * 76)
    print(f'Data directory:   {data_dir}')
    print(f'Output directory: {output_dir}')
    print(f'Receptor length:  {args.receptor_residues} residues')
    print(f'Output:           PNG + TIFF, {DPI} dpi')
    print('Palette:          Okabe–Ito colorblind-friendly')
    print('Legend:           inside, upper-right, with reserved empty headroom')
    print('=' * 76)

    data = load_all_data(data_dir)
    plot_hbonds(data, output_dir)
    plot_rmsd(data, output_dir)
    plot_rg(data, output_dir)
    plot_rmsf_receptor(data, output_dir)
    plot_rmsf_toxin(data, output_dir)
    plot_rmsf_complex(data, output_dir, args.receptor_residues)
    plot_summary(data, output_dir, args.receptor_residues)

    print('\nFigures generated successfully:')
    for file in sorted(output_dir.iterdir()):
        if file.suffix.lower() in {'.png', '.tiff'}:
            print(' ', file.name)


if __name__ == '__main__':
    main()

```

## 14. Molecular Dynamics Simulation Commands for OaHLEα19
**Gromacs v.2023**

### `minim.sh`

```bash
#!/bin/bash

#SBATCH --job-name=Oa19_min
#SBATCH --partition=main
#SBATCH --constraint=gromacs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=min-%j.out
#SBATCH --error=min-%j.err

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"

if [[ "$(hostname)" == bellatrix* ]]; then
    source /home/node03/anaconda3/etc/profile.d/conda.sh
    conda activate gmx2023
else
    source /usr/local/gromacs/bin/GMXRC
fi

echo "======================================"
echo "Etapa: minimização"
echo "Nó: $(hostname)"
echo "Início: $(date)"
echo "GROMACS: $(command -v gmx)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-N/A}"
echo "SLURM_JOB_GPUS: ${SLURM_JOB_GPUS:-N/A}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-N/A}"
gmx --version
if command -v nvidia-smi >/dev/null 2>&1; then
    echo "GPUs visíveis:"
    nvidia-smi -L || true
fi
echo "======================================"

gmx grompp \
    -f step4.0_minimization.mdp \
    -o step4.0_minimization.tpr \
    -c step3_input.gro \
    -r step3_input.gro \
    -p topol.top \
    -n index.ndx


gmx mdrun \
    -v \
    -ntomp "$SLURM_CPUS_PER_TASK" \
    -deffnm step4.0_minimization \
    -nb gpu

test -s step4.0_minimization.gro

echo "SUCESSO: minimização concluída"
echo "Fim: $(date)"
```

### `step4.0_minimization.mdp`

```ini
define                  = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator              = steep
emtol                   = 1000.0
nsteps                  = 5000
nstlist                 = 10
cutoff-scheme           = Verlet
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
constraints             = h-bonds
constraint_algorithm    = LINCS
```

### `equil.sh`

```bash
#!/bin/bash

#SBATCH --job-name=Oa16_eq
#SBATCH --partition=main
#SBATCH --constraint=gromacs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=eq-%j.out
#SBATCH --error=eq-%j.err

set -eo pipefail

cd "$SLURM_SUBMIT_DIR"

if [[ "$(hostname)" == bellatrix* ]]; then
    source /home/node03/anaconda3/etc/profile.d/conda.sh
    conda activate gmx2023
else
    source /usr/local/gromacs/bin/GMXRC
fi

echo "======================================"
echo "Etapa: equilibração"
echo "Nó: $(hostname)"
echo "Início: $(date)"
echo "GROMACS: $(command -v gmx)"
echo "SLURM_JOB_ID: ${SLURM_JOB_ID:-N/A}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-N/A}"
echo "SLURM_JOB_GPUS: ${SLURM_JOB_GPUS:-N/A}"
echo "CUDA_VISIBLE_DEVICES: ${CUDA_VISIBLE_DEVICES:-N/A}"
gmx --version
if command -v nvidia-smi >/dev/null 2>&1; then
    echo "GPUs visíveis:"
    nvidia-smi -L || true
fi
echo "======================================"

if [[ ! -s step4.0_minimization.gro ]]; then
    echo "ERRO: step4.0_minimization.gro não encontrado."
    exit 1
fi

gmx grompp \
    -f step4.1_equilibration.mdp \
    -o step4.1_equilibration.tpr \
    -c step4.0_minimization.gro \
    -r step3_input.gro \
    -p topol.top \
    -n index.ndx

gmx mdrun \
    -v \
    -ntomp "$SLURM_CPUS_PER_TASK" \
    -deffnm step4.1_equilibration \
    -nb gpu

test -s step4.1_equilibration.gro

echo "SUCESSO: equilibração concluída"
echo "Fim: $(date)"
```

### `step4.1_equilibration.mdp`

```ini
define                  = -DPOSRES -DPOSRES_FC_BB=400.0 -DPOSRES_FC_SC=40.0
integrator              = md
dt                      = 0.001
nsteps                  = 125000
nstxout-compressed      = 5000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
coulombtype             = PME
rcoulomb                = 1.2
;
tcoupl                  = v-rescale
tc_grps                 = SOLU SOLV
tau_t                   = 1.0 1.0
ref_t                   = 309.65 309.65
;
constraints             = h-bonds
constraint_algorithm    = LINCS
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SOLU SOLV
;
gen-vel                 = yes
gen-temp                = 309.65
gen-seed                = -1
```

### `prod.sh`

```bash
#!/bin/bash

#SBATCH --job-name=Oa19_prod
#SBATCH --partition=main
#SBATCH --constraint=gromacs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=14-00:00:00
#SBATCH --output=prod-%j.out
#SBATCH --error=prod-%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

if [[ "$(hostname)" == bellatrix* ]]; then
    source /home/node03/anaconda3/etc/profile.d/conda.sh
    conda activate gmx2023
else
    source /usr/local/gromacs/bin/GMXRC
fi

echo "======================================"
echo "Nó: $(hostname)"
echo "Início: $(date)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-não definido}"
echo "======================================"

if [[ ! -s step4.1_equilibration.gro ]]; then
    echo "ERRO: step4.1_equilibration.gro não encontrado."
    exit 1
fi

if [[ ! -s step5_production.tpr ]]; then

    gmx grompp \
        -f step5_production.mdp \
        -o step5_production.tpr \
        -c step4.1_equilibration.gro \
        -p topol.top \
        -n index.ndx

fi

if [[ -s step5_production.cpt ]]; then

    echo "Checkpoint encontrado."
    echo "Retomando produção..."

    gmx mdrun \
        -v \
        -ntomp "$SLURM_CPUS_PER_TASK" \
        -deffnm step5_production \
        -cpi step5_production.cpt \
        -nb gpu

else

    echo "Iniciando nova produção..."

    gmx mdrun \
        -v \
        -ntomp "$SLURM_CPUS_PER_TASK" \
        -deffnm step5_production \
        -nb gpu

fi

echo "Fim: $(date)"
```

### `step5_production.mdp`

```ini
integrator              = md
dt                      = 0.002
nsteps                  = 100000000
nstxout-compressed      = 50000
nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstcalcenergy           = 100
nstenergy               = 1000
nstlog                  = 1000
;
cutoff-scheme           = Verlet
nstlist                 = 20
vdwtype                 = Cut-off
vdw-modifier            = Force-switch
rvdw_switch             = 1.0
rvdw                    = 1.2
rlist                   = 1.2
rcoulomb                = 1.2
coulombtype             = PME
;
tcoupl                  = v-rescale
tc_grps                 = SOLU SOLV
tau_t                   = 1.0 1.0
ref_t                   = 309.65 309.65
;
pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 5.0
compressibility         = 4.5e-5
ref_p                   = 1.0
;
constraints             = h-bonds
constraint_algorithm    = LINCS
continuation            = yes
;
nstcomm                 = 100
comm_mode               = linear
comm_grps               = SOLU SOLV
;
```

### `submit_all.sh`

```bash
#!/bin/bash

J1=$(sbatch --parsable minim.sh)

J2=$(sbatch --parsable \
    --dependency=afterok:$J1 \
    equil.sh)

J3=$(sbatch --parsable \
    --dependency=afterok:$J2 \
    prod.sh)

echo "Minimização: $J1"
echo "Equilibração: $J2"
echo "Produção:     $J3"

squeue -u $USER
```

### Post-Molecular Dynamics Analysis: plot_postMD_OaHLEα19.py

```python
#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

DEFAULT_DATA_DIR = Path('postMD_OaHLEa19/data')
DEFAULT_OUTPUT_DIR = Path('postMD_OaHLEa19/figures_replotted_v3')
DEFAULT_RECEPTOR_RESIDUES = 545
DPI = 300

# Okabe-Ito colorblind-friendly palette
COLORS = {
    'blue': '#0072B2',
    'orange': '#E69F00',
    'green': '#009E73',
    'vermillion': '#D55E00',
    'purple': '#CC79A7',
    'skyblue': '#56B4E9',
    'yellow': '#F0E442',
    'black': '#000000',
}

SERIES_STYLE = {
    'Complex': {'color': COLORS['blue'], 'linestyle': '-', 'linewidth': 1.35},
    'Octopamine receptor': {'color': COLORS['vermillion'], 'linestyle': '-', 'linewidth': 1.35},
    'OaHLEa19': {'color': COLORS['green'], 'linestyle': '-', 'linewidth': 1.35},
}


def read_xvg(path: Path) -> np.ndarray:
    rows = []
    with path.open('r', encoding='utf-8', errors='ignore') as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(('#', '@', '&')):
                continue
            try:
                values = [float(x) for x in line.split()]
            except ValueError:
                continue
            if len(values) >= 2:
                rows.append(values)
    if not rows:
        raise RuntimeError(f'No numerical data found in: {path}')
    return np.asarray(rows, dtype=float)


def require_file(data_dir: Path, filename: str) -> Path:
    path = data_dir / filename
    if not path.is_file():
        raise FileNotFoundError(f'Required file not found: {path}')
    return path


def load_all_data(data_dir: Path) -> dict[str, np.ndarray]:
    names = {
        'hbonds': 'hbonds_receptor_toxin.xvg',
        'rmsd_complex': 'rmsd_complex.xvg',
        'rmsd_receptor': 'rmsd_receptor.xvg',
        'rmsd_toxin': 'rmsd_toxin.xvg',
        'rg_complex': 'rg_complex.xvg',
        'rg_receptor': 'rg_receptor.xvg',
        'rg_toxin': 'rg_toxin.xvg',
        'rmsf_complex': 'rmsf_complex_CA.xvg',
        'rmsf_receptor': 'rmsf_receptor_CA.xvg',
        'rmsf_toxin': 'rmsf_toxin_CA.xvg',
    }
    return {key: read_xvg(require_file(data_dir, filename)) for key, filename in names.items()}


def time_ns(data: np.ndarray) -> np.ndarray:
    return data[:, 0] / 1000.0


def angstrom(values: np.ndarray) -> np.ndarray:
    return values * 10.0


def configure_matplotlib() -> None:
    plt.rcParams.update({
        'font.size': 12,
        'axes.labelsize': 14,
        'axes.titlesize': 15,
        'legend.fontsize': 11,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'axes.linewidth': 1.0,
        'figure.dpi': 150,
        'savefig.dpi': DPI,
    })


def save_figure(fig: plt.Figure, output_base: Path) -> None:
    output_base.parent.mkdir(parents=True, exist_ok=True)

    # Raster formats
    fig.savefig(
        output_base.with_suffix('.png'),
        dpi=DPI,
        bbox_inches='tight'
    )

    fig.savefig(
        output_base.with_suffix('.tiff'),
        dpi=DPI,
        bbox_inches='tight'
    )

    # Vector format
    fig.savefig(
        output_base.with_suffix('.svg'),
        format='svg',
        bbox_inches='tight'
    )

    plt.close(fig)


def legend_inside_upper_right(ax: plt.Axes) -> None:
    ax.legend(
        loc='upper right',
        bbox_to_anchor=(0.985, 0.985),
        frameon=True,
        facecolor='white',
        edgecolor='0.75',
        framealpha=0.94,
        fancybox=False,
        borderpad=0.6,
    )


def add_legend_headroom(ax: plt.Axes, y_arrays, fraction: float = 0.28, force_zero: bool = False) -> None:
    arrays = [np.asarray(y, dtype=float) for y in y_arrays]
    finite = [a[np.isfinite(a)] for a in arrays if a.size > 0]
    if not finite:
        return
    values = np.concatenate(finite)
    ymin = float(np.min(values))
    ymax = float(np.max(values))

    if force_zero:
        ymin_plot = 0.0
    else:
        data_range = ymax - ymin
        if data_range <= 0:
            data_range = max(abs(ymax), 1.0)
        ymin_plot = ymin - 0.04 * data_range

    data_range = ymax - ymin_plot
    if data_range <= 0:
        data_range = max(abs(ymax), 1.0)

    ax.set_ylim(ymin_plot, ymax + fraction * data_range)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.10, 1.055, label, transform=ax.transAxes,
            fontsize=15, fontweight='bold', va='top', ha='left')


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(direction='out', length=4, width=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def plot_hbonds(data, output_dir):
    d = data['hbonds']
    y = d[:, 1]
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    ax.plot(time_ns(d), y, color=COLORS['purple'], linestyle='-', linewidth=1.25,
            label='Receptor–OaHLEa19')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Number of intermolecular hydrogen bonds')
    ax.set_title('Intermolecular Hydrogen Bonds')
    add_legend_headroom(ax, [y], fraction=0.22, force_zero=True)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'intermolecular_hydrogen_bonds')


def plot_rmsd(data, output_dir):
    fig, ax = plt.subplots(figsize=(12.5, 7.5))
    plotted_y = []
    for label, d in [
        ('Complex', data['rmsd_complex']),
        ('Octopamine receptor', data['rmsd_receptor']),
        ('OaHLEa19', data['rmsd_toxin']),
    ]:
        y = angstrom(d[:, 1])
        plotted_y.append(y)
        ax.plot(time_ns(d), y, label=label, **SERIES_STYLE[label])
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('RMSD')
    add_legend_headroom(ax, plotted_y, fraction=0.34)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'RMSD_complex_octopamine_receptor_OaHLEa19')


def plot_rg(data, output_dir):
    fig, ax = plt.subplots(figsize=(12.5, 7.5))
    plotted_y = []
    for label, d in [
        ('Complex', data['rg_complex']),
        ('Octopamine receptor', data['rg_receptor']),
        ('OaHLEa19', data['rg_toxin']),
    ]:
        y = angstrom(d[:, 1])
        plotted_y.append(y)
        ax.plot(time_ns(d), y, label=label, **SERIES_STYLE[label])
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Radius of gyration (Å)')
    ax.set_title('Radius of Gyration')
    add_legend_headroom(ax, plotted_y, fraction=0.34)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'radius_of_gyration_complex_octopamine_receptor_OaHLEa19')


def plot_rmsf_receptor(data, output_dir):
    d = data['rmsf_receptor']
    x = np.arange(1, len(d) + 1)
    y = angstrom(d[:, 1])
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    ax.plot(x, y, label='Octopamine receptor', **SERIES_STYLE['Octopamine receptor'])
    ax.set_xlabel('Residue index')
    ax.set_ylabel('Cα RMSF (Å)')
    ax.set_title('Cα RMSF – Octopamine Receptor')
    add_legend_headroom(ax, [y], fraction=0.26, force_zero=True)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'RMSF_octopamine_receptor')


def plot_rmsf_toxin(data, output_dir):
    d = data['rmsf_toxin']
    x = np.arange(1, len(d) + 1)
    y = angstrom(d[:, 1])
    fig, ax = plt.subplots(figsize=(11.5, 7.0))
    ax.plot(x, y, label='OaHLEa19', **SERIES_STYLE['OaHLEa19'])
    ax.set_xlabel('Residue index')
    ax.set_ylabel('Cα RMSF (Å)')
    ax.set_title('Cα RMSF – OaHLEa19')
    add_legend_headroom(ax, [y], fraction=0.26, force_zero=True)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'RMSF_OaHLEa19')


def plot_rmsf_complex(data, output_dir, receptor_residues):
    d = data['rmsf_complex']
    x = np.arange(1, len(d) + 1)
    y = angstrom(d[:, 1])
    if receptor_residues >= len(d):
        raise ValueError(f'Boundary {receptor_residues} is outside the complex RMSF range ({len(d)} residues).')
    fig, ax = plt.subplots(figsize=(12.0, 7.0))
    ax.plot(x, y, label='Complex', **SERIES_STYLE['Complex'])
    ax.axvline(receptor_residues + 0.5, color=COLORS['yellow'], linestyle='--',
               linewidth=2.0, label='Receptor–toxin boundary')
    ax.set_xlabel('Residue index')
    ax.set_ylabel('Cα RMSF (Å)')
    ax.set_title('Cα RMSF – Receptor–OaHLEa19 Complex')
    add_legend_headroom(ax, [y], fraction=0.30, force_zero=True)
    legend_inside_upper_right(ax)
    style_axes(ax)
    fig.tight_layout()
    save_figure(fig, output_dir / 'RMSF_complex')


def plot_summary(data, output_dir, receptor_residues):
    fig, axes = plt.subplots(2, 2, figsize=(20.0, 14.5))
    ax_a, ax_b, ax_c, ax_d = axes.flatten()

    d = data['hbonds']
    y_hb = d[:, 1]
    ax_a.plot(time_ns(d), y_hb, color=COLORS['purple'], linestyle='-', linewidth=1.20,
              label='Receptor–OaHLEa19')
    ax_a.set_xlabel('Time (ns)')
    ax_a.set_ylabel('Number of intermolecular hydrogen bonds')
    ax_a.set_title('Intermolecular Hydrogen Bonds')
    add_legend_headroom(ax_a, [y_hb], fraction=0.24, force_zero=True)
    legend_inside_upper_right(ax_a)
    style_axes(ax_a)
    add_panel_label(ax_a, 'A')

    rmsd_y = []
    for label, d in [
        ('Complex', data['rmsd_complex']),
        ('Octopamine receptor', data['rmsd_receptor']),
        ('OaHLEa19', data['rmsd_toxin']),
    ]:
        y = angstrom(d[:, 1])
        rmsd_y.append(y)
        ax_b.plot(time_ns(d), y, label=label, **SERIES_STYLE[label])
    ax_b.set_xlabel('Time (ns)')
    ax_b.set_ylabel('RMSD (Å)')
    ax_b.set_title('RMSD')
    add_legend_headroom(ax_b, rmsd_y, fraction=0.38)
    legend_inside_upper_right(ax_b)
    style_axes(ax_b)
    add_panel_label(ax_b, 'B')

    rg_y = []
    for label, d in [
        ('Complex', data['rg_complex']),
        ('Octopamine receptor', data['rg_receptor']),
        ('OaHLEa19', data['rg_toxin']),
    ]:
        y = angstrom(d[:, 1])
        rg_y.append(y)
        ax_c.plot(time_ns(d), y, label=label, **SERIES_STYLE[label])
    ax_c.set_xlabel('Time (ns)')
    ax_c.set_ylabel('Radius of gyration (Å)')
    ax_c.set_title('Radius of Gyration')
    add_legend_headroom(ax_c, rg_y, fraction=0.38)
    legend_inside_upper_right(ax_c)
    style_axes(ax_c)
    add_panel_label(ax_c, 'C')

    d = data['rmsf_complex']
    x = np.arange(1, len(d) + 1)
    y_rmsf = angstrom(d[:, 1])
    ax_d.plot(x, y_rmsf, label='Complex', **SERIES_STYLE['Complex'])
    ax_d.axvline(receptor_residues + 0.5, color=COLORS['yellow'], linestyle='--',
                 linewidth=2.0, label='Receptor–toxin boundary')
    ax_d.set_xlabel('Residue index')
    ax_d.set_ylabel('Cα RMSF (Å)')
    ax_d.set_title('Cα RMSF – Receptor–OaHLEa19 Complex')
    add_legend_headroom(ax_d, [y_rmsf], fraction=0.34, force_zero=True)
    legend_inside_upper_right(ax_d)
    style_axes(ax_d)
    add_panel_label(ax_d, 'D')

    fig.suptitle('Post-MD Structural Analysis of the Octopamine Receptor–OaHLEa19 Complex',
                 fontsize=18, y=0.985)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.07, top=0.92,
                        wspace=0.24, hspace=0.30)
    save_figure(fig, output_dir / 'postMD_summary_4panel')


def main():
    parser = argparse.ArgumentParser(description='Replot OaHLEa19 post-MD XVG data.')
    parser.add_argument('--data-dir', type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument('--receptor-residues', type=int, default=DEFAULT_RECEPTOR_RESIDUES)
    args = parser.parse_args()

    configure_matplotlib()
    data_dir = args.data_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print('=' * 76)
    print('OaHLEa19 post-MD figure generator — v3')
    print('=' * 76)
    print(f'Data directory:   {data_dir}')
    print(f'Output directory: {output_dir}')
    print(f'Receptor length:  {args.receptor_residues} residues')
    print(f'Output:           PNG + TIFF, {DPI} dpi')
    print('Palette:          Okabe–Ito colorblind-friendly')
    print('Legend:           inside, upper-right, with reserved empty headroom')
    print('=' * 76)

    data = load_all_data(data_dir)
    plot_hbonds(data, output_dir)
    plot_rmsd(data, output_dir)
    plot_rg(data, output_dir)
    plot_rmsf_receptor(data, output_dir)
    plot_rmsf_toxin(data, output_dir)
    plot_rmsf_complex(data, output_dir, args.receptor_residues)
    plot_summary(data, output_dir, args.receptor_residues)

    print('\nFigures generated successfully:')
    for file in sorted(output_dir.iterdir()):
        if file.suffix.lower() in {'.png', '.tiff'}:
            print(' ', file.name)


if __name__ == '__main__':
    main()
```
