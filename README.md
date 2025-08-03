# Snakemake-Repeat-Annotation
Workflow to run repeat annotation.

Runs:
* `RepeatMasker`
* `ModDotPlot`
* `HumAS-SD` or `HumAS-HMMER`

## Getting Started
```bash
git clone https://github.com/logsdon-lab/Snakemake-Repeat-Annotation.git --recursive
```

## Usage
```bash
snakemake -np --configfile config/config.yaml
```

## Configuration
```yaml
# Input files
samples:
  # Sample name
  name:
    # Fasta (gzip or ungzipped)
    fa: "test/input/HG00096_cens.fa.gz"
    # Bedfile
    # bed: ...

# Workflow dirs
output_dir: "results"
log_dir: "logs"
benchmark_dir: "benchmarks"

repeatmasker:
  species: "human"
  engine: "rmblast"
  threads: 12
  mem: 20GB

moddotplot:
  mem: 20GB
  window: 5000
  ident_thr: 0.7

humas_annot:
  mode: "sd" # or "hmmer"
  threads: 12
  mem: 20GB
  hmm_profile: "data/AS-HORs-hmmer3.4-071024.hmm.gz"
```

To only run select workflows, just comment/omit the unwanted sections in the config:

```yaml
# Input files
samples:
  # Sample name
  name:
    # Fasta
    fa: "test/input/HG00096_cens.fa.gz"
    # Bedfile
    # bed: ...

# Workflow dirs
output_dir: "results"
log_dir: "logs"
benchmark_dir: "benchmarks"

repeatmasker:
  species: "human"
  engine: "rmblast"
  threads: 12
  mem: 20GB
```
