# cattleHolPanSV

Workflows for building a Holstein-centred cattle **pangenome** and calling
structural variants from it.

The pipeline covers de novo HiFi assembly, Minigraph-Cactus graph construction,
SV and SNV calling across long- and short-read platforms, benchmarking of the
resulting callsets against each other, and SV-based GWAS.

The imputation of these SVs into large genotyped cohorts lives in a companion
repository, [cattlePanSVimp](https://github.com/xyxss/cattlePanSVimp).

> **Scope.** This repository is code only. It was written for a Slurm HPC cluster
> and assumes Singularity for the containerised tools. Reproducing it end to end
> requires the sequence data listed under [Data availability](#data-availability)
> and substantial compute — the Cactus step alone is provisioned for 48 cores and
> 380 GB of memory.

---

## Repository layout

```
scripts/
  1.calling_pipeline/     Assembly, pangenome construction, SV and SNV calling
  2.comparison_AI_org/    Benchmarking, refactored into single-purpose scripts
  2.comparison_raw/       The original benchmarking scripts, as run
  3.SV_GWAS/              SV-based association analysis
config.sh.example         Site paths - copy to config.sh and edit
```

### About the two comparison directories

`2.comparison_raw/` holds the scripts as they were actually executed for the
manuscript. They are long, interleave several analyses per file, and contain
commented-out exploratory work.

`2.comparison_AI_org/` is the same logic reorganised into smaller
single-purpose scripts, easier to read and reuse. It was restructured with LLM
assistance and has **not** been re-run end to end to confirm it reproduces the
original outputs bit for bit.

**If you are reproducing published numbers, use `2.comparison_raw/`.** Use
`2.comparison_AI_org/` to understand the method or to adapt a step.

---

## Setup

### 1. Site configuration

All cluster-specific paths live in one git-ignored file:

```bash
cp config.sh.example config.sh
$EDITOR config.sh
source config.sh
```

Every script asserts the variables it needs at the top and exits with the
missing variable's name if you haven't sourced it.

### 2. Software

Most tools come from conda environments (`clair3`, `pangenie`, `svtools`, `cactus`)
or Singularity images. Point `CONDA_BASE` and `SOFTWARE_DIR` at your installs.
Versions used for the published run:

| Tool | Version | Stage |
|---|---|---|
| hifiasm | — | HiFi assembly |
| Minigraph-Cactus | cactus-bin 2.7.1 | Pangenome graph |
| DeepVariant | 1.6.1 (Docker/Singularity) | Short-read and hybrid SNVs |
| Clair3 | — | HiFi SNVs |
| PanGenie | — | Graph-based SV genotyping |
| SV callers | sniffles, cuteSV, pbsv, svim, svim-asm, SVision, PAV | SV discovery |
| smoove / svtools / CNVnator | — | Short-read SV and CNV |
| Truvari, RTG Tools | — | Benchmarking |
| gfastats, panacus | — | Graph statistics |
| SRA Toolkit | 3.1.0 | Public data retrieval |

### 3. Reference genome

**ARS-UCD2.0**, as `${REF_DIR}/ARS_UCD_v2.0.fa` (and `ARS_UCD_v2.0.chr.fa`, the
chromosome-only subset).

> **A naming caveat:** in the Cactus seqfile the reference haplotype is *labelled*
> `bosTau9`, and that label appears in graph paths and some VCF sample columns.
> It is only a label — the assembly is ARS-UCD2.0. By the usual UCSC convention
> `bosTau9` means ARS-UCD1.2, which is **not** what is used here.

---

## Pipeline

### Stage 1 — Assembly and variant calling (`1.calling_pipeline/`)

| Step | Script | What it does |
|---|---|---|
| 0 | `0.assembliesDownload.md` | Retrieve public *Bos taurus* assemblies for breed diversity |
| 1 | `1.asseblyAndSta.sh` | hifiasm assembly from HiFi (±Hi-C); QC with BUSCO, QUAST, Inspector |
| 2 | `2.Pangenime.minigraph-cactus.sh` | Build the pangenome graph; emits GBZ, GFA, VCF |
| 3.1 | `3.1.long-readsSVcallingAssembly-based.sh` | Assembly-based SVs (svim-asm, PAV) |
| 3.2 | `3.2.long-readsSVcallingReadsMapping-based.sh` | Read-mapping SVs (sniffles, cuteSV, pbsv, svim, SVision) |
| 4.1 | `4.1.HiFi.SNVcallingClair3.sh` | HiFi SNVs and indels via Clair3 |
| 4.2 | `4.2.WGS.SNVcallingDeepV.sh` | Illumina SNVs via DeepVariant, joint-called with GLnexus |
| 4.3 | `4.3.Hybrid.SNVcallingDeepV.sh` | Hybrid HiFi+Illumina SNVs via DeepVariant |
| 5 | `5.WGS.SVcalling.sh` | Short-read SVs (smoove/Lumpy, svtools) and CNVs (CNVnator) |
| 6 | `6.WGS.PangenieSVgenotype.sh` | Genotype graph SVs in short-read samples with PanGenie |

### Stage 2 — Benchmarking (`2.comparison_raw/`, `2.comparison_AI_org/`)

Compares pangenomes and SV callers against one another:

- **Graph annotation and statistics** — annotate graph VCFs with topology, then
  compute node/edge/base counts and per-haplotype coverage (`gfastats`, `panacus`)
- **SV classification** — label variants Holstein-specific, bovine-multi, or
  shared; summarise counts and lengths by type and allele frequency
- **Caller overlap** — reciprocal `bedtools intersect` across callers and platforms
- **Truvari benchmarking** — precision, recall, F1, and genotype concordance,
  stratified by SV size, SV type, and repeat class
- **SNV concordance** — `rtg vcfeval` across HiFi, Illumina, hybrid, PanGenie, and
  RNA-seq callsets, stratified by repeat class
- **Assembly QC aggregation** — collate BUSCO, QUAST, and Inspector metrics

### Stage 3 — SV-based GWAS (`3.SV_GWAS/`)

| Script | What it does |
|---|---|
| `2.1.emmax-cdcb.sh` | EMMAX GWAS on SV-only, SNP-only, and combined marker sets, sharing one kinship matrix |
| `2.2.downsample_compareSNPandSV.sh` | Downsamples SNPs to match SV marker count, so SNP and SV association power are compared at equal marker density |

---

## Reproducibility notes

**Unseeded sampling.** `3.SV_GWAS/2.2.downsample_compareSNPandSV.sh` downsamples
SNPs with `shuf -n 61249`, which reseeds from the OS on every run. Re-running it
selects a different SNP set and shifts the comparison. To make it reproducible,
replace `shuf` with a seeded shuffle:

```bash
SAMPLING_SEED=1
seeded_shuf () {
    shuf --random-source=<(openssl enc -aes-256-ctr -pass "pass:${SAMPLING_SEED}" \
        -nosalt </dev/zero 2>/dev/null) "$@"
}
```

(The same idiom is used in `cattlePanSVimp`. It is deterministic for a given GNU
coreutils version.)

**Known issue — mismatched assembly coordinates.** In `1.calling_pipeline/5.WGS.SVcalling.sh`
(and the same block in `4.2`, `4.3`, and `6.*`), `smoove call --exclude` is given
`Bos_taurus.ARS-UCD1.2.dna.toplevel.genomic_gaps.txt` while `--fasta` is
ARS-UCD2.0. The two assemblies have different coordinates, so the excluded
intervals do not correspond to actual ARS-UCD2.0 gaps. This has not been changed
here because doing so would alter published results, but it should be corrected
before the pipeline is reused.

**Version pinning.** Exact versions are recorded in [Setup](#2-software) where
known. Several tools were run from module or conda environments without a pinned
version; record yours when reproducing.

**Line endings.** All files are LF, enforced by `.gitattributes`. This matters:
CRLF endings make `bash` fail to parse `function name (){` and leave stray `\r`
characters in AWK fields.

You can check every script parses before submitting anything:

```bash
find scripts -name '*.sh' -exec bash -n {} \;
```

**Loop structure repairs.** Every script now parses (`bash -n` clean). Three in
`2.comparison_raw/` carried stray `done` lines, from being run interactively
section by section rather than top to bottom:

- `1.6.sv.com.sh` and `1.6.sv.com-jerhap.sh` each had one stray `done` at the
  shell level, plus one *inside* an `sbatch --wrap="..."` payload.
- `1.51.intersect.ab_rb_sr-pangenie.sh` parsed fine but had a stray `done`
  inside a `--wrap` payload, so the submitted job would have failed at runtime
  even though the wrapper looked healthy. Its `if [ -f 1.truvari/  ]` guard was
  also testing an incomplete path, so the `rm -rf` meant to clear a stale
  Truvari output directory never ran; the path now matches the `rm` on the
  following line, as it does in the sibling scripts.

Structurally identical blocks elsewhere in the same files were used to confirm
which `done` was the stray in each case.

---

## Data availability

This repository holds **code only** — no sequence data, genotypes, or phenotypes,
and none should be added (see `CONTRIBUTING.md`).

- Public assemblies are retrieved by accession in `0.assembliesDownload.md`.
- HiFi, Hi-C, and Illumina sequence data are deposited under the accessions
  listed in the manuscript.
- CDCB-derived genotypes and phenotypes used in stage 3 are **restricted-access**
  and cannot be redistributed. Access requires an agreement with the Council on
  Dairy Cattle Breeding.

---

## Citation

See `CITATION.cff`. If you use this pipeline, please cite the accompanying
manuscript once it is published.

## License

MIT — see `LICENSE`.

## Contact

Maintainer contact is listed in `CITATION.cff`. Questions and issue reports are
welcome; see `CONTRIBUTING.md`.
