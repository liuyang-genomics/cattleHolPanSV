# cattleHolPanSV

# P1.🧬 Structural Variant and SNV Calling Workflow (Annotated)

This document integrates all shell scripts provided into a unified and annotated markdown workflow based on the methods described in the manuscript. Each section provides detailed explanations of the corresponding analysis step, the rationale for method selection, and command-line execution used in the project.

---

## 1. Public Assembly Download

**Script**: `0.assembliesDownload.md`

Download publicly available Bos taurus genome assemblies from NCBI or Zenodo to ensure representative breed diversity for pangenome construction.

```bash
# Example command (GCA ID)
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/947/034/695/GCA_947034695.1_seqoccin.Bt.char.v1.0/GCA_947034695.1_seqoccin.Bt.char.v1.0_genomic.fna.gz
```

These assemblies serve as input references for Minigraph-Cactus and are included in the M13H or breed-specific pangenome graphs.

---

## 2. PacBio HiFi Assembly and Quality Assessment

**Script**: `1.asseblyAndSta.sh`

Generate high-quality haplotype-resolved de novo genome assemblies using PacBio HiFi data with or without Hi-C integration. Evaluate and report contiguity, completeness, and accuracy metrics.

Main tools and steps:

- **Hifiasm** for diploid/primary assembly
    
- **fastp** for QC and filtering (length ≥ 1 kb)
    
- **gfatools** to convert `.gfa` to `.fasta`
    
- **minimap2** for alignment to reference genome
    
- **BUSCO**, **QUAST**, **Inspector** for quality metrics
    

```bash
hifiasm -t $nthreads -r 3 1.clean_data/${sample}.hifi_reads.fastq.gz -o 2.assembly/${sample}
```

---

## 3. Pangenome Graph Construction (Minigraph-Cactus)

**Script**: `2.Pangenime.minigraph-cactus.sh`

Build variation-aware pangenome graphs by integrating multiple breed assemblies using Minigraph-Cactus (MC) pipeline. This allows for comprehensive SV discovery and graph-based genotyping.

```bash
cactus-pangenome \
  --binariesMode local \
  ${cactusVcf_name}.js ${cactusVcf_name}-seqfile.txt \
  --reference bosTau9 --giraffe clip filter --vcf \
  --permissiveContigFilter --haplo \
  --consCores 48 --consMemory 380G \
  --outDir $cactusVcf_name
```

The resulting GBZ, GFA, and VCF formats are used for variant extraction and downstream pangenome-based genotyping.

---

## 4. Long-read SV Detection

### (a) Assembly-based SV Calling

**Script**: `3.1.long-readsSVcallingAssembly-based.sh`

Align assembled contigs back to the reference genome and detect SVs based on contig-reference alignment.

Tools:

- `svim-asm`: supports haploid and diploid mode
    
- `PAV`: SV discovery using phased contigs
    

```bash
svim-asm haploid --min_sv_size 50 --max_sv_size 200000 ...
singularity run -B $work_dir pav_latest.sif -c $nthreads
```

### (b) Read Mapping-based SV Calling

**Script**: `3.2.long-readsSVcallingReadsMapping-based.sh`

Align HiFi reads to the reference using `pbmm2` or `minimap2`, followed by SV calling using multiple tools.

Tools:

- `sniffles`, `cuteSV`, `pbsv`, `svim`, `SVision`
    

```bash
pbmm2 align $ref_fa $sample.fastq.gz --preset CCS -j $nthreads --sort > sample.bam
sniffles -i sample.bam -v sample.vcf --minsupport 3 --minsvlen 50
```

Each SV callset is indexed and prepared for benchmarking and integration.

---

## 5. SNV Calling

### (a) PacBio HiFi Variant Calling

**Script**: `4.1.HiFi.SNVcallingClair3.sh`

Use `Clair3` to call SNPs and INDELs from HiFi BAMs with a pre-trained HiFi model.

```bash
run_clair3.sh \
  --bam_fn=sample.bam --ref_fn=$ref_fa --model_path=$model_path \
  --platform=hifi --output=output_dir --threads=$nthreads --gvcf
```

### (b) Illumina WGS Variant Calling

**Script**: `4.2.WGS.SNVcallingDeepV.sh`

Align short reads with `bwa-mem2`, then call variants using DeepVariant (v1.6.1). Post-process with GLnexus for joint genotyping.

```bash
singularity run docker://google/deepvariant:1.6.1 \
  --model_type=WGS --ref=$ref_fa --reads=sample.bam --output_vcf=sample.vcf.gz
```

### (c) Hybrid Variant Calling

**Script**: `4.3.Hybrid.SNVcallingDeepV.sh`

Merge HiFi and Illumina alignments, then call variants with DeepVariant using the hybrid model.

```bash
run_deepvariant \
  --model_type=HYBRID_PACBIO_ILLUMINA --reads=sample.hybri.bam
```

---

## 6. WGS-based SV Calling (Lumpy + CNVnator)

**Script**: `5.WGS.SVcalling.sh`

Call structural variants from Illumina short-read BAMs using `smoove` and `svtools`, then assess CNVs with `CNVnator`.

```bash
smoove call --outdir outdir --name sample --fasta $ref_fa --genotype sample.bam
cnvnator -root sample.root -tree sample.bam -chrom 1-29 X Y
```

The result is a merged and pruned SV callset with CNV annotations.

---

## 7. PanGenie-based Genotyping from Short Reads

**Script**: `6.WGS.PangenieSVgenotype.sh`

Genotype SVs against the Minigraph-Cactus-generated pangenome using PanGenie’s k-mer based genotyping model.

```bash
PanGenie -f index_dir -i <(zcat sample_1.fq.gz sample_2.fq.gz) -s sample -j $nthreads -t $nthreads -o output
```

Steps include:

- Filtering and preprocessing the VCF
    
- Genotyping single samples and merging callsets across samples
    

> ✅ _All paths and environment variables are derived from the provided scripts. Replace template values (e.g., `$sample`, `$work_path`) for actual use. Refer to your SLURM environment and module system when scheduling jobs._

---

# P2 Comparison between pangenomes and SV callers (for raw code)
*We have a AI organized version for good understood (Please see 2.comparison_AI_org floder)*

---

## 8. Pangenome Annotation and Statistics

**Scripts**: `1.0.pangenieAnno.sh`, `1.1.mc-stats.v2.sh`, `1.1.mc-stats.v3-jer-hapPri.sh`

Annotate graph-based SVs with topological features from GFA, and compute statistics on pangenome structure and SV types.

Tools:

- `annotate_vcf.py` for GFA-aware annotation
    
- `bcftools` for VCF filtering
    
- `gfastats`, `panacus` for graph stats
    

Example commands:

```bash
python3 annotate_vcf.py -vcf input.vcf -gfa input.gfa -o annotated.vcf
bcftools view -i 'F_MISSING<0.2' annotated.vcf.gz -Oz -o filtered.vcf.gz
```

Pangenome statistics include:

- node/edge/base counts
    
- per-haplotype coverage
    
- ordering across haplotypes
    

---

## 9. VCF Transformation, Merging, and Classification

**Scripts**: `1.2.pav.vcf.comp.sh`, `1.4.hol.unique.v2.sh`

- Extract BED representations of filtered SVs for different populations or pangenomes
    
- Annotate variants as "Holstein-specific", "Bovine-multi", or "Shared"
    
- Summarize SV counts and lengths by type, population, and allele frequency
    

```bash
zcat input.vcf.gz | awk -f class.id.awk > output.class
awk -f class.stat.awk output.class > output.stat
```

This enables comparison across pangenomes and identification of population-specific variants.

---

## 10. SV Overlap Across Callers and Platforms

**Scripts**: `1.5.intersect.ab_rb_sr.sh`, `1.5.intersect.ab_rb_sr-jerhap.sh`

Construct BED files from pangenome VCFs and various SV callers (Sniffles, cuteSV, pbsv, PanGenie, etc.), then compute reciprocal overlaps using `bedtools intersect`.

```bash
bedtools intersect -a panSV.bed -b callerSV.bed -f 0.9 -r -wo > overlap.bed
```

This provides pairwise SV sharing metrics, enabling benchmarking and caller agreement evaluation. Scripts also quantify:

- Total SV counts per tool
    
- Length-distribution binning
    
- Shared vs non-shared SV proportions
    

---

## 11. Truvari Benchmarking of SV Callsets

**Script**: `1.6.sv.com.sh`

Benchmark each SV callset against a pangenome-derived VCF using `truvari` across multiple categories:

- SV size (e.g., 50–200 bp, 200–1k, etc.)
    
- SV type (INS, DEL)
    
- Repeat category (e.g., LINE, SINE, LTR)
    

```bash
truvari bench -b pan.vcf.gz -c tool.vcf.gz -r 2000 -C 2000 --passonly --sizemin 50 --sizemax 1000000 -o output_dir
```

Each comparison exports precision, recall, F1 score, and genotype concordance to `.tsv` summaries.

---

## 12. Global Assembly Quality Checks

**Script**: `0.0.assemly.all.check.sh`

Aggregate BUSCO, QUAST, and Inspector statistics across all assemblies:

```bash
cat *.busco-hifiasm-bp/summary | grep 'C:'
cat *.inspector-hifiasm-bp/summary_statistics | grep 'QV'
```

These are compiled to evaluate coverage, mapping quality, and completeness of all HiFi assemblies used in pangenome construction.

---

## 13. SNV Accuracy Benchmarking via RTG vcfeval

**Script**: `1.7.snv.comp.sh`

Compare SNV callsets from different platforms (HiFi, Illumina, hybrid, PanGenie, RNA-seq) using `rtg vcfeval` across reference regions and repeat categories.

Workflow:

- Prepare per-sample VCFs for each method (e.g., `snv_hifi`, `snv_deepv`, `pangenie-sr`, etc.)
    
- Normalize and sort VCFs using `bcftools`
    
- Run `rtg vcfeval` for pairwise comparisons
    

```bash
rtg vcfeval -b baseline.vcf.gz -c call.vcf.gz -t ref.sdf -o output_dir -T 4
```

- Additional stratification by repeat class (e.g., LINE, SINE, LTR) using `-e $ref_rm/rm.LTR.bed`
    
- Final summary statistics exported to `all.summary.table` with fields:
    
    - `Precision`, `Recall`, `F1`, `GT Concordance`
        

This module provides SNV-level accuracy metrics across tools and sequencing platforms.

---

> ✅ _These additional sections (8–13) provide comprehensive benchmarking, comparative statistics, and cross-platform SV/SNV consistency evaluations aligned with the described manuscript pipeline._



# P3 Pangenome SV-based GWAS 

---
## 14. GWAS of SVs and SNPs Using EMMAX

**Script**: `2.5.emmax-cdcb.sh`

Perform genome-wide association studies (GWAS) using `EMMAX` for structural variants (SVs), SNPs, and combined variant datasets based on pangenome genotyping.

Workflow steps:

- Generate genotype files using `plink2` from pangenome VCF (SV, SNP, combined)
    
- Filter by `MAF > 0.05` and `missing < 10%`
    
- Compute kinship matrix using thinned SNPs (1 SNP per 10 kb)
    
- Run association tests for multiple traits with the same kinship matrix
    

Example commands:

```
plink2 --vcf input.vcf.gz --make-bed --maf 0.05 --mind 0.1 --out output_prefix
./emmax -v -d 10 -t genotypes.tped -p phenotypes.txt -k kinship_matrix -o results
```

- Separate models are run for SVs (`.sv.bfile`), SNPs (`.gsnp.bfile`), and full variants (`.var.bfile`)
    
- Results are post-processed to extract significant hits (`p < genome-wide threshold`)
    
- Position files are generated to convert marker IDs to coordinates for downstream annotation
    

This pipeline enables direct comparison of association signals between structural and sequence-level variants for complex traits in cattle.


---

### 15. Downsampled SNP vs SV GWAS Comparison

**Script**: `2.2.downsample_compareSNPandSV.sh`

Perform GWAS comparisons between SNP-only and SNP+SV datasets using EMMAX, after random downsampling to balance marker counts.

**Workflow steps:**

- **Generate `.bfile` using plink2 from SNP VCF.**
    
- **Randomly downsample SNPs to match SV counts:**
    

```bash
shuf -n 61249 input.bim | cut -f2 > snplist
plink2 --bfile input --extract snplist --make-bed --out downsampled
```

- **Recode to transposed PLINK format for EMMAX.**
    
- **Run GWAS on both full and downsampled sets:**
    

```bash
./emmax-intel64 -v -d 10 -t data.tped -p trait.txt -k kinship.kinf -o output_prefix
```

This allows for direct comparison of GWAS power between SNP-only, SV-only, and SNP+SV datasets under controlled marker count conditions.




# 📚 References and Citation

If you use this pipeline or its components in your research, please cite:
xxx

## 👥 Authors and Contact

**Pipeline Lead:**  
Dr. Yang Liu (杨柳)  
Postdoctoral Researcher, University of Maryland / USDA-ARS

**Collaborators:**

- Dr. Li Ma, University of Maryland
    
- Dr. George Liu, USDA-ARS
    

For questions or contributions, please open an issue on GitHub or contact:  
📧 **Email:** [yangqism@gmail.com]


## 💡 Contributing

We welcome pull requests and issue reports!  
Please follow the [Contributor Guidelines](CONTRIBUTING.md) if you wish to contribute new features, improvements, or fixes.
