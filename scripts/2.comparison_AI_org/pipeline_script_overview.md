# 🛠️ Pipeline Script Overview with Functions and Key Code Snippets

This document lists all organized and cleaned pipeline scripts, describing their roles and showing representative code snippets that illustrate their key operations.

---

## 1️⃣ VCF Annotation and Filtering

- **`annotate_graph_vcf.sh`**  
  Annotates graph-based VCFs using GFA files and outputs compressed, sorted `.anno.vcf.gz` files.

  ```bash
  python3 annotate_vcf.py -vcf input.vcf -gfa input.gfa -o annotated.vcf
  ```

- **`pangenie_remove50.awk`**  
  AWK script to add `SVTYPE`, `END`, and `SVLEN` fields to VCF records and remove small variants (≤50bp).

  ```awk
  if (len > 50) {
      $8 = "SVTYPE=" type ";END=" ($2 + len) ";SVLEN=" len ";" $8
      print
  }
  ```

- **`filter_vcf_add_info.sh`**  
  Runs sorting, indexing, and filtering on VCFs, using the above AWK script, preparing them for Truvari evaluation.

  ```bash
  bcftools sort input.vcf.gz -Oz -o sorted.vcf.gz
  tabix -p vcf sorted.vcf.gz
  awk -f pangenie_remove50.awk sorted.vcf > filtered.vcf
  ```

---

## 2️⃣ Structural Variant Intersections and Statistics

- **`intersect_sv_callers.sh`**  
  Performs `bedtools intersect` between pangenome SVs and tool SVs, then runs downstream summary AWK scripts.

  ```bash
  bedtools intersect -a panSV.bed -b toolSV.bed -f 0.9 -r -wo > intersect.bed
  ```

- **`pan_stats.awk`**  
  Counts SV occurrences and cumulative lengths across types, chromosomes, and length bins.

  ```awk
  array["count", bed, chr, type, lenCate, extp]++
  array["length", bed, chr, type, lenCate, extp] += length
  ```

- **`pan_share_stats.awk`**  
  Compares shared SVs between pangenome and tool outputs, including exact match statistics for high-resolution overlap.

  ```awk
  if (uniq_pan[$7] == 1) { bed = "pan"; each_count_point(...) }
  if (uniq_tool[$14] == 1) { bed = "tool"; each_count_point(...) }
  ```

---

## 3️⃣ Truvari Benchmarking

- **`run_truvari_benchmark.sh`**  
  Runs Truvari benchmark for each sample and SV tool, stratified by repeat regions (LINE, SINE, LTR, etc.).

  ```bash
  truvari bench -b baseline.vcf.gz -c callset.vcf.gz -o output_dir --passonly --sizemin 50 --sizemax 1000000
  ```

- **`summarize_truvari_results.sh`**  
  Aggregates Truvari summary results (precision, recall, F1, TP, FP, FN) into a single combined table.

  ```bash
  awk -F '\t' 'NR==2 {print $2, $3, $4, $5, $6, $7}' summary.tsv
  ```

---

## 4️⃣ RTG SNV Comparison

- **`compare_snv_rtg.sh`**  
  Runs RTG `vcfeval` for pairwise SNV callset comparisons, stratified by genome-wide and repeat regions.

  ```bash
  rtg vcfeval -b baseline.vcf.gz -c callset.vcf.gz -t ref.sdf -o output_dir -T 4
  ```

- **`summarize_rtg_results.sh`**  
  Collects RTG summary results into a combined TSV file, allowing cross-method performance comparison.

  ```bash
  awk '/^Overall/ {getline; print $2, $3, $4, $5, $6, $7}' summary.txt
  ```

---

## 5️⃣ Assembly Quality and Summary

- **`summarize_assembly_metrics.sh`**  
  Extracts BUSCO completeness, QUAST N50, and gfastats total lengths from multiple assemblies.

  ```bash
  grep "C:" short_summary.json
  grep "Contig N50" report.txt
  grep "Total scaffold length" report.txt
  ```

- **`merge_assembly_summaries.sh`**  
  Merges all individual summaries into a single combined table, preparing data for plotting or reporting.

  ```bash
  paste busco.tsv quast.tsv gfastats.tsv | awk '{print $1, $2, $4, $6, $8, $10}'
  ```

---

## 6️⃣ Miscellaneous

- **`generate_sample_lists.sh`**  
  Generates `.sample` lists for Holstein, Jersey, Public, and combined cohorts.

  ```bash
  awk '$1 ~ /pattern/ {print $1}' seqfile.txt > sample.list
  ```

- **`filter_and_normalize_vcf.sh`**  
  Applies `bcftools` filters (e.g., missingness) and normalizes VCFs to biallelic format.

  ```bash
  bcftools view -i 'F_MISSING<0.2' input.vcf.gz | bcftools norm -m- -Oz -o output.vcf.gz
  ```

---
