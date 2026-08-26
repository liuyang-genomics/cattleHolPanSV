#!/bin/bash
# generate_sample_lists.sh
# Generate sample lists for Holstein, Jersey, Public, and All categories

set -euo pipefail

seqfile="../minigraph-cactus/allbovinePan-2024-07-03/allbovinePan-2024-07-03-seqfile.txt"

awk '$1 ~ /sample/ {gsub(/\..*/, "", $1); print $1}' "$seqfile" | sort -u > hol.sample
awk '$1 ~ /Jersey/ {gsub(/\..*/, "", $1); print $1}' "$seqfile" | sort -u > jer.sample
awk '$1 ~ /hifiasm/ || $1 ~ /GCA/ {gsub(/\..*/, "", $1); print $1}' "$seqfile" | sort -u > pub.sample

cat hol.sample jer.sample pub.sample | sort -u > all.sample
