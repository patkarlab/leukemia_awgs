#!/bin/bash
# Count reads linking two panel intervals via supplementary alignments.
# Usage: pair_rescue.sh <bam> <chrA> <startA> <endA> <chrB> <startB> <endB>
BAM=$1; CA=$2; SA_=$3; EA=$4; CB=$5; SB=$6; EB=$7

samtools view "$BAM" "${CA}:${SA_}-${EA}" \
  | awk -v c="$CB" -v lo="$SB" -v hi="$EB" '
      $5>=50 && index($0,"SA:Z:"c",")>0 {
        s=substr($0,index($0,"SA:Z:"));
        split(s,a,",");
        if (a[2]+0>=lo && a[2]+0<=hi) print $1, $4, a[2], $5
      }' | sort -u
