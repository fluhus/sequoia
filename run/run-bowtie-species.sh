#!/usr/bin/bash
#SBATCH -c 4
#SBATCH --mem 1g
#SBATCH -p free

#aSBATCH -A katrine_lab

set -e

if [ -z "$1" ]; then
  echo "no input dir"
  exit 1
fi
if [ -z "$2" ]; then
  echo "no ref dir"
  exit 1
fi
if [ -z "$3" ]; then
  echo "no output dir"
  exit 1
fi

indir="$1"
refdir="$2"
outdir="$3"

i=$SLURM_ARRAY_TASK_ID
nrefs=$(find "$refdir" -name \*.fasta | wc -l)
isample=$(((i-1) / nrefs + 1))
iref=$(((i-1) % nrefs + 1))

r=$(find "$refdir" -name \*.fasta | sort | sed 's/.fasta//' | line $iref)
rb=$(basename "$r")
s=$(ls -1S "$indir" | grep .1.fastq.zst | sed 's/.1.fastq.zst//' | line $isample)
outfile=$outdir/${rb}__$s.sam.zst

echo "Sample:    ($isample) $s"
echo "Reference: ($iref) $rb"
echo "Output:    $outfile"

bowtie2 -p 4 --no-head --no-unal -t --very-fast \
  -x $r \
  -1 $indir/$s.1.fastq.zst \
  -2 $indir/$s.2.fastq.zst \
  -S >(zstd > $outfile)

sleep 1
zstd --test $outfile

echo Done
