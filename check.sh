#!/bin/bash

if [ $# -ne 1 ]; then
  echo "Usage: ./test.sh <benchmark>"
  echo "Example: ./test.sh doitgen"
  exit
fi

BM=$1

cd benchmarks/

if [ ! -f $BM.c ]; then
  echo "File $BM.c not found. Aborting ..."
  exit
fi

TRIALS=3

export TARGET=$BM
export OMP_PROC_BIND=True
export OMP_NUM_THREADS=20

make base
for run in `seq $TRIALS`; do
  ./$BM-base.exe
done

make seq
for run in `seq $TRIALS`; do
  msg=`./$BM-seq.exe`
  echo $msg | sed -e 's/CHECK.*$//g'
  res=`echo "$msg" | sed -e 's/^.*CHECK //g'`
done
echo $res

make par
for run in `seq $TRIALS`; do
  msg=`./$BM-par.exe`
  echo $msg | sed -e 's/CHECK.*$//g'
  res=`echo "$msg" | sed -e 's/^.*CHECK //g'`
done
echo $res
