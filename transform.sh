#!/usr/bin/env bash

if [ $# -ne 1 ]; then
  echo "Usage: ./transform.sh <benchmark>"
  echo "Example: ./transform.sh doitgen"
  exit
fi

BM=$1

./transform.py ./benchmarks/$BM.c.in ./benchmarks/$BM.c
