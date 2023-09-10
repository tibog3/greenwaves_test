#!/bin/bash

MAT_SIZE_ARG=""
while getopts "s:" opt; do
  case "$opt" in
    s)
      MAT_SIZE_ARG="$OPTARG"
      ;;
    *)
      echo "Usage: $0 [-s MAT_SIZE]" >&2
      exit 1
      ;;
  esac
done

# Exécute le make avec MAT_SIZE
make clean all run platform=gvsoc MAT_SIZE="${MAT_SIZE_ARG:-64}"