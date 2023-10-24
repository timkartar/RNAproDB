#!/bin/bash

# Script to download files from RCSB http file download services.
# Use the -h switch to get help on usage.

if ! command -v curl &> /dev/null
then
    echo "'curl' could not be found. You need to install 'curl' for this script to work."
    exit 1
fi

PROGNAME=$0
BASE_URL="https://files.rcsb.org/download"

usage() {
  cat << EOF >&2
Usage: $PROGNAME -f <file> [-o <dir>] [-c] [-p]

 -f <file>: the input file containing a comma-separated list of PDB ids
 -o  <dir>: the output dir, default: current dir
 -m  <int>: maximum number of assemblies to try and download
EOF
  exit 1
}

download() {
    url="$BASE_URL/$1"
    out=$2/$1
    curl -s -f $url -o $out || echo "$1"
}

listfile=""
outdir="."
num_assemblies="4"
while getopts f:o:m: o
do
  case $o in
    (f) listfile=$OPTARG;;
    (o) outdir=$OPTARG;;
    (m) num_assemblies=$OPTARG;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"

if [ "$listfile" == "" ]
then
  echo "Parameter -f must be provided"
  exit 1
fi
contents=$(cat $listfile)

# see https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash#tab-top
IFS=',' read -ra tokens <<< "$contents"

for token in "${tokens[@]}"; do
    #for i in `seq $num_assemblies`; do
    #    download ${token}.pdb${i}.gz $outdir
    #done
    
    s=${outdir}/${token}
    #echo $s
    if ls $s* &>/dev/null
    then
        #echo "$s* exists."
        continue
    else
        download ${token}-assembly1.cif.gz $outdir
    fi
done









