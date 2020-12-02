#!/bin/sh

set -e

. $(dirname $(dirname $0))/functions.sh

# --- Install software dependencies ------------------------------------------

if [ ! -x "$(command -v hmmsearch)" ]; then
    log Installing executable dependencies with aptitude
    apt update
    apt install -y hmmer
fi

log Installing Python dependencies with pip
pip install -U coverage tqdm pyhmmer

# --- Install data dependencies ----------------------------------------------

#mkdir -p ci/cache
#mkdir -p build/lib/gecco/data/hmms
#
#if [ "$CI_SERVER" = "yes" ]; then
#    QUIET="-q"
#else
#    QUIET=""
#fi
#
#for ini_file in gecco/hmmer/*.ini; do
#  url=$(grep "url" $ini_file | cut -d'=' -f2 | sed 's/ //g')
#  hmm=$(grep "id" $ini_file | cut -d'=' -f2 | sed 's/ //g')
#  version=$(grep "version" $ini_file | cut -d'=' -f2 | sed 's/ //g')
#
#  hmm_file=${ini_file%.ini}.hmm
#  cache_base="ci/cache/${hmm}.${version}.hmm"
#
#  if ! [ -e "$cache_base" ]; then
#    log Downloading $hmm v$version
#    wget "$url" $QUIET -O- | gunzip > "$cache_base"
#  else
#    log Using cached $hmm v$version
#  fi
#
#  mkdir -p $(dirname "build/lib/$hmm_file")
#  cp "$cache_base" "build/lib/$hmm_file"
#  hmmpress "build/lib/$hmm_file"
#done
