#!/bin/bash

JOBDIR=$1
ISDATA=$2
DATASET=$3
NUM=$4
STORAGE_ELEMENT=$5
STORAGE_PATH=$6
USER_REMOTE_DIR=$7

PSET=makingBacon_Data.py
CRABCFG=crab.cfg.data
if [ ${ISDATA} -eq 0 ]
  then
    PSET=makingBacon_MC.py
    CRABCFG=crab.cfg.mc
fi

mkdir $JOBDIR
cp $PSET $JOBDIR

sed -e "s|XX-DATASET-XX|${DATASET}|g" \
    -e "s|XX-NUM-XX|${NUM}|g" \
    -e "s|XX-STORAGE-ELEMENT-XX|${STORAGE_ELEMENT}|g" \
    -e "s|XX-STORAGE-PATH-XX|${STORAGE_PATH}|g" \
    -e "s|XX-USER-REMOTE-DIR-XX|${USER_REMOTE_DIR}|g" \
    < $CRABCFG \
    > $JOBDIR/crab.cfg
