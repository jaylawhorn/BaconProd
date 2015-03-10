#!/bin/bash

DIR=/eos/cms/store/user/ksung/production/phys14test/Phys14-PU20bx25_ZJetsToNuNu_HT-400to600_Tune4C

eos ls -l ${DIR} > files.tmp

for n in `cat files.tmp | awk -F'_' '{print $2}' | uniq -d`
do
  NDUPS=`grep Output_${n}_ files.tmp | grep -c root`
  if [ $NDUPS -eq 2 ]; then
    f1=$( grep Output_${n}_ files.tmp | head -n1 | awk '{print $5}' )
    f2=$( grep Output_${n}_ files.tmp | tail -n1 | awk '{print $5}' )
    if [ $f1 -eq $f2 ]; then
      filename=$( grep Output_${n}_ files.tmp | head -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}
    fi
    if [ $f1 -gt $f2 ]; then
      filename=$( grep Output_${n}_ files.tmp | tail -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}
    fi
    if [ $f1 -lt $f2 ]; then
      filename=$( grep Output_${n}_ files.tmp | head -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}
    fi
  fi
  if [ $NDUPS -eq 3 ]; then
    f1=$( grep Output_${n}_ files.tmp | head -n1 | awk '{print $5}' )
    f2=$( grep Output_${n}_ files.tmp | head -n2 | tail -n1 | awk '{print $5}' )
    f3=$( grep Output_${n}_ files.tmp | tail -n1 | awk '{print $5}' )
    if [ $f1 -eq $f2 ] && [ $f1 -eq $f3 ]; then
      filename=$( grep Output_${n}_ files.tmp | head -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}     
      filename=$( grep Output_${n}_ files.tmp | tail -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}    
    fi  
    if [ $f1 -gt $f2 ] && [ $f1 -gt $f3 ]; then
      filename=$( grep Output_${n}_ files.tmp | head -n2 | tail -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}     
      filename=$( grep Output_${n}_ files.tmp | tail -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}    
    fi  
    if [ $f2 -gt $f1 ] && [ $f2 -gt $f3 ]; then
      filename=$( grep Output_${n}_ files.tmp | head -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}     
      filename=$( grep Output_${n}_ files.tmp | tail -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}    
    fi  
    if [ $f3 -gt $f1 ] && [ $f3 -gt $f2 ]; then
      filename=$( grep Output_${n}_ files.tmp | head -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}     
      filename=$( grep Output_${n}_ files.tmp | head -n2 | tail -n1 | awk '{print $9}' )
      eos rm ${DIR}/${filename}    
    fi  
  fi  
done

rm -f files.tmp
