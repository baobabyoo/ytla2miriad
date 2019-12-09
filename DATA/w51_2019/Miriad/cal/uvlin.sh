#!/bin/bash

infiles='W51XX_lsb W51XX_usb W51YY_lsb W51YY_usb'
fields_ytla=$(seq 0 1 49)

for infile in $infiles
   do

   for field in $fields_ytla
      do

      inname=$infile.$field.miriad
   
      # flagging
      uvflag vis=$inname select='ant(2)' flagval=flag
   
      # uvlin / lsb: 92.88109 GHz / usb: 95.11891 GHz
      # Increment: 2.240000  Restfreq=Freq
      outnam=$inname'_ch0.'$field'.miriad'
      rm -rf $outname
      echo $outname
      uvlin vis=$inname chans=100,500 mode=chan0 out=$infile'_ch0.'$field'.miriad'

      done
   done
