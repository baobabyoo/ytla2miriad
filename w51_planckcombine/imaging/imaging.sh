#!/bin/bash

aca_cell=10
aca_imsize=256
acatp_niters=20000
acatp_cutoff=0.3
acatp_method='clean'


rm -rf acatp.map
rm -rf acatp.beam
invert "vis=./allvis/*" options=systemp,double,mosaic,mfs \
       map=acatp.map beam=acatp.beam cell=$aca_cell imsize=$aca_imsize robust=2.0


rm -rf acatp.model
if [ $acatp_method == 'clean' ]
then
    mossdi map=acatp.map beam=acatp.beam out=acatp.model gain=0.1 \
           niters=$acatp_niters cutoff=$acatp_cutoff
else
    mosmem map=acatp.map beam=acatp.beam out=acatp.model rmsfac=1.5
fi

rm -rf acatp.clean
rm -rf acatp.residual
restor map=acatp.map beam=acatp.beam model=acatp.model \
       mode=clean out=acatp.clean
restor map=acatp.map beam=acatp.beam model=acatp.model \
       mode=residual out=acatp.residual

fits in=acatp.map op=xyout out=acatp.map.fits
fits in=acatp.clean op=xyout out=acatp.clean.fits
fits in=acatp.residual op=xyout out=acatp.residual.fits
fits in=acatp.beam op=xyout out=acatp.beam.fits
fits in=acatp.model op=xyout out=acatp.model.fits

