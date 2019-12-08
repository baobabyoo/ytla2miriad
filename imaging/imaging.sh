#/bin/bash

path='../Miriad_converter/'

cp -r $path'ytla_usb.miriad' ./w51_usb
cp -r $path'ytla_lsb.miriad' ./w51_lsb

uvflag vis='w51_lsb' select='ant(2)' flagval=flag
uvflag vis='w51_usb' select='ant(2)' flagval=flag
#uvflag vis='w51_lsb' select='ant(5)' flagval=flag
#uvflag vis='w51_usb' select='ant(5)' flagval=flag


rm -rf ./*.temp
invert vis='w51_usb,w51_lsb' map=map.temp beam=beam.temp options=mosaic,mfs,systemp robust=2.0 cell=20 imsize=64 \
       line=channel,400,100,1,1

mossdi map=map.temp beam=beam.temp out=model.temp niters=1000 cutoff=1.0

restor map=map.temp beam=beam.temp model=model.temp out=clean.temp mode=clean
restor map=map.temp beam=beam.temp model=model.temp out=residual.temp mode=residual

rm -rf ./*.fits
fits in=clean.temp op=xyout out=clean.fits
fits in=residual.temp op=xyout out=residual.fits
fits in=model.temp op=xyout out=model.fits
fits in=beam.temp op=xyout out=beam.fits
