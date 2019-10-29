# ytla2miriad
Converting YTLA HDF5 data to Miriad format.
Modified from the idl2miriad task of the MIR IDL package (https://github.com/qi-molecules/sma-mir)
The source code ytla2miriad is under the Miriad_converter directory.
It depends on the external library files under the idl_sav directory.


History: 

2019.Oct.27 (Baobab) : Can run. But absolute flux is factor ~2 lower than expected, 
                       which may be related to convention about XX and YY (related to I).
                       Presently only permit outputing one source, one sideband (and 1 spectral window per sideband)
                       at a time. The velocity (vlsr) information is incorrect, which will be fixed when needed.
                       


Example of using the code:
IDL> .compile ytla2miriad.pro
IDL> filename = 'W51-ab.ytla7X.mrgh5'
IDL> load_ytla, ytla, filename, /verbose
IDL> ytla2miriad, ytla, dir='ytla_lsb.miriad', sideband='lsb'


Example of Miriad imaging command (under linux command lines)
# Created dirty images
invert vis=ytla_lsb.miriad,ytla_usb.miriad options=systemp,mosaic,mfs robust=2.0 map=ytla_map beam=ytla_beam cell=10 imsize=128
# mosaic clean
mossdi map=ytla_map beam=ytla_beam out=ytla_model niters=100 cutoff=1.0
# create cleaned image
restor map=ytla_map beam=ytla_beam model=ytla_model out=ytla_clean mode=clean
# create residual image
restor map=ytla_map beam=ytla_beam model=ytla_model out=ytla_clean mode=residual
# create fits output
fits in=ytla_clean op=xyout out=ytla_clean.fits
fits in=ytla_residual op=xyout out=ytla_residual.fits