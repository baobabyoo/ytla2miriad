#!/bin/bash

##### Parameters #########################################

# README -------------------------------------------------
#
#
# Latest update: 2019 Nov. 10 by Baobab Liu, adapted
#                from the script to combine ALMA, ACA, and TP data.
#
# Compatible with outputs of CASA 5.4 and Miriad-carma 4.3.8
# For combining spectral line cubes.
#
#
# It would be very much appreciated if you can cite
# https://ui.adsabs.harvard.edu/abs/2013ApJ...770...44L/abstract
# when using this script
#
#
# Flow ---------------------------------------------------
#
# 1. Convert input data from FITS to Miriad format 
#
# 2. Correct headers if necessary
#
# 3. Generate TP visibilities at ACA pointings
#
# 4. Jointly image ACA and TP visibilities
#
# 5. Re-generate ACA+TP visibilities at 12m pointings
#
# --------------------------------------------------------


# flow control -------------------------------------------
if_fitstomiriad='yes'
if_setheaders='yes'
if_tpdeconvolve='yes'
if_tp2vis='yes'
if_tprewt='yes'
if_duplicateACATP='yes'
if_acaim='yes'
if_finalim='no'
if_cleanup='yes'
if_reimmerge='no'
# --------------------------------------------------------


# global information -------------------------------------
linerestfreq='92.88109' # in GHz unit

obsfreq='92.88109' # in GHz unit

# set the starting channel and number of channels
ch_start="1"
numch="1"

# set the dimension for the YTLA images
aca_cell=10
aca_imsize=256

# parameters for deconvolving the TP image
cutoff_tp=1.0
niters_tp=100

# --------------------------------------------------------

# YTLA data information
path_ytla='../../DATA/w51_2019/Miriad/ch0/'
name_ytla='W51XX_lsb_ch0'
# fields_ytla='0 1 2 3 4 5 6 9 10 11 12 13 14 24 26 28 29' #$(seq 0 1 32)
fields_ytla=$(seq 0 1 49)
pbfwhm_ytla='635.05'

# TP data information
path_tp='../../../planck_img/'
name_tp='PlanckHFI_Type3_tan_100GHz'
pbfwhm_tp='579.06'
crval3_tp=167.4129
cdelt3_tp=-0.3175914
crpix3_tp=1
tsys_tp=6500.0
tp_unit_convert="yes"
tp_brightness_unit="Jy/beam"
tp_conv_f=2072.0  # place holder
tp_add_f=0.0
uvmax=0.4

# - - Defining global variables - - - - - - - - - - - - -#
ch='channel,'$numch',1,1,1'
chout='channel,'$numch',1,1,1'

# ACA+TP imaging parameter
acatp_niters=10000
acatp_cutoff=0.3
acatp_method='clean'
outvis_end='Xl'


##########################################################



##### Converting FITS data to Miriad format ##############

if [ $if_fitstomiriad == 'yes' ] 
then


   # Planck
   echo '########## Importing Planck data ##########'
   outname=$name_tp'.image.miriad'
   rm -rf $outname
   fits in=$path_tp$name_tp'.fits' \
	op=xyin \
	out=$outname

   # YTLA
   echo '########## Importing YTLA data ##########'
   for field_id in $fields_ytla
     do

        outname=$name_ytla'.'$field_id'.uv.miriad'
        rm -rf $outname
        cp -r $path_ytla$name_ytla'.'$field_id'.miriad' './'$outname
     done


fi

##########################################################



##### Reset headers to allow Miriad processing ###########

if [ $if_setheaders == 'yes' ]
then

   # YTLA data
   for field_id in $fields_ytla
     do
        pb="gaus("$pbfwhm_ytla")"
        puthd in=$name_ytla'.'$field_id'.uv.miriad'/telescop \
              value='single' type=a
        puthd in=$name_ytla'.'$field_id'.uv.miriad'/pbtype \
              value=$pb type=a
        puthd in=$name_ytla'.'$field_id'.uv.miriad'/systemp \
	      value=100.0 type=r
        puthd in=$name_ytla'.'$field_id'.uv.miriad'/jyperk \
              value=100.0 type=r
#        puthd in=$name_ytla'.'$field_id'.uv.miriad'/ctype3 value='VELO-LSR' type=ascii
#        puthd in=$name_ytla'.'$field_id'.uv.miriad'/cunit3 value='km/s    ' type=ascii
#        puthd in=$name_ytla'.'$field_id'.uv.miriad'/crval3 value=$crval3_tp type=double
#        puthd in=$name_ytla'.'$field_id'.uv.miriad'/cdelt3 value=$cdelt3_tp type=double
#        puthd in=$name_ytla'.'$field_id'.uv.miriad'/crpix3 value=$crpix3_tp type=double
	 
     done


   # TP
   puthd in=$name_tp'.image.miriad'/bmaj value=$pbfwhm_tp,arcsec type=double
   puthd in=$name_tp'.image.miriad'/bmin value=$pbfwhm_tp,arcsec type=double
   puthd in=$name_tp'.image.miriad'/bpa  value=0,degree type=double
#   puthd in=$name_tp'.image.miriad'/ctype3 value='VELO-LSR' type=ascii
#   puthd in=$name_tp'.image.miriad'/cunit3 value='km/s    ' type=ascii
#   puthd in=$name_tp'.image.miriad'/crval3 value=$crval3_tp type=double
#   puthd in=$name_tp'.image.miriad'/cdelt3 value=$cdelt3_tp type=double
#   puthd in=$name_tp'.image.miriad'/crpix3 value=$crpix3_tp type=double

   puthd in=$name_tp'.image.miriad'/ctype3 value='FREQ' type=ascii
   puthd in=$name_tp'.image.miriad'/cdelt3 value=1 type=int
   puthd in=$name_tp'.image.miriad'/crpix3 value=1 type=double
   puthd in=$name_tp'.image.miriad'/crval3 value=$obsfreq type=double
   puthd in=$name_tp'.image.miriad'/naxis3 value=1 type=int


   rm -rf single_input.miriad
   if [ $tp_unit_convert == 'yes' ]
   then
     maths exp="(($name_tp.image.miriad)*$tp_conv_f)" \
	   out=single_input.miriad options=unmask
   else
     cp -r $name_tp'.image.miriad' single_input.miriad
   fi

   puthd in=single_input.miriad/bmaj value=$pbfwhm_tp,arcsec type=double
   puthd in=single_input.miriad/bmin value=$pbfwhm_tp,arcsec type=double
   puthd in=single_input.miriad/bpa  value=0,degree type=double
   puthd in=single_input.miriad/bunit value='Jy/beam' type=ascii

fi

##########################################################



##### Deconvolve TP map ##################################

if [ $if_tpdeconvolve == 'yes' ]
then

   # Generate the TP Gaussian Beam
   rm -rf tp_beam
   imgen out=tp_beam imsize=$aca_imsize cell=$aca_cell \
         object=gaussian \
         spar=1,0,0,$pbfwhm_tp,$pbfwhm_tp,0


   for field_id in $fields_ytla
     do
        # Creat template ACA maps for regriding TP maps
        rm -rf single_input.aca_$field_id.temp.miriad
	rm -rf temp.beam
	echo $ch
        invert vis=$name_ytla'.'$field_id'.uv.miriad'   \
               imsize=$aca_imsize cell=$aca_cell options=double \
               map=single_input.aca_$field_id.temp.miriad beam=temp.beam line=$ch

	# Regrid TP maps
        rm -rf single_input.aca_$field_id.regrid.miriad       
        regrid in=single_input.miriad tin=single_input.aca_$field_id.temp.miriad \
               out=single_input.aca_$field_id.regrid.miriad \
	       project=sin

        # Deconvolve the TP Map
	rm -rf single_input.aca_$field_id.deconv.miriad
        clean map=single_input.aca_$field_id.regrid.miriad beam=tp_beam \
	      out=single_input.aca_$field_id.deconv.miriad \
	      niters=$niters_tp cutoff=$cutoff_tp gain=0.05

	# Restore the deconvolved TP map for a sanity check
	rm -rf single_input.aca_$field_id.restor.miriad
        restor map=single_input.aca_$field_id.regrid.miriad beam=tp_beam \
	       model=single_input.aca_$field_id.deconv.miriad \
               mode=clean out=single_input.aca_$field_id.restor.miriad

        rm -rf single_input.aca_$field_id.residual.miriad
        restor map=single_input.aca_$field_id.regrid.miriad beam=tp_beam \
               model=single_input.aca_$field_id.deconv.miriad \
               mode=residual out=single_input.aca_$field_id.residual.miriad

	# Apply the ACA primary beam to TP clean models
	rm -rf temp1
	rm -rf single_input.aca_$field_id.demos.miriad
        demos map=single_input.aca_$field_id.deconv.miriad vis=$name_ytla'.'$field_id'.uv.miriad' \
	      out=temp
	mv temp1 single_input.aca_$field_id.demos.miriad

     done


     # clean up
     rm -rf temp.beam
     rm -rf tp_beam
     rm -rf single_input.aca_*.temp.miriad

fi

##########################################################



##### Convolve TP map to visibility ######################

if [ $if_tp2vis == 'yes' ]
then

   rm -rf uv_random.miriad
   uvrandom npts=$npts freq=$linerestfreq inttime=10 uvmax=$uvmax nchan=1 \
            gauss=true out=uv_random.miriad

   
   for field_id in $fields_ytla
     do

	rm -rf single_input.aca_$field_id.uvmodel.miriad

	# faking headers to avoid crashes.
#        puthd in=single_input.aca_$field_id.demos.miriad/ctype3 value='VELO-LSR' type=ascii
#        puthd in=single_input.aca_$field_id.demos.miriad/cunit3 value='km/s    ' type=ascii
#        puthd in=single_input.aca_$field_id.demos.miriad/crval3 value=$crval3_tp type=double
#        puthd in=single_input.aca_$field_id.demos.miriad/cdelt3 value=$cdelt3_tp type=double
#        puthd in=single_input.aca_$field_id.demos.miriad/crpix3 value=$crpix3_tp type=double
#	puthd in=single_input.aca_$field_id.demos.miriad/restfreq value=$crval3_tp type=double
        puthd in=single_input.aca_$field_id.demos.miriad/ctype3 value='FREQ' type=ascii
        puthd in=single_input.aca_$field_id.demos.miriad/cdelt3 value=2.240000 type=double
        puthd in=single_input.aca_$field_id.demos.miriad/crpix3 value=1 type=double
        puthd in=single_input.aca_$field_id.demos.miriad/crval3 value=$obsfreq type=double
	puthd in=single_input.aca_$field_id.demos.miriad/restfreq value=$obsfreq type=double
        puthd in=single_input.aca_$field_id.demos.miriad/naxis3 value=1 type=int
        puthd in=single_input.aca_$field_id.demos.miriad/bmaj value=$pbfwhm_tp,arcsec type=double
        puthd in=single_input.aca_$field_id.demos.miriad/bmin value=$pbfwhm_tp,arcsec type=double
        puthd in=single_input.aca_$field_id.demos.miriad/bpa  value=0,degree type=double
#        puthd in=single_input.aca_$field_id.demos.miriad/bunit value='Jy/beam' type=ascii


        uvmodel vis=uv_random.miriad model=single_input.aca_$field_id.demos.miriad \
                options=replace,imhead \
                out=single_input.aca_$field_id.uvmodel.miriad "select=uvrange(0,$uvmax)" \


        rm -rf temp
        uvputhd vis=single_input.aca_$field_id.uvmodel.miriad hdvar='telescop' \
		varval='TP  ' type=a out=temp
        rm -rf single_input.aca_$field_id.uvmodel.miriad
        mv temp single_input.aca_$field_id.uvmodel.miriad

        pb="gaus("$pbfwhm_ytla")"
        puthd in=single_input.aca_$field_id.uvmodel.miriad/telescop \
                 value='single' type=a
        puthd in=single_input.aca_$field_id.uvmodel.miriad/pbtype \
                 value=$pb type=a



     done

fi

##########################################################



##### Manually reweight the TP visibilities ##############
if  [ $if_tprewt == 'yes' ]
then

   echo '##### Reweighting TP visibility assuming Tsys ='$tsys_tp' Kelvin'

   for field_id in $fields_ytla
     do

	 outname=single_input.aca_$field_id.uvmodel.rewt.miriad
         rm -rf $outname
         uvputhd vis=single_input.aca_$field_id.uvmodel.miriad hdvar=systemp type=r length=1 \
		 varval=$tsys_tp out=$outname
	 puthd in=$outname/jyperk value=1.0 type=r

	 pb="gaus("$pbfwhm_ytla")"
	 echo $pb '##################'
         puthd in=$outname/telescop \
               value='single' type=a
         puthd in=$outname/pbtype \
               value=$pb type=a


     done

fi

##########################################################




##### Make a copy of relevant files for imaging ##########
if [ $if_duplicateACATP == 'yes' ]
then

   rm -rf intermediate_vis
   mkdir intermediate_vis

   # YTLA
#   uvlin vis=$path_ytla/w51_lsb_epoch2000 out=./intermediate_vis/w51_lsb_epoch2000 \
#	  chans='100,400' mode=chan0
#   uvlin vis=$path_ytla/w51_usb_epoch2000 out=./intermediate_vis/w51_usb_epoch2000 \
#          chans='100,400' mode=chan0
#   cp -r ./ytla_vis/Miriad/ch0/w51_lsb_epoch2000.ch0.miriad ./intermediate_vis/
#   uvflag vis=./ytla_vis/Miriad/ch0/w51_lsb_epoch2000.ch0.miriad select='ant(2)' flagval=flag
#   uvaver vis=./ytla_vis/Miriad/ch0/w51_lsb_epoch2000.ch0.miriad \
#	      out=./intermediate_vis/w51_lsb_epoch2000.ch0.miriad \
#	      stokes='ii' options=nocal,nopass,nopol

   for field_id in $fields_ytla
   do
      cp -r $name_ytla'.'$field_id'.uv.miriad' ./intermediate_vis/INT_$field_id.$outvis_end
   done

   # TP
   for field_id in $fields_ytla
   do
      cp -r 'single_input.aca_'$field_id'.uvmodel.rewt.miriad' ./intermediate_vis/sd_$field_id.$outvis_end
   done

   combinedvis=$name_ytla'.combined.miriad'
   rm -rf $combinedvis
   uvaver "vis=./intermediate_vis/*" options=nocal,nopass,nopol out=./$combinedvis

fi
##########################################################



##### Imaging ACA and TP visibilities together ###########
if [ $if_acaim == 'yes' ]
then


   rm -rf acatp.map
   rm -rf acatp.beam
   invert vis=./$combinedvis options=systemp,double,mosaic,mfs \
         map=acatp.map beam=acatp.beam cell=$aca_cell imsize=$aca_imsize robust=2.0



   rm -rf acatp.model
   if [ $acatp_method == 'clean' ]
   then
       mossdi map=acatp.map beam=acatp.beam out=acatp.model gain=0.1 \
	      niters=$acatp_niters cutoff=$acatp_cutoff \
	      options=positive cutoff=$acatp_cutoff
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

fi
##########################################################



##### Removing meta data #################################
if [ $if_cleanup == 'yes' ]
then
   echo '##### Removing meta data #############'
   rm -rf ./*uvmodel*
   rm -rf ./*regrid*
   rm -rf ./*temp*
   rm -rf ./single_input.miriad
   rm -rf ./single_input*.restor.miriad
   rm -rf ./single_input*.residual.miriad
   rm -rf ./temp.*
   rm -rf ./*.temp
   rm -rf ./uv_random.miriad
   # rm -rf ./acatp*
   # rm -rf ./intermediate_vis
   rm -rf ./*demos*
   rm -rf ./*deconv*
fi
##########################################################



##### Re-immerge #########################################

if [ $if_reimmerge == 'yes' ]
then


   # Gegrid TP maps
   rm -rf GBT.regrid.miriad
   regrid in=./GBT/W51_feathered_img_v3.miriad tin=acatp.clean \
          out=GBT.regrid.miriad \
          project=sin
   puthd in=GBT.regrid.miriad/bmaj value=9.0,arcsec type=d
   puthd in=GBT.regrid.miriad/bmin value=9.0,arcsec type=d
   puthd in=GBT.regrid.miriad/bpa value=0.0,degree type=d
   puthd in=GBT.regrid.miriad/naxis value=2 type=i
   puthd in=GBT.regrid.miriad/ctype3 value='VELO-LSR' type=ascii
   puthd in=GBT.regrid.miriad/cunit3 value='km/s    ' type=ascii
   puthd in=GBT.regrid.miriad/crval3 value=$crval3_tp type=double
   puthd in=GBT.regrid.miriad/cdelt3 value=$cdelt3_tp type=double
   puthd in=GBT.regrid.miriad/crpix3 value=$crpix3_tp type=double
   puthd in=GBT.regrid.miriad/restfreq value=100 type=double
   puthd in=acatp.clean/restfreq value=100 type=double

   pb="gaus("9.0")"
   puthd in=GBT.regrid.miriad/telescop \
               value='single' type=a
   puthd in=GBT.regrid.miriad/pbtype \
               value=$pb type=a

   puthd in=acatp.clean/bmaj value=269.389,arcsec type=d
   puthd in=acatp.clean/bmin value=214.045,arcsec type=d
   puthd in=acatp.clean/bpa value=38.5,degree type=d   
   pb="gaus("200.0")"
   puthd in=acatp.clean/telescop \
               value='single' type=a
   puthd in=acatp.clean/pbtype \
               value=$pb type=a


   rm -rf temp
   regrid in=acatp.clean tin=GBT.regrid.miriad project=sin out=temp axes=1,2

   rm -rf combined.clean.reimmerge
   immerge in=GBT.regrid.miriad,acatp.clean factor=1.0 \
	   out=combined.clean.reimmerge

   rm -rf combined.clean.reimmerge.fits
   fits in=combined.clean.reimmerge op=xyout out=combined.clean.reimmerge.fits

fi

##########################################################
