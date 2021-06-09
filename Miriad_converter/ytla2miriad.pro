pro load_ytla, ytla, filename, verbose=verbose
; ################################################################
;
; Loading the raw YTLA HDF5 data to IDL
;
; Version:
;    v0. : Developed by Kyle Lin on 2019.Oct.15
;    v1. : Modified by Baobab on 2019.Oct.18
;
; Input  : 
;     ytla     : A variable name to handle the loaded YTLA data.
;     filename : [str] name of the input YTLA HDF5 data

; Example :
;   .compile ytla2miriad.pro
;   filename = 'W51-ab.ytla7X.mrgh5'
;   load_ytla, ytla, filename, /verbose
;
;   .RESET_SESSION to clean up the memory
;
; ################################################################


;; Sanity check --------------------------------------------------
if ( N_ELEMENTS(filename) EQ 0 ) then begin
   print, 'QUIT! Please specify input filename'
   return
endif else begin
   if KEYWORD_SET(verbose) then begin
      print, '# file summary:'
      h5_list, filename
      print, ' '
   endif
endelse
;; ---------------------------------------------------------------


;; Loading YTLA HDF5 data ----------------------------------------

fid = h5f_open(filename)

;; loading the LO freq
lomhz = h5a_read(  h5a_open_name(fid, 'LO')  )
ytla  = create_struct('LO', lomhz)

bwmhz = h5a_read(  h5a_open_name(fid, 'bw')  )
ytla  = create_struct(ytla, 'BW', bwmhz)

;; loading the visibilities
cstruct = h5d_read(  h5d_open(fid, 'cross')  )
cross   = dcomplex(cstruct.r, cstruct.i)
ytla    = create_struct(ytla, 'cross', cross)
if KEYWORD_SET(verbose) then begin
   print, '# array size of /cross:'
   s = size(cross)
   print, '## ndims, nunits, nch, nbl, nsb, type, size'
   print, s
   print, ''
endif

;; loading the variances
vstruct = h5d_read(  h5d_open(fid, 'variance')  )
cvar    = dcomplex(vstruct.r, vstruct.i)
ytla    = create_struct(ytla, 'variance', cvar)
ytla    = create_struct(ytla, 'nbsl', ( size(ytla.cross) )[3] ) ;number of baselines
ytla    = create_struct(ytla, 'npol', 1 ) ;number of polarizations
if KEYWORD_SET(verbose) then begin
  print, '# array size of /variance:'
  s = size(cvar)
  print, '## ndims, nunits, nch, nbl, nsb, type, size'
  print, s
  print, ''
endif
weight = 2. / (vstruct.r + vstruct.i)
ytla   = create_struct(ytla, 'weight', weight)

;; loading the flags
flag = h5d_read(  h5d_open(fid, 'flag')  )
ytla = create_struct(ytla, 'flag', flag)
if KEYWORD_SET(verbose) then begin
   print, '# array size of /flag:'
   s = size(flag)
   print, '## ndims, nunits, nch, nbl, nsb, type, size'
   print, s
   print, ''
endif

;; loading the baseline lengths
blmeter = h5d_read(  h5d_open(fid, 'blmeter')  )
ytla    = create_struct(ytla, 'blmeter', blmeter)
if KEYWORD_SET(verbose) then begin
   print, '# array size of /blmeter:'
   s = size(blmeter)
   print, '## ndims, nunits, nxy, nbl, type, size'
   print, s
   print, ''
endif

;; loading observing date
obsdate = h5d_read(  h5d_open(fid, 'obsdate')  )
ytla    = create_struct(ytla, 'obsdate', obsdate)
if KEYWORD_SET(verbose) then begin
   print, '# array size of /obsdate:'
   s = size(obsdate)
   print, '## ndims, nint, type, size'
   print, s
   print, ''
endif

;; loading time stamp and convert to Juliand date
time    = h5d_read(  h5d_open(fid, 'time')  )
ytla    = create_struct(ytla, 'time', time)
ytla         = create_struct(ytla, 'nint', (size(ytla.time))[1]) ; number of integrations
if KEYWORD_SET(verbose) then begin
   print, '# array size of /time:'
   s = size(time)
   print, '## ndims, nint, type, size'
   print, s
   print, ''
endif

nunit = n_elements(time)
jd    = dblarr(nunit)
for k = 0, nunit-1 do begin
   jd[k] = systime(elapsed=time[k], /julian, /utc)
endfor
ytla  = create_struct(ytla, 'jd', jd)


;; loading target
target  = h5d_read(  h5d_open(fid, 'target')  )
ytla    = create_struct(ytla, 'target', target)
if KEYWORD_SET(verbose) then begin
   print, '# array size of /target:'
   s = size(target)
   print, '## ndims, nint, type, size'
   print, s
   print, ''
endif


;; loading pointings
pointing    = h5d_read(  h5d_open(fid, 'pointing')  )
pnt_header  = h5d_read(  h5d_open(fid, 'pnt_header')  )
ytla         = create_struct(ytla, 'pointing', pointing)
ytla         = create_struct(ytla, 'pnt_header', pnt_header)
if KEYWORD_SET(verbose) then begin
   print, '# pointing header'
   print, pnt_header
   print, '# array size of /pointing:'
   s = size(pointing)
   print, '## ndims, nint, type, size'
   print, s
   print, ''
endif

;; loading the pointing offset and reference coordinates
pid  = h5d_open(fid, 'poffset')
poff = h5d_read(pid)
rc   = h5a_read(  h5a_open_name(pid, 'ref_coord')  )
ytla = create_struct(ytla, 'poffset', poff, 'ref_coord', rc, 'note_poff', $ 
                     'pointing offset dRA,dDEC in arcmin',                $
                     'note_rc', 'ref_coord RA,DEC in deg J2000')


; help, ytla

h5f_close, fid
;; ---------------------------------------------------------------


end



function ra_radiantohms, ra_radian
; ################################################################
;
; Converting R.A. from radian to hh:mm:ss
;
; Input  : 
;     ra_radian : [double] R.A. in radian
;
; Return :
;     ra_h     : [double] hh
;     ra_m     : [double] mm
;     ra_s     : [double] ss
;
; ################################################################

  ra_tmp = ( ra_radian / (2d*!PI) ) * 24d
  ra_h = floor(ra_tmp)
  ra_tmp = (ra_tmp - ra_h) * 60d
  ra_m = floor(ra_tmp)
  ra_tmp = (ra_tmp - ra_m) * 60d
  ra_s = ra_tmp

  return, [ra_h, ra_m, ra_s]

end


function dec_radiantodms, dec_radian
; ################################################################
;
; Converting Decl. from radian to dd:mm:ss
;
; Input  : 
;     dec_radian : [double] Decl. in radian
;
; Return :
;     dec_d     : [double] hh
;     dec_m     : [double] mm
;     dec_s     : [double] ss
;
; ################################################################

  dec_tmp = double( dec_radian / (2d*!PI) * 360d)
  if (dec_tmp ge 0) then begin
    dec_d = floor(dec_tmp)
    dec_tmp = double((dec_tmp - dec_d) * 60.)
    dec_m = floor(dec_tmp)
    dec_tmp = double((dec_tmp - dec_m) * 60.)
    dec_s = dec_tmp
  endif else begin
    dec_d = ceil(dec_tmp)
    dec_tmp = double((dec_d - dec_tmp) * 60.)
    dec_m = floor(dec_tmp)
    dec_tmp = double((dec_tmp - dec_m) * 60.)
    dec_s = dec_tmp
  endelse

  return, [dec_d, dec_m, dec_s]

end



PRO JULDATE, DATE, JD, PROMPT = prompt
;+                                                                  
; NAME:
;     JULDATE
; PURPOSE:                                   
;     Convert from calendar to Reduced Julian Date
;
; EXPLANATION:
;     Julian Day Number is a count of days elapsed since Greenwich mean noon 
;     on 1 January 4713 B.C.  The Julian Date is the Julian day number
;     followed by the fraction of the day elapsed since the preceding noon. 
;
;     This procedure duplicates the functionality of the JULDATE() function in
;     in the standard IDL distribution, but also allows interactive input and
;     gives output as Reduced Julian date (=JD - 2400000.)  
;     (Also note that prior to V5.1 there was a bug in JULDATE() that gave 
;     answers offset by 0.5 days.)
;
; CALLING SEQUENCE:
;     JULDATE, /PROMPT           ;Prompt for calendar Date, print Julian Date
;               or
;     JULDATE, date, jd      
;
; INPUT:
;     DATE -  3 to 6-element vector containing year,month (1-12),day, and 
;              optionally hour, minute, and second all specified as numbers
;              (Universal Time).   Year should be supplied with all digits.
;              Years B.C should be entered as negative numbers (and note that
;              Year 0 did not exist).  If Hour, minute or seconds are not 
;              supplied, they will default to 0. 
;
;  OUTPUT:
;       JD - Reduced Julian date, double precision scalar.  To convert to
;               Julian Date, add 2400000.   JULDATE will print the value of
;               JD at the terminal if less than 2 parameters are supplied, or 
;               if the /PROMPT keyword is set
;      
;  OPTIONAL INPUT KEYWORD:
;       /PROMPT - If this keyword is set and non-zero, then JULDATE will prompt
;               for the calendar date at the terminal.
;
;  RESTRICTIONS:
;       The procedure HELIO_JD can be used after JULDATE, if a heliocentric
;       Julian date is required.
;
;  EXAMPLE:
;       A date of 25-DEC-1981 06:25 UT may be expressed as either
;
;       IDL> juldate, [1981, 12, 25, 6, 25], jd       
;       IDL> juldate, [1981, 12, 25.2673611], jd 
;
;       In either case, one should obtain a Reduced Julian date of 
;       JD = 44963.7673611
;
;  PROCEDURE USED:
;       GETOPT()
;  REVISION HISTORY
;       Adapted from IUE RDAF (S. Parsons)                      8-31-87
;       Algorithm from Sky and Telescope April 1981   
;       Added /PROMPT keyword, W. Landsman    September 1992
;       Converted to IDL V5.0   W. Landsman   September 1997
;       Make negative years correspond to B.C. (no year 0), work for year 1582
;       Disallow 2 digit years.    W. Landsman    March 2000
;-
 On_error,2

 if ( N_params() EQ 0 ) and (not keyword_set( PROMPT ) ) then begin
     print,'Syntax - JULDATE, date, jd          or JULDATE, /PROMPT'
     print, $
     '  date - 3-6 element vector containing [year,month,day,hour,minute,sec]'
     print,'  jd - output reduced julian date (double precision)'
     return
 endif

 if ( N_elements(date) EQ 0 ) then begin

    opt = ''
    rd: read,' Enter Year,Month,Day,Hour, Minute, Seconds (All Numeric): ',opt
    date = getopt( opt, 'F' )

 endif

 case N_elements(date) of

    6:
    5: date = [ date, 0.0d]
    4: date = [ date, 0.0d,0.0d]
    3: date = [ date, 0.0d, 0.0d,0.0d]
    else: message,'Illegal DATE Vector - must have a least 3 elements'

  endcase

 iy = floor( date[0] )
 if iy lt 0 then iy = iy +1  else $
    if iy EQ 0 then message,'ERROR - There is no year 0'
 im = fix( date[1] )
 date = double(date)
 day = date[2] + ( date[3] + date[4]/60.0d + date[5]/3600.0d) / 24.0d
;
 if ( im LT 3 ) then begin   ;If month is Jan or Feb, don't include leap day

     iy= iy-1 & im = im+12

 end

 a = long(iy/100)
 ry = float(iy)

 jd = floor(ry*0.25d) + 365.0d*(ry -1860.d) + fix(30.6001d*(im+1.)) + $
      day  - 105.5d

;Gregorian Calendar starts on Oct. 15, 1582 (= RJD -100830.5)
 if jd GT -100830.5 then jd = jd + 2 - a + floor(a/4)

 if N_params() LT 2 or keyword_set( PROMPT) then begin
    yr = fix( date[0] )
    print, FORM='(A,I4,A,I3,A,F9.5)',$
       ' Year ',yr,'    Month', fix(date[1] ),'    Day', day
    print, FORM='(A,F15.5)',' Reduced Julian Date:',JD
 endif

 return
 end                                  ; juldate



PRO EPO2JUL, epoch, code, julian
; copied from epo2jul in miriad library ephem.o

  if (CODE eq ' ') then begin
     julian = (EPOCH gt 1984)
  endif else begin
     julian = ((CODE eq 'J') or (CODE eq 'j'))
  endelse

  if (julian) then begin
     julian = 365.25d0       *(EPOCH-2000.0) + 2451545.0d0
  endif else begin
     julian = 365.242198781d0 * (EPOCH-1900.0) + 2415020.31352d0
  endelse

end



PRO JULLST, jday, longitude, lst

;  Converted from Jullst(jday,longitude,lst) in MIRIAD /subs/ephem.for
;
;c  Reference: Explanatory Supplement to the Astronomical Almanac, p50-52.
;c  Accuracy appears to be 0.01 sec of time.
;c
;c  Input:
;c    jday       Julian day of interest.
;c    longitude      Observatory longitude (radians). East of Greenwich is
;c               positive.
;c  Output:
;c    lst        Local mean sidereal time (radians), in the range [0,2*pi].
;c
        dpi = 3.14159265358979323846d0

        T = double(floor(jday - 1.0) + 0.5)
        UT = double(jday - T)
        T = double((T - 2451545.0d0) / 36525.0d0)
;c
        GMST = double(24110.54841d0 +(8640184.812866d0 + (0.093104d0 - 6.2d-6*T)*T)*T)
        GMST = double(GMST / (3600.*24.) + UT * (1.002737909350795d0 + (5.9006d-11 - 5.9d-15*T)*T))

        lst = double(2*dpi*((GMST+longitude/(2*dpi))-floor(GMST+longitude/(2*dpi))))
        if (lst lt 0) then lst = lst + 2*dpi

END



FUNCTION MOBLIQ,jday

;c* Mobliqu -- Mean obliquity of the ecliptic
;c
;c  Return the mean obliquity of the ecliptic.
;c
;c  Input:
;c    jday       Julian day.
;c  Output:
;c    mobliq     Mean obliquity of the ecliptic, in radians.
;c
;c  Reference:
;c    Explanatory Supplement ... page 114.

        dpi = double(3.14159265358979323846)

;c
;c  Centuries from J2000
;c
        T = double((jday - 2451545.0d0) / 36525.0d0)
;c
;c  Mean obliquity.
;c
        mobliq = double(84381.448d0 - (46.8150d0+(0.00059d0-0.001813d0*T)*T)*T)
        mobliq = double(dpi/(180.0*3600.0) * mobliq)

        return,mobliq
END



PRO PRECESS, jday1, ra1, dec1, jday2, ra2, dec2
; copied from epo2jul in miriad library ephem.o

;c  A simple precession routine, to precess from one set of mean
;c  equatorial coordinates (RA,DEC), to another at a different epoch.
;c  This is accurate to order 0.3 arcsec over 50 years.
;c
;c  Reference:
;c    Explanatory Supplement to the Astronomical Almanac, 1993. p 105-106.
;c
;c  NOTE: This does not take account of atmospheric refraction,
;c  nutation, aberration nor gravitational deflection.
;c
;c  Input:
;c    jday1      Julian day of the known epoch.
;c    ra1,dec1   RA,DEC at the jday1 epoch (radians).
;c    jday2      Julian day of the new epoch.
;c  Output:
;c    ra2,dec2   Precessed coordinates (radians).

        dpi = 3.14159265358979323846d0

        T = double((jday1 - 2451545.0d0)/36525.0)
        M = dpi/180.0 * (1.2812323d0 + (0.0003879d0 + 0.0000101d0*T)*T)*T
        N = dpi/180.0 * (0.5567530d0 - (0.0001185d0 + 0.0000116d0*T)*T)*T
        rm = ra1 - 0.5*(M + N*sin(ra1)*tan(dec1))
        dm = dec1 - 0.5*N*cos(rm)
;c
;c  J2000 coordinates.
;c
        r0 = ra1 - M - N*sin(rm)*tan(dm)
        d0 = dec1 - N*cos(rm)
;c
;c  Coordinates of the other epoch.
;c
        T = double((jday2 - 2451545.0d0)/36525.0)
        M = dpi/180 * (1.2812323d0 + (0.0003879d0 + 0.0000101d0*T)*T)*T
        N = dpi/180 * (0.5567530d0 - (0.0001185d0 + 0.0000116d0*T)*T)*T
        rm = r0 + 0.5*(M + N*sin(r0)*tan(d0))
        dm = d0 - 0.5*N*cos(rm)
;c
        ra2 = r0 + M + N*sin(rm)*tan(dm)
        dec2 = d0 + N*cos(rm)

end



PRO Nutate, jday, rmean, dmean, rtrue, dtrue
; copied from epo2jul in miriad library ephem.o

;c
;c  Convert between mean and true equatorial coordinates, by
;c  accounting for nutation.
;c
;c  Input:
;c    jday       Julian day.
;c    rmean,dmean Mean (RA,DEC) at jday.
;c  Output:
;c    rtrue,dtrue True (RA,DEC) at jday.
;c
;c  Nutation parameters.
;c
        nuts,jday,dpsi,deps

;c  True obliquity.

        eps = mobliq(jday) + deps

;c  Various parameters.
        sineps = sin(eps)
        coseps = cos(eps)
        sinra  = sin(rmean)
        cosra  = cos(rmean)
        tandec = tan(dmean)

        rtrue = rmean + (coseps + sineps*sinra*tandec)*dpsi - cosra*tandec*deps
        dtrue = dmean + sineps*cosra*dpsi + sinra*deps

end



PRO ABERRATE,jday,ra,dec,raap,dapp,libfile=libfile, smamiriad=smamiriad

; common global

;c* Aberrate -- Convert RA,DEC from true to geocentric apparent coords.

;c
;c  Account for the effect of annual aberration, to convert
;c  from a true (RA,DEC) to a geocentric apparent (RA,DEC).
;c
;c  Input:
;c    jday       Julian date.
;c    ra,dec     True (RA,DEC).
;c  Output:
;c    rapp,dapp  Geocentric apparent (RA,DEC).

        cmks = double(299792458.0d0)

        pos=make_array(3,/double)
        vel=make_array(3,/double)

        if not keyword_set(libfile) then begin
           if (strpos(!VERSION.ARCH,'86') ge 0) then begin
              if (strpos(!VERSION.ARCH,'64') ge 0) then begin
                 if keyword_set(smamiriad) then libfile=e.idl_sav+'libidlsmamiriad.so' else libfile=e.idl_sav+'libidl64mir.so'
              endif else begin
                 libfile=e.idl_sav+'libidlmir.so'
              endelse
           endif else begin
              print, 'QUIT !!! NO MIRIAD LIBRARY AVAILABLE FOR ',!VERSION.ARCH
              return
           endelse
        endif
        result=CALL_EXTERNAL(libfile, $
          'idl_vearth',double(jday),pos,vel)

        sinra = double(sin(ra))
        cosra = double(cos(ra))
        sindec = double(sin(dec))
        cosdec = double(cos(dec))

        rapp = double(ra +  (-vel(0)*sinra + vel(1)*cosra)/(0.001*cmks*cosdec))
        dapp = double(dec + (-vel(0)*cosra*sindec - vel(1)*sinra*sindec + vel(2)*cosdec)/(0.001*cmks))

END



function get_jyperK, reff, rant
; ################################################################
;
; Obtaining the Jansky per Kelvin conversion factor.
;
; Input:
;   reff      : [double] aperature efficiency, dimensionless
;   rant      : [double] antenna radius in meter
;
; ################################################################

   jyperK = float((2 * 1.38e3 / !PI) / (reff * rant^2))
   return, jyperK

end



function get_mirpol, polar
; ################################################################
;
; Obtaining the Miriad polarization code array mirpol,
; based on the input polar index.
;
; Input  :
;   polar     : [int] 0: non-polarization data / 1: circular / 2: linear
;
; Return :
;   mirpol    : [int array] Miriad polarization codes
;
; ################################################################

  if ((polar gt 2) or (polar lt 1)) then polar = 0

  if (polar eq 0) then begin
    print,"The MIR data are considered as non-polarization data, will use XX as default."
    mirpol = long(-5)
  endif
  if (polar eq 1) then begin
    print,"The polarization codes in MIR data are considered as circularly polarization."
    ; 0 -> unknown, 1 -> rr, 2-> rl, 3-> lr, 4-> ll
    mirpol = [long(0),long(-1),long(-3),long(-4),long(-2)]
  endif
  if (polar eq 2) then begin
    print,"The polarization codes in MIR data are considered as linearly polarization."
    ; 0 -> unknown, 1 -> xx, 2-> xy, 3-> yx, 4-> yy
    mirpol = [long(0),long(-5),long(-7),long(-8),long(-6)]
  endif

  return, mirpol

end



function get_veltype, vtype
; ################################################################
;
; Yielding the veltype [string] header
;
; Input  :
;   vtype   : [int] an index from 0 to 4
;
; ################################################################

  case vtype of
    0   : veltype='VELO-LSR'
    1   : veltype='VELO-HEL'
    2   : veltype='VELO-HEL'
    3   : veltype='VELO-OBS'
    4   : veltype='YTLA_TMP'
  else: begin
        print,'Velocity reference system UNrecognized!!!'
        veltype='VELO-UNK'
        endelse
  endcase

  return, veltype

end



function ytla_countpnts, ytla, mosaic_flag
; ################################################################
;
; Counting total number of distinct pointings
;
; Input  : 
;     ytla     : A variable name to handle the loaded YTLA data.
;
; Return :
;     npts     : [int] total number of dinstinct pointings
;
; ################################################################

  ; Checking source positions/pointings
  ra_list  = reform( ytla.poffset[0, *] )
  dec_list = reform( ytla.poffset[1, *] )

  distinct_ras  = ra_list(  uniq( ra_list,  sort(ra_list)   ) )
  distinct_decs = dec_list( uniq( dec_list, sort(dec_list)  ) )
  nra           = n_elements(distinct_ras)
  ndec          = n_elements(distinct_decs)

  ra_index  = make_array(n_elements(ra_list)  ,/int)
  dec_index = make_array(n_elements(dec_list) ,/int)

  for i = 0, (nra-1) do begin
     matchidx = where( ra_list eq distinct_ras[i] )
     ra_index[matchidx] = i
  endfor

  for i = 0, (ndec - 1) do begin
     matchidx = where( dec_list eq distinct_decs[i] )
     dec_index[matchidx] = i
  endfor
  off_index = ra_index * ndec + dec_index

  tmp_off_list = uniq(off_index)
  noff         = n_elements( tmp_off_list )

  if ( noff gt 1 ) then begin
     mosaic_flag = 1
     print, "*** Dataset with mosaic fields ***"
  endif else begin
     mosaic_flag = 0
     print, "*** Dataset with a single source/pointing ***"
  endelse

  print, "The dataset has ", noff, " positions/pointings."
  return, noff

end



function ytla_getbaselines, nants
; ################################################################
;
; Yielding the ytla_baselines array:
;    each element is a 2-elements array to register ant1 and ant2 
;    for each baseline.
;
; Input  :
;   nants    : [int] number of antennae
;
; Return :
;   ytla_baselines: [int array] as described at the beginning.
;
; ################################################################
   
   n_baselines = nants * (nants-1) / 2
   ytla_baselines = make_array(n_baselines, 2 ,/int)

   idx = 0

   for i = 0, nants-1, 1 do begin
      for j = 0, nants-1, 1 do begin

        if ( j gt i ) then begin
           ytla_baselines[idx, 0] = i+1
           ytla_baselines[idx, 1] = j+1
           idx = idx + 1
        endif

      endfor
   endfor

   return, ytla_baselines

end



function ytla_getobscooord, ytla, srcra, srcdec
; ################################################################
;
; NOT IMPLEMENTED! place holder
; Yielding the apparent geocentric coordinate systems.
;
; Input  :
;    ytla     : A variable name to handle the loaded YTLA data.
;               This is passed from the load_ytla task.
;    srcra    : [double] ra
;    srcdec   : [double] dec
;
; Output :
;    obscoord : [2 element double array] ra and dec in the 0th and 1th elements
;
; ################################################################

    ; *** MIR IDL source code
    ;epoch  = float(in[pil[0]].epoch)
    ;epo2jul,epoch,' ',eporef

    ;getdate,c.ref_time(in(pil(0)).iref_time), 0.0, truedate
    ;juldate, truedate, timeref
    ;timeref=timeref+2400000

    ;precess,eporef,srcra,srcdec,timeref,obsra,obsdec
    ;nutate,timeref,obsra,obsdec,r0,d0
    ;aberrate,timeref,r0,d0,obsra,obsdec, libfile=libfile, smamiriad=smamiriad

    obscoord = make_array(2, /double)
    obscoord[0] = srcra
    obscoord[1] = srcdec

    return, obscoord

end



function ytla_getinttime, ytla, integration
; ################################################################
;
; Obtaining integration time information for a certian ingegration.
;
; Input  :
;    ytla         : A variable name to handle the loaded YTLA data.
;    integration  : [int] the ID of the integration.
;
; Output :
;    inttime      : [double] integration time in units of seconds
;
; #################################################################

   inttime = ytla.pointing[integration, 3] - ytla.pointing[integration, 2]

   return, inttime

end



function ytla_getjd, ytla, integration
; #################################################################
;
; obtaining observing date/time and convert it to Julian date
;
; Input:
;   ytla        : A variable name to handle the loaded YTLA data
;   integraion  : [int] the ID of the integration
;
; Output:
;   jd          : the true Juliand date of the observing unit
;                 not the reduced Juliand date (do not add 2400000. to this)
;
; #################################################################

   epoch = (ytla.pointing[integration, 2] + ytla.pointing[integration, 3]) / 2.d
   jd = systime(elapsed=epoch, /julian, /utc)

   return, jd
end



function ytla_getdate, ytla, integration
; ################################################################
;
; Obtaining observing date information for a certian ingegration.
;
; Input  :
;    ytla         : A variable name to handle the loaded YTLA data.
;    integration  : [int] the ID of the integration.
;
; Output :
;     true        :  4-element vector containing year,month (1-12),day, and 
;                    hour, minute (Universal Time).   
;                    Year should be supplied with all digits.
;                    Years B.C should be entered as negative numbers (and note that
;                    Year 0 did not exist).  If Hour, minute or seconds are not 
;                    supplied, they will default to 0. 
;
; #################################################################

   indate = (ytla.obsdate[integration])
   yr     = double( '20' + strmid(indate, 0, 2) )
   mon    = double( strmid(indate, 2, 2) )
   day    = double( strmid(indate, 4, 2) )
   hr     = 0d  ; this is a place holder
   truedate = [yr, mon, day, hr]

   return, truedate

end



pro ytla_wttoTsys, ytla, integration, wbtsys, wbflag, delta_nu, delta_t, jyperK, sideband, ytla_baselines
; ################################################################
;
; Converting YTLA (baseline based) weight to (baseline based) Tsys
; for a specific integration
;
; Input  :
;    ytla         : A variable name to handle the loaded YTLA data.
;    integration  : [int] the ID of the integration.
;    wt       : [double] YTLA baseline based weight
;    delta_nu : [double] bandwidth in units of Hz
;    delta_t  : [double] integration time in unit of second
;    jyperK   : [double] the conversion factor from Jy to Kelvin,
;                        which can be given by function get_jyperK.
;    sideband : [str] 'usb' or 'lsb'. No default.  (need to include an option for double sideband)
;    ytla_baselines: [int array] 
;
;
; Update :
;    wbtsys     : [double array] baseline based system temperature
;    wbflag : [int array   ] tsys flag
;
;
; ################################################################

  ; allowing only one sideband and one spectral window per sideband
  nspec = 1

  if (sideband eq 'lsb') then sidebandid = 0 ; lsb
  if (sideband eq 'usb') then sidebandid = 1 ; usb


  for j = 0, ( ytla.nbsl - 1 ) do begin

    ant1 = ytla_baselines[j, 0]
    ant2 = ytla_baselines[j, 1]

    wt = sqrt(  total( ytla.weight[integration, 10:510, j, sidebandid] )  )
    Var = (1d / wt)
    Tsysbl = sqrt(   ( Var * 500d * delta_nu * delta_t)  ) / jyperK


    ; initialize narrow band (spectral window) baseline-based tsys array
    ; presently only one spectral window is allowed
    for k = 0, ( nspec - 1 ) do begin

      wbtsys[(ant1-1), (ant2-1), k] = Tsysbl
      wbtsys[(ant2-1), (ant1-1), k] = Tsysbl
      wbflag[(ant1-1), (ant2-1), k] = 1

    endfor


  endfor

end



pro ytla_getwtsys, ytla, nants, wbtsys, wbflag, wtsys, wsolflag, ant_wflag
; ################################################################
;
; Solving the antenna-based wideband tsys array
;
; Input  : 
;     ytla     : A variable name to handle the loaded YTLA data.
;    nants     : [int] total number of antennae
;   wbtsys     : baseline based Tsys table
;   wbflag     : baseline based Tsys flag
;
; Update :
;      wtsys    : Antenna based Tsys table
;   wsolflag    : Tsys channel flag
;  ant_wflag    : Antenna flag
;
; ################################################################

  ; presently only permits outputing one sideband
  nsb = 1

  ; assuming no duplication of antennae in the list
  inflag = make_array(nants,/int,value=0)


  for k=0,( nsb - 1 ) do begin

    tosb = 0 ; allowing only one sideband

    resultx=max(total(wbflag[*,*,tosb],1), wbflagx)
    resulty=max(total(wbflag[*,*,tosb],2), wbflagy)

    if (resultx ge resulty) then begin
      refant=wbflagx
      antl=where(wbflag[*,refant,tosb] gt 0)
    endif else begin
      refant=wbflagy
      antl=where(wbflag[refant,*,tosb] gt 0)
    endelse

    antn=n_elements(antl)

    wsolflag[tosb]=0
    if (antn ge 2) then begin
       tosol=1
       p=0
       q=p+1
       wsolflag[tosb]=1
    endif else begin
       tosol=0
    endelse

    while (tosol eq 1) do begin
      an1=antl(p)
      an2=antl(q)
      if (wbflag[an1,an2,tosb] eq 1) then begin
         tref=wbtsys[refant,an1,tosb]*wbtsys[refant,an2,tosb]/wbtsys[an1,an2,tosb]
         wsolflag[tosb]=1
         tosol = 0
      endif else begin
        if (q lt antn-1) then begin
          q=q+1
        endif else begin
          if (p lt antn-2) then begin
            p=p+1
            q=p+1
          endif else begin
            tosol = 0
          endelse
        endelse
      endelse
    endwhile


    ; wide band tsys
    if (wsolflag[tosb] eq 1) then begin
      wtsys[*,tosb] =  wbtsys[refant,*,tosb]^2/tref
      wtsys[refant,tosb] = tref
      ant_wflag[*,tosb] = 1L


      zerolist = where((wtsys[*,tosb] eq 0) and (inflag[*] eq 1), zerocnt)
      if (zerocnt gt 0) then begin
        for iz = 0, zerocnt -1 do begin
          tmpxlist = where((wbtsys[*,zerolist[iz],tosb]) gt 0, tmpxcnt)
          tmpylist = where((wbtsys[zerolist[iz],*,tosb]) gt 0, tmpycnt)
          if (tmpxcnt eq 0 and tmpycnt eq 0) then begin
            print, " ### NOTE: NO TSYS DATA FOR ANT ",(zerolist[iz]+1), " AT INT ",int_list[i], $
                   "  SIDEBAND ", c.sb[all_sb[k]], " ###"
            print,'Arbitrary system temperature of 400 K is inserted for this antenna'
            print,'and data on associated baselines will be flagged bad'
            wtsys[zerolist[iz],tosb] = 400.0
            ant_wflag[zerolist[iz],tosb] = 0
          endif else begin
            zsol=0
            if (tmpxcnt gt 0 and zsol eq 0) then begin
              itmp = 0
              while (zsol eq 0 and itmp lt tmpxcnt) do begin
                if (ant_wflag[tmpxlist[itmp],tosb] eq 1) then begin
                  wtsys[zerolist[iz],tosb] =  wbtsys[tmpxlist[itmp],zerolist[iz],tosb]^2/wtsys[tmpxlist[itmp],tosb]
                  zsol=1
                endif
                itmp = itmp + 1
              endwhile
            endif

            if (tmpycnt gt 0 and zsol eq 0) then begin
              itmp = 0
              while (zsol eq 0 and itmp lt tmpycnt) do begin
                if (ant_wflag[tmpylist[itmp],tosb] eq 1) then begin
                  wtsys[zerolist[iz],tosb] =  wbtsys[zerolist[iz],tmpylist[itmp],tosb]^2/wtsys[tmpylist[itmp],tosb]
                  zsol=1
                endif
                itmp = itmp + 1
              endwhile
            endif

            if (zsol eq 0) then begin
              print, " ### NOTE: NOT SOLVABLE FOR ANT ",(zerolist[iz]+1), " AT INT ",int_list[i], $
                     "  SIDEBAND ", " ###"
              print,'Arbitrary system temperature of 400 K is inserted for this antenna'
              print,'and data on associated baselines will be flagged bad'
              wtsys[zerolist[iz],tosb] = 400.0
              ant_wflag[zerolist[iz],tosb] = 0
            endif
          endelse
        endfor
      endif

    endif else begin
      print,'No solution for wideband ',' at this integration ','.'
      print,'Arbitrary system temperature of 400 K is inserted for all antennas'
      print,'and data will be flagged bad'
      wtsys[*,tosb] = 400.0
    endelse


  endfor



end



pro ytla_gettsys, ytla, nants, nbtsys, nbflag, tsys, solflag, ant_flag
; ################################################################
;
; Solving the antenna-based tsys array
;
; Input  : 
;     ytla     : A variable name to handle the loaded YTLA data.
;    nants     : [int] total number of antennae
;   nbtsys     : baseline based Tsys table
;   nbflag     : baseline based Tsys flag
;
; Update :
;       tsys    : Antenna based Tsys table
;    solflag    : Tsys channel flag
;   ant_flag    : Antenna flag
;
; ################################################################

  ; presently only permits one spectral window for each sideband
  nspect = 1

  ; assuming no duplication of antennae in the list
  inflag = make_array(nants,/int,value=0)


  for k=0,(nspect-1) do begin


    resultx=max(total(nbflag[*,*,k],1),nbflagx)
    resulty=max(total(nbflag[*,*,k],2),nbflagy)

    if (resultx ge resulty) then begin
      refant=nbflagx
      antl=where(nbflag[*,refant,k] gt 0)
    endif else begin
      refant=nbflagy
      antl=where(nbflag[refant,*,k] gt 0)
    endelse

    antn=n_elements(antl)

    solflag[k]=0
    if (antn ge 2) then begin
       tosol=1
       p=0
       q=p+1
       solflag[k]=1
    endif else begin
       tosol=0
    endelse

    while (tosol eq 1) do begin
      an1=antl(p)
      an2=antl(q)
      if (nbflag[an1,an2,k] eq 1) then begin
         tref=nbtsys[refant,an1,k]*nbtsys[refant,an2,k]/nbtsys[an1,an2,k]
         solflag[k]=1
         tosol = 0
      endif else begin
        if (q lt antn-1) then begin
          q=q+1
        endif else begin
          if (p lt antn-2) then begin
            p=p+1
            q=p+1
          endif else begin
            tosol = 0
          endelse
        endelse
      endelse
    endwhile

    ; narrow band tsys
    if (solflag[k] eq 1) then begin
      ; channel tsys
      tsys[*,k] =  nbtsys[refant,*,k]^2/tref
      tsys[refant,k] = tref
      ant_flag[*,k] = 1L

      zerolist = where((tsys[*,k] eq 0) and (inflag[*] eq 1), zerocnt)
      if (zerocnt gt 0) then begin
        for iz = 0, zerocnt -1 do begin
          tmpxlist = where((nbtsys[*,zerolist[iz],k]) gt 0, tmpxcnt)
          tmpylist = where((nbtsys[zerolist[iz],*,k]) gt 0, tmpycnt)
          if (tmpxcnt eq 0 and tmpycnt eq 0) then begin
            tsys[zerolist[iz],k] = 400.0
            ant_flag[zerolist[iz],k] = 0
          endif else begin
            zsol=0
            if (tmpxcnt gt 0 and zsol eq 0) then begin
              itmp = 0
              while (zsol eq 0 and itmp lt tmpxcnt) do begin
                if (ant_flag[tmpxlist[itmp],k] eq 1) then begin
                  tsys[zerolist[iz],k] =  nbtsys[tmpxlist[itmp],zerolist[iz],k]^2/tsys[tmpxlist[itmp],k]
                  zsol=1
                endif
                itmp = itmp + 1
              endwhile
            endif

            if (tmpycnt gt 0 and zsol eq 0) then begin
              itmp = 0
              while (zsol eq 0 and itmp lt tmpycnt) do begin
                if (ant_flag[tmpylist[itmp],k] eq 1) then begin
                  tsys[zerolist[iz],k] =  nbtsys[zerolist[iz],tmpylist[itmp],k]^2/tsys[tmpylist[itmp],k]
                  zsol=1
                endif
                itmp = itmp + 1
              endwhile
            endif

            if (zsol eq 0) then begin
              tsys[zerolist[iz],k] = 400.0
              ant_flag[zerolist[iz],k] = 0
            endif
          endelse
        endfor
      endif

    endif else begin
      print,'No solution for narrowband ',k,' at integration ',int_list[i],'.'
      print,'Arbitrary system temperature of 400 K is inserted for all antennas'
      print,'and data will be flagged bad'
      ; channel tsys
      tsys[*,k] = 400.0
    endelse

  endfor



end



function ytla_getpreamble, ytla, intid, blid, LO, ytla_baselines
; ################################################################
;
; Obtaining the visibility header 'Preamble' for Miriad
; Preamble is a four element array.
; The zeroth and first elements are u and v coordinates in units of nanoseconds;
; the second element is time in units of Julian date; the last element
; is the baseline number.
;
; baseline number is defined as bl = 256 * A1 + A2, 
; where A1 and A2 are teh index of the first and second antennae,
; respectively. Need to ensure that A1 < A2.
;
; Input   :
;     intid     : [int] index of integration
;     blid      : [int] index of baseline
;     LO        : [double] Observing frequency in units of [GHz]
; ytla_baselines: [int array]
;
; ################################################################

  cvel_mks       = 299792458d
  kilowavelength_mks = (cvel_mks * 1e3) / ( LO * 1e9 )

  preamble = make_array(4,/double)

  ; visibility header/data
  preamble[0]   = -( ( ytla.blmeter[intid,0,blid] / kilowavelength_mks ) * 1e3 ) / LO
  preamble[1]   = -( ( ytla.blmeter[intid,1,blid] / kilowavelength_mks ) * 1e3 ) / LO
  preamble[2]   = ytla.jd[intid]

  ant1 = ytla_baselines[blid, 0]
  ant2 = ytla_baselines[blid, 1]
  preamble[3] = 2048*ant1 + ant2 + 65536

  return, preamble

end



function ytla_getcomplexvis, ytla, integration, blid, sideband, nschan, flags
; ################################################################
;
; Returning complex visibility in a format that Miriad like.
;
; Input  : 
;     ytla        : A variable name to handle the loaded YTLA data.
;     integration : [int] integration id
;     blid        : [int] baseline index
;     sideband    : [str] 'usb' or 'lsb'. No default.  (need to include an option for double sideband)
;     nschan      : [int] number of spectral channels
;
; Return :
;     complexvis  : [complex array] real and imaginary parts
;
; Updates:
;     flags	  : [long] flags of data 1 = good; 0 = bad?
;
; ################################################################

  flags = make_array(nschan, /long, value=1)

  if (sideband eq 'lsb') then sidebandid = 0 ; lsb
  if (sideband eq 'usb') then sidebandid = 1 ; usb

  real = reform(  real_part( ytla.cross[integration, 0:(nschan-1), blid, sidebandid] )  )
  img  = reform(  imaginary( ytla.cross[integration, 0:(nschan-1), blid, sidebandid] )  )

   infid = where( finite(real) eq 0 )
     real[infid] = 0d
     img[infid]  = 0d
     flags[infid] = 0
   infid = where( finite(img) eq 0 )
     real[infid] = 0d
     img[infid]  = 0d
     flags[infid] = 0

   maskid = where( ytla.flag[integration, 0:(nschan-1), blid, sidebandid] gt 0)
   flags[maskid] = 0

   complexvis = complex(real, img)

  return, complexvis

end



pro ytla2miriad, ytla,                                                      $
                 dir=dir, source=source, sideband=sideband,                 $
                ; band=band, 
                 edge_trim=edge_trim,polar=polar,verbose=verbose, $
                 oldweight=oldweight,   $
                 libpath=libpath, $
                 ; libfile=libfile, 
                 state=state, $
                 smamiriad=smamiriad,  $
                 ; trimflux=trimflux,
                 pntsplit=pntsplit, $
                 intsplit=intsplit, $
                 auto=auto
; ################################################################
;
; Exporting YTLA data to Miriad file
;
; Dependence:
;    Depending on the external library files in the ../idl_sav directory.
;    They were copied overed from the MIR IDL software package.
;
;
; Version:
;    v0. : Baobab Modified from idl2Miriad on 2019.Oct.19
;
;
; Input  :
;    ytla     : A variable name to handle the loaded YTLA data.
;               This is passed from the load_ytla task.
;
;    dir     : [str] the output directory
;
;    source  : [str] output target source name. Default: 'TEMP'
;
;    sideband : [str] 'usb' or 'lsb'. No default.  (need to include an option for double sideband)
;
;
;    libpath : [str] the path where the external library for exporting data locates.
;                    Default: '../idl_sav/'
;
; Keyword :
;    pntsplit : If set, split individual pointings to separated Miriad files.
;               Files are indexed with integer numbers.
;
;    intsplit : If set, will split individual integrations into two.
;               This is to cheat Miriad. Miriad does not operate correctly if
;               a file only has a single integration.
;               This keyword option can only be activated when /pntsplit is activated.
;
;
; Example :
;
;   .RESET_SESSION to clean up the memory
;
; ################################################################


; common global
; common data_set


;; Physical constants
  cvel_mks = 299792458d

; ***** check verbose keyword *****

if (not keyword_set(verbose)) then verbose = 0


;; Setting path to the library file ------------------------------
  if not KEYWORD_SET(libpath) then begin
    tmp = getenv('IDLLIB')
    if (tmp eq '') then begin
      libpath = '../idl_sav/'
    endif else begin
      libpath = tmp
    endelse
  endif

  if (strpos(!VERSION.ARCH,'86') ge 0) then begin
     if (strpos(!VERSION.ARCH,'64') ge 0) then begin
        if KEYWORD_SET(smamiriad) then begin 
           libfile=libpath+'libidlsmamiriad.so' 
        endif else begin 
           libfile=libpath+'libidl64mir.so'
        endelse
     endif else begin
        libfile=libpath+'libidlmir.so'
     endelse
  endif else begin
     print, 'QUIT !!! No Miriad Library available for ',!VERSION.ARCH
     return
  endelse


; ***** preparations ------------------------------------
  if (not KEYWORD_SET(dir)) then begin
     print,"Please name the output directory"
     return
  endif

  if (file_test(dir) eq 1) then begin
    linux_command = 'rm -rf ' + dir
    spawn, linux_command
  endif

  if (not KEYWORD_SET(source)) then begin
    source = 'TEMP'
    print, 'Setting a tentatitive target source name: TEMP'
  endif

  if (not KEYWORD_SET(sideband)) then begin
    print, 'Please select upper sideband (usb) or lower sideband (lsb).'
    return
  endif else begin
     if ( (sideband ne 'usb') and (sideband ne 'lsb')  ) then begin
        print, sideband
        print, 'Please select upper sideband (usb) or lower sideband (lsb).'
        return
     endif
  endelse

  if ( KEYWORD_SET(pntsplit) ) then begin
    print, "************************************************************************"
    print, "***** Keywork pntsplit is set                                      *****"
    print, "***** Will separate individual pointings to separated Miriad files *****"
    if ( KEYWORD_SET(intsplit) ) then begin
      print, "***** Keywork intsplit is set                                      *****"
      print, "***** Will separate individual pointings to separated Miriad files *****"
    endif
    print, "************************************************************************"
  endif


; ***** observatory/telescope specific information ------

  telescop = 'YTLA'
  version  = 'YTLA 1.0'

  ; currently using SMA's longitude and latitude  YTLA: Needs this meta data ; place holder
  YTLAlatitude = 19.82420526391d0
  YTLAlongitude = (360. - (155.+(28.+37.20394/60.)/60.0))/360.*2.*!PI
  YTLAlatitude = YTLAlatitude * 2 * !PI / 360.0

  ; ***** check polar keyword *****

    if (not keyword_set(polar)) then polar = 0
    mirpol = get_mirpol(polar)

  ; ***** obtaining LO info *****

    ; evaluating LO frequencies in units of GHz
      YTLA_nch = 1024 ; number of channels in each sideband
      YTLALO = ( 84000d + ytla.LO )
      YTLABW = ytla.BW
      if ( sideband eq 'usb' ) then begin
         RF = YTLALO + YTLABW * findgen(YTLA_nch) / 1024d
      endif else begin
         RF = YTLALO - YTLABW * findgen(YTLA_nch) / 1024d
      endelse
      YTLALO = YTLALO / 1000d
      RF     = RF     / 1000d
      LO     = mean(RF)


  ; ***** Antenna/baseline related information *****

    ; number of antennae
    nants = long(7)

    ; assign baselines: an array with n_baselines elements
    ytla_baselines = ytla_getbaselines(nants)

    ; radius of antennae  YTLA: Needs this meta data
    rant = (1.182/2)

    ; fitted beam width = 64.6 * lambda / D (in deg) 
    pbfwhm = 64.6 * (cvel_mks / (LO * 1.e9 * rant * 2d) ) * 3600.	; fwhm im arcsec

    ; get the assumed aperture efficiency
    ;reff = 0.5 will give approximately jyperK = 5000, which relates Tsys=200 to SEFD=1MJy
    reff = 0.5

    ; get jyperK information
    jyperK = get_jyperK(reff, rant)

  ; ***** observed source information *****

    epoch = 2000d

    ; Checking source positions/pointings
    mosaic_flag = 0
    noff        = ytla_countpnts(ytla, mosaic_flag)

    ; Coordintaes of reference pointings
      ; converting from degree to radian units (following the syntax of idl2miriad)
      srcra  = ytla.ref_coord[0] * ( !pi / 180d )
      srcdec = ytla.ref_coord[1] * ( !pi / 180d )

      srcra_hms  = ra_radiantohms( srcra )
      srcdec_dms = dec_radiantodms( srcdec )

      print, ""
      print, "Referencing Position :"
      print, "RA= ", srcra_hms[0],":",srcra_hms[1],":", srcra_hms[2], ", (in radian) ", srcra
      print, "Dec =", srcdec_dms[0],":",srcdec_dms[1],":",srcdec_dms[2], ", (in radian) ", srcdec

  ; ***** Obtaining ra and dec in geocentric coordinate system ***** (may not be important; not implemented yet)

    obscoord = ytla_getobscooord(ytla, srcra, srcdec)


  ; ***** observational setup information *****

    ; *** polarization ***
    ; for YTLA, reading one polarization at a time
    all_pol = [0]
    npol    = long(1)

    ; *** velocity system ***
    ; for YTLA, this is not yet important (2019.Oct.19). Now focusing on continuum
    vtype   = 0 ; Assign a tentative velocity type to bypass many header issues (use vtype=0 for SMA data)
    veltype = get_veltype(vtype)


    ; *** spectral setup ***

      ; currently only permit one total number of sideband (and wideband; either USB or LSB)
      nsb     = 1
      nspect  = 1   ;  currently permits one spectral window per sideband
      sbchunk = 1   ;  currently permits one spectral window per sideband
      sbchan  = YTLA_nch   ; a dummy variable to fit to the syntax of MIR IDL

      speccode = make_array(nspect,/int)     ; an index for wide spectral band
      restfreq = make_array(nspect,/double)  ; rest frequency
      sdf      = make_array(nspect,/double)  ; Frequency resolution
      sfreq    = make_array(nspect,/double)  ; Sky frequency (approximated with rest frequency)
      nschan   = make_array(nspect,/long)    ; Number of spectral channels (fixed to 1024 for YTLA)
      ischan   = make_array(nspect,/long)
      nbw      = make_array(nspect,/float)
      spvdop   = make_array(nspect,/double)  ; doppler velocity

       sbcount = 0
      ; "Rest Freq" in GHz unit
      restfreq[sbcount] = LO

      ; "Freq resolution" in GHz unit
      ;sdf[sbcount] = ( YTLABW / 1024d ) / 1000d
      if ( sideband eq 'usb' ) then begin
         sdf[sbcount] = ( YTLABW / 1024d ) / 1000d
      endif else begin
         sdf[sbcount] = -( YTLABW / 1024d ) / 1000d
      endelse

      ; "Sky Freq for central channel" (tentative, not accurate), in GHz unit
      sfreq[sbcount] = LO
      fsky           = LO

      nschan[sbcount] = 1024

      case vtype of
      0   : begin
              vabs = (1- (fsky/restfreq[sbcount])^2)/(1+(restfreq[sbcount])^2)*cvel_mks/1000.0 ; absolute vel in km/s
              ; sp[psl[j[0]]].vel VLSR recorded
              ; Deriving doppler velocity from the observer to the LSR (Here not implemented yet)
              ; spvdop[sbcount] = (vabs - sp[psl[j[0]]].vel)
              spvdop[sbcount] = vabs ; this line is a place holder in case of header conflicts.
               ; Taking the AVERAGE of Doppler velocity derived over all narrow bands
                ; (to remove the relatively small numerical errors)
                veldop = float(total(spvdop)/nspect)
                vsource =float(0.0)
            end
      else: begin
                print,"CAUTION: calculation of vtype velocity not implemented!!!"
                vsource=0.0
                veldop=0.0
            endelse
      endcase


      ; some more modification for all narrow windows

        ; Shifting the frequency from the window center (MIR/IDL) to the window edge (center of 1st channel) (MIRIAD)
        ; For a chuck of 128ch, for example, MIR has sfreq at the end ch64, use 63.5 * sdf to shift to the center of 1ch.
        ; For a cont. chan of a ch. MIR has sfreq at the center, no need to shift, or a shift of 0.
        sfreq = sfreq - sdf * ( (nschan/2.0) - 0.5 )

        ; narrow band band width
        nbw   = (nschan * sdf)

        ; total number of spectral channels                     
        numchan = long(total(nschan))

        ; first channel-index of a spectral window
        ischan[0]=1
        for sbcount=1, (nspect-1) do begin
           ischan[sbcount] = total(nschan[0:sbcount-1])+1
        endfor


      ; *** wideband spectral setup info ***
        nwide  = 1 ; currently only allow one wide spectral window
        wfreq  = make_array(nwide,/float)
        wwidth = make_array(nwide,/float)

        wfreq[0]  = sfreq
        wwidth[0] = double(nbw)
        wfreq     = wfreq + 0.5 * wwidth

;    print,"speccode",speccode
    print,"restfreq [GHz]",restfreq
    print,"sfreq [GHz]",sfreq
    print,"frequency resolution [MHz]",sdf*1000d
    print,"nschan",nschan
;    print,"ischan",ischan
;    print,"wfreq [GHz]",wfreq
;    print,"wwidth [GHz]",wwidth

  ; *** All header information prepared ***
  one=long(1)

    PRINT,"--- READY TO CYCLE THROUGH ALL INTEGRATIONS ---"
    unit = 0
    for i = 0, ( ytla.nint - 1 ), 1 do begin

        dra  = float(ytla.poffset[0, i])  ; ra  offset in arcmin units
        ddec = float(ytla.poffset[1, i])  ; dec offset in arcmin units
        dra =  dra  * float(2. * !PI / (60. * 360.))
        ddec = ddec * float(2. * !PI / (60. * 360.))

        ; ***** setting output file handle unit *****
        if ( KEYWORD_SET(pntsplit) ) then begin
          mosaic_flag = 0
          dir  = source + '_' + sideband + '.' + strtrim( string(i), 1) + '.miriad'
          if (verbose) then print, ' - - - - -Outputting data into a new file- - - - - '
          if (verbose) then print, dir
        endif

        if ( mosaic_flag eq 0 ) then begin
        ; outputting single-source file
           srcra_hd  = srcra  + dra / cos(srcdec)
           srcdec_hd = srcdec + ddec
           obscoord = ytla_getobscooord(ytla, srcra_hd, srcdec_hd)
        endif else begin
        ; outputting multi-source file
           srcra_hd  = srcra
           srcdec_hd = srcdec
        endelse

        ; control whether or not a new file is opened for output
        if_newfile = 'no'
        if (                           $    
            ( KEYWORD_SET(pntsplit) )  $
           ) then begin
          if_newfile = 'yes'
          unit = 0
        endif

        if ( i eq 0 ) then begin
          if_newfile = 'yes'
          unit = 0
        endif

      if ( if_newfile eq 'yes' ) then begin
        ;print,"--- Creating/Opening new data directory ---"

          result=CALL_EXTERNAL(libfile, $
                           'idl_uvopen',dir,unit)

        ;print,"--- Writing brief history entry ---"
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_hisappend',unit)
          result=CALL_EXTERNAL(libfile, $
                           'idl_hiswrite',unit,'IDL/MIR-to-MIRIAD, Version: Beta')
          result=CALL_EXTERNAL(libfile, $
                           'idl_hiswrite',unit,'based on library version xx-xx-xx')
          result=CALL_EXTERNAL(libfile, $
                           'idl_hiswrite',unit,'Target Source: '+source)
          result=CALL_EXTERNAL(libfile, $
                           'idl_hiswrite',unit,'LO Frequency: '+string(LO)+' GHz')
          result=CALL_EXTERNAL(libfile, $
                           'idl_hisclose',unit)

        ;print,"--- Inserting observer header info ---"
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'a','observer', 'Baobab')
  
        ;print,"--- Inserting observatory specific header info ---"
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'a','telescop', telescop)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'a','version', version)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','latitud',YTLAlatitude,one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','longitu',YTLAlongitude,one)
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'i','nants',nants,one)
       ;    result=CALL_EXTERNAL(libfile, $
       ;                    'idl_uvputvr',unit,'r','jyperk',jyperk,one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'r','pbfwhm',pbfwhm,one)

        ;print,"--- Inserting source header info ---"
  
          if (strlen(source) gt 8) then begin
             print,'NOTICE!!! MIRIAD only allows source name with a maximum of 8 characters'
             print,'Source name ',source,' in the header will be truncated into ',strmid(source,0,8)
             source=strmid(source,0,8)
          endif
          source = strupcase(source)
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'a','source',source)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','ra',srcra_hd,one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','dec',srcdec_hd,one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'r','epoch',float(epoch),one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','obsra',obscoord[0],one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','obsdec',obscoord[1],one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'a','purpose','S')
  
          ; put primary beam type header
          pb_string = 'gaus(' + string(pbfwhm, format='(f5.1)' ) + ')'
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'a','pbtype', pb_string)
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'a','veltype',veltype)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'r','veldop',veldop,one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'r','vsource',vsource,one)
  
        ;print,"--- Inserting observing setup header info ---"
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','lo1',LO,one)
  
        ;print,"---         wide band channel header info ---"
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'i','nwide',nwide,one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'r','wfreq',wfreq,nwide)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'r','wwidth',abs(wwidth),nwide)
  
        ;print,"---       narrow band channel header info ---"
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'i','numchan',numchan,one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'i','nspect',nspect,one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','restfreq',restfreq,nspect)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','sfreq',sfreq,nspect)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','sdf',sdf,nspect)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'i','nschan',nschan,nspect)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'i','ischan',ischan,nspect)
  
        ;print,"---              polarization header info ---"
  
  
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'i','npol',npol,one)
          result=CALL_EXTERNAL(libfile, $
                           'idl_wrhdi',unit,'npol',npol)

        ;print,"---              fake antpos/corr header info ---"
  
          fakecormode = long(1)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'i','cormode',fakecormode,one)
          fakecorfin  = float([100.,200.,300.,400.])
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'r','corfin',fakecorfin,4)
          fakecorbw  = float([20.,20.])
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'r','corbw',fakecorbw,2)
          fakeantpos = make_array(nants*3,/double,value=100.0)
          result=CALL_EXTERNAL(libfile, $
                           'idl_uvputvr',unit,'d','antpos',fakeantpos,nants*3)

        if_newfile = 'no'
      endif ; finish opening new file and inserting header



    ; ###########################################################


         if (verbose) then print,"---  initialize antenna-based tsys matrices ---"
  
  
         ;  wideband channel tsys
         wtsys    = make_array(nants, nspect,/float,value=0.0)
         wsolflag = make_array(nspect,/float,value=0.0)
  
         ;  wideband antenna tsys flag
         ant_wflag = make_array(nants, nspect,/int,value=0)
  
         ;  temporary wideband tsys matrices
         wbtsys  = make_array(nants,nants,1,value=0.0)
         wbflag  = make_array(nants,nants,1,/int,value=0)
  
         ;  narrowband channel tsys
         tsys = make_array(nants,nspect,/float,value=0.0)
         solflag = make_array(nspect,/float,value=0.0)
  
         ;  narrowband antenna tsys flag
         ant_flag = make_array(nants,nspect,/int,value=0)
  
         ;  temporary narrowband matrices
         nbtsys  = make_array(nants,nants,nspect,/float,value=0.0)
         nbflag  = make_array(nants,nants,nspect,/int,value=0)
  
  
         if (verbose) then print,"---          integration time header info ---"
  
           ; obtaining integration time information for integration i
           inttime = float( ytla_getinttime(ytla, i) )

         intrepeats = 1
         if ( KEYWORD_SET(intsplit) ) then begin
           inttime = inttime / 2d
           inttime = float(inttime)
           intrepeats = 2
         endif

         for intrep = 0, intrepeats-1, 1 do begin ; # This for loop permits separating one integration into 
                                             ; # intrepeats times
               
               result=CALL_EXTERNAL(libfile, $
                                   'idl_uvputvr',unit,'r','inttime',inttime,one)
      
      	       jd       = ytla.jd[i]  
               jd_interval = ( double(inttime) / 86400d ) / double(intrepeats)
               if (intrepeats eq 2) then begin
                 jd = double(ytla.jd[i]) - ( 0.5d - double(intrep) )  * jd_interval
               endif
               ut       = double( ((jd-0.5) - floor(jd-0.5)) * 2. * !PI)
               result=CALL_EXTERNAL(libfile, $
                                   'idl_uvputvr',unit,'d','ut',ut, one)
      
               jullst,jd,YTLAlongitude,lst
               result=CALL_EXTERNAL(libfile, $
                                   'idl_uvputvr',unit,'d','lst',lst, one)
      
             if (mosaic_flag eq 1) then begin      
               result=CALL_EXTERNAL(libfile, $
                              'idl_uvputvr',unit,'r','dra',dra, one)
      
               result=CALL_EXTERNAL(libfile, $
                              'idl_uvputvr',unit,'r','ddec',ddec, one)
             endif
      
             if (verbose) then print,'  - cycling through all baselines to get baseline-based tsys  -'
             delta_nu = abs(sdf[0]) * 1e9
             ytla_wttoTsys, ytla, i, wbtsys, wbflag, delta_nu, inttime, jyperK, sideband, ytla_baselines
             ytla_wttoTsys, ytla, i, nbtsys, nbflag, delta_nu, inttime, jyperK, sideband, ytla_baselines
      
             ; Solving the antenna-based tsys values
      
               ; wideband
               ytla_getwtsys, ytla, nants, wbtsys, wbflag, wtsys, wsolflag, ant_wflag
               wtsys2 = reform( wtsys, nants*1 )
               ; narrowband
                ;print, nbtsys
               ytla_gettsys, ytla, nants, nbtsys, nbflag, tsys, solflag, ant_flag
               tsys2 = reform(tsys, nants*nspect)
                ;print, tsys2
    
      
             ; Exporting Tsys
             result=CALL_EXTERNAL(libfile, $
                                  'idl_uvputvr',unit,'r','wsystemp',wtsys2, long(n_elements(wtsys2)))
      
             result=CALL_EXTERNAL(libfile, $
                                  'idl_uvputvr',unit,'r','systemp',tsys2, long(n_elements(tsys2)))
    
    
    
           if (verbose) then print,'  - cycling through all baselines again to store header/data -'
           for j = 0, ( ytla.nbsl - 1 ) do begin
             for jpol = 0, ( ytla.npol -1 ) do begin
    
                   ; exporting polarization code
                   thispol = mirpol[all_pol[jpol]]
                   result=CALL_EXTERNAL(libfile, $
                                        'idl_uvputvr',unit,'i','pol',thispol, one)
    
                   ; setting visibility header
                   preamble = ytla_getpreamble(ytla, i, j, LO, ytla_baselines)

                   if (intrepeats eq 2) then begin ; correct time stamp when splitting one integration into two
                     preamble[2] =  double(ytla.jd[i]) - ( 0.5d - double(intrep) )  * jd_interval
                   endif
 
                   wdata=make_array(nwide,/complex)
                   wflags = make_array(nwide,/long, value=0)
    
                   data=make_array(numchan,/complex)
                   flags = make_array(numchan,/long, value=0)
    
                   chunkcount = 0
                   chancount  = 0
    
                   ; *** wide bands ***
                   if (verbose) then print,'  - continuum channel -'
    
                   for k = 0, ( nsb - 1 ) do begin
    
                     
                     if (sideband eq 'lsb') then sidebandid = 0 ; lsb
                     if (sideband eq 'usb') then sidebandid = 1 ; usb
    
                      wdata[0] = ytla.cross[i,0,j,sidebandid]
    
                   endfor
    
                   ; *** narrow bands ***
                   if (verbose) then print,'  - spectral window -'
    
                   chancount = 0
                   for k = 0, ( nspect - 1 ) do begin
    
                     data = ytla_getcomplexvis( ytla, i, j, sideband, nschan[k], flags)
                     chancount = chancount + nschan[k]
    
                   endfor
    
    

                   result=CALL_EXTERNAL(libfile, $
                                        'idl_uvputvr',unit,'r','jyperk', jyperK, one)

                   result=CALL_EXTERNAL(libfile, $
                                        'idl_uvwwrite',unit,wdata,wflags,nwide)
    
                   result=CALL_EXTERNAL(libfile, $
                                        'idl_uvwrite',unit,preamble,data,flags,numchan)
    
    
             endfor ; loop for polarization
           endfor ; loop for baseline
          endfor  ; loop for intrepeats    

          if ( KEYWORD_SET(pntsplit) ) then begin
            if (verbose) then print,"--- Closing up data directory ---"
            result=CALL_EXTERNAL(libfile, $
                             'idl_uvclose',unit)
          endif

    endfor ; loop for integration


  if ( not KEYWORD_SET(pntsplit) ) then begin
    if (verbose) then print,"--- Closing up data directory ---"
    result=CALL_EXTERNAL(libfile, $
                     'idl_uvclose',unit)
  endif


end




pro run_convert, name
    ;;; given the name, load all associated data into IDL
    ;;; and convert them with ytla2miriad.pro
    ;;; (X and Y; LSB and USB; 4 miriad folder in total)
    ;;;
    ;;; Input:
    ;;;	    name	    a string to identify the data
    ;;;			    e.g. W51-ab
    ;;;
    ;;;	Output:
    ;;;	    <name>.X.lsb.mir
    ;;;	    <name>.X.usb.mir
    ;;;	    <name>.Y.lsb.mir
    ;;;	    <name>.Y.usb.mir
    ;;;

    if (n_elements(name) eq 0) then stop, 'need to provide a name.'
    bname = file_basename(name)
    dname = file_dirname(name)


    loc = name.indexof('.ytla')
    if (loc gt -1) then begin		; strip trailing characters
	name = name.remove(loc, -1)
    endif
    loc2 = bname.indexof('.ytla')
    if (loc2 gt -1) then begin		; strip trailing characters
	bname = bname.remove(loc2, -1)
    endif

    pols = ['X', 'Y']
    nsb  = 2

    foreach p, pols do begin
	pname = '.ytla7' + p

	; ytla input HDF name
	fname = name + pname + '.mrgh5'
	if (file_test(fname)) then begin
	    print, fname, ' --> exists'
	    load_ytla, ytla, fname

	    for s = 0, nsb-1 do begin
		if (s eq 0) then sb = 'lsb'
		if (s eq 1) then sb = 'usb'
		sname = '.' + sb

		; miriad output dir name
		oname = bname + '.' + p + sname + '.vis'
		print, oname

		ytla2miriad, ytla, dir=oname, sideband=sb
	    endfor

	endif else begin
	    print, fname, ' --> not found'
	endelse
    endforeach

end



