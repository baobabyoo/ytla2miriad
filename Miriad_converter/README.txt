1. The two HDF5 files contain the YTLA data on W51 from UT 190925,190926,191003,191004

2. 'load_ytla.pro' is my first attempt to read HDF5 into IDL
    usage:
        .r load_ytla.pro
        load_ytla, data [, /verbose]

    the results are stored in a IDL structure 'data' with the following content:
    ** Structure <3122968>, 10 tags, length=103236856, data length=103236856, refs=1:
       LO              DOUBLE           10000.000
       CROSS           DCOMPLEX  Array[50, 1024, 21, 2]
       VARIANCE        DCOMPLEX  Array[50, 1024, 21, 2]
       WEIGHT          DOUBLE    Array[50, 1024, 21, 2]
       FLAG            LONG64    Array[50, 1024, 21, 2]
       BLMETER         DOUBLE    Array[50, 2, 21]
       POFFSET         DOUBLE    Array[2, 50]
       REF_COORD       DOUBLE    Array[2]
       NOTE_POFF       STRING    'pointing offset dRA,dDEC in arcmin'
       NOTE_RC         STRING    'ref_coord RA,DEC in deg J2000'

    data.LO: the second LO frequency in MHz
    data.cross: complex visibilities for 50 tracking units, 1024 channels, 21 baselines and 2 sidebands(0 for lsb, 1 for usb)
    data.variance: the corresponding estimates of the data variance in real and imag parts
    data.weight: defined as 2/(data.variance.real + data.variance.imag)
    data.flag: corresponding data flag (0 = good; >0 = bad)
    data.blmeter: baseline vectors in X/Y in meters for each tracking units; X/Y corresponds to RA/DEC.
    data.poffset: dRA/dDEC in arcmin from the reference coordinates for each tracking unit
    data.ref_coord: reference coordinates in degree J2000

3. conversion from channels to RF:
    RF = 84000 + LO +/- 2240 * ch/1024
    (+ for USB, - for LSB)


Edit: 2019/Oct/15


