pro load_ytla, ytla, verbose=verbose
; examples of HDF5 operation

;print, ytla

file = 'W51-ab.ytla7X.mrgh5'
if keyword_set(verbose) then begin
    print, '# file summary:'
    h5_list, file
    print , ''
endif



fid = h5f_open(file)

;; loading the LO freq
lomhz = h5a_read(h5a_open_name(fid, 'LO'))
ytla = create_struct('LO', lomhz)


;; loading the visibilities
cstruct = h5d_read(h5d_open(fid, 'cross'))
cross = dcomplex(cstruct.r, cstruct.i)
ytla = create_struct(ytla, 'cross', cross)
if keyword_set(verbose) then begin
    print, '# array size of /cross:'
    s = size(cross)
    print, '## ndims, nunits, nch, nbl, nsb, type, size'
    print, s
    print, ''
endif


;; loading the variances
vstruct = h5d_read(h5d_open(fid, 'variance'))
cvar = dcomplex(vstruct.r, vstruct.i)
ytla = create_struct(ytla, 'variance', cvar)
if keyword_set(verbose) then begin
    print, '# array size of /variance:'
    s = size(cvar)
    print, '## ndims, nunits, nch, nbl, nsb, type, size'
    print, s
    print, ''
endif

weight = 2. / (vstruct.r + vstruct.i)
ytla = create_struct(ytla, 'weight', weight)


;; loading the flags
flag = h5d_read(h5d_open(fid, 'flag'))
ytla = create_struct(ytla, 'flag', flag)
if keyword_set(verbose) then begin
    print, '# array size of /flag:'
    s = size(flag)
    print, '## ndims, nunits, nch, nbl, nsb, type, size'
    print, s
    print, ''
endif


;; loading the baseline lengths
blmeter = h5d_read(h5d_open(fid, 'blmeter'))
ytla = create_struct(ytla, 'blmeter', blmeter)
if keyword_set(verbose) then begin
    print, '# array size of /blmeter:'
    s = size(blmeter)
    print, '## ndims, nunits, nxy, nbl, type, size'
    print, s
    print, ''
endif

;; loading the pointing offset and reference coordinates
pid = h5d_open(fid, 'poffset')
poff = h5d_read(pid)
rc = h5a_read(h5a_open_name(pid, 'ref_coord'))
ytla = create_struct(ytla, 'poffset', poff, 'ref_coord', rc, 'note_poff', 'pointing offset dRA,dDEC in arcmin', 'note_rc', 'ref_coord RA,DEC in deg J2000')

help, ytla


h5f_close, fid

end

