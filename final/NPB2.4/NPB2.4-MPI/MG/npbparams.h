c NPROCS = 2 CLASS = C
c  
c  
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the class of the NPB
c  in this directory. Do not modify it by hand.
c  
        integer nprocs_compiled
        parameter (nprocs_compiled = 2)
        integer nx_default, ny_default, nz_default
        parameter (nx_default=512, ny_default=512, nz_default=512)
        integer nit_default, lm, lt_default
        parameter (nit_default=20, lm = 9, lt_default=9)
        integer debug_default
        parameter (debug_default=0)
        integer ndim1, ndim2, ndim3
        parameter (ndim1 = 9, ndim2 = 9, ndim3 = 8)
        logical  convertdouble
        parameter (convertdouble = .false.)
        character*11 compiletime
        parameter (compiletime='04 May 2012')
        character*3 npbversion
        parameter (npbversion='2.4')
        character*6 cs1
        parameter (cs1='mpif77')
        character*6 cs2
        parameter (cs2='mpif77')
        character*22 cs3
        parameter (cs3='-L/usr/local/lib -lmpi')
        character*20 cs4
        parameter (cs4='-I/usr/local/include')
        character*4 cs5
        parameter (cs5='-O3 ')
        character*6 cs6
        parameter (cs6='(none)')
        character*6 cs7
        parameter (cs7='randi8')
