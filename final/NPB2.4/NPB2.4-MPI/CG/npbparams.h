c NPROCS = 2 CLASS = C
c  
c  
c  This file is generated automatically by the setparams utility.
c  It sets the number of processors and the class of the NPB
c  in this directory. Do not modify it by hand.
c  
        integer            na, nonzer, niter
        double precision   shift, rcond
        parameter(  na=150000,
     >              nonzer=15,
     >              niter=75,
     >              shift=110.,
     >              rcond=1.0d-1 )

c number of nodes for which this version is compiled
        integer    nnodes_compiled
        parameter( nnodes_compiled = 2)
        integer    num_proc_cols, num_proc_rows
        parameter( num_proc_cols=2, num_proc_rows=1 )
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
