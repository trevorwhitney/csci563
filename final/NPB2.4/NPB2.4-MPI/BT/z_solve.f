c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     Performs line solves in Z direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix, 
c     and then performing back substitution to solve for the unknow
c     vectors of each line.  
c     
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer c, kstart, stage,
     >     first, last, recv_id, error, r_status(MPI_STATUS_SIZE),
     >     isize,jsize,ksize,send_id

      kstart = 0

c---------------------------------------------------------------------
c     in our terminology stage is the number of the cell in the y-direction
c     i.e. stage = 1 means the start of the line stage=ncells means end
c---------------------------------------------------------------------
      do stage = 1,ncells
         c = slice(3,stage)
         isize = cell_size(1,c) - 1
         jsize = cell_size(2,c) - 1
         ksize = cell_size(3,c) - 1
c---------------------------------------------------------------------
c     set last-cell flag
c---------------------------------------------------------------------
         if (stage .eq. ncells) then
            last = 1
         else
            last = 0
         endif

         if (stage .eq. 1) then
c---------------------------------------------------------------------
c     This is the first cell, so solve without receiving data
c---------------------------------------------------------------------
            first = 1
            call lhsz(c)
            call z_solve_cell(first,last,c)
         else
c---------------------------------------------------------------------
c     Not the first cell of this line, so receive info from
c     processor working on preceeding cell
c---------------------------------------------------------------------
            first = 0
            call z_receive_solve_info(recv_id,c)
c---------------------------------------------------------------------
c     overlap computations and communications
c---------------------------------------------------------------------
            call lhsz(c)
c---------------------------------------------------------------------
c     wait for completion
c---------------------------------------------------------------------
            call mpi_wait(send_id,r_status,error)
            call mpi_wait(recv_id,r_status,error)
c---------------------------------------------------------------------
c     install C'(kstart+1) and rhs'(kstart+1) to be used in this cell
c---------------------------------------------------------------------
            call z_unpack_solve_info(c)
            call z_solve_cell(first,last,c)
         endif

         if (last .eq. 0) call z_send_solve_info(send_id,c)
      enddo

c---------------------------------------------------------------------
c     now perform backsubstitution in reverse direction
c---------------------------------------------------------------------
      do stage = ncells, 1, -1
         c = slice(3,stage)
         first = 0
         last = 0
         if (stage .eq. 1) first = 1
         if (stage .eq. ncells) then
            last = 1
c---------------------------------------------------------------------
c     last cell, so perform back substitute without waiting
c---------------------------------------------------------------------
            call z_backsubstitute(first, last,c)
         else
            call z_receive_backsub_info(recv_id,c)
            call mpi_wait(send_id,r_status,error)
            call mpi_wait(recv_id,r_status,error)
            call z_unpack_backsub_info(c)
            call z_backsubstitute(first,last,c)
         endif
         if (first .eq. 0) call z_send_backsub_info(send_id,c)
      enddo
10000 format (6(2x,e8.2))

      return
      end
      
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      
      subroutine z_unpack_solve_info(c)
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     unpack C'(-1) and rhs'(-1) for
c     all i and j
c---------------------------------------------------------------------

      include 'header.h'

      integer i,j,m,n,ptr,c,kstart 

      kstart = 0
      ptr = 0
      do i=0,IMAX-1
         do j=0,JMAX-1
            do m=1,BLOCK_SIZE
               do n=1,BLOCK_SIZE
                  lhs(m,n,cc,i,j,kstart-1,c) = out_buffer(ptr+n)
               enddo
               ptr = ptr+BLOCK_SIZE
            enddo
            do n=1,BLOCK_SIZE
               rhs(n,i,j,kstart-1,c) = out_buffer(ptr+n)
            enddo
            ptr = ptr+BLOCK_SIZE
         enddo
      enddo

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      
      subroutine z_send_solve_info(send_id,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     pack up and send C'(kend) and rhs'(kend) for
c     all i and j
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer i,j,m,n,ksize,ptr,c,ip,jp
      integer error,send_id,buffer_size

      ksize = cell_size(3,c)-1
      ip = cell_coord(1,c) - 1
      jp = cell_coord(2,c) - 1
      buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*
     >     (BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)

c---------------------------------------------------------------------
c     pack up buffer
c---------------------------------------------------------------------
      ptr = 0
      do i=0,IMAX-1
         do j=0,JMAX-1
            do m=1,BLOCK_SIZE
               do n=1,BLOCK_SIZE
                  in_buffer(ptr+n) = lhs(m,n,cc,i,j,ksize,c)
               enddo
               ptr = ptr+BLOCK_SIZE
            enddo
            do n=1,BLOCK_SIZE
               in_buffer(ptr+n) = rhs(n,i,j,ksize,c)
            enddo
            ptr = ptr+BLOCK_SIZE
         enddo
      enddo

c---------------------------------------------------------------------
c     send buffer 
c---------------------------------------------------------------------
      call mpi_isend(in_buffer, buffer_size,
     >     dp_type, successor(3),
     >     BOTTOM+ip+jp*NCELLS, comm_solve,
     >     send_id,error)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_send_backsub_info(send_id,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     pack up and send U(jstart) for all i and j
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer i,j,n,ptr,c,kstart,ip,jp
      integer error,send_id,buffer_size

c---------------------------------------------------------------------
c     Send element 0 to previous processor
c---------------------------------------------------------------------
      kstart = 0
      ip = cell_coord(1,c)-1
      jp = cell_coord(2,c)-1
      buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      ptr = 0
      do i=0,IMAX-1
         do j=0,JMAX-1
            do n=1,BLOCK_SIZE
               in_buffer(ptr+n) = rhs(n,i,j,kstart,c)
            enddo
            ptr = ptr+BLOCK_SIZE
         enddo
      enddo

      call mpi_isend(in_buffer, buffer_size,
     >     dp_type, predecessor(3), 
     >     TOP+ip+jp*NCELLS, comm_solve, 
     >     send_id,error)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_unpack_backsub_info(c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     unpack U(ksize) for all i and j
c---------------------------------------------------------------------

      include 'header.h'

      integer i,j,n,ptr,c

      ptr = 0
      do i=0,IMAX-1
         do j=0,JMAX-1
            do n=1,BLOCK_SIZE
               backsub_info(i,j,n,c) = out_buffer(ptr+n)
            enddo
            ptr = ptr+BLOCK_SIZE
         enddo
      enddo

      return
      end


c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_receive_backsub_info(recv_id,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     post mpi receives
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer error,recv_id,ip,jp,c,buffer_size
      ip = cell_coord(1,c) - 1
      jp = cell_coord(2,c) - 1
      buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      call mpi_irecv(out_buffer, buffer_size,
     >     dp_type, successor(3), 
     >     TOP+ip+jp*NCELLS, comm_solve, 
     >     recv_id, error)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_receive_solve_info(recv_id,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     post mpi receives 
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer ip,jp,recv_id,error,c,buffer_size
      ip = cell_coord(1,c) - 1
      jp = cell_coord(2,c) - 1
      buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*
     >     (BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)
      call mpi_irecv(out_buffer, buffer_size,
     >     dp_type, predecessor(3), 
     >     BOTTOM+ip+jp*NCELLS, comm_solve,
     >     recv_id, error)

      return
      end
      
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_backsubstitute(first, last, c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     back solve: if last cell, then generate U(ksize)=rhs(ksize)
c     else assume U(ksize) is loaded in un pack backsub_info
c     so just use it
c     after call u(kstart) will be sent to next cell
c---------------------------------------------------------------------

      include 'header.h'

      integer first, last, c, i, k
      integer m,n,j,jsize,isize,ksize,kstart
      
      kstart = 0
      isize = cell_size(1,c)-end(1,c)-1      
      jsize = cell_size(2,c)-end(2,c)-1
      ksize = cell_size(3,c)-1
      if (last .eq. 0) then
         do j=start(2,c),jsize
            do i=start(1,c),isize
c---------------------------------------------------------------------
c     U(jsize) uses info from previous cell if not last cell
c---------------------------------------------------------------------
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rhs(m,i,j,ksize,c) = rhs(m,i,j,ksize,c) 
     >                    - lhs(m,n,cc,i,j,ksize,c)*
     >                    backsub_info(i,j,n,c)
                  enddo
               enddo
            enddo
         enddo
      endif
      do k=ksize-1,kstart,-1
         do j=start(2,c),jsize
            do i=start(1,c),isize
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rhs(m,i,j,k,c) = rhs(m,i,j,k,c) 
     >                    - lhs(m,n,cc,i,j,k,c)*rhs(n,i,j,k+1,c)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine z_solve_cell(first,last,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(KMAX) and rhs'(KMAX) will be sent to next cell.
c---------------------------------------------------------------------

      include 'header.h'

      integer first,last,c
      integer i,j,k,isize,ksize,jsize,kstart

      kstart = 0
      isize = cell_size(1,c)-end(1,c)-1
      jsize = cell_size(2,c)-end(2,c)-1
      ksize = cell_size(3,c)-1

c---------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c---------------------------------------------------------------------
      if (first .eq. 1) then 
         do j=start(2,c),jsize 
            do i=start(1,c),isize

c---------------------------------------------------------------------
c     multiply c(i,j,kstart) by b_inverse and copy back to c
c     multiply rhs(kstart) by b_inverse(kstart) and copy to rhs
c---------------------------------------------------------------------
               call binvcrhs( lhs(1,1,bb,i,j,kstart,c),
     >                        lhs(1,1,cc,i,j,kstart,c),
     >                        rhs(1,i,j,kstart,c) )

            enddo
         enddo
      endif

c---------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c---------------------------------------------------------------------
      do k=kstart+first,ksize-last
         do j=start(2,c),jsize
            do i=start(1,c),isize

c---------------------------------------------------------------------
c     subtract A*lhs_vector(k-1) from lhs_vector(k)
c     
c     rhs(k) = rhs(k) - A*rhs(k-1)
c---------------------------------------------------------------------
               call matvec_sub(lhs(1,1,aa,i,j,k,c),
     >                         rhs(1,i,j,k-1,c),rhs(1,i,j,k,c))

c---------------------------------------------------------------------
c     B(k) = B(k) - C(k-1)*A(k)
c     call matmul_sub(aa,i,j,k,c,cc,i,j,k-1,c,bb,i,j,k,c)
c---------------------------------------------------------------------
               call matmul_sub(lhs(1,1,aa,i,j,k,c),
     >                         lhs(1,1,cc,i,j,k-1,c),
     >                         lhs(1,1,bb,i,j,k,c))

c---------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,j,1) by b_inverse(i,j,1) and copy to rhs
c---------------------------------------------------------------------
               call binvcrhs( lhs(1,1,bb,i,j,k,c),
     >                        lhs(1,1,cc,i,j,k,c),
     >                        rhs(1,i,j,k,c) )

            enddo
         enddo
      enddo

c---------------------------------------------------------------------
c     Now finish up special cases for last cell
c---------------------------------------------------------------------
      if (last .eq. 1) then
         do j=start(2,c),jsize
            do i=start(1,c),isize

c---------------------------------------------------------------------
c     rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
c---------------------------------------------------------------------
               call matvec_sub(lhs(1,1,aa,i,j,ksize,c),
     >                         rhs(1,i,j,ksize-1,c),rhs(1,i,j,ksize,c))

c---------------------------------------------------------------------
c     B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
c     call matmul_sub(aa,i,j,ksize,c,
c     $              cc,i,j,ksize-1,c,bb,i,j,ksize,c)
c---------------------------------------------------------------------
               call matmul_sub(lhs(1,1,aa,i,j,ksize,c),
     >                         lhs(1,1,cc,i,j,ksize-1,c),
     >                         lhs(1,1,bb,i,j,ksize,c))

c---------------------------------------------------------------------
c     multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
c---------------------------------------------------------------------
               call binvrhs( lhs(1,1,bb,i,j,ksize,c),
     >                       rhs(1,i,j,ksize,c) )

            enddo
         enddo
      endif
 1000 format (6(2x,e8.2))

      return
      end
      





