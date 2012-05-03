c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine y_solve

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     Performs line solves in Y direction by first factoring
c     the block-tridiagonal matrix into an upper triangular matrix, 
c     and then performing back substitution to solve for the unknow
c     vectors of each line.  
c     
c     Make sure we treat elements zero to cell_size in the direction
c     of the sweep.
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer 
     >     c, jstart, stage,
     >     first, last, recv_id, error, r_status(MPI_STATUS_SIZE),
     >     isize,jsize,ksize,send_id

      jstart = 0

c---------------------------------------------------------------------
c     in our terminology stage is the number of the cell in the y-direction
c     i.e. stage = 1 means the start of the line stage=ncells means end
c---------------------------------------------------------------------
      do stage = 1,ncells
         c = slice(2,stage)
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
            call lhsy(c)
            call y_solve_cell(first,last,c)
         else
c---------------------------------------------------------------------
c     Not the first cell of this line, so receive info from
c     processor working on preceeding cell
c---------------------------------------------------------------------
            first = 0
            call y_receive_solve_info(recv_id,c)
c---------------------------------------------------------------------
c     overlap computations and communications
c---------------------------------------------------------------------
            call lhsy(c)
c---------------------------------------------------------------------
c     wait for completion
c---------------------------------------------------------------------
            call mpi_wait(send_id,r_status,error)
            call mpi_wait(recv_id,r_status,error)
c---------------------------------------------------------------------
c     install C'(jstart+1) and rhs'(jstart+1) to be used in this cell
c---------------------------------------------------------------------
            call y_unpack_solve_info(c)
            call y_solve_cell(first,last,c)
         endif

         if (last .eq. 0) call y_send_solve_info(send_id,c)
      enddo

c---------------------------------------------------------------------
c     now perform backsubstitution in reverse direction
c---------------------------------------------------------------------
      do stage = ncells, 1, -1
         c = slice(2,stage)
         first = 0
         last = 0
         if (stage .eq. 1) first = 1
         if (stage .eq. ncells) then
            last = 1
c---------------------------------------------------------------------
c     last cell, so perform back substitute without waiting
c---------------------------------------------------------------------
            call y_backsubstitute(first, last,c)
         else
            call y_receive_backsub_info(recv_id,c)
            call mpi_wait(send_id,r_status,error)
            call mpi_wait(recv_id,r_status,error)
            call y_unpack_backsub_info(c)
            call y_backsubstitute(first,last,c)
         endif
         if (first .eq. 0) call y_send_backsub_info(send_id,c)
      enddo
10000 format (6(2x,e8.2))

      return
      end
      
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      
      subroutine y_unpack_solve_info(c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     unpack C'(-1) and rhs'(-1) for
c     all i and k
c---------------------------------------------------------------------

      include 'header.h'

      integer i,k,m,n,ptr,c,jstart 

      jstart = 0
      ptr = 0
      do i=0,IMAX-1
         do k=0,KMAX-1
            do m=1,BLOCK_SIZE
               do n=1,BLOCK_SIZE
                  lhs(m,n,cc,i,jstart-1,k,c) = out_buffer(ptr+n)
               enddo
               ptr = ptr+BLOCK_SIZE
            enddo
            do n=1,BLOCK_SIZE
               rhs(n,i,jstart-1,k,c) = out_buffer(ptr+n)
            enddo
            ptr = ptr+BLOCK_SIZE
         enddo
      enddo

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------
      
      subroutine y_send_solve_info(send_id,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     pack up and send C'(jend) and rhs'(jend) for
c     all i and k
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer i,k,m,n,jsize,ptr,c,ip,kp
      integer error,send_id,buffer_size 

      jsize = cell_size(2,c)-1
      ip = cell_coord(1,c) - 1
      kp = cell_coord(3,c) - 1
      buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*
     >     (BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)

c---------------------------------------------------------------------
c     pack up buffer
c---------------------------------------------------------------------
      ptr = 0
      do i=0,IMAX-1
         do k=0,KMAX-1
            do m=1,BLOCK_SIZE
               do n=1,BLOCK_SIZE
                  in_buffer(ptr+n) = lhs(m,n,cc,i,jsize,k,c)
               enddo
               ptr = ptr+BLOCK_SIZE
            enddo
            do n=1,BLOCK_SIZE
               in_buffer(ptr+n) = rhs(n,i,jsize,k,c)
            enddo
            ptr = ptr+BLOCK_SIZE
         enddo
      enddo

c---------------------------------------------------------------------
c     send buffer 
c---------------------------------------------------------------------
      call mpi_isend(in_buffer, buffer_size,
     >     dp_type, successor(2),
     >     SOUTH+ip+kp*NCELLS, comm_solve,
     >     send_id,error)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine y_send_backsub_info(send_id,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     pack up and send U(jstart) for all i and k
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer i,k,n,ptr,c,jstart,ip,kp
      integer error,send_id,buffer_size

c---------------------------------------------------------------------
c     Send element 0 to previous processor
c---------------------------------------------------------------------
      jstart = 0
      ip = cell_coord(1,c)-1
      kp = cell_coord(3,c)-1
      buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      ptr = 0
      do i=0,IMAX-1
         do k=0,KMAX-1
            do n=1,BLOCK_SIZE
               in_buffer(ptr+n) = rhs(n,i,jstart,k,c)
            enddo
            ptr = ptr+BLOCK_SIZE
         enddo
      enddo
      call mpi_isend(in_buffer, buffer_size,
     >     dp_type, predecessor(2), 
     >     NORTH+ip+kp*NCELLS, comm_solve, 
     >     send_id,error)

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine y_unpack_backsub_info(c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     unpack U(jsize) for all i and k
c---------------------------------------------------------------------

      include 'header.h'

      integer i,k,n,ptr,c 

      ptr = 0
      do i=0,IMAX-1
         do k=0,KMAX-1
            do n=1,BLOCK_SIZE
               backsub_info(i,k,n,c) = out_buffer(ptr+n)
            enddo
            ptr = ptr+BLOCK_SIZE
         enddo
      enddo

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine y_receive_backsub_info(recv_id,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     post mpi receives
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer error,recv_id,ip,kp,c,buffer_size
      ip = cell_coord(1,c) - 1
      kp = cell_coord(3,c) - 1
      buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*BLOCK_SIZE
      call mpi_irecv(out_buffer, buffer_size,
     >     dp_type, successor(2), 
     >     NORTH+ip+kp*NCELLS, comm_solve, 
     >     recv_id, error)
      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine y_receive_solve_info(recv_id,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     post mpi receives 
c---------------------------------------------------------------------

      include 'header.h'
      include 'mpinpb.h'

      integer ip,kp,recv_id,error,c,buffer_size
      ip = cell_coord(1,c) - 1
      kp = cell_coord(3,c) - 1
      buffer_size=MAX_CELL_DIM*MAX_CELL_DIM*
     >     (BLOCK_SIZE*BLOCK_SIZE + BLOCK_SIZE)
      call mpi_irecv(out_buffer, buffer_size, 
     >     dp_type, predecessor(2), 
     >     SOUTH+ip+kp*NCELLS,  comm_solve, 
     >     recv_id, error)

      return
      end
      
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine y_backsubstitute(first, last, c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     back solve: if last cell, then generate U(jsize)=rhs(jsize)
c     else assume U(jsize) is loaded in un pack backsub_info
c     so just use it
c     after call u(jstart) will be sent to next cell
c---------------------------------------------------------------------

      include 'header.h'

      integer first, last, c, i, k
      integer m,n,j,jsize,isize,ksize,jstart
      
      jstart = 0
      isize = cell_size(1,c)-end(1,c)-1      
      jsize = cell_size(2,c)-1
      ksize = cell_size(3,c)-end(3,c)-1
      if (last .eq. 0) then
         do k=start(3,c),ksize
            do i=start(1,c),isize
c---------------------------------------------------------------------
c     U(jsize) uses info from previous cell if not last cell
c---------------------------------------------------------------------
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rhs(m,i,jsize,k,c) = rhs(m,i,jsize,k,c) 
     >                    - lhs(m,n,cc,i,jsize,k,c)*
     >                    backsub_info(i,k,n,c)
                  enddo
               enddo
            enddo
         enddo
      endif
      do k=start(3,c),ksize
         do j=jsize-1,jstart,-1
            do i=start(1,c),isize
               do m=1,BLOCK_SIZE
                  do n=1,BLOCK_SIZE
                     rhs(m,i,j,k,c) = rhs(m,i,j,k,c) 
     >                    - lhs(m,n,cc,i,j,k,c)*rhs(n,i,j+1,k,c)
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end

c---------------------------------------------------------------------
c---------------------------------------------------------------------

      subroutine y_solve_cell(first,last,c)

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c     performs guaussian elimination on this cell.
c     
c     assumes that unpacking routines for non-first cells 
c     preload C' and rhs' from previous cell.
c     
c     assumed send happens outside this routine, but that
c     c'(JMAX) and rhs'(JMAX) will be sent to next cell
c---------------------------------------------------------------------

      include 'header.h'

      integer first,last,c
      integer i,j,k,isize,ksize,jsize,jstart

      jstart = 0
      isize = cell_size(1,c)-end(1,c)-1
      jsize = cell_size(2,c)-1
      ksize = cell_size(3,c)-end(3,c)-1

c---------------------------------------------------------------------
c     outer most do loops - sweeping in i direction
c---------------------------------------------------------------------
      if (first .eq. 1) then 
         do k=start(3,c),ksize 
            do i=start(1,c),isize

c---------------------------------------------------------------------
c     multiply c(i,jstart,k) by b_inverse and copy back to c
c     multiply rhs(jstart) by b_inverse(jstart) and copy to rhs
c---------------------------------------------------------------------
               call binvcrhs( lhs(1,1,bb,i,jstart,k,c),
     >                        lhs(1,1,cc,i,jstart,k,c),
     >                        rhs(1,i,jstart,k,c) )

            enddo
         enddo
      endif

c---------------------------------------------------------------------
c     begin inner most do loop
c     do all the elements of the cell unless last 
c---------------------------------------------------------------------
      do k=start(3,c),ksize 
         do j=jstart+first,jsize-last
            do i=start(1,c),isize

c---------------------------------------------------------------------
c     subtract A*lhs_vector(j-1) from lhs_vector(j)
c     
c     rhs(j) = rhs(j) - A*rhs(j-1)
c---------------------------------------------------------------------
               call matvec_sub(lhs(1,1,aa,i,j,k,c),
     >                         rhs(1,i,j-1,k,c),rhs(1,i,j,k,c))

c---------------------------------------------------------------------
c     B(j) = B(j) - C(j-1)*A(j)
c---------------------------------------------------------------------
               call matmul_sub(lhs(1,1,aa,i,j,k,c),
     >                         lhs(1,1,cc,i,j-1,k,c),
     >                         lhs(1,1,bb,i,j,k,c))

c---------------------------------------------------------------------
c     multiply c(i,j,k) by b_inverse and copy back to c
c     multiply rhs(i,1,k) by b_inverse(i,1,k) and copy to rhs
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
         do k=start(3,c),ksize 
            do i=start(1,c),isize

c---------------------------------------------------------------------
c     rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
c---------------------------------------------------------------------
               call matvec_sub(lhs(1,1,aa,i,jsize,k,c),
     >                         rhs(1,i,jsize-1,k,c),rhs(1,i,jsize,k,c))

c---------------------------------------------------------------------
c     B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
c     call matmul_sub(aa,i,jsize,k,c,
c     $              cc,i,jsize-1,k,c,bb,i,jsize,k,c)
c---------------------------------------------------------------------
               call matmul_sub(lhs(1,1,aa,i,jsize,k,c),
     >                         lhs(1,1,cc,i,jsize-1,k,c),
     >                         lhs(1,1,bb,i,jsize,k,c))

c---------------------------------------------------------------------
c     multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
c---------------------------------------------------------------------
               call binvrhs( lhs(1,1,bb,i,jsize,k,c),
     >                       rhs(1,i,jsize,k,c) )

            enddo
         enddo
      endif
 1000 format (6(2x,e8.2))

      return
      end
      


