      program pdsdrv7
c
c     Parallel version of EXAMPLES/SYM/dsdrv7: all modes with
c     lambda <= LMAX of K = diag(i**2) via shift-invert (OP = K**-1)
c     and the Ritz value threshold set with dsetrt. K is distributed
c     by blocks of NLOC rows. Run with any number of processes.
c
      include 'mpif.h'
      integer          nloc, nev, ncv, lworkl
      parameter        (nloc=500, nev=40, ncv=60,
     &                  lworkl=ncv*(ncv+8))
      Double precision lmax
      parameter        (lmax=900.0D+0)
      integer          iparam(11), ipntr(11), ido, info, j, nin, nconv,
     &                 rank, ierr
      logical          select(ncv)
      Double precision resid(nloc), v(nloc,ncv), workd(3*nloc),
     &                 workl(lworkl), d(nev), z(nloc,nev), sigma, tol
c
      call MPI_INIT (ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
c
      tol = 1.0D-10
      ido = 0
      info = 0
      iparam(1) = 1
      iparam(3) = 1000
      iparam(7) = 1
c
      call dsetrt (1.0D+0/lmax)
c
 10   continue
      call pdsaupd ( MPI_COMM_WORLD, ido, 'I', nloc, 'LA', nev, tol,
     &               resid, ncv, v, nloc, iparam, ipntr, workd, workl,
     &               lworkl, info )
      if (ido .eq. -1 .or. ido .eq. 1) then
         do 20 j = 1, nloc
            workd(ipntr(2)+j-1) = workd(ipntr(1)+j-1)
     &                            / dble(rank*nloc+j)**2
 20      continue
         go to 10
      end if
      call dclrrt ()
      if (info .lt. 0) call MPI_ABORT (MPI_COMM_WORLD, 1, ierr)
c
      call pdseupd ( MPI_COMM_WORLD, .false., 'A', select, d, z, nloc,
     &               sigma, 'I', nloc, 'LA', nev, tol, resid, ncv, v,
     &               nloc, iparam, ipntr, workd, workl, lworkl, info )
      if (info .ne. 0) call MPI_ABORT (MPI_COMM_WORLD, 2, ierr)
      nconv = iparam(5)
c
      nin = 0
      do 30 j = 1, nconv
         if (1.0D+0/d(j) .le. lmax*(1.0D+0+1.0D-6)) nin = nin + 1
 30   continue
      if (rank .eq. 0) print *, 'converged:', nconv,
     &                          ' modes with lambda <= lmax:', nin
      if (nin .ne. 30 .or. nconv .ge. nev)
     &   call MPI_ABORT (MPI_COMM_WORLD, 3, ierr)
      call MPI_FINALIZE (ierr)
      end
