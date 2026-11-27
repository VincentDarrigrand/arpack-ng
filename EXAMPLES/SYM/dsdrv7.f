      program dsdrv7
c
c     Ritz value threshold: all modes below a cutoff, without knowing
c     how many there are.
c
c     K x = lambda x with K = diag(i**2), i = 1..N. Shift-invert at
c     sigma = 0 gives OP = K**-1 whose Ritz values are 1/lambda, so
c     the modes with lambda <= LMAX are the LARGEST Ritz values of OP
c     (>= 1/LMAX). dsetrt(1/LMAX) makes dsaupd stop as soon as all of
c     them have converged; NEV is only an upper bound.
c
      integer          n, nev, ncv, lworkl, ldv
      parameter        (n=2000, nev=40, ncv=60, ldv=n,
     &                  lworkl=ncv*(ncv+8))
      Double precision lmax
      parameter        (lmax=900.0D+0)
      integer          iparam(11), ipntr(11), ido, info, j, nin, nconv
      logical          select(ncv)
      Double precision resid(n), v(ldv,ncv), workd(3*n),
     &                 workl(lworkl), d(nev), z(n,nev), sigma, tol
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
      call dsaupd ( ido, 'I', n, 'LA', nev, tol, resid, ncv, v, ldv,
     &              iparam, ipntr, workd, workl, lworkl, info )
      if (ido .eq. -1 .or. ido .eq. 1) then
         do 20 j = 1, n
            workd(ipntr(2)+j-1) = workd(ipntr(1)+j-1) / dble(j)**2
 20      continue
         go to 10
      end if
      call dclrrt ()
      if (info .lt. 0) stop 1
c
      call dseupd ( .false., 'A', select, d, z, ldv, sigma, 'I', n,
     &              'LA', nev, tol, resid, ncv, v, ldv, iparam, ipntr,
     &              workd, workl, lworkl, info )
      if (info .ne. 0) stop 1
      nconv = iparam(5)
c
      nin = 0
      do 30 j = 1, nconv
         if (1.0D+0/d(j) .le. lmax*(1.0D+0+1.0D-6)) nin = nin + 1
 30   continue
      print *, 'converged:', nconv, ' modes with lambda <= lmax:', nin
c
c     30 modes (i = 1..30) lie in range; none may be missing and the
c     solver must have stopped before computing all NEV
      if (nin .ne. 30 .or. nconv .ge. nev) stop 1
c
c     Same cutoff with WHICH='SA' and OP = K: stop once all
c     eigenvalues <= LMAX have converged.
c
      ido = 0
      info = 0
      iparam(1) = 1
      iparam(3) = 5000
      iparam(7) = 1
      call dsetrt (lmax)
 40   continue
      call dsaupd ( ido, 'I', n, 'SA', nev, tol, resid, ncv, v, ldv,
     &              iparam, ipntr, workd, workl, lworkl, info )
      if (ido .eq. -1 .or. ido .eq. 1) then
         do 50 j = 1, n
            workd(ipntr(2)+j-1) = workd(ipntr(1)+j-1) * dble(j)**2
 50      continue
         go to 40
      end if
      call dclrrt ()
      if (info .lt. 0) stop 1
      call dseupd ( .false., 'A', select, d, z, ldv, sigma, 'I', n,
     &              'SA', nev, tol, resid, ncv, v, ldv, iparam, ipntr,
     &              workd, workl, lworkl, info )
      if (info .ne. 0) stop 1
      nconv = iparam(5)
      nin = 0
      do 60 j = 1, nconv
         if (d(j) .le. lmax*(1.0D+0+1.0D-6)) nin = nin + 1
 60   continue
      print *, 'SA converged:', nconv, ' lambda <= lmax:', nin
      if (nin .ne. 30 .or. nconv .ge. nev) stop 1
      end
