c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: ssetrt
c
c\Description:
c  Optional early termination of ssaupd/pssaupd for WHICH = 'LA' or
c  'SA'. After ssetrt(THRESH), the next run(s) stop as soon as every
c  wanted Ritz value on the wanted side of THRESH has converged and at
c  least one converged Ritz value lies beyond it:
c     'LA': all wanted Ritz values >= THRESH converged, one < THRESH;
c     'SA': all wanted Ritz values <= THRESH converged, one > THRESH.
c  NEV then acts as an upper bound on the number of eigenvalues and
c  IPARAM(5) returns how many were found. This is useful with
c  shift-invert, to obtain "all modes up to a cutoff" without knowing
c  how many there are. By default (and after sclrrt) the threshold is
c  off and the solvers behave as before. ssetrt/sclrrt are the single
c  precision counterparts (separate setting).
c
c  The setting is process-global state (like the debug/stat commons) and
c  is not thread safe. It is ignored for other values of WHICH.
c
c  ssetrt(THRESH)     enable with threshold THRESH
c  sclrrt()           disable
c  ssrthk(...)        test used by ssaup2/pssaup2
c
c\EndDoc
c-----------------------------------------------------------------------
      subroutine ssetrt (thresh)
      Real thresh
      logical on
      Real thr
      save on, thr
      data on / .false. /, thr / 0.0E+0 /
      logical lon
      Real lthr
      on = .true.
      thr = thresh
      return
      entry sclrrt ()
      on = .false.
      return
      entry sgetrt (lon, lthr)
      lon = on
      lthr = thr
      return
      end
c
      subroutine ssrthk (which, nev, ritz, bounds, tol, nconv, hit)
      character  which*2
      integer    nev, nconv
      logical    hit
      Real tol, ritz(nev), bounds(nev)
      integer    j
      logical    on, wanted, beyond
      Real thr, eps23, temp
      Real slamch
      external   sgetrt, slamch
c
      hit = .false.
      call sgetrt (on, thr)
      if (.not. on .or. nconv .le. 0) return
      if (which .ne. 'LA' .and. which .ne. 'SA') return
c
      eps23 = slamch('Epsilon-Machine')**(2.0E+0 / 3.0E+0)
      hit = .true.
      beyond = .false.
      do j = 1, nev
         temp = max( eps23, abs(ritz(j)) )
         if (which .eq. 'LA') then
            wanted = ritz(j) .ge. thr
         else
            wanted = ritz(j) .le. thr
         end if
         if (wanted) then
            if (bounds(j) .gt. tol*temp) hit = .false.
         else if (bounds(j) .le. tol*temp) then
            beyond = .true.
         end if
      end do
      hit = hit .and. beyond
      return
      end
