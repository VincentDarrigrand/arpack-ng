c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsetrt
c
c\Description:
c  Optional early termination of dsaupd/pdsaupd for WHICH = 'LA' or
c  'SA'. After dsetrt(THRESH), the next run(s) stop as soon as every
c  wanted Ritz value on the wanted side of THRESH has converged and at
c  least one converged Ritz value lies beyond it:
c     'LA': all wanted Ritz values >= THRESH converged, one < THRESH;
c     'SA': all wanted Ritz values <= THRESH converged, one > THRESH.
c  NEV then acts as an upper bound on the number of eigenvalues and
c  IPARAM(5) returns how many were found. This is useful with
c  shift-invert, to obtain "all modes up to a cutoff" without knowing
c  how many there are. By default (and after dclrrt) the threshold is
c  off and the solvers behave as before. ssetrt/sclrrt are the single
c  precision counterparts (separate setting).
c
c  The setting is process-global state (like the debug/stat commons) and
c  is not thread safe. It is ignored for other values of WHICH.
c
c  dsetrt(THRESH)     enable with threshold THRESH
c  dclrrt()           disable
c  dsrthk(...)        test used by dsaup2/pdsaup2
c
c\EndDoc
c-----------------------------------------------------------------------
      subroutine dsetrt (thresh)
      Double precision thresh
      logical on
      Double precision thr
      save on, thr
      data on / .false. /, thr / 0.0D+0 /
      logical lon
      Double precision lthr
      on = .true.
      thr = thresh
      return
      entry dclrrt ()
      on = .false.
      return
      entry dgetrt (lon, lthr)
      lon = on
      lthr = thr
      return
      end
c
      subroutine dsrthk (which, nev, ritz, bounds, tol, nconv, hit)
      character  which*2
      integer    nev, nconv
      logical    hit
      Double precision tol, ritz(nev), bounds(nev)
      integer    j
      logical    on, wanted, beyond
      Double precision thr, eps23, temp
      Double precision dlamch
      external   dgetrt, dlamch
c
      hit = .false.
      call dgetrt (on, thr)
      if (.not. on .or. nconv .le. 0) return
      if (which .ne. 'LA' .and. which .ne. 'SA') return
c
      eps23 = dlamch('Epsilon-Machine')**(2.0D+0 / 3.0D+0)
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
