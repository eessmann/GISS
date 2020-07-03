program CW263
    !
    !   Stream function wave theory code:
    !   * automatic selection of order
    !   * uniform current
    !   * returns proportion of limiting height
    !   j.r.chaplin@soton.ac.uk
    !
    !   Subroutine cw260 solves the wave.  Kinematics are then available
    !   through subtroutine kmts as illustrated below in the main program.
    !   Stores results in CW263.PSI in a format compatible with that used
    !      by CW6.FOR (24/3/99)
    !
    character ans*1
    !
    10 write(*, '(a)') ' Water depth (m)   = '
    read*, d
    write(*, '(a)') ' Period      (s)   = '
    read*, t
    write(*, '(a)') ' Wave height (m)   = '
    read*, h
    write(*, '(a)') ' Current     (m/s) = '
    read*, u
    nverb = 1
    !
    call cw260(d, t, h, u, nverb, n, el)
    !
    20 write(*, '(/a)')&
            ' (H)orizontal, (V)ertical, (S)urface, (N)ew wave, (Q)uit : '
    read(*, '(a)') ans
    goto (20, 21, 22, 23, 10, 29) (index('HhVvSsNnQq', ans) + 3) / 2
    !
    !   Horizontal
    !
    21 write(*, '(a)') ' y (m) = '
    read*, yy
    npt = 21
    write(*, 3)
    do i = 1, npt
        xx = el * (i - 1) / float(npt - 1)
        tt = 0.0
        call kmts(n, xx, yy, tt, uu, vv, ut, vt, du, dv, etah)
        ans = ' '
        if (yy.gt.etah) ans = '*'
        write(*, '(f9.3,8f8.3,1x,a)') xx, yy, uu, vv, ut, vt, du, dv, etah, ans
    end do
    goto 20
    !
    !   Vertical
    !
    22 write(*, '(a)') ' x/L = '
    read*, xl
    xx = xl * el
    yy = 0.0
    tt = 0.0
    call kmts(n, xx, yy, tt, uu, vv, ut, vt, du, dv, etah)
    npt = 21
    write(*, 3)
    do i = 1, npt
        yy = (d + etah) * (npt - i) / float(npt - 1) - d
        call kmts(n, xx, yy, tt, uu, vv, ut, vt, du, dv, etah)
        write(*, '(f9.3,8f8.3,1x,1a)') xx, yy, uu, vv, ut, vt, du, dv, etah
    end do
    goto 20
    !
    !   Free surface
    !
    23 npt = 21
    write(*, 3)
    do i = 1, npt
        xx = el * (i - 1) / float(npt - 1)
        yy = h
        tt = 0.0
        call kmts(n, xx, yy, tt, uu, vv, ut, vt, du, dv, etah)
        write(*, '(f9.3,8f8.3)') xx, etah, uu, vv, ut, vt, du, dv, etah
    end do
    goto 20
    !
    29 stop
    3 format('      x       y       u       v      ut   ', &
            '   vt      du      dv      eta')
end
!
!
!
subroutine cw260(zd, zt, zh, zu, nverb, nfun, zel)
    !
    !  Input:  zd=depth; zt=period; zh=height; u=current; nverb=verbosity
    !  Output: nfun=order; zel=wavelength
    !
    implicit double precision (a-h, o-z)
    parameter (nmax = 25)
    real zd, zt, zh, zu, zel
    double precision k
    character itl*79, datim*22
    common /one/ d, t, h, u, k&
            /two/ eta(nmax), c(nmax), amp(0:nmax)
    !
    pi = 4 * atan(1.0)
    !
    d = zd
    t = zt
    hw = zh
    u = zu
    !
    !  Get required solution order
    !
    call wavecel(t, d, u, tr, cel)
    call limit(hw, d, tr, rat, 1)
    dl0 = d / (9.81 * t * t / (2 * pi))
    a = 0.86 / sqrt(dl0)
    b = 7 + 2.2 * log(dl0)
    cc = 2.7 - 3 * log(dl0)
    nve = nint((a + b * rat + cc * rat**2) / 2) * 2
    nve = nve + 2
    hb = hw / rat
    !
    k = 2 * pi / (cel * t)
    !
    !  Start with height=<Hb/2
    !    First step up orders; then heights ...
    !
    if (nverb/=0) write(*, '(/2a)')&
            '      d       T       H       U    order  iter  ', &
            'rms error  code      L'
    nfun = 6
    h = min(hw, 0.5 * hb)
    do i = 0, nmax
        amp(i) = 0.0
    end do
    amp(1) = h / 2
    12 continue
    call cw261(nfun, iter, fsumsq, ifail)
    if (ifail==0) then
        el = 2 * pi / k
        if (nverb/=0) write(*, '(a,4f8.3,i5,i7,1p1e12.3,0p,i5,f10.3)')&
                ' ', d, t, h, u, nfun, iter, fsumsq, ifail, el
    else
        if (nverb/=0) write(*, '(a,4f8.3,i5,i7,1p1e12.3,0p,i5)')&
                ' ', d, t, h, u, nfun, iter, fsumsq, ifail
        stop
    endif
    if (nfun<nve) then
        nfun = nfun + 2
        goto 12
    endif
    !
    if (hw>0.5 * hb) then
        fac = 1.1
        11    hm = h
        if (h * fac>hw) then
            h = hw
            ilast = 1
        else
            h = h * fac
            ilast = 0
        endif
        do i = 1, nfun - 1
            amp(i) = (h / hm) * amp(i)
        end do
        call cw261(nfun, iter, fsumsq, ifail)
        if (ifail==0) then
            el = 2 * pi / k
            if (nverb/=0)&
                    write(*, '(a,4f8.3,i5,i7,1p1e12.3,0p,i5,f10.3)')&
                            ' ', d, t, h, u, nfun, iter, fsumsq, ifail, el
        else
            if (nverb/=0) write(*, '(a,4f8.3,i5,i7,1p1e12.3,0p,i5)')&
                    ' ', d, t, h, u, nfun, iter, fsumsq, ifail
            stop
        endif
        if (ilast/=1) then
            fac = 0.995 * fac
            goto 11
        endif
    endif
    zel = 2 * pi / k
    !
    !$$$      open(12,file='cw263.psi')
    !$$$      itl= 'Solved by CW263'
    !$$$      datim= 'Space for date & time:'
    !$$$      ver= 5.01
    !$$$      el0= 9.81*zt**2/(2*pi)
    !$$$      hl0= zh/el0
    !$$$      dl0= zd/el0
    !$$$      wl0= zel/el0
    !$$$      write(12,101) itl,datim,ver,hl0,dl0,wl0,zh,zd,zt,nfun-1,zu
    !$$$      write(12,102) (eta(i),c(i+1),i=1,nfun-1),eta(nfun)
    !$$$      write(12,107) (amp(i)/zh,i=0,nfun-1)
    !$$$      close(12)
    !
    101 format(a/a, f10.2/3f16.10/3f16.10, i5, f16.10)
    102 format(1p2e25.16)
    107 format(1p1e25.16)
    !
    return
end
!
!
!
subroutine cw261(nfun, iter, fsumsq, ifail)
    implicit double precision (a-h, o-z)
    parameter (nmax = 25)
    double precision k
    dimension x(nmax), f(nmax), etas(2 * nmax), bb(0:nmax)
    common /one/ d, t, h, u, k&
            /two/ eta(nmax), c(nmax), amp(0:nmax)
    !
    pi = 4 * atan(1.0)
    !
    do i = 1, nfun - 2
        th = (i - 1) * pi / (nfun - 1)
        x(i) = 0.0
        do j = 1, nfun - 1
            x(i) = x(i) + cos(j * th) * amp(j)
        end do
    end do
    x(nfun - 1) = k
    jverb = 0
    call gaf(nfun, nfun - 1, x, f, fsumsq, jverb, iter, ifail)
    if (ifail==1) return
    k = x(nfun - 1)
    do i = 1, nfun
        etas(i) = eta(i)
    end do
    do i = 1, nfun - 1
        etas(nfun + i) = eta(nfun - i)
    end do
    call four(etas, 2 * nfun - 2, amp, bb, nfun - 1)
    amp(nfun) = 0.0
    bb(nfun) = 0.0
    return
end
!
!
!
subroutine lsfun(nfun, x, f, jac, sq, ifail)
    !
    !  Computes function errors and Jacobian
    !  x(1:nfun-1) = unknowns,  f(1:nfun) = functions
    !
    implicit double precision (a-h, o-z)
    parameter (nmax = 25)
    double precision k, jac, ks
    dimension x(nmax), f(nmax), jac(nmax, nmax)
    dimension dex(nmax, nmax), th(nmax), s(nmax), &
            chcs(nmax, nmax), chsn(nmax, nmax), shcs(nmax, nmax), shsn(nmax, nmax), &
            fu(nmax, nmax), df(nmax, nmax), a(nmax), da(nmax, nmax), &
            dax(nmax, nmax), dcx(nmax, nmax), dfx(nmax, nmax), &
            us(nmax), vs(nmax), dux(nmax, nmax), dvx(nmax, nmax)
    common /one/ d, t, h, u, k&
            /two/ eta(nmax), c(nmax), amp(0:nmax)
    !
    pi = 4 * atan(1.0d0)
    g = 9.81d0
    ifail = 0
    om = 2 * pi / t
    !
    !  Get surface elevations
    !
    do i = 1, nfun - 2
        eta(i) = x(i)
    end do
    eta(nfun) = eta(1) - h
    sm = 0.0
    do i = 2, nfun - 2
        sm = sm + eta(i)
    end do
    eta(nfun - 1) = -sm - (eta(1) + eta(nfun)) / 2
    do i = 2, nfun
        if (eta(i)>eta(i - 1) + 5 * h / nfun) then
            ifail = 1
            !        write(*,'(7f11.3)') (eta(m),m=1,nfun-1)
            return
        endif
    end do
    k = x(nfun - 1)
    !
    !  Get dex(i,j) = d(eta(i))/d(eta(j)), 1 <= j <= nfun-2
    !
    do i = 1, nfun
        dex(i, nfun - 1) = 0.0d0
        do j = 1, nfun - 2
            if (i<nfun - 1) then
                if (i==j) then
                    dex(i, j) = 1.0d0
                else
                    dex(i, j) = 0.0d0
                endif
            else if (i==nfun - 1) then
                dex(i, j) = -1.0d0
            else if (i==nfun) then
                if (j==1) then
                    dex(i, j) = 1.0d0
                else
                    dex(i, j) = 0.0d0
                endif
            endif
        end do
    end do
    !
    !  Set some useful functions
    !
    do i = 1, nfun
        th(i) = pi * (i - 1) / float(nfun - 1)
        s(i) = d + eta(i)
        ks = k * s(i)
        do n = 1, nfun
            ch = cosh(n * ks)
            sh = sinh(n * ks)
            cs = cos(n * th(i))
            sn = sin(n * th(i))
            chcs(i, n) = ch * cs
            shsn(i, n) = sh * sn
            chsn(i, n) = ch * sn
            shcs(i, n) = sh * cs
        end do
    end do
    !
    !  Get normalised stream function coefficients a(j)
    !
    do i = 1, nfun
        fu(i, 1) = k * s(i)
        df(i, 1) = k
        do j = 2, nfun
            n = j - 1
            fu(i, j) = shcs(i, n)
            df(i, j) = chcs(i, n) * n * k
        end do
    end do
    call trans2(nfun, fu, df, a, da)
    !
    !  Get dax(i,j) = d(a(i))/d(eta(j)), 1 <=j <= nfun-2
    !
    do i = 1, nfun
        do j = 1, nfun - 2
            dij = 0.0
            do n = 1, nfun
                dij = dij + da(i, n) * dex(n, j)
            end do
            dax(i, j) = dij
        end do
    end do
    !
    !   Get dax(i,nfun-1) = d(a(i))/dk
    !
    do i = 1, nfun
        sm = 0.0
        do m = 1, nfun
            sm = sm + da(i, m) * s(m) / k
        end do
        dax(i, nfun - 1) = sm
    end do
    !
    !  De-normalise to get correct reverse mean flow
    !
    a1k2 = a(1) * k**2
    r = (u * k - om) / a1k2
    do i = 1, nfun
        c(i) = a(i) * r
        do j = 1, nfun - 1
            drx = -(u * k - om) * dax(1, j) * (k / a1k2)**2
            if (j==nfun - 1) then
                drx = drx + u / a1k2 - 2 * a(1) * k * (u * k - om) / a1k2**2
            endif
            dcx(i, j) = r * dax(i, j) + a(i) * drx
        end do
    end do
    !
    !  Get surface velocities and functions
    !
    do i = 1, nfun
        su = k * c(1)
        sv = 0.0d0
        do n = 1, nfun - 1
            su = su + n * k * c(n + 1) * chcs(i, n)
            sv = sv + n * k * c(n + 1) * shsn(i, n)
        end do
        us(i) = su
        vs(i) = sv
    end do
    sm = 0.0d0
    do i = 1, nfun
        f(i) = (us(i)**2 + vs(i)**2) / (2 * g) + eta(i)
        sm = sm + f(i)
    end do
    sm = sm / float(nfun)
    sq = 0.0d0
    do i = 1, nfun
        f(i) = f(i) - sm
        sq = sq + f(i)**2
    end do
    sq = sqrt(sq / nfun) / h
    !
    !  Get d(u(i))/d(eta(j)) and d(u(i))/dk
    !
    do i = 1, nfun
        do j = 1, nfun - 2
            su = k * dcx(1, j)
            sv = 0.0d0
            do n = 1, nfun - 1
                su = su + n * k * dcx(n + 1, j) * chcs(i, n) + (n * k)**2 * c(n + 1) * shcs(i, n) * dex(i, j)
                sv = sv + n * k * dcx(n + 1, j) * shsn(i, n) + (n * k)**2 * c(n + 1) * chsn(i, n) * dex(i, j)
            end do
            dux(i, j) = su
            dvx(i, j) = sv
        end do
        j = nfun - 1
        su = c(1) + k * dcx(1, j)
        sv = 0.0d0
        do n = 1, nfun - 1
            su = su + (c(n + 1) + k * dcx(n + 1, j)) * n * chcs(i, n) + &
                    n**2 * k * s(i) * c(n + 1) * shcs(i, n)
            sv = sv + (c(n + 1) + k * dcx(n + 1, j)) * n * shsn(i, n) + &
                    n**2 * k * s(i) * c(n + 1) * chsn(i, n)
        end do
        dux(i, j) = su
        dvx(i, j) = sv
    end do
    !
    !  Get derivatives of functions
    !
    do i = 1, nfun
        do j = 1, nfun - 1
            dfx(i, j) = (us(i) * dux(i, j) + vs(i) * dvx(i, j)) / g + dex(i, j)
        end do
    end do
    !
    do j = 1, nfun - 1
        sm = 0.0
        do i = 1, nfun
            sm = sm + dfx(i, j)
        end do
        sm = sm / nfun
        do i = 1, nfun
            jac(i, j) = dfx(i, j) - sm
        end do
    end do
    !
    return
end
!
!
!
subroutine trans2(n, f, df, a, da)
    !
    !  Gramm-Schmidt orthogonalisation with derivatives
    !
    implicit double precision (a-h, o-z)
    parameter (nmax = 25)
    dimension f(nmax, nmax), df(nmax, nmax), a(nmax), da(nmax, nmax)
    dimension c(nmax, nmax), b(nmax), cr(nmax, nmax), g(nmax, nmax), &
            dly(nmax), dc(nmax, nmax), br(nmax, nmax), db(nmax)
    do i = 1, n
        do j = 1, n
            c(i, j) = 0.0d0
        end do
    end do
    !
    st = 0.0d0
    do i = 1, n
        st = st + f(i, 1)**2
    end do
    st = dsqrt(st)
    !
    do i = 1, n
        g(i, 1) = f(i, 1) / st
    end do
    c(1, 1) = st
    do k = 2, n
        km = k - 1
        do j = 1, km
            st = 0.0d0
            do i = 1, n
                st = st + g(i, j) * f(i, k)
            end do
            c(j, k) = st
        end do
        c(k, k) = 1.0d0
        do i = 1, n
            st = 0.0d0
            do j = 1, km
                st = st + c(j, k) * g(i, j)
            end do
            g(i, k) = f(i, k) - st
        end do
        st = 0.0d0
        do i = 1, n
            st = st + g(i, k)**2
        end do
        st = dsqrt(st)
        do i = 1, n
            g(i, k) = g(i, k) / st
        end do
        c(k, k) = st
        do j = 1, km
            st = 0.0d0
            do i = 1, n
                st = st + g(i, j) * g(i, k)
            end do
            db(j) = st
        end do
        do i = 1, n
            st = 0.0d0
            do j = 1, km
                st = st + db(j) * g(i, j)
            end do
            g(i, k) = g(i, k) - st
        end do
    end do
    do j = 1, n
        st = 0.0d0
        do i = 1, n
            st = st + g(i, j)
        end do
        b(j) = st
    end do
    st = 0.0d0
    do i = 1, n
        sb = 0.0d0
        do j = 1, n
            sb = sb + b(j) * g(i, j)
        end do
        st = st + (sb - 1.0d0)**2
    end do
    call trinv(n, c, cr)
    do i = 1, n
        st = 0.0d0
        do j = 1, n
            st = st + cr(i, j) * b(j)
        end do
        a(i) = st
    end do
    st = 0.0d0
    do i = 1, n
        sb = 0.0d0
        do j = 1, n
            sb = sb + a(j) * f(i, j)
        end do
        st = st + (sb - 1.0d0)**2
    end do
    st = dsqrt(st / float(n))
    do k = 1, n
        do i = 1, n
            do j = 1, n
                dc(i, j) = 0.0d0
            end do
        end do
        dly(1) = 2.0d0 * f(k, 1) * df(k, 1)
        dc(1, 1) = 0.5d0 * dly(1) / c(1, 1)
        do j = 2, n
            dc(1, j) = (f(k, j) * df(k, 1) + df(k, j) * f(k, 1)) / c(1, 1)&
                    - 0.5d0 * dly(1) * c(1, j) / c(1, 1)**2
        end do
        do l = 2, n
            lm = l - 1
            st = 0.0d0
            do m = 1, lm
                st = st + c(m, l) * dc(m, l)
            end do
            dly(l) = 2.0d0 * (f(k, l) * df(k, l) - st)
            dc(l, l) = 0.5d0 * dly(l) / c(l, l)
            if (l==n) goto 84
            lp = l + 1
            do m = lp, n
                st = 0.0d0
                do j = 1, lm
                    st = st + c(j, m) * dc(j, l) + dc(j, m) * c(j, l)
                end do
                dc(l, m) = (f(k, m) * df(k, l) + df(k, m) * f(k, l)) / c(l, l)&
                        - st / c(l, l) - 0.5d0 * dly(l) * c(l, m) / c(l, l)**2
            end do
            84 continue
        end do
        do l = 1, n
            do m = 1, n
                st = 0.0d0
                do j = 1, n
                    st = st + g(l, j) * dc(j, m)
                end do
                br(l, m) = -st
            end do
        end do
        do m = 1, n
            br(k, m) = br(k, m) + df(k, m)
        end do
        do j = 1, n
            sb = 0.0d0
            do i = 1, n
                st = 0.0d0
                do l = 1, n
                    st = st + br(i, l) * cr(l, j)
                end do
                sb = sb + st
            end do
            db(j) = sb
        end do
        do i = 1, n
            st = db(i)
            do j = 1, n
                st = st - dc(i, j) * a(j)
            end do
            dly(i) = st
        end do
        do i = 1, n
            st = 0.0d0
            do j = 1, n
                st = st + cr(i, j) * dly(j)
            end do
            da(i, k) = st
        end do
    end do
    return
end
!
!
!
subroutine trinv(n, c, w)
    !
    !  3-diagonal system solver
    !
    implicit double precision (a-h, o-z)
    parameter (nmax = 25)
    dimension c(nmax, nmax), w(nmax, nmax)
    do i = 1, n
        do j = 1, n
            w(i, j) = 0.0d0
        end do
    end do
    do jj = 1, n
        j = n + 1 - jj
        w(j, j) = 1.0d0
        if (j==1) goto 12
        k = j - 1
        do ii = 1, k
            i = k - ii + 1
            sum = 0.0d0
            do l = 1, ii
                m = n + 2 - jj - l
                sum = sum + c(i, m) * w(m, j) / c(m, m)
            end do
            w(i, j) = -sum
        end do
        12 continue
    end do
    do i = 1, n
        sum = c(i, i)
        do j = 1, n
            w(i, j) = w(i, j) / sum
        end do
    end do
    return
end
!
!
!
subroutine wavecel(ta, d, u, tr, c)
    !
    !  Linear theory C by series approximation for waves on a current u
    !    d  = still water depth
    !    ta = absolute period (fixed reference frame)
    !    tr = relative period (reference frame moving with the current)
    !    c  = celerity (reference frame moving with the current)
    !
    implicit double precision (a-h, o-z)
    pi = 4 * atan(1.0d0)
    g = 9.81
    sigma = 2 * pi / ta
    y = sigma * sigma * d / g
    a = 1.0 / (1.0 + y * (0.6667 + y * (0.3556 + y * (0.1608 + y * (0.06321 + &
            y * (0.02744 + y * 0.01))))))
    c = dsqrt((d * g) / (y + a))
    if (abs(u)<1.0d-6) then
        tr = ta
        return
    else
        el = c * ta
        iter = 0
        10     tr = el / (el / ta - u)
        elp = (g * tr**2 / (2 * pi)) * tanh(2 * pi * d / el)
        del = elp - el
        el = el + del * 0.4
        if (abs(del / el)>1.0d-6) then
            iter = iter + 1
            if (iter==100) then
                write(*, '(a)') ' WAVECEL error'
                stop
            endif
            goto 10
        endif
        tr = el / (el / ta - u)
        c = el / tr
        return
    endif
end
!
!
!
subroutine slpds(a, b, n, cc)
    !
    !  Linear system solution
    !
    implicit double precision (a-h, o-z)
    parameter (nmax = 25)
    dimension a(nmax, nmax), b(nmax), cc(nmax)
    n1 = n - 1
    do k = 1, n1
        c = a(k, k)
        k1 = k + 1
        if (abs(c)<1.0d-10) then
            write(*, '(a,i5)') ' Matrix error 1: ', k
            stop
        endif
        do j = k1, n
            a(k, j) = a(k, j) / c
        end do
        b(k) = b(k) / c
        do i = k1, n
            c = a(i, k)
            do j = k1, n
                a(i, j) = a(i, j) - c * a(k, j)
            end do
            b(i) = b(i) - c * b(k)
        end do
    end do
    if (abs(a(n, n))<1.0d-10) then
        write(*, '(a,i5)') ' Matrix error 2: ', n
        stop
    endif
    b(n) = b(n) / a(n, n)
    do l = 1, n1
        k = n - l
        k1 = k + 1
        do j = k1, n
            b(k) = b(k) - a(k, j) * b(j)
        end do
    end do
    do i = 1, n
        cc(i) = b(i)
    end do
    return
end
!
!
!
subroutine gaf(nf, nv, xc, fvecc, fsumsq, nverb, iter, ifail)
    !
    !   Non-linear system error minimisation
    !   ifail =  0: OK
    !            1: looks hopeless
    !            2: poor convergence
    !
    implicit double precision (a-h, o-z)
    parameter (nmax = 25)
    dimension xc(nmax), fvecc(nmax)
    dimension fjacc(nmax, nmax), &
            aa(nmax, nmax), bb(nmax), cc(nmax), xcm(nmax)
    !
    iter = 0
    icalls = 0
    orf = 2.0 / nf
    fsumsm = 100.0
    !
    15 do i = 1, nv
        xcm(i) = xc(i)
    end do
    19 call lsfun(nf, xc, fvecc, fjacc, fsumsq, ifl)
    icalls = icalls + 1
    if ((fsumsq>fsumsm.and.iter>1).or.ifl/=0) then
        if (orf<0.05) then
            ifail = 1
            return
        endif
        orf = 0.8 * orf
        iter = max(iter - 1, 0)
        do i = 1, nv
            xc(i) = xcm(i)
        end do
        fsumsm = 100.0
        goto 19
    endif
    call monit(nf, fvecc, icalls, nverb)
    !
    do i = 1, nv
        do j = 1, nv
            ra = 0.0d0
            do l = 1, nf
                ra = ra + fjacc(l, i) * fjacc(l, j)
            end do
            aa(i, j) = ra
        end do
    end do
    !
    do i = 1, nv
        ra = 0.0d0
        do l = 1, nf
            ra = ra + fvecc(l) * fjacc(l, i)
        end do
        bb(i) = -ra
    end do
    !
    call slpds(aa, bb, nv, cc)
    !
    dxmax = 0.0d0
    do i = 1, nv
        dxmax = max(dxmax, abs(cc(i)))
        xc(i) = xc(i) + orf * cc(i)
    end do
    !
    iter = iter + 1
    fsumsm = fsumsq
    orf = min(1.0d0, orf * 1.1)
    !
    if (iter>=50.and.fsumsq<1.0d-4) then
        ifail = 2
        return
    endif
    if (iter>=50) then
        ifail = 1
        return
    endif
    if (fsumsq>1.0d-6) goto 15
    ifail = 0
    return
end
!
!
!
subroutine monit(nfun, f, icalls, nverb)
    implicit double precision (a-h, o-z)
    double precision k
    parameter (nmax = 25)
    dimension f(nmax)
    common /one/ d, t, h, u, k&
            /two/ eta(nmax), c(nmax), amp(0:nmax)
    !
    !  Outputs monitoring data
    !
    if (nverb==0) return
    sm = 0.0d0
    do i = 1, nfun
        sm = sm + f(i)**2
    end do
    sq = sqrt(sm / nfun) / h
    write(*, '(/i11,1p1e11.3)') icalls, sq
    write(*, '(1p7e11.3)') (eta(i), i = 1, nfun)
    return
end
!
!
!
subroutine limit(h, d, t, rat, nverb)
    double precision h, d, t, rat
    !
    !  Estimates H/(limiting height for d and t) = rat
    !
    dimension dl0(18), hl0(18)
    data dl0/2, 0.578, 0.440, 0.356, 0.293, 0.243, 0.201, 0.166, 0.1359, &
            0.1100, 0.0876, 0.0686, 0.0524, 0.0390, 0.0277, 0.01879, &
            0.01168, 0.00638/
    data hl0/0.1682, 0.1665, 0.1613, 0.1531, 0.1423, 0.1298, 0.1159, &
            0.1017, 0.0873, 0.0735, 0.0605, 0.0487, 0.0380, 0.0289, &
            0.0208, 0.01440, 0.00911, 0.00501/
    pi = 4 * atan(1.0)
    el0 = 9.81 * t**2 / (2 * pi)
    ha = h / el0
    da = d / el0
    if (da>dl0(1)) then
        rat = ha / hl0(1)
    else if (da<dl0(18)) then
        rat = ha / (0.8 * da)
    else
        do i = 2, 18
            if (dl0(i)<da) goto 11
        end do
        11    x1 = log(dl0(i))
        x2 = log(dl0(i - 1))
        y1 = log(hl0(i))
        y2 = log(hl0(i - 1))
        r = (log(da) - x1) / (x2 - x1)
        hb = exp(y1 + r * (y2 - y1))
        rat = ha / hb
    endif
    !
    !      if (nverb.ne.0.or.rat.gt.1.0) then
    !         write(*,'(a,f5.3)') ' H/Hb              = ',rat
    !      endif
    if (rat>1.0) stop
    !
    return
end
!
!
!
subroutine four(f, n, a, b, nb)
    !
    !  Fourier analysis
    !
    implicit double precision (a-h, o-z)
    dimension a(0:nb), b(0:nb), f(n)
    pi = 4 * atan(1.0d0)
    rn = 2.0d0 / n
    t = 2 * pi / n
    c = cos(t)
    s = sin(t)
    vk = 0.0d0
    vl = -1.0d0
    do k = 0, nb
        t = c * vk
        ck = t - vl
        vl = vk
        vk = ck + t
        t = ck + ck
        ul = 0.0d0
        um = f(n)
        do mm = 3, n
            m = n + 2 - mm
            u0 = ul
            ul = um
            um = t * ul - u0 + f(m)
        end do
        a(k) = (ck * um - ul + f(1)) * rn
        b(k) = s * vl * um * rn
    end do
    a(0) = a(0) * 0.5d0
    if (2 * nb/=n) return
    a(nb) = a(nb) * 0.5d0
    b(nb) = 0.0d0
    return
end
!
!
!
subroutine kmts(nfun, xx, yy, tt, uu, vv, ut, vt, du, dv, etah)
    !
    !   Computes
    !        horizontal and vertical velocity components (u,v)
    !        horizontal and vertical local acceleration components (ut,vt)
    !        horizontal and vertical total acceleration components (du,dv)
    !        water surface elevation (etah)
    !   at t=tt; x=xx; y=yy.
    !   If yy>eta kinematics are returned at the free surface
    !
    implicit double precision (a-h, o-z)
    real xx, yy, tt, uu, vv, ut, vt, du, dv, etah
    parameter (nmax = 25)
    double precision k, ks
    common /one/ d, t, h, u, k&
            /two/ eta(nmax), c(nmax), amp(0:nmax)
    pi = 4 * atan(1.0d0)
    om = 2 * pi / t
    theta = k * xx - om * tt
    !
    etah = 0.0
    do i = 1, nfun - 1
        etah = etah + cos(i * theta) * amp(i)
    end do
    !
    ks = k * (d + min(yy, etah))
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    s4 = 0.0
    !
    do i = 1, nfun - 1
        ip = i + 1
        ch = cosh(i * ks)
        sh = sinh(i * ks)
        cs = cos(i * theta)
        sn = sin(i * theta)
        s1 = s1 + i * ch * cs * c(ip)
        s2 = s2 + i * sh * sn * c(ip)
        s3 = s3 + i * i * ch * sn * c(ip)
        s4 = s4 + i * i * sh * cs * c(ip)
    end do
    uu = u + k * s1
    vv = k * s2
    ut = k * om * s3
    vt = -k * om * s4
    ux = -k * k * s3
    vx = k * k * s4
    uy = vx
    vy = -ux
    du = ut + uu * ux + vv * uy
    dv = vt + uu * vx + vv * vy
    !
    return
end
