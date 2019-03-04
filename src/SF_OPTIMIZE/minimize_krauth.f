C========+=========+=========+=========+=========+=========+=========+=$
C COMMENT: This is a most reliable conjugent gradient routine! It has
C          served us well for many years, and is capable to cope with
C          a very large number of variables. Unfortunately, we don't
C          know who wrote this routine (original name: 'va10a'), and 
C          we find it very obscure.
C          Don't worry, it works just fine.
Cnoprint=+=========+=========+=========+=========+=========+=========+=$
      subroutine minimize_krauth_(funct, n, x, f, g, h, w, dfn, xm,
     $  hh, eps, mode, maxfn, iprint, iexit,itn)
      implicit double precision (a-h,o-z)
      dimension x(*), g(*), h(*), w(*), xm(*)
      external funct
      data zero, half, one, two /0.0d0, 0.5d0, 1.0d0, 2.0d0/
      if (iprint .ne. 0) write (6,1000)
 1000 format (' entry into minimize')
      np = n + 1
      n1 = n - 1
      nn=(n*np)/2
      is = n
      iu = n
      iv = n + n
      ib = iv + n
      idiff = 1
      iexit = 0
      if (mode .eq. 3) go to 15
      if (mode .eq. 2) go to 10
      ij = nn + 1
      do 5 i = 1, n
      do 6 j = 1, i
      ij = ij - 1
   6  h(ij) = zero
   5  h(ij) = one
      go to 15
  10  continue
      ij = 1
      do 11 i = 2, n
      z = h(ij)
      if (z .le. zero) return
      ij = ij + 1
      i1 = ij
      do 11 j = i, n
      zz = h(ij)
      h(ij) = h(ij) / z
      jk = ij
      ik = i1
      do 12 k = i, j
      jk = jk + np - k
      h(jk) = h(jk) - h(ik) * zz
      ik = ik + 1
  12  continue
      ij = ij + 1
  11  continue
      if (h(ij) .le. zero) return
  15  continue
      ij = np
      dmin = h(1)
      do 16 i = 2, n
      if (h(ij) .ge. dmin) go to 16
      dmin = h(ij)
  16  ij = ij + np - i
      if (dmin .le. zero) return
      z = f
      itn = 0
      call funct (n, x, f)
      ifn = 1
      df = dfn
      if (dfn .eq. zero) df = f - z
      if (dfn .lt. zero) df = abs (df * f)
      if (df .le. zero) df = one
  17  continue
      do 19 i = 1, n
      w(i) = x(i)
  19  continue
      link = 1
      if (idiff - 1) 100, 100, 110
  18  continue
      if (ifn .ge. maxfn) go to 90
  20  continue
      if (iprint .eq. 0) go to 21
      if (mod (itn, iprint) .ne. 0) go to 21
       write (6,1001) itn, ifn
1001  format (1x,'itn = ',i5,' ifn = ',i5)
      write (6,1002) f
1002  format (1x,'f = ',e15.7)
      if (iprint .lt. 0) go to 21
      write (6,1003) (x(i), i = 1, n)
***
***
1003  format (1x,'x = ',4e15.7 / (5x, 4e15.7))
      write (6,1004) (g(i), i = 1, n)
1004  format (1x,'g = ',4e15.7 / (5x, 4e15.7))
**
***
  21  continue
      itn = itn + 1
      w(1) = -g(1)
      do 22 i = 2, n
      ij = i
      i1 = i - 1
      z = -g(i)
      do 23 j = 1, i1
      z = z - h(ij) * w(j)
      ij = ij + n - j
  23  continue
  22  w(i) = z
      w(is+n) = w(n) / h(nn)
      ij = nn
      do 25 i = 1, n1
      ij = ij - 1
      z = zero
      do 26 j = 1, i
      z = z + h(ij) * w(is+np-j)
      ij = ij - 1
  26  continue
  25  w(is+n-i) = w(n-i) / h(ij) - z
      z = zero
      gs0 = zero
      do 29 i = 1, n
      if (z * xm(i) .ge. abs (w(is+i))) go to 28
      z = abs (w(is+i)) / xm(i)
  28  gs0 = gs0 + g(i) * w(is+i)
  29  continue
      aeps = eps / z
      iexit = 2
      if (gs0 .ge. zero) go to 92
      alpha = -two * df / gs0
      if (alpha .gt. one) alpha = one
      ff = f
      tot = zero
      int = 0
      iexit = 1
  30  continue
      if (ifn .ge. maxfn) go to 90
      do 31 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  31  continue
      call funct (n, w, f1)
      ifn = ifn + 1
      if (f1 .ge. f) go to 40
      f2 = f
      tot = tot + alpha
  32  continue
      do 33 i = 1, n
      x(i) = w(i)
  33  continue
      f = f1
      if (int - 1) 35, 49, 50
  35  continue
      if (ifn .ge. maxfn) go to 90
      do 34 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  34  continue
      call funct (n, w, f1)
      ifn = ifn + 1
      if (f1 .ge. f) go to 50
      if ((f1 + f2 .ge. f + f) .and.
     $  (7.0d0 * f1 + 5.0d0 * f2 .gt. 12.0d0 * f)) int = 2
      tot = tot + alpha
      alpha = two * alpha
      go to 32
  40  continue
      if (alpha .lt. aeps) go to 92
      if (ifn .ge. maxfn) go to 90
      alpha = half * alpha
      do 41 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  41  continue
      call funct (n, w, f2)
      ifn = ifn + 1
      if (f2 .ge. f) go to 45
      tot = tot + alpha
      f = f2
      do 42 i = 1, n
      x(i) = w(i)
  42  continue
      go to 49
  45  continue
      z = 0.1d0
      if (f1 + f .gt. f2 + f2)
     $  z = one + half * (f - f1) / (f + f1 - f2 - f2)
      if (z .lt. 0.1d0) z = 0.1d0
      alpha = z * alpha
      int = 1
      go to 30
  49  continue
      if (tot .lt. aeps) go to 92
  50  continue
      alpha = tot
      do 56 i = 1, n
      w(i) = x(i)
      w(ib+i) = g(i)
  56  continue
      link = 2
      if (idiff - 1) 100, 100, 110
  54  continue
      if (ifn .ge. maxfn) go to 90
      gys = zero
      do 55 i = 1, n
      w(i) = w(ib+i)
      gys = gys + g(i) * w(is+i)
  55  continue
      df = ff - f
      dgs = gys - gs0
      if (dgs .le. zero) go to 20
      link = 1
      if (dgs + alpha * gs0 .gt. zero) go to 52
      do 51 i = 1, n
      w(iu + i) = g(i) - w(i)
  51  continue
      sig = one / (alpha * dgs)
      go to 70
  52  continue
      zz = alpha / (dgs - alpha * gs0)
      z = dgs * zz - one
      do 53 i = 1, n
      w(iu+i) = z * w(i) + g(i)
  53  continue
      sig = one / (zz * dgs * dgs)
      go to 70
  60  continue
      link = 2
      do 61 i = 1, n
      w(iu+i) = w(i)
  61  continue
      if (dgs + alpha * gs0 .gt. zero) go to 62
      sig = one / gs0
      go to 70
  62  continue
      sig = -zz
  70  continue
      w(iv+1) = w(iu+1)
      do 71 i = 2, n
      ij = i
      i1 = i - 1
      z = w(iu+i)
      do 72 j = 1, i1
      z = z - h(ij) * w(iv+j)
      ij = ij + n - j
  72  continue
      w(iv+i) = z
  71  continue
      ij = 1
      do 75 i = 1, n
      z = h(ij) + sig * w(iv+i) * w(iv+i)
      if (z .le. zero) z = dmin
      if (z .lt. dmin) dmin = z
      h(ij) = z
      w(ib+i) = w(iv+i) * sig / z
      sig = sig - w(ib+i) * w(ib+i) * z
      ij = ij + np - i
  75  continue
      ij = 1
      do 80 i = 1, n1
      ij = ij + 1
      i1 = i + 1
      do 80 j = i1, n
      w(iu+j) = w(iu+j) - h(ij) * w(iv+i)
      h(ij) = h(ij) + w(ib+i) * w(iu+j)
      ij = ij + 1
  80  continue
      go to (60, 20), link
  90  continue
      iexit = 3
      go to 94
  92  continue
      if (idiff .eq. 2) go to 94
      idiff = 2
      go to 17
  94  continue
      if (iprint .eq. 0) return
      write (6,1005) itn, ifn, iexit
1005  format (1x,'itn = ',i5, ' ifn = ',i5,' iexit = ',i5)
      write (6,1002) f
      write (6,1003) (x(i), i = 1, n)
      write (6,1004) (g(i), i = 1, n)
      return
 100  continue
      do 101 i = 1, n
      z = hh * xm(i)
      w(i) = w(i) + z
      call funct (n, w, f1)
      g(i) = (f1 - f) / z
      w(i) = w(i) - z
 101  continue
      ifn = ifn + n
      go to (18, 54), link
 110  continue
      do 111 i = 1, n
      z = hh * xm(i)
      w(i) = w(i) + z
      call funct (n, w, f1)
      w(i) = w(i) - z - z
      call funct (n, w, f2)
      g(i) = (f1 - f2) / (two * z)
      w(i) = w(i) + z
 111  continue
      ifn = ifn + n + n
      go to (18, 54), link
      end 
