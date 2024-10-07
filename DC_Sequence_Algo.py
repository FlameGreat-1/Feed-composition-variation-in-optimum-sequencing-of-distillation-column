import math
import numpy as np

def feedbp(kk, x, xf, fd, ncomp, spheat, at, ap, thi, kmax, tol):
    
  xfl = [xf[j] for j in range(ncomp)]
  xfv = [0.0] * ncomp

  at1 = at[:, 0]                     # Antoine equation constants
  at2 = at[:, 1] 
  at3 = at[:, 2] 
  at4 = at[:, 3] 
  at5 = at[:, 4] 
  at6 = at[:, 5] 

  ap1 = ap[:, 0]                       # Poynting correction constants
  ap2 = ap[:, 1] 
  ap3 = ap[:, 2] 
  ap4 = ap[:, 3] 
  aps = ap[:, 4] 

  # Iterate to find bubble point temperature
  t = thi
  k = 0
  while k < kmax:
    dfoft = 0.0
    sum_j = 0.0
    for j in range(ncomp):
      ek = math.exp(at1[j] / (t ** 2) + at6[j] + 
                        ap1[j] * math.log(x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0) + 
                        ap3[j] / (x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0))
      sum_j += xfl[j] * ek
      dfoft += xfl[j] * (-2 * at1[j] * ek / (t ** 3))

    foft = sum_j - 1.0
    tfcalc = t - foft / dfoft
    if abs(tfcalc - t) < tol:
      break

    t = tfcalc
    if t < 0.0:
      t = thi - 100.0
    k += 1

  # Calculate vapor mole fractions and enthalpy
    
  tfb = tfcalc
  sumxfl = 1.0
  sumxfv = 0.0
  for l in range(ncomp):
    ek = math.exp(at1[l] / (tfb ** 2) + at6[l] + 
                      ap1[l] * math.log(x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0) + 
                      ap3[l] / (x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0))
    xfv[l] = xfl[l] * ek
    sumxfv += xfv[l]

  for i in range(ncomp):
    xfv[i] /= sumxfv

  hf = 0.0
  for i in range(ncomp):
    f = fd * xf[i]
    hfl = f * spheat[i] * (tfb * 5.0 / 9.0 - 298.0)
    hf += hfl

  return hf, tfb

def ghfug(qf, kk, x, fd, xf, xdhk, xblk, ncomp, lk, hk, t, at, ap):

  alpha = relvol(kk, x, ncomp, t, hk, at, ap)
  dg = np.zeros((2, 2))               
  f = np.zeros(ncomp) 
  # Set recovery fractions
  recdlk = 0.998
  recbhk = 0.995

  # Calculate feed flow rates for each component
  for i in range(ncomp):
    f[i] = fd * xf[i]

  # Perform Newton-Raphson iteration to find distillate and bottoms compositions
  kmax = 50
  delta = 0.00005
  for k in range(kmax):
    xd, xb, d, b8, a2 = hg(kk, x, recdlk, recbhk, ncomp, lk, hk, xf, fd, at, ap, t)  # Calculate distillate and bottoms compositions using Hengstebeck-Geddes method
    func = [xd[hk] - xdhk, xb[lk] - xblk]

    # Estimate partial derivatives
    s = [recdlk, recbhk]
    g = func.copy()
    for i in range(2):
      s[i] += delta
      xd, xb, d, b8, a2 = hg(kk, recdlk, recbhk, ncomp, lk, hk, xf, fd, at, ap, t)
      func = [xd[hk] - xdhk, xb[lk] - xblk]
      for j in range(2):
        dg[j][i] = (func[j] - g[j]) / delta
      s[i] -= delta

    # Update recovery fractions
    dm = dg[0][0] * dg[1][1] - dg[0][1] * dg[1][0]
    dels1 = (g[1] * dg[0][1] - g[0] * dg[1][1]) / dm
    dels2 = (g[0] * dg[1][0] - g[1] * dg[0][0]) / dm
    s[0] += dels1
    s[1] += dels2
    if s[0] > 1.0:
      s[0] = 0.9999
    if s[1] > 1.0:
      s[1] = 0.9999

  # Calculate minimum number of plates using Fenske's equation
  nminw = a2 
  nmin = fenske(nminw, d[lk], d[hk], b8[lk], b8[hk], alpha[lk], kk) 

  # Calculate minimum reflux ratio and Underwood parameter using Underwood's equation
  phi, rmin = undwd(qf, kk, x,  xd, xf, t, hk,  ncomp, at, ap)

  # Find number of theoretical plates at specified reflux ratio using Gilliland's correlation
  nt, rf = glld(x, rmin, nmin, kk)

  # Find feed plate location using Kirkbride's equation
  m = kirbde(sum(b8), sum(d), xd[hk], xb[lk], xf[lk], xf[hk], nt)

  return nmin, phi, rmin, nt, m, xd, xb, d, b8, rf

def hg(kk, x, recdlk, recbhk, ncomp, lk, hk, xf, fd, at, ap, t):

  f = [0.0] * ncomp
  d = [0.0] * ncomp
  b8 = [0.0] * ncomp
  xd = [0.0] * ncomp
  xb = [0.0] * ncomp
    
  alpha = relvol(kk, x, ncomp, t, hk, at, ap)
  # Calculate flow rates of light and heavy key components in distillate and bottoms
  for i in range(ncomp):
    f[i] = xf[i] * fd
      
  d[lk] = recdlk * f[lk]
  b8[lk] = f[lk] - d[lk]
  d[hk] = (1.0 - recbhk) * f[hk]
  b8[hk] = f[hk] - d[hk]

  # Calculate Hengstebeck-Geddes correlation constant A1 and A2
  a1 = math.log(d[hk] / b8[hk])
  a2 = math.log(d[lk] / d[hk] * b8[hk] / b8[lk]) / math.log(alpha[lk])
  # Calculate flow rates and mole fractions for other components
  for i in range(ncomp):
    if i == lk or i == hk:
      continue
    b8[i] = f[i] / (1.0 + math.exp(a1 + a2 * math.log(alpha[i])))
    d[i] = f[i] - b8[i]

  # Calculate total flow rates and mole fractions
  sumd = sum(d)
  sumb = sum(b8)
    
  for i in range(ncomp):
    xd[i] = d[i] / sumd
    xb[i] = b8[i] / sumb

  return xd, xb, d, b8, a2

def fenske(nmin1, dlk, dhk, blk, bhk, alfalk, kk):
  nmin = math.log((dlk / dhk) * (bhk / blk)) / math.log(alfalk)
  if nmin1 == nmin:
    print(True)
  else:
    print(False)
  return nmin

def undwd(qf, kk, x,  xd, xf, t, hk,  ncomp, at, ap):

  alpha = relvol(kk, x, ncomp, t, hk, at, ap)
  # Set initial limits for PHI
  phill = alpha[hk]
  phiul = alpha[hk - 1]

  # Find PHI using bisection method
  for i in range(ncomp):
    phi = 0.5 * (phill + phiul)
    fphi = qf - 1.0
    for j in range(ncomp):
      fphi += alpha[j] * xf[j] / (alpha[j] - phi)
    if fphi < 0:
      phill = phi
    else:
      phiul = phi

  # Calculate minimum reflux ratio (RMIN)
  rmin = 1.0
  for i in range(ncomp):
    rmin += alpha[i] * xd[i] / (alpha[i] - phi)

  return phi, rmin

def glld(x, rmin, nmin, kk):

  # Calculate actual reflux ratio and Gilliland correlation parameters
  rf = x[kk] * rmin
  xxx = (rf - rmin) / (rf + 1.0)

  # yy = 1.0 - math.exp((1.0 + 54.4 + xxx) * (xxx - 1.0) / ((11.0 + 117.2 + xxx) * math.sqrt(xxx)))  
  yy = 0.75 * (1.0 - (xxx ** 0.5668))  # Alternative calculation
  
  # Calculate number of theoretical trays
  nt = (yy + nmin) / (1.0 - yy)
  print("nt", nt)

  return nt, rf

def kirbde(sumb, sumd, xdhk, xblk, xflk, xfhk, nt):

  s = 0.206 * math.log((sumb / sumd) * (xfhk / xflk) * ((xblk / xdhk) ** 2))
  m = (nt * (10 ** s)) / (1.0 + (10 ** s))

  return m

def relvol(kk, x, ncomp, t, hk, at, ap):

  at1 = at[:, 0]                     # Antoine equation constants
  at2 = at[:, 1] 
  at3 = at[:, 2] 
  at4 = at[:, 3] 
  at5 = at[:, 4] 
  at6 = at[:, 5] 

  ap1 = ap[:, 0]                       # Poynting correction constants
  ap2 = ap[:, 1] 
  ap3 = ap[:, 2] 
  ap4 = ap[:, 3] 
  aps = ap[:, 4] 

  alpha = [0.0] * ncomp
  # Calculate equilibrium constants for each component
  eqk = [0.0] * ncomp
  for i in range(ncomp):
    eqk[i] = math.exp(at1[i] / (t ** 2) + at6[i] + 
                      ap1[i] * math.log(x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0) + 
                      ap3[i] / (x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0))

  # Calculate relative volatilities with respect to the heavy key
  for ii in range(ncomp):
    alpha[ii] = eqk[ii] / eqk[hk]

  return alpha

def distil(kk, x, xd, spheat, d, hvap, sumd, cv, at, ap, tol, tlo, kmax, ncomp, thi):

  cv1 = cv[:, 0] 
  cv2 = cv[:, 1] 
  cv3 = cv[:, 2] 
  cv4 = cv[:, 3] 

    
  xl = [0.0] * 15                       # Liquid mole fractions
  xv = [0.0] * 15                       # Vapor mole fractions

  # Calculate dew point temperature
  xv = xd.copy()
  td = dewpt(kk, x, tol, tlo, xv, kmax, at, ap, ncomp)
  tod = td

  # Calculate bubble point temperature
  xl = xv.copy()
  tb = bubbpt(kk, x, xl, thi, tol, ncomp, kmax, at, ap)
  tob = tb

  # Calculate enthalpy of distillate stream
  hd = 0.0
  for i in range(ncomp):
    hdl = d[i] * spheat[i] * (tob * 5.0 / 9.0 - 298.0)
    hd += hdl

  # Calculate enthalpy of top vapor stream
  hdv = 0.0
  for i in range(ncomp):
    hdt = (hvap[i] + cv1[i] * ((tod * 5.0 / 9.0) - 298.0) + 
           cv2[i] / 2.0 * (((tod * 5.0 / 9.0) ** 2) - (298.0 ** 2)) +
           cv3[i] / 3.0 * (((tod * 5.0 / 9.0) ** 3) - (298.0 ** 3)) +
           cv4[i] / 4.0 * (((tod * 5.0 / 9.0) ** 4) - (298.0 ** 4)))
    hdv += xd[i] * sumd * hdt

  return hd, hdv, tod, tob

def bottms(kk, x, xb, spheat, b8, thi, tol, ncomp, kmax, at, ap ):

  # Calculate bubble point temperature
  xl = xb.copy()
  tb = bubbpt(kk, x, xl, thi, tol, ncomp, kmax, at, ap)
  tbb = tb

  hb = 0.0
  for j in range(ncomp):
    hbl = b8[j] * spheat[j] * ((tbb * 5.0 / 9.0) - 298.0)
    hb += hbl 
      
  return hb, tbb

def dewpt(kk, x, tol, tlo, xv, kmax, at, ap, ncomp):
    
  at1 = at[:, 0]                     # Antoine equation constants
  at2 = at[:, 1] 
  at3 = at[:, 2] 
  at4 = at[:, 3] 
  at5 = at[:, 4] 
  at6 = at[:, 5] 

  ap1 = ap[:, 0]                       # Poynting correction constants
  ap2 = ap[:, 1] 
  ap3 = ap[:, 2] 
  ap4 = ap[:, 3] 
  aps = ap[:, 4] 

  xl = [0.0] * ncomp
    
  # Initialize temperature and iteration counter
  t = 479.0
  k = 0

  # Iterate to find dew point temperature
  while k < kmax:
    dfoft = 0.0
    sum_k = 0.0

    # Calculate K values and summations
    for ki in range(ncomp):
      ek = math.exp(at1[ki] / (t ** 2) + at6[ki] + ap1[ki] * math.log(x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0) + ap3[ki] / (x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0))
      sum_k += xv[ki] / ek
      dfoft -= (xv[ki] / (ek ** 2)) * (-2 * at1[ki] * ek / (t ** 3))

    foft = sum_k - 1.0
    tcalc = t - foft / dfoft

    # Check for convergence
    if abs(tcalc - t) < tol:
      break

    t = tcalc

    # Temperature limit check
    if t < 0.0:
      t = tlo + 120.0

    k += 1

  # Calculate liquid mole fractions
  dwptt = tcalc
  td = dwptt
  sumxv = 1.0
  sumx = 0.0

  for i in range(ncomp):
    ek = math.exp(at1[i] / (td ** 2) + at6[i] + ap1[i] * math.log(x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0) + ap3[i] / (x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0))
    xl[i] = xv[i] / ek
    sumx += xl[i]

  for l in range(ncomp):
    xl[l] /= sumx

  return td

def bubbpt(kk, x, xl, thi, tol, ncomp, kmax, at, ap):
    
  at1 = at[:, 0]                     # Antoine equation constants
  at2 = at[:, 1] 
  at3 = at[:, 2] 
  at4 = at[:, 3] 
  at5 = at[:, 4] 
  at6 = at[:, 5] 

  ap1 = ap[:, 0]                       # Poynting correction constants
  ap2 = ap[:, 1] 
  ap3 = ap[:, 2] 
  ap4 = ap[:, 3] 
  aps = ap[:, 4] 

  xv = [0.0] * ncomp

  # Initialize temperature and iteration counter
  t = thi
  k = 0

  # Iterate to find bubble point temperature
  while k < kmax:
    dfoft = 0.0
    sum_j = 0.0

    # Calculate K values and summations
    for j in range(ncomp):
      ek = math.exp(at1[j] / (t ** 2) + at6[j] + ap1[j] * math.log(x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0) + ap3[j] / (x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0))
      sum_j += xl[j] * ek
      dfoft += xl[j] * (-2 * at1[j] * ek / (t ** 3.0))

    foft = sum_j - 1.0
    tcalc = t - foft / dfoft

    # Check for convergence
    if abs(tcalc - t) < tol:
      break

    t = tcalc

    # Temperature limit check
    if t < 0.0:
      t = thi - 100.0

    k += 1

  # Calculate vapor mole fractions
  bubptt = tcalc
  tb = bubptt
  sumxl = 1.0
  sumy = 0.0

  for l in range(ncomp):
    ek = math.exp(at1[l] / (tb ** 2) + at6[l] + ap1[l] * math.log(x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0) + ap3[l] / (x[kk + 1 + 3 * (ncomp - 1)] * 14.7 / 760.0))
    xv[l] = xl[l] * ek
    sumy += xv[l]

  for i in range(ncomp):
    xv[i] /= sumy

  return tb

def costfn(x, fd, xf, ncomp, twin, twrise, ucond, htcwat, plife, ureb, ucswat, trayef, spheat, hvap, cv, at, ap, thi, kmax, tol, tlo):
    # Initialize cost variables
    sumac = 0.0
    sumsstm = 0.0
    sumwat = 0.0
    sumcap = 0.0
    sumop = 0.0
    
    for kk in range(ncomp - 1):
        print(kk, "iteration")
        if kk == 0:
            fd = fd  # Assume 'fd' is the total feed flow rate
            lk = 0
            hk = 1
            xdhk = 0.0010
            xblk = 0.0001
            x[kk + 1 + 3 + (ncomp - 1)] = 4851.0  # Set pressure for the first column
        else: 
            lk = 0
            hk = 1
            # Set recovery fractions and pressure for subsequent columns (based on your specific problem)
            if kk == 1:
                xdhk = 0.0004
                xblk = 0.0001
                x[kk + 1 + 3 * (ncomp - 1)] = 4169.0  
            elif kk == 2:
                xblk = 0.0001
                xdhk = 0.0008
                x[kk + 1 + 3 * (ncomp - 1)] = 3659.0 
            elif kk == 3:
                xdhk = 0.0008
                xblk = 0.0010
                x[kk + 1 + 3 * (ncomp - 1)] = 2639.0 
    
        qf = 0.8  # Assume saturated liquid feed
        
        hf, tfb = feedbp(kk, x, xf, fd, ncomp, spheat, at, ap, thi, kmax, tol)
        nmin, phi, rmin, nt, m, xd, xb, d, b8, rf = ghfug(qf, kk, x, fd, xf, xdhk, xblk, ncomp, lk, hk, tfb, at, ap)
        hd, hdv, tod, tob = distil(kk, x, xd, spheat, d, hvap, sum(d), cv, at, ap, tol, tlo, kmax, ncomp, thi)
        hb, tbb = bottms(kk, x, xb, spheat, b8, thi, tol, ncomp, kmax, at, ap)

        fd = sum(b8[1:])
        xf = xb[1:]
        ncomp = ncomp - 1
        
        # Condenser calculations
        twout = twin + twrise
        dutytp = (rf + 1.0) * abs(hdv - hd)
        acond = dutytp / (ucond * x[kk + 2 * (ncomp - 1)])
        
        # Reboiler calculations
        tsteam = (tbb * 5.0 / 9.0) + 25.0
        stmvht = 0.521294e+05 - 0.142381e+02 * tsteam - 0.445721e-01 * (tsteam ** 2)
        dutybt = hd + hb + dutytp - hf
        areb = dutybt / (ureb * x[kk + ncomp - 1])
        
        # Cost of utilities
        watcst = (8500.0 * ucswat * dutytp) / (twrise * htcwat)
        sumwat += watcst
        
        ucsstm = 0.708252e+01 + 0.228804e-01 * tsteam
        stmcost = (8500.0 * ucsstm * dutybt) / stmvht 
        sumsstm += stmcost
        
        # Column design
        hc = 0.5092 * (nt / trayef) + 3.048
        vapr = 0.761 + ((1.0 / (x[kk + 1 + 3 * (ncomp - 1)] / 760.00)) ** 0.5)
        dd = ((28.0 / (22.0 * vapr)) * (rf + 1.0) * 22.2 * ((np.abs(tod) * 5.0 / 9.0 / 273.0) * (1.0 / (x[kk + 1 + 3 * (ncomp - 1)] / 760.0)) * (1.0 / 3600.0)) ** 0.5)
        dc = dd * (sum(d) ** 0.5)
        
        print("dd, vapr, rf, tod", dd, vapr, rf, tod)
        fprime = (0.778 - 0.000082 * hc) * math.log(3.281 * dc) + 0.9199 * math.sqrt(hc) - 1.433
        print("hc, dc", hc, dc)
        wcstin = math.exp((1.33 * fprime) - 0.541) * 1.0e+03 
        print("fprime", fprime)
        
        # Pressure correction factor
        if x[kk + 1 + 3 * (ncomp - 1)] <= 3040.0:
            cp = 1.05
        elif x[kk + 1 + 3 * (ncomp - 1)] <= 6080.0:
            cp = 1.05
        elif x[kk + 1 + 3 * (ncomp - 1)] <= 11400.0:
            cp = 1.15 
        elif x[kk + 1 + 3 * (ncomp - 1)] <= 15200.0:
            cp = 1.20
        elif x[kk + 1 + 3 * (ncomp - 1)] <= 30400.0:
            cp = 1.60 
        else:
            cp = 2.50 
        
        fm = 1.0  # Material type factor (carbon steel)
        colcst = wcstin * (cp * fm) 
        print("wcstin, cp, fm", wcstin, cp, fm)
        
        ft = 1.8  # Tray type factor (bubble cap)
        fm = 0.0  # Material type factor (carbon steel)
        
        trcost = (1.0 + ft + fm) * nt * (0.030 + 0.038 * (dc ** 2.0)) * 1000.0 
        cincst = 8000.0  # Instrumentation cost
        cmncst = 0.02 * (trcost + colcst)  # Maintenance cost of column
        
        # Pressure correction factor for heat exchanger
        if x[kk + 1 + 3 * (ncomp - 1)] <= 7600.0:
            cphx = 0.0 
        elif x[kk + 1 + 3 * (ncomp - 1)] <= 15200.0:
            cphx = 0.10 
        elif x[kk + 1 + 3 * (ncomp - 1)] <= 30400.0:
            cphx = 0.40 
        else:
            cphx = 0.55
        
        ft = 0.35  # Exchanger type factor (kettle)
        fm = 0.0  # Material type factor (carbon steel)
        
        concst = (cphx + 1.0 + fm + ft) * (0.73 + 0.30 * (acond ** 0.65)) * 1000.0 
        rebcost = (cphx + 1.0 + fm + ft) * (0.73 + 0.30 * (areb ** 0.65)) * 1000.0 
        hxcost = concst + rebcost 
        
        opcost = stmcost + watcst + 0.02 * hxcost 
        sumop += opcost 
        capcst = (colcst + trcost + cmncst + cincst + hxcost) / plife 
        print(capcst)
        tacst = opcost + (colcst + trcost + cmncst + cincst + hxcost) / plife 
        
        print("________________________++++++++here we go +++++++++___________________")
        print(opcost, colcst, trcost, cmncst, cincst, hxcost, plife)
        print("________________________++++++++here we go+++++++++___________________")
        
        tacost = tacst 
        sumac += tacst 
        sumcap += capcst
    
        # Print final design information
        print("+" * 19 + " THE FINAL DESIGN " + "+" * 19)
        print(f"FOR KK= {kk:12}   DIAMETER = {dc:9.4f} M")
        print(f"    HEIGHT= {hc:9.4f} M   NUMBER OF PLATES= {nt:10.4f}")

        print("+" * 28 + " COST INFORMATION " + "+" * 28)
        print(f"UCSSTM= {ucsstm}")
        print(f"UCSWAT= {ucswat}")
        print(f"PLIFE= {plife}")
    
        pcstm = (stmcost / tacost) * 100.0
        pcwat = (watcst * 100.0 / tacost)
        pccap = (capcst * 100.0 / tacost)
        print(f"kk= {kk:12}   PCSTM= {pcstm:6.1f} PER CENT")
        print(f"    PCWAT= {pcwat:6.1f} PER CENT")
        print(f"    PCCAP= {pccap:6.1f} PER CENT")
        print(f"FOR kk= {kk:12}   STEAM COST= {stmcost:13.4f}") 
        print(f"    COOLING WATER COST= {watcst:12.4f}   EMPTY COLUMN COST= {colcst:14.4f}")
        print(f"FOR kk= {kk:12}   TRAY COST= {trcost:12.4f}")
        print(f"    MAINTENANCE COST= {cmncst:12.4f},'CAPCST=', {capcst:13.4f}")
        print(f"FOR kk= {kk:12}   CONDENSER COST= {concst:14.4f}")
        print(f"    REBOILER COST= {rebcost:14.4f}") 
        print(f"FOR kk= {kk:12}   HEAT EXCHANGER COST= {hxcost:13.4f}") 
        print(f"    OPERATING COST= {opcost:14.4f}")
        print(f"    TOTAL ANNUAL COST= {tacost:14.4f}") 
        print(f"FOR kk= {kk:12}   CONDENSER DUTY= {dutytp:18.3f}")
        print(f"    REBOILER DUTY= {dutybt:18.3f}") 
        print(f"    CONDENSER AREA= {acond:16.3f}")
        print(f"    REBOILER AREA= {areb:16.3f}")
        
        # Print overall cost percentages
        stmpc = (sumsstm * 100.0) / sumac
        watpc = (sumwat * 100.0) / sumac
        opcpc = (sumop * 100.0) / sumac 
        cappc = (sumcap * 100.0) / sumac 
        print(f"STMPC= {stmpc:6.1f} PER CENT")
        print(f"WATPC= {watpc:6.1f} PER CENT")
        print(f"OPCPC= {opcpc:6.1f} PER CENT")
        print(f"CAPPC= {cappc:6.1f} PER CENT")

    return sumac, sumcap, sumsstm, sumwat, sumop

# Constants and data
at = np.array([
    [-970688.5625, 0.0, 0.0, 0.0, 0.0, 7.15059],
    [-1166846.0, 0.0, 0.0, 0.0, 0.0, 7.72668],
    [-1280557.0, 0.0, 0.0, 0.0, 0.0, 7.94986],
    [-1481583.0, 0.0, 0.0, 0.0, 0.0, 7.58071],
    [-1524891.0, 0.0, 0.0, 0.0, 0.0, 7.33129]
])

ap = np.array([
    [-0.76984, 0.0, 6.90244, 0.0, 0.0],
    [-0.92213, 0.0, 0.0, 0.0, 0.0],
    [-0.96455, 0.0, 0.0, 0.0, 0.0],
    [-0.93159, 0.0, 0.0, 0.0, 0.0],
    [-0.89143, 0.0, 0.0, 0.0, 0.0]
])



cv = np.array([
    [-4.0444, 3.0480E-01, -1.5721E-04, 3.1736E-08],
    [-7.9131, 4.1617E-01, -2.3006E-04, 4.9907E-08],
    [3.9565, 3.7149E-01, -1.8338E-04, 3.5002E-08],
    [-9.5166, 5.2042E-01, -2.9714E-04, 6.6403E-08],
    [6.7742, 4.5427E-01, -2.2462E-04, 4.2287E-08]
])

xf = [0.05, 0.15, 0.25, 0.20, 0.35]
spheat = [106.3447, 133.5589, 133.5589, 159.0984, 159.0984]
hvap = [15100.0, 18800.0, 21100.0, 24600.0, 26500.0]
ncomp = 5
ncomp, fd, tlo, thi, kmax = (5, 907.2, 460.0, 760.0, 20)
tol, plife, ureb, ucond, twrise, twin = [0.00001, 10.0, 2093.40, 2512.08, 20.0, 305.0]
x = np.arange(7, 101)
htcwat, ucswat, ucsstm, trayef = (75.3624, 0.000286, 0.02967, 0.8)


if __name__ == "__main__":
    sumac, sumcap, sumsstm, sumwat, sumop = costfn(x, fd, xf, ncomp, twin, twrise, ucond, htcwat, plife, ureb, ucswat, trayef, spheat, hvap, cv, at, ap, thi, kmax, tol, tlo)
    print("sumac, sumcap, sumsstm, sumwat, sumop", sumac, sumcap, sumsstm, sumwat, sumop)

