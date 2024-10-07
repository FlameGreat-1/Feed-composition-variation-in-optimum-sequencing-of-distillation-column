import sys
import math
import numpy as np
from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QTextEdit, QScrollArea, QTabWidget, QTableWidget, QTableWidgetItem, QHeaderView, QSplitter)
from PyQt5.QtGui import QFont, QPalette, QColor
from PyQt5.QtCore import Qt


# Include all your original functions here
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
      xd, xb, d, b8, a2 = hg(kk, x, recdlk, recbhk, ncomp, lk, hk, xf, fd, at, ap, t)
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


def undwd(qf, kk, x, xd, xf, t, hk, ncomp, at, ap):
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

def bottms(kk, x, xb, spheat, b8, thi, tol, ncomp, kmax, at, ap):
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
    sumac = sumsstm = sumwat = sumcap = sumop = 0.0
    iteration_results = []
    
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

        tacost = tacst 
        sumac += tacst 
        sumcap += capcst


        pcstm = (stmcost / tacost) * 100.0
        pcwat = (watcst * 100.0 / tacost)
        pccap = (capcst * 100.0 / tacost)


        # Store detailed results for each iteration
        iteration_result = {
            'kk': kk,
            'diameter': dc,
            'height': hc,
            'num_plates': nt,
            'ucsstm': ucsstm,
            'ucswat': ucswat,
            'plife': plife,
            'pcstm': pcstm,
            'pcwat': pcwat,
            'pccap': pccap,
            'steam_cost': stmcost,
            'cooling_water_cost': watcst,
            'empty_column_cost': colcst,
            'tray_cost': trcost,
            'maintenance_cost': cmncst,
            'capcst': capcst,
            'condenser_cost': concst,
            'reboiler_cost': rebcost,
            'heat_exchanger_cost': hxcost,
            'operating_cost': opcost,
            'total_annual_cost': tacost,
            'condenser_duty': dutytp,
            'reboiler_duty': dutybt,
            'condenser_area': acond,
            'reboiler_area': areb,
            'dd': dd,
            'vapr': vapr,
            'rf': rf,
            'tod': tod,
            'tob': tob,
            'fprime': fprime,
            'wcstin': wcstin,
            'cp': cp,
            'fm': fm,

        }
        iteration_results.append(iteration_result)

    return sumac, sumcap, sumsstm, sumwat, sumop, iteration_results



                    # GRAPHICAL USER INTERFACE (GUI) SECTION


class DistillationGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Feed Composition Variation In Optimum Sequencing Of Distillation Column")
        self.setGeometry(100, 100, 1200, 800)

        self.setStyleSheet("""
            QMainWindow {
                background-color: #f0f0f0;
            }
            QLabel {
                font-size: 14px;
                font-weight: bold;
            }
            QLineEdit, QPushButton {
                font-size: 13px;
                padding: 5px;
            }
            QPushButton {
                background-color: #4CAF50;
                color: white;
                border: none;
                border-radius: 4px;
            }
            QPushButton:hover {
                background-color: #45a049;
            }
            QTableWidget {
                gridline-color: #d0d0d0;
            }
            QHeaderView::section {
                background-color: #e0e0e0;
                font-weight: bold;
                border: none;
                border-right: 1px solid #d0d0d0;
                border-bottom: 1px solid #d0d0d0;
            }
        """)

        # Main widget and layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout()
        main_widget.setLayout(main_layout)

        # Custom title label
        title_label = QLabel("FEED COMPOSITION VARIATION IN OPTIMUM SEQUENCING OF DISTILLATION COLUMN")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("""
            font-size: 18px;
            font-weight: bold;
            color: black;
            padding: 10px;
            background-color: #f0f0f0;
            border-bottom: 1px solid #d0d0d0;
        """)

        # Title label to the main layout
        main_layout.addWidget(title_label)


        # Input section
        input_widget = QWidget()
        input_layout = QHBoxLayout()
        input_widget.setLayout(input_layout)

        self.fd_input = QLineEdit()
        self.fd_input.setPlaceholderText("Enter feed flow rate (e.g., 907.2)")
        self.xf_input = QLineEdit()
        self.xf_input.setPlaceholderText("Enter feed composition (e.g., 0.05,0.15,0.25,0.20,0.35)")
        calculate_button = QPushButton("Compute")
        calculate_button.clicked.connect(self.run_calculation)

        input_layout.addWidget(QLabel("Feed flow rate (fd):"))
        input_layout.addWidget(self.fd_input)
        input_layout.addWidget(QLabel("Feed composition (xf):"))
        input_layout.addWidget(self.xf_input)
        input_layout.addWidget(calculate_button)

        main_layout.addWidget(input_widget)

        # Output section
        self.output_tabs = QTabWidget()
        main_layout.addWidget(self.output_tabs)

        
        for i in range(6):  # 4 iterations + 1 summary + 1 additional results
            if i < 4:
                tab = self.create_iteration_tab()
                self.output_tabs.addTab(tab, f"Column {i}")
            elif i == 4:
                tab = QTextEdit()
                tab.setReadOnly(True)
                self.output_tabs.addTab(tab, "Cost Summary")
            else:
                tab = self.create_additional_results_tab()
                self.output_tabs.addTab(tab, "Distillation Column Design Information")

    def create_iteration_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        splitter = QSplitter(Qt.Vertical)

        design_table = self.create_table(["Parameter", "Value"], 3)
        cost_table = self.create_table(["Parameter", "Value"], 17)
        additional_table = self.create_table(["Parameter", "Value"], 4)

        splitter.addWidget(design_table)
        splitter.addWidget(cost_table)
        splitter.addWidget(additional_table)

        layout.addWidget(splitter)
        return tab

    def create_additional_results_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        tab.setLayout(layout)

        table = self.create_table(["Iteration", "Parameter", "Value 1", "Value 2", "Value 3"], 24)
        layout.addWidget(table)

        return tab

    def create_table(self, headers, rows):
        table = QTableWidget(rows, len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        table.verticalHeader().setVisible(False)
        return table

    def populate_table(self, table, data):
        for row, (key, value) in enumerate(data):
            table.setItem(row, 0, QTableWidgetItem(key))
            table.setItem(row, 1, QTableWidgetItem(value))

    def run_calculation(self):
        try:
            fd = float(self.fd_input.text())
            xf = [float(x) for x in self.xf_input.text().split(',')]

            # Constants and data (you may need to adjust these)
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

            spheat = [106.3447, 133.5589, 133.5589, 159.0984, 159.0984]
            hvap = [15100.0, 18800.0, 21100.0, 24600.0, 26500.0]
            ncomp = 5
            tlo, thi, kmax = 460.0, 760.0, 20
            tol, plife, ureb, ucond, twrise, twin = 0.00001, 10.0, 2093.40, 2512.08, 20.0, 305.0
            x = np.arange(7, 101)
            htcwat, ucswat, trayef = 75.3624, 0.000286, 0.8

            # Run the calculation
            sumac, sumcap, sumsstm, sumwat, sumop, iteration_results = costfn(
                x, fd, xf, ncomp, twin, twrise, ucond, htcwat, plife, ureb, ucswat, trayef, spheat, hvap, cv, at, ap, thi, kmax, tol, tlo)

            # Update the GUI with the results
            self.update_output(sumac, sumcap, sumsstm, sumwat, sumop, iteration_results)
        except ValueError as e:
            # Handle invalid input
            error_message = QMessageBox()
            error_message.setIcon(QMessageBox.Critical)
            error_message.setText("Invalid input")
            error_message.setInformativeText(str(e))
            error_message.setWindowTitle("Error")
            error_message.exec_()
    

    def update_output(self, sumac, sumcap, sumsstm, sumwat, sumop, iteration_results):
       for i, result in enumerate(iteration_results):
           if i >= 4:  # We only have 4 iteration tabs
              break
           tab = self.output_tabs.widget(i)
           design_table, cost_table, additional_table = tab.findChildren(QTableWidget)
           
           design_data = [
            ("DIAMETER", f"{result['diameter']:.4f} M"),
            ("HEIGHT", f"{result['height']:.4f} M"),
            ("NUMBER OF PLATES", f"{result['num_plates']:.4f}")
        ]
           self.populate_table(design_table, design_data)
           
           cost_data = [
            ("UCSSTM", f"{result['ucsstm']:.6f}"),
            ("UCSWAT", f"{result['ucswat']:.6f}"),
            ("PLIFE", f"{result['plife']:.1f}"),
            ("PCSTM", f"{result['pcstm']:.1f} PER CENT"),
            ("PCWAT", f"{result['pcwat']:.1f} PER CENT"),
            ("PCCAP", f"{result['pccap']:.1f} PER CENT"),
            ("STEAM COST", f"{result['steam_cost']:.4f}"),
            ("COOLING WATER COST", f"{result['cooling_water_cost']:.4f}"),
            ("EMPTY COLUMN COST", f"{result['empty_column_cost']:.4f}"),
            ("TRAY COST", f"{result['tray_cost']:.4f}"),
            ("MAINTENANCE COST", f"{result['maintenance_cost']:.4f}"),
            ("CAPCST", f"{result['capcst']:.4f}"),
            ("CONDENSER COST", f"{result['condenser_cost']:.4f}"),
            ("REBOILER COST", f"{result['reboiler_cost']:.4f}"),
            ("HEAT EXCHANGER COST", f"{result['heat_exchanger_cost']:.4f}"),
            ("OPERATING COST", f"{result['operating_cost']:.4f}"),
            ("TOTAL ANNUAL COST", f"{result['total_annual_cost']:.4f}")
        ]
           
           self.populate_table(cost_table, cost_data)
           
           additional_data = [
            ("CONDENSER DUTY", f"{result['condenser_duty']:.3f}"),
            ("REBOILER DUTY", f"{result['reboiler_duty']:.3f}"),
            ("CONDENSER AREA", f"{result['condenser_area']:.3f}"),
            ("REBOILER AREA", f"{result['reboiler_area']:.3f}")
        ]
           
           self.populate_table(additional_table, additional_data)
           
           
           summary_table = self.create_table(["Parameter", "Value"], 5)
           summary_data = [
              ("Total Annual Cost (sumac)", f"{sumac:.4e}"),
              ("Capital Cost (sumcap)", f"{sumcap:.4e}"),
              ("Steam Cost (sumsstm)", f"{sumsstm:.4e}"),
              ("Water Cost (sumwat)", f"{sumwat:.4e}"),
              ("Operating Cost (sumop)", f"{sumop:.4e}")
           ]
           
           self.populate_table(summary_table, summary_data)
              
           summary_tab = self.output_tabs.widget(4)
           if isinstance(summary_tab, QTextEdit):
            self.output_tabs.removeTab(4)
            self.output_tabs.insertTab(4, summary_table, " Cost Summary")


          
           
           additional_results_table = self.output_tabs.widget(5).findChild(QTableWidget)
           row = 0
           for i, result in enumerate(iteration_results):
              additional_results_table.setItem(row, 0, QTableWidgetItem(f"{i}"))
              additional_results_table.setItem(row, 1, QTableWidgetItem("nt"))
              additional_results_table.setItem(row, 2, QTableWidgetItem(f"{result['num_plates']:.4f}"))
              row += 1

              additional_results_table.setItem(row, 0, QTableWidgetItem(f"{i}"))
              additional_results_table.setItem(row, 1, QTableWidgetItem("dd, vapr, rf, tod"))
              additional_results_table.setItem(row, 2, QTableWidgetItem(f"{result['dd']:.4f}"))
              additional_results_table.setItem(row, 3, QTableWidgetItem(f"{result['vapr']:.4f}"))
              additional_results_table.setItem(row, 4, QTableWidgetItem(f"{result['rf']:.4f}"))
              row += 1

              additional_results_table.setItem(row, 0, QTableWidgetItem(f"{i}"))
              additional_results_table.setItem(row, 1, QTableWidgetItem("tod"))
              additional_results_table.setItem(row, 2, QTableWidgetItem(f"{result['tod']:.4f}"))
              row += 1

              additional_results_table.setItem(row, 0, QTableWidgetItem(f"{i}"))
              additional_results_table.setItem(row, 1, QTableWidgetItem("hc, dc"))
              additional_results_table.setItem(row, 2, QTableWidgetItem(f"{result['height']:.4f}"))
              additional_results_table.setItem(row, 3, QTableWidgetItem(f"{result['diameter']:.4f}"))
              row += 1

              additional_results_table.setItem(row, 0, QTableWidgetItem(f"{i}"))
              additional_results_table.setItem(row, 1, QTableWidgetItem("fprime"))
              additional_results_table.setItem(row, 2, QTableWidgetItem(f"{result['fprime']:.4f}"))
              row += 1

              additional_results_table.setItem(row, 0, QTableWidgetItem(f"{i}"))
              additional_results_table.setItem(row, 1, QTableWidgetItem("wcstin, cp, fm"))
              additional_results_table.setItem(row, 2, QTableWidgetItem(f"{result['wcstin']:.4f}"))
              additional_results_table.setItem(row, 3, QTableWidgetItem(f"{result['cp']:.4f}"))
              additional_results_table.setItem(row, 4, QTableWidgetItem(f"{result['fm']:.4f}"))
              row += 1
    

def populate_table(self, table, data):
    for row, (key, value) in enumerate(data):
        table.setItem(row, 0, QTableWidgetItem(key))
        table.setItem(row, 1, QTableWidgetItem(value))


if __name__ == '__main__':
    app = QApplication(sys.argv)
    gui = DistillationGUI()
    gui.show()
    sys.exit(app.exec_())


    