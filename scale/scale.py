# -*- coding: utf-8 -*-
import ROOT, time, os, sys, ConfigParser, cPickle
import numpy as np
from atlas import Atlas
from scipy.interpolate import UnivariateSpline

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html#scipy.interpolate.UnivariateSpline
class BSpline:  # class # class # class # class # class # class # class # class # class # class # class # class # class
  def Reset(self,V,E,dE):
    self.V, self.E, self.dE, self.pars = V, E, dE, []
  def __call__(self, x, p):
    if [p[0],p[1]] != self.pars:
      self.pars, X, Y, W = [p[0],p[1]], [], [], []
      for i in range(len(self.V)):
        e = p[0] + p[1]*self.V[i]
        X.append(e);  Y.append(self.E[i]-e); W.append(1./self.dE[i])
      self.bspline = UnivariateSpline(X, Y, W)
    return self.bspline(x[0])


class Scale(Atlas): # class # class # class # class # class # class # class # class # class # class # class # class # class
  Colors = [ROOT.kRed+2, ROOT.kGreen+3, ROOT.kBlue+2]
  Styles = [20,  22,  24];    Sizes  = [1.25, 1.0, 1.0]
  emin, emax = 100., 10000.  # just for a case when no *.cfg file is present

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def __init__(self, outfile, cfg_file, action):
    self.outfile = outfile
    Atlas.__init__(self)
    self.InitParameters(cfg_file)
    if action == 'calibration':
      self.cv   =  ROOT.TCanvas('cv','HPGe calibration', 2, 2, 1002, 1002)
      self.PADS = [ROOT.TPad('PAD5', 'Scale',      0.01, 0.01, 0.49, 0.36, 0, 1),
                   ROOT.TPad('PAD6', 'Resolution', 0.51, 0.01, 0.99, 0.36, 0, 1),
                   ROOT.TPad('PAD0', 'Spectrum',   0.01, 0.66, 0.99, 0.99, 0, 1),
                   ROOT.TPad('PAD1', 'Fit 1',      0.01, 0.37, 0.25, 0.65, 0, 1),
                   ROOT.TPad('PAD2', 'Fit 2',      0.25, 0.37, 0.50, 0.65, 0, 1),
                   ROOT.TPad('PAD3', 'Fit 3',      0.50, 0.37, 0.75, 0.65, 0, 1),
                   ROOT.TPad('PAD4', 'Fit 4',      0.75, 0.37, 0.99, 0.65, 0, 1)]
      for pad in range(7):
        self.cv.cd(); self.PADS[pad].Draw()
        if pad > 2:   self.PADS[pad].SetLogy()

      self.SP = ROOT.TSpectrum()
      # Peaks Fitting Section
      self.asps = ROOT.TF1('asps', LineShape(), 0.0, 1.0, 7);    self.asps.SetLineColor(ROOT.kBlue+2)
      self.asps.SetParNames('Amp','E_{0}, keV', '#sigma_{R}, keV', '#sigma_{L}/#sigma_{R}, keV', 'Compton', 'Background', '#kappa')
      self.p_limits = ((1.0, 1.e+7), (50.0, 1.e+4), (0.5, 10.00),  (1.0, 5.0), (0.0,100.0), (0.0, 1.e+6), (1.0, 10.0) )
    else:
      self.cv   =  ROOT.TCanvas('cv','HPGe calibration', 2, 2, 1002, 502)
      self.PADS = [ROOT.TPad('PAD5', 'Scale',      0.01, 0.01, 0.49, 0.99, 0, 1),
                   ROOT.TPad('PAD6', 'Resolution', 0.51, 0.01, 0.99, 0.99, 0, 1)]
      for pad in self.PADS:  self.cv.cd(); pad.Draw()
    # Energy Scale Section
    self.NLG      = ROOT.TMultiGraph(); self.NLG.SetTitle('E_{FIT} - E_{REF} [keV]')
    self.lisc     = ROOT.TF1('lisc', '[0] + [1]*x',                self.emin, self.emax) # linear scale calibration
    self.bspline  = BSpline()
    self.SplineB  = ROOT.TF1('SplineB', self.bspline, self.emin, self.emax, 2)
    self.SplineB.SetLineColor(self.Colors[1]); self.SplineB.SetLineWidth(2); self.SplineB.SetNpx(250);
    # Energy Resolution Section
    self.ERG      = ROOT.TMultiGraph(); self.ERG.SetTitle('#sigma_{E} / E [%]')
    self.eres     = ROOT.TF1('eres', ResolutionModel(), -self.emax, self.emax, 3) # combined resolution model [%]
    self.eres.SetParameters(1.0,  0.24, 0.1e-6, 100.0); self.eres.SetLineColor(self.Colors[0])
    self.eres.SetParLimits(0, 0.5, 10.0); self.eres.SetParLimits(1, 0.05, 0.5); self.eres.SetParLimits(2, 0.05, 0.5)

    self.pulr     = ROOT.TF1('pulr', '100.0*([0]/abs(x)+[1])',  -self.emax, self.emax) # pulser resolution [%]
    self.pulr.SetParameters(1.0, 0.0); self.pulr.SetLineColor(self.Colors[1])
    self.InitGraphics()
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def __del__(self):
    Atlas.__del__(self);
    try: 
      self.cv.cd(); self.cv.Clear()
    except AttributeError:
      pass  

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Do(self, utb, ute, pb5, hps, ptype):
    self.UTB, self.UTE, self.hps, self.PB5, self.ptype = utb, ute, hps, pb5, ptype
    self.nbins = self.hps.GetNbinsX()
    self.hps.GetXaxis().SetTitle(' E_{#gamma} [keV]')
    self.hps.SetBins(self.nbins, self.zero, self.zero + self.gain * self.nbins)
    self.findPeaks()
    self.identifyPeaks()
    self.ShowScale(0)
    if len(self.ScalePeaks)>2:
      for itr in range(1, self.nitr+1):
        self.fitPeaks(L = self.fitL, R = self.fitR)
        self.Show_Energy_Scale(self.Do_Energy_Scale(self.nitr == itr))
        self.hps.SetBins(self.nbins, self.zero, self.zero + self.gain * self.nbins)
      self.ShowScale(1)
      self.PeaksTable() # there p3 = SL/SR is replaced by p3 = p3*SR = SL
      self.Show_Energy_Resolution(self.Do_Energy_Resolution())
      self.Save_Calibration()
    else:
      raw_input('Scale calibration is impossible!')
    return self.zero, self.gain
      #      for i in range(7): self.PADS[i].SaveAs('PAD%1d.eps' % i)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def findPeaks(self):
    rebin, self.nFoundPeaks, self.FoundPeaks = 2, 0, []
    self.PADS[2].cd(); self.PADS[2].Clear(); self.PADS[2].SetGrid(); self.PADS[2].SetLogy()
    H = self.hps.Rebin(rebin,"H")
    B = self.SP.Background(H)
    self.SP.Search(H, 2.0, 'nodraw', self.amin) # the parameters are sigma, options, threshold
    self.hps.DrawCopy('HIST'); self.cv.Modified(); self.cv.Update()

    N  =  self.SP.GetNPeaks()
    LX = [self.SP.GetPositionX()[n] for n in range(N)]

    for peak in range(N):
      xbin = 1 + int((LX[peak]-self.zero)/self.gain/rebin)
      YB      = B.GetBinContent(xbin)
      YC, dYC = H.GetBinContent(xbin),    H.GetBinError(xbin)
      YL, dYL = H.GetBinContent(xbin-15), H.GetBinError(xbin-15)
      YR, dYR = H.GetBinContent(xbin+15), H.GetBinError(xbin+15)
      if (YL+2*dYL)<(YC-2*dYC)>(YR+2*dYR):
        self.nFoundPeaks+=1;  self.FoundPeaks.append({'E':LX[peak], 'A':(YC-YB)/rebin, 'B':YB/rebin})
    H.Delete();   B.Delete()
    del H, B, LX
#    for n in range(self.nFoundPeaks):  print self.FoundPeaks[n]


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def identifyPeaks(self):
    self.ScalePeaks, self.OtherPeaks = [], []
#    self.PB5 = []
    if self.PB5:
      self.PulsePeaks = [{'name':'P%02d %5.3f V' % (self.PB5.index(v), v), 'V':v, 'E':self.zero_p + self.gain_p*v, 'dE':self.eerr_p} for v in self.PB5]
      self.PB5 = True
    else:
      self.PulsePeaks = []
      self.PB5 = False
    unkn = 0
    while len(self.FoundPeaks):
      name, x, A, B = False, self.FoundPeaks[0]['E'], self.FoundPeaks[0]['A'], self.FoundPeaks[0]['B']
      if self.PB5:
        for line in self.PulsePeaks: # first try to find pulser peaks
          if abs(x - line['E']) < self.erec:
            name, i  = line['name'], self.PulsePeaks.index(line)
            self.PulsePeaks[i]['X'], self.PulsePeaks[i]['A'], self.PulsePeaks[i]['B'] = x, A, B
      if not name:                   # second try to find real peaks
        for nucl,lines in self.atlas.iteritems():
          for line in lines:
            if abs(x - line['W']) < self.erec:
              name = "%5s(%4s)" % (nucl, lines.index(line))
              V = {'name':name, 'E':line['W'], 'dE':line['dW'] + self.tbpa, 'X':x, 'A':A, 'B':B}
              if name in self.scale: self.ScalePeaks.append(V)
              else:                  self.OtherPeaks.append(V)
              break
            else: continue
          if name: break
      if not name and A>10.:
        print 'Unknown line: % 8.2f keV (Amp=%4d)' % (x, A)
        if x > 2000.:
          name = "Unknown(%2d)" % unkn
          V = {'name':name, 'E':x, 'dE':2*self.gain, 'X':x, 'A':A, 'B':B}
          self.OtherPeaks.append(V)
          unkn +=1

      self.FoundPeaks.pop(0)

    self.ScalePeaks.sort(key = lambda p: p['E'])
    self.OtherPeaks.sort(key = lambda p: p['E'])
    if self.PB5:
      self.PulsePeaks = [el for el in self.PulsePeaks if el.has_key('A')]
      self.PulsePeaks.sort(key = lambda p: p['E'])
    if len(self.PulsePeaks) < 4 : self.PB5=False


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def fitPeaks(self, L = 2.0, R = 2.0):
    for p in self.ScalePeaks + self.PulsePeaks + self.OtherPeaks:
      PB5, ALREADY = p['name'][-1]=='V', p.has_key('shape')
      if   ALREADY: # if peak was fitted already:   Amplitude,  Position, Sigma_R, Sigma_L/Sigma_R, Compton, Background, kappa
        V = p['shape']
        if PB5: p0,p1,p2,p3,p4,p5,p6 = V['p0'], V['p1'], V['p2'], 1.0,     100.000, V['p5'], 10.0
        else:   p0,p1,p2,p3,p4,p5,p6 = V['p0'], V['p1'], V['p2'], V['p3'], V['p4'], V['p5'], V['p6']
      else:   # if a peak was not fitted yet:   Amplitude,  Position, Sigma_R, Sigma_L/Sigma_R, Compton, Background, kappa
        p0,p1,p2,p3,p4,p5,p6  = p['A'], p['X'], 1.0, 1.0, 0.00,  p['B'], 1.5
      self.asps.SetParameters(p0,p1,p2,p3,p4,p5,p6)
      self.asps.SetRange(p1 - L*p2*p3, p1 + R*p2)
      for np in xrange(7): self.asps.SetParLimits(np, self.p_limits[np][0], self.p_limits[np][1])
      if PB5:
        self.asps.FixParameter(3, 1.0)
        self.asps.FixParameter(4, 100.0)
      elif self.exp0:
        self.asps.FixParameter(6, self.kapp.Eval(p1) )
      for attempt in xrange(3):
        Result = self.hps.Fit('asps','RSQN')
        if not Result.Status(): break # print Result.Print()
      if PB5:
        self.asps.ReleaseParameter(3)
        self.asps.ReleaseParameter(4)
      elif self.exp0:
        self.asps.ReleaseParameter(6)
      P = fitParameters(self.asps)
      A_OK = abs(P['dp0']/P['p0']) < 0.50
      E_OK = abs(P['dp1']/P['p1']) < 0.005
      S_OK = abs(P['dp2']/P['p2']) < 0.50
      OK = A_OK and E_OK and S_OK and not Result.Status()
      if OK or not ALREADY:
        p['shape'] = P
      else:
        if   p in self.ScalePeaks: self.ScalePeaks.remove(p)
        elif p in self.PulsePeaks: self.PulsePeaks.remove(p)
        elif p in self.OtherPeaks: self.OtherPeaks.remove(p)
      if   self.names[0] in p['name']: self.Show_Peak(p, L, R, 3)
      elif self.names[1] in p['name']: self.Show_Peak(p, L, R, 4)
      elif self.names[2] in p['name']: self.Show_Peak(p, L, R, 5)
      elif self.names[3] in p['name']: self.Show_Peak(p, L, R, 6)
      self.cv.Modified(); self.cv.Update()
    return


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Show_Peak(self, p, VLEVO, VPRAVO, npad):
    C, SR, SL = p['shape']['p1'], p['shape']['p2'], p['shape']['p2']*p['shape']['p3']
    L, R, N, M = C-SL*(VLEVO+1), C+SR*(VPRAVO+1), '', ''
    for c in p['name'][0:5]:
      if   c.isdigit(): N+=c
      elif c.isspace(): pass
      else:             M+=c
    name  = '^{%s}%s' % (N,M)
    self.PADS[npad].cd(); self.PADS[npad].Clear(); self.PADS[npad].SetGrid()
    H = self.hps.DrawCopy();  H.GetXaxis().SetRangeUser(L, R)
    H.SetNameTitle('Spectrum', '%5s (%.3f keV) #chi^{2}/ndf = %.1f/%d' % (name, p['E'], p['shape']['Chi2'], p['shape']['NDF']))
    self.asps.DrawCopy('SAME')
    self.cv.Modified(); self.cv.Update()
#    self.PADS[npad].SaveAs('%.0fkeV.png' % p['E'])

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Save_Peak(self, p, VLEVO, VPRAVO):
    C, SR, SL = p['shape']['p1'], p['shape']['p2'], p['shape']['p2']*p['shape']['p3']
    L, R, N, M = C-SL*(VLEVO+1), C+SR*(VPRAVO+1), '', ''
    for c in p['name'][0:5]:
      if   c.isdigit(): N+=c
      elif c.isspace(): pass
      else:             M+=c
    name  = '^{%s}%s' % (N,M)
    ct = ROOT.TCanvas('cv','HPGe calibration', 2, 2, 1002, 1002)
    ct.cd(); ct.SetGrid(); ct.SetLogy(); ct.GetYaxis.SetMoreLogLabels()
    H = self.hps.DrawCopy();  H.GetXaxis().SetRangeUser(L, R)
    H.SetNameTitle('Spectrum', '%5s (%.3f keV) #chi^{2}/ndf = %.1f/%s' % (name, p['E'], p['shape']['Chi2'], p['shape']['NDF']))
    self.asps.DrawCopy('SAME')
    ct.Modified(); ct.Update()
    ct.SaveAs('%.0fkeV.pdf' % p['E'])


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Do_Energy_Scale(self, fin):
    K = 'shape'
# 1) LINEAR SCALE CORRECTION WITH CALIBRATION GAMMA LINES
    self.S0.Set(len(self.ScalePeaks))
    for pk in self.ScalePeaks:
      n = self.ScalePeaks.index(pk)
      self.S0.SetPoint(     n, pk[ 'E'],  pk[K]['p1'] - pk['E'])
      self.S0.SetPointError(n, pk['dE'], (pk['dE']**2 + pk[K]['dp1']**2)**0.5)
    self.linear_scale_fit_result = self.S0.Fit('lisc','QRNS')
    k0,  k1 = self.lisc.GetParameter(0), self.lisc.GetParameter(1)
    self.zero = (self.zero-k0) / (1.+k1)
    self.gain =  self.gain     / (1.+k1)

# 2) JUST SHOW OTHER GAMMA LINES
#    self.S2.Set(len(self.OtherPeaks))
#    for pk in self.OtherPeaks:
#      n = self.OtherPeaks.index(pk)
#      self.S2.SetPoint(     n, pk[ 'E'], pk[K][ 'p1'] - pk['E'])
#      self.S2.SetPointError(n, pk['dE'], pk[K]['dp1'])

# 3) pulser LINES
    self.spline = True
    if self.PB5:
      V,E,dE = [],[],[]
      self.S1.Set(len(self.PulsePeaks))
      for pk in self.PulsePeaks:
        n = self.PulsePeaks.index(pk)
        dY = (pk['dE']**2 + pk[K]['dp1']**2)**0.5
        self.S1.SetPoint(n,  pk['E'], pk[K]['p1'] - pk['E'])
        self.S1.SetPointError(n, 0.0, dY)
        V.append(pk['V']);  E.append(pk[K]['p1']);  dE.append(dY)
      if fin:
        W = [V,E,dE]
        with open('pulser.table','w') as fp:  cPickle.dump(W, fp, -1)
    elif fin:
      self.S1.Set(0)
      try:
        with open('pulser.table','rb') as fp:  V,E,dE = cPickle.load(fp)
      except IOError:
        self.spline = False
        print 'pulser.table does not exist!'
    else:
        self.spline = False
    if self.spline:
      self.bspline.Reset(V,E,dE)
      self.SplineB.SetParameters(self.zero_p, self.gain_p)
      self.pulser_correction_fit_result = self.S0.Fit('SplineB','QRNS')
      if self.pulser_correction_fit_result.Prob() < 0.00001:
        self.spline = False
      else:
        x, e = np.fromiter(E, 'float64'), np.ndarray(len(E), 'float64')
        self.pulser_correction_fit_result.GetConfidenceIntervals(len(x), 1, 1, x, e, 0.683, False)
        P = fitParameters(self.SplineB)
        self.zero_p, self.gain_p = P['p0'], P['p1']
        if not fin:
          for pk in self.PulsePeaks:
            pk['E']  = P['p0'] + P['p1']*pk['V']
#        pk['dE'] = e[self.PulsePeaks.index(pk)]
    if self.spline:
      return '#chi^{2}/NDF = %5.1f/%3d' % (self.pulser_correction_fit_result.Chi2(), self.pulser_correction_fit_result.Ndf())
    else:
      return '#chi^{2}/NDF = %5.1f/%3d' % (self.linear_scale_fit_result.Chi2(),      self.linear_scale_fit_result.Ndf())

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Show_Energy_Scale(self, quality):
    xmin, xmax, ymin, ymax = MultiGraphLimits(self.NLG)
    self.PADS[0].cd();  self.PADS[0].Clear();  self.PADS[0].SetGrid()
    self.NLG.Draw('AP');
    self.NLG.GetXaxis().SetLimits(   xmin, xmax); self.NLG.GetXaxis().SetTitle(' E_{#gamma} [keV]')
    self.NLG.GetYaxis().SetRangeUser(ymin, ymax); self.NLG.GetYaxis().SetDecimals()
    self.Lg1.Clear()
    self.Lg1.AddEntry(self.S0, 'reference lines', 'lpe')
    self.Lg1.AddEntry(self.S1, 'pulser lines',    'lpe')
    self.Lg1.AddEntry(self.S2, 'other lines',     'lpe')
    self.Lg1.AddEntry(self.SplineB, quality, 'lpe')
    self.Lg1.Draw('SAME')
    if self.spline:  self.SplineB.Draw('SAME')
    self.cv.Modified();    self.cv.Update()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Do_Energy_Resolution(self):
    K = 'shape'

# 1. Resolution by pulser (noise level)
    if self.PB5:
      self.RP.Set(len(self.PulsePeaks))
      for p in self.PulsePeaks:
        n = self.PulsePeaks.index(p);
        self.RP.SetPoint(n, p['E'], 100*p[K]['p2']/p['E']); self.RP.SetPointError(n, p['dE'], 100*p[K]['dp2']/p['E'])
      self.pulser_resolution_fit_result = self.RP.Fit('pulr','QRNS')
      P = fitParameters(self.pulr); PP = (P['p0'], P['dp0'], 1e+6*P['p1'], 1e+6*P['dp1'], P['Chi2'], P['NDF'])
      print ' ╔ Resolution by pulser: ═══╤══════════════════════════════╤══════════════════════╗'
      print ' ║ Noise: %5.3f ± %5.3f keV │ Slope:  %6.3f ± %5.3f ppm   │ χ²/NDF:   %5.1f/%3d  ║' % PP
      print ' ╚══════════════════════════╧══════════════════════════════╧══════════════════════╝\n'
    else: self.RP.Set(0)

# 2. Resolution by isotopes for calibration
    N = len(self.ScalePeaks); self.RC.Set(2*N)
    for p in self.ScalePeaks:
      n = self.ScalePeaks.index(p)
      self.RC.SetPoint(N-n-1, -p['E'], 100*p[K]['p3']/p['E']); self.RC.SetPointError(N-n-1, p['dE'], 100*p[K]['dp3']/p['E'])
      self.RC.SetPoint(N+n,    p['E'], 100*p[K]['p2']/p['E']); self.RC.SetPointError(N+n,   p['dE'], 100*p[K]['dp2']/p['E'])

# 3. Resolution for other isotopes
    OP = [p for p in self.OtherPeaks if p not in self.ScalePeaks]; N = len(OP); self.RO.Set(2*N)
    for p in OP:
      n = OP.index(p);
      self.RO.SetPoint(N-n-1, -p['E'], 100*p[K]['p3']/p['E']); self.RO.SetPointError(N-n-1, p['dE'], 100*p[K]['dp3']/p['E'])
      self.RO.SetPoint(N+n,    p['E'], 100*p[K]['p2']/p['E']); self.RO.SetPointError(N+n,   p['dE'], 100*p[K]['dp2']/p['E'])

    """
    if self.rep3: self.eres.FixParameter(3, self.rep3)
    else:         self.eres.SetParLimits(3, -1500.00, 500.)
    self.real_resolution_fit_result = self.RC.Fit('eres','QRNS')
    if self.rep3: self.eres.ReleaseParameter(3)
    P = fitParameters(self.eres); quality = '#chi^{2}/NDF = %5.1f/%3d' % (P['Chi2'], P['NDF'])
    PP1 = (P['p0'], P['dp0'], 1e+6*P['p2'], 1e+6*P['dp2'], P['Chi2'], P['NDF'])
    PP2 = (P['p1'], P['dp1'],      P['p3'],      P['dp3'], self.real_resolution_fit_result.Prob())
    print ' ╔ Resolution by isotopes: ═╤══════════════════════════════╤══════════════════════╗'
    print ' ║ Noise: %5.3f ± %5.3f keV │ Trapping:  %5.3f ± %5.3f ppm │ χ²/NDF:   %5.1f/%3d  ║' % PP1
    print ' ║ Fano fact: %5.3f ± %5.3f │ Threshold: %5.0f ± %5.0f keV │ Probability:  %5.3f  ║' % PP2
    print ' ╚══════════════════════════╧══════════════════════════════╧══════════════════════╝\n'
    """
    self.real_resolution_fit_result = self.RC.Fit('eres','QRNS')
    P = fitParameters(self.eres)
    PP1 = (P['p0'], P['dp0'], P['Chi2'], P['NDF'])
    PP2 = (P['p1'], P['dp1'], P['p2'],   P['dp2'], self.real_resolution_fit_result.Prob())
    print ' ╔ Resolution by isotopes: ═╤══════════════════════════════╤══════════════════════╗'
    print ' ║ Noise: %5.3f ± %5.3f keV │ ---------------------------- │ χ²/NDF:   %5.1f/%3d  ║' % PP1
    print ' ║ Fano fact: %5.3f ± %5.3f │ Fano left: %5.3f ± %5.3f keV │ Probability:  %5.3f  ║' % PP2
    print ' ╚══════════════════════════╧══════════════════════════════╧══════════════════════╝\n'


    return '#chi^{2}/NDF = %5.1f/%3d' % (P['Chi2'], P['NDF'])

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Show_Energy_Resolution(self, quality):
    xmin, xmax, ymin, ymax = MultiGraphLimits(self.ERG)
    self.Lg2.Clear()
    self.Lg2.AddEntry(self.RC, '#sigma_{L}, #sigma_{R} reference #gamma', 'lpe')
    self.Lg2.AddEntry(self.RO, '#sigma_{L}, #sigma_{R} observed #gamma', 'lpe')
    self.Lg2.AddEntry(self.RP, '#sigma  pulser ', 'lpe')
    self.Lg2.AddEntry(self.eres, quality, 'lpe')
    self.PADS[1].cd(); self.PADS[1].Clear(); self.PADS[1].SetGrid()
    self.ERG.Draw('AP'); self.eres.DrawCopy('SAME'); self.Lg2.Draw('SAME')
    if self.PB5:         self.pulr.DrawCopy('SAME')
    self.ERG.GetXaxis().SetLimits( -xmax, xmax); self.ERG.GetXaxis().SetTitle(' E_{#gamma} [keV]');
    self.ERG.GetYaxis().SetRangeUser(0.0, ymax); self.ERG.GetYaxis().SetDecimals()
    self.cv.Modified(); self.cv.Update()
    return


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def PeaksTable(self):
    K='shape'
    for p in self.ScalePeaks + self.OtherPeaks + self.PulsePeaks:
      p[K]['p3' ] = p[K]['p3']*p[K]['p2']
      p[K]['dp3'] = ((p[K]['p2']*p[K]['dp3'])**2 + (p[K]['p3']*p[K]['dp2'])**2)**0.5
    if self.PB5:
      print ' ╔ pulser: ══╤═════════╤═════════╤════════╤════════╤═════════╤════════╤═══════════╗'
      print ' ║ Line Name │  Counts │ Eγ,keV  │ σR,keV │ σL,keV │ Compton │ Backgr │   χ²/NDF  ║'
      print ' ╟───────────┼─────────┼─────────┼────────┼────────┼─────────┼────────┼───────────╢'
      for p in self.PulsePeaks:
        PP = (p['name'], p[K]['p0']/self.gain, p[K]['p1'], p[K]['p2'], p[K]['p3'], p[K]['p5'], p[K]['Chi2'], p[K]['NDF'])
        print ' ║%11s│ %7.0f │ %7.2f │ %6.4f │ %6.4f │   ---   │ %6.0f │ %5.1f/%3d ║' % PP
      print ' ╚═══════════╧═════════╧═════════╧════════╧════════╧═════════╧════════╧═══════════╝\n'
    print ' ╔ Calibration Isotopes: ════════╤════════╤════════╤═════════╤════════╤═══════════╗'
    print ' ║ Line Name │  Counts │ Eγ,keV  │ σR,keV │ σL,keV │ exptail │ Backgr │   χ²/NDF  ║'
    print ' ╟───────────┼─────────┼─────────┼────────┼────────┼─────────┼────────┼───────────╢'
    for p in self.ScalePeaks:
      PP = (p['name'], p[K]['p0']/self.gain, p[K]['p1'], p[K]['p2'], p[K]['p3'], p[K]['p6'], p[K]['p5'], p[K]['Chi2'], p[K]['NDF'])
      print ' ║%11s│ %7.0f │ %7.2f │ %6.4f │ %06.4f │ %7.3f │ %6.0f │ %5.1f/%3d ║' % PP
    print ' ╚═══════════╧═════════╧═════════╧════════╧════════╧═════════╧════════╧═══════════╝\n'
    print ' ╔ Other Isotopes: ════╤═════════╤════════╤════════╤═════════╤════════╤═══════════╗'
    print ' ║ Line Name │  Counts │ Eγ,keV  │ σR,keV │ σL,keV │ exptail │ Backgr │   χ²/NDF  ║'
    print ' ╟───────────┼─────────┼─────────┼────────┼────────┼─────────┼────────┼───────────╢'
    for p in self.OtherPeaks:
      PP = (p['name'], p[K]['p0']/self.gain, p[K]['p1'], p[K]['p2'], p[K]['p3'], p[K]['p6'], p[K]['p5'], p[K]['Chi2'], p[K]['NDF'])
      print ' ║%11s│ %7.0f │ %7.2f │ %6.4f │ %06.4f │ %7.3f │ %6.0f │ %5.1f/%3d ║' % PP
    print ' ╚═══════════╧═════════╧═════════╧════════╧════════╧═════════╧════════╧═══════════╝\n'

    G = 1./(self.gain*self.LiveT)
#    for p in self.OtherPeaks:
#      if 'In116' in p['name']:
#        filename = p['name'] + '.dat'
#        t, dt = 0.5*(self.UTB + self.UTE), 0.5*(self.UTE - self.UTB)
#        outs = "%10d  %5d  %6.4f  %6.4f\n" % (t, dt, G*p[K]['p0'], G*p[K]['dp0'])
#        with open(filename,'a') as f:  f.write(outs)

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def ShowScale(self,n):
    status = ['═══ Preset',' Corrected']
    print ' ╔%00000009s energy scale: ═══════════════╤═══════════════════════════════════════╗' % status[n]
    if   n==0:
      print ' ║ γ-lines:                               │ pulser:                               ║'
    elif n==1:
      sa = '%5.1f / %1d'   % (self.linear_scale_fit_result.Chi2(),       self.linear_scale_fit_result.Ndf()     )
      if self.spline:
        sb = '%5.2f / %1d' % (self.pulser_correction_fit_result.Chi2(),  self.pulser_correction_fit_result.Ndf())
      else: 
        sb = 'none/nan'
      print ' ║ linear fit:         χ²/NDF = %9s │ pulser fit:        χ²/NDF = %9s ║' % (sa, sb)
    print ' ║ Zero = %5.3f keV, Gain = %6.4f keV/ch │ Zero = %5.2f keV, Gain = %6.2f keV/V ║' % (self.zero, self.gain, self.zero_p, self.gain_p)
    print ' ╚════════════════════════════════════════╧═══════════════════════════════════════╝\n'



# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def InitParameters(self, cfg_file):
    self.scale = ['Cs137(   0)', 'Co60 (   0)', 'Co60 (   1)', 'Tl208(  2)']
    self.names = self.scale
    cfg = ConfigParser.ConfigParser(); cfg.read(cfg_file)

    S = 'scale'
    if cfg.has_section(S):
      self.zero  = cfg.getfloat(S, 'zero');    self.gain = cfg.getfloat(S, 'gain')
      self.emin  = cfg.getfloat(S, 'emin');    self.emax = cfg.getfloat(S, 'emax')
      self.fitL  = cfg.getfloat(S, 'fitL');    self.fitR = cfg.getfloat(S, 'fitR')
      self.nitr  = cfg.getint(  S, 'nitr');    self.amin = cfg.getfloat(S, 'amin')
      self.erec  = cfg.getfloat(S, 'erec');    self.tbpa = cfg.getfloat(S, 'tbpa')
      self.expt  = cfg.get(S, 'expt')
      self.exp0  = cfg.getfloat(S, 'exp0');    self.exp1 = cfg.getfloat(S, 'exp1')
      self.kapp  = ROOT.TF1('kapp', self.expt, self.emin, self.emax)
      self.kapp.SetParameters(self.exp0, self.exp1)
      self.scale = cfg.get(S, 'escale').split(', ')
      self.names = cfg.get(S, 'toshow').split(', ')[0:4]
    else:
      print 'Can not read configuration for scale!'

    S = 'pulser'
    if cfg.has_section(S):
      self.zero_p = cfg.getfloat(S, 'zero_p'); self.gain_p = cfg.getfloat(S, 'gain_p'); self.eerr_p = cfg.getfloat(S, 'eerr_p')
    else:
      print 'Can not read configuration for pulser!'


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def InitGraphics(self):
    # Energy Scale Section
    self.S0 = ROOT.TGraphErrors();          self.S1 = ROOT.TGraphErrors();          self.S2 = ROOT.TGraphErrors()
    self.S0.SetMarkerColor(self.Colors[0]); self.S1.SetMarkerColor(self.Colors[1]); self.S2.SetMarkerColor(self.Colors[2])
    self.S0.SetMarkerStyle(self.Styles[0]); self.S1.SetMarkerStyle(self.Styles[1]); self.S2.SetMarkerStyle(self.Styles[2])
    self.S0.SetMarkerSize(self.Sizes[0]);   self.S1.SetMarkerSize(self.Sizes[1]);   self.S2.SetMarkerSize(self.Sizes[2])
    self.NLG.Add(self.S0);                  self.NLG.Add(self.S1);                  self.NLG.Add(self.S2)
    self.Lg1= ROOT.TLegend(0.62, 0.75, 0.99, 1.00, '', 'brNDC')
    # Energy Resolution Section
    self.RC = ROOT.TGraphErrors();          self.RP = ROOT.TGraphErrors();          self.RO = ROOT.TGraphErrors();
    self.RC.SetMarkerColor(self.Colors[0]); self.RP.SetMarkerColor(self.Colors[2]); self.RO.SetMarkerColor(self.Colors[1])
    self.RC.SetMarkerStyle(self.Styles[0]); self.RP.SetMarkerStyle(self.Styles[1]); self.RO.SetMarkerStyle(self.Styles[2])
    self.RC.SetMarkerSize(self.Sizes[0]);   self.RP.SetMarkerSize(self.Sizes[1]);   self.RO.SetMarkerSize(self.Sizes[2])
    self.ERG.Add(self.RC);                  self.ERG.Add(self.RP);                  self.ERG.Add(self.RO)
    self.Lg2= ROOT.TLegend(0.62, 0.75, 0.99, 1.00, '', 'brNDC')


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Save_Calibration(self):
    W = {}
    W['infos'] = [self.UTB, self.UTE, self.ptype, self.PB5]
    W['scale'] = [self.zero, self.gain, self.zero_p, self.gain_p]
    W['peaks'] = [self.ScalePeaks, self.PulsePeaks, self.OtherPeaks]
    if self.PB5:
      W['scl_Q'] = [self.pulser_correction_fit_result.Chi2(), self.pulser_correction_fit_result.Ndf()]
    else:
      W['scl_Q'] = [0.0, 0.0]
    W['res_Q'] = [self.real_resolution_fit_result.Chi2(),   self.real_resolution_fit_result.Ndf()]
    with open(self.outfile,'a') as fp:  cPickle.dump(W, fp, -1)
    print 'Calibration results have been saved to %s' % self.outfile

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Get_Calibration(self, t = 0, w = 1000.0, pause = False):
    data = []
    OD = {'T':[], 'dT':[], 'Z':[], 'G':[], 'N':[], 'dW':[], 'C':[], 'dC':[], 'R':[], 'dR':[], 'L':[], 'dL':[], 'X':[]}
    with open(self.outfile,'rb') as fp:
      while True:
        try:   data.append(cPickle.load(fp))
        except EOFError: break
    for rec in data:
      self.UTB,    self.UTE,  self.ptype,  self.PB5     = rec['infos'][0], rec['infos'][1], rec['infos'][2], rec['infos'][3]
      if t==0 or (self.UTB - 120) < t < (self.UTE + 120):
        self.zero,   self.gain, self.zero_p, self.gain_p  = rec['scale'][0], rec['scale'][1], rec['scale'][2], rec['scale'][3]
        self.ScalePeaks, self.PulsePeaks, self.OtherPeaks = rec['peaks'][0], rec['peaks'][1], rec['peaks'][2]
        allp = [int(el['shape']['p1']) for el in self.ScalePeaks + self.PulsePeaks + self.OtherPeaks]
        self.OtherPeaks = []
        if len(self.PulsePeaks) < 3: self.PB5 = False
        self.Show_Energy_Scale(self.Do_Energy_Scale(True))
        self.Show_Energy_Resolution(self.Do_Energy_Resolution())
        x = np.ndarray(2, 'float64');   x[0], x[1] =  w, -w;   er = np.ndarray(2, 'float64')
        self.linear_scale_fit_result.GetConfidenceIntervals(1, 1, 1, x, er, 0.683, False)
        dW = er[0]
        if self.spline:
          self.pulser_correction_fit_result.GetConfidenceIntervals(1, 1, 1, x, er, 0.683, False)
          C,  dC = self.SplineB.Eval(w), er[0]
        else: C, dC = 0.0, 0.0
        self.real_resolution_fit_result.GetConfidenceIntervals(2, 1, 1, x, er, 0.683, False)
        SR, dSR = 0.01*w*self.eres.Eval( w), 0.01*w*er[0]
        SL, dSL = 0.01*w*self.eres.Eval(-w), 0.01*w*er[1]
        """
        https://root.cern.ch/root/html/ROOT__Fit__FitResult.html
        GetConfidenceIntervals(unsigned int n, unsigned int stride1, unsigned int stride2, const double* x, double* ci, double cl, bool norm = true)
        get confidence intervals for an array of n points x.
        stride1 indicates the stride in the coordinate space while stride2 the stride in dimension space.
        For 1-dim points : stride1=1, stride2=1
        for multi-dim points arranged as (x0,x1,...,xN,y0,....yN)          stride1=1      stride2=n
        for multi-dim points arraged  as (x0,y0,..,x1,y1,...,xN,yN,..)     stride1=ndim,  stride2=1
        the confidence interval are returned in the array ci
        cl is the desired confidedence interval value
        norm is a flag to control if the intervals need to be normalized to the chi2/ndf value
        By default the intervals are corrected using the chi2/ndf value of the fit if a chi2 fit is performed
        """
        OD['T' ].append( 0.5 * (self.UTB + self.UTE) )
        OD['dT'].append( 0.5 * (self.UTE - self.UTB) )
        OD['Z' ].append(self.zero) ;  OD['G' ].append(self.gain);
        OD['N' ].append(self.ptype);  OD['dW'].append(dW)
        OD['C' ].append(C)         ;  OD['dC'].append(dC)
        OD['R' ].append(SR)        ;  OD['dR'].append(dSR)
        OD['L' ].append(SL)        ;  OD['dL'].append(dSL)
        if t!=0.0:                    OD['X' ].append(allp)
        if pause: raw_input('Press <Enter> to proceed')
    return(OD)

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Check_Resolution_Model(self):
    from numpy import asarray
    ROOT.gStyle.SetOptFit(1111)
    
    def Weighted_Averages(D):
      S0, S1, S2 = 0.0, 0.0, 0.0
      for [v, s] in D:
        w = 1.0 if not s else 1.0/s**2
        S0 += v*w; S1 += v*v*w; S2 += w
      return S0/S2, (S1/S2 - (S0/S2)**2)**0.5
      
    d = [] 
    with open(self.outfile,'rb') as fp:
      while True:
        try:   d.append(cPickle.load(fp))
        except EOFError: break
    
    hs = ROOT.TH1I('hs','', 20, 0., 20.)
    fs = ROOT.TF1( 'fs', '[0]*ROOT::Math::chisquared_pdf(x,[1])', 0., 20.)
    hr = ROOT.TH1I('hr','', 20, 0., 40.)
    fr = ROOT.TF1( 'fr', '[0]*ROOT::Math::chisquared_pdf(x,[1])', 0., 40.)
    D = {}
    ScPks = {}
    n, ndfs, ndfr = 0.0, 0.0, 0.0
    for p in d[0]['peaks'][0]: 
      D[p['E']] = {'R':[], 'L':[], 'K':[]}
      ScPks[p['E']] = {'name':p['name'], 'Chi2':[], 'Ndf':[]}
    for rec in d:
      hs.Fill(rec['scl_Q'][0]); ndfs += rec['scl_Q'][1]
      hr.Fill(rec['res_Q'][0]); ndfr += rec['res_Q'][1]
      n   += 1.0
      ScalePeaks = rec['peaks'][0]
      for p in ScalePeaks:
        D[p['E']]['R'].append([p['shape']['p2'], p['shape']['dp2']])
        D[p['E']]['L'].append([p['shape']['p3'], p['shape']['dp3']])
        D[p['E']]['K'].append([p['shape']['p6'], p['shape']['dp6']])
        ScPks[p['E']]['Chi2'].append(p['shape']['Chi2'])
        ScPks[p['E']]['Ndf' ].append(p['shape']['NDF' ])
    print '%d points found' % n
    ndfs = ndfs / n; ndfr = ndfr / n
    hs.SetNameTitle('hs', 'scale fit: #chi^{2} distribution (NDF = %.2f)' % ndfs)
    hs.SetLineWidth(2);    hs.GetXaxis().SetTitle('#chi^{2}')
    fs.SetParameters(n, ndfs); fs.SetLineColor(ROOT.kRed)
    hr.SetNameTitle('hr', 'resolution fit: #chi^{2} distribution (NDF = %.2f)' % ndfr)
    hr.SetLineWidth(2);    hr.GetXaxis().SetTitle('#chi^{2}')
    fr.SetParameters(2*n, ndfr); fr.SetLineColor(ROOT.kRed)

    n = 0
    gR, gL, gK = ROOT.TGraphErrors(), ROOT.TGraphErrors(), ROOT.TGraphErrors()
    for k,v in D.iteritems():
      R, dR = Weighted_Averages(v['R']); gR.SetPoint(n, k, R); gR.SetPointError(n, 0.0, dR)
      L, dL = Weighted_Averages(v['L']); gL.SetPoint(n, k, L); gL.SetPointError(n, 0.0, dL)
      K, dK = Weighted_Averages(v['K']); gK.SetPoint(n, k, K); gK.SetPointError(n, 0.0, dK)
      n += 1
    gR.SetMarkerStyle(20); gR.SetMarkerColor(ROOT.kRed);  gR.SetLineColor(ROOT.kRed)
    gL.SetMarkerStyle(21); gL.SetMarkerColor(ROOT.kBlue); gL.SetLineColor(ROOT.kBlue)
    gK.SetMarkerStyle(20); gK.SetMarkerColor(ROOT.kRed);  gK.SetLineColor(ROOT.kRed)
    fR = ROOT.TF1('fR', 'sqrt([0]**2 + [1]*2.96e-3*x)')
    fL = ROOT.TF1('fL', 'sqrt([0]**2 + [1]*2.96e-3*x)')
    gK.SetNameTitle('', 'f(x) = ' + self.expt)
    self.cv.Clear();  self.cv.Divide(2, 2) 
    self.cv.cd(1);    self.cv.GetPad(1).SetGrid() 
    gL.Draw('AP');    gL.Fit('fL');  gL.GetYaxis().SetTitle('#sigma_{L} [keV]')
    gR.Draw('PSAME'); gR.Fit('fR');  gR.GetYaxis().SetTitle('#sigma_{R} [keV]'); gR.GetXaxis().SetTitle('keV')
    self.cv.cd(2);    self.cv.GetPad(2).SetGrid() 
    gK.Draw('AP');                   gK.GetYaxis().SetTitle('#kappa [a.u.]');  gK.GetXaxis().SetTitle('keV')
    if self.exp0 == 0.0:  gK.Fit('kapp')
    else:                 self.kapp.Draw('SAME')
    self.cv.cd(3);    self.cv.GetPad(3).SetGrid();    hs.Draw('HIST'); fs.Draw('SAME')
    self.cv.cd(4);    self.cv.GetPad(4).SetGrid();    hr.Draw('HIST'); fr.Draw('SAME')
    self.cv.Modified(); self.cv.Update()
    
    n, plt = 0, {'C':[], 'H':[], 'F':[]}
    for e,t in ScPks.iteritems():
      plt['C'].append(ROOT.TCanvas(str(e), ScPks[e]['name'], n*10+10, n*10+10, 400, 400))
      plt['H'].append(ROOT.TH1I('H%02d' % n, '', 100, 0., 200.))
      plt['F'].append(ROOT.TF1( 'F%02d' % n, '[0]*ROOT::Math::chisquared_pdf(x,[1])', 0., 200.))
      np, ndf = len(t['Chi2']), 0.0
      for i in range(np):
        plt['H'][n].Fill(t['Chi2'][i])
        ndf += t['Ndf'][i]
      ndf = ndf/np
      plt['C'][n].cd(); plt['C'][n].SetGridx(); plt['C'][n].SetGridy()
      plt['H'][n].SetNameTitle('H%02d' % n, 'E = %.3f keV, NDF = %.1f' % (e,ndf))
      plt['H'][n].Draw('HIST')
      plt['F'][n].SetParameters(2*np, ndf)
      plt['F'][n].SetLineColor(ROOT.kRed)
      plt['F'][n].Draw('SAME')
      plt['C'][n].Modified(); plt['C'][n].Update()
      n += 1
    raw_input('Have a look, then press <Enter> to exit.')

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Check_Scale(self, energy, prompt):
    from numpy import asarray
    E = self.Get_Calibration(t=0.0, w=energy, pause = prompt)
    for k,v in E.iteritems():  E[k] = asarray(v)
    mchan = (energy - E['Z'].mean()) / E['G'].mean()
    NE, TF = len(E['T']), '#splitline{%b %d}{%H:%M}%F1970-01-01 00:00:00'
    Wl, Wc, Dc = [], [], []
    for n in range(NE):
      W = E['Z'][n] + E['G'][n]*mchan
      Wl.append(W)
      Wc.append(W - E['C'][n])
      Dc.append((E['dW'][n]**2 + E['dC'][n]**2)**0.5)
    Wl = asarray(Wl);    Wc = asarray(Wc);     Dc = asarray(Dc)
    e1 = ROOT.TGraphErrors(NE, E['T'], Wl,     E['dT'], E['dW'])
    e1.SetMarkerStyle(20); e1.SetMarkerColor(ROOT.kRed);  e1.SetLineColor(ROOT.kRed)
    e2 = ROOT.TGraphErrors(NE, E['T'], Wc,     E['dT'], Dc)
    e2.SetMarkerStyle(21); e2.SetMarkerColor(ROOT.kBlue); e2.SetLineColor(ROOT.kBlue)
    e3 = ROOT.TGraphErrors(NE, E['T'], E['R'], E['dT'], E['dR'])
    e3.SetMarkerStyle(20); e3.SetMarkerColor(ROOT.kRed);  e3.SetLineColor(ROOT.kRed); e3.SetLineWidth(3)
    e4 = ROOT.TGraphErrors(NE, E['T'], E['L'], E['dT'], E['dL'])
    e4.SetMarkerStyle(21); e4.SetMarkerColor(ROOT.kBlue); e4.SetLineColor(ROOT.kBlue)
    C = ROOT.TF1('C','[0]')
    e1.Fit('C','Q'); p1, dp1 = C.GetParameter(0), C.GetParError(0)
    e2.Fit('C','Q'); p2, dp2 = C.GetParameter(0), C.GetParError(0)
    e3.Fit('C','Q'); p3, dp3 = C.GetParameter(0), C.GetParError(0)
    e4.Fit('C','Q'); p4, dp4 = C.GetParameter(0), C.GetParError(0)
    D = p2-p1; dD = (dp1**2 + dp2**2)**0.5
    R = p4/p3; dR = R*((dp4/p4)**2 + (dp3/p3)**2)**0.5
    T1 = 'CH. %d energy equivalent. Average pulser correction: %5.3f #pm %5.3f keV.' % (mchan, D, dD)
    g1 = ROOT.TMultiGraph('g1','scale');      g1.Add(e1); g1.Add(e2); g1.SetTitle(T1)
    l1 = ROOT.TLegend(0.75, 0.75, 0.98, 0.95, '', 'brNDC');
    l1.AddEntry(e1, 'linear scale',    'lpe')
    l1.AddEntry(e2, 'corrected scale', 'lpe')
    T2 = 'Resolution for %5.0f keV #gamma-rays. Average #sigma_{L}/#sigma_{R} = %5.3f #pm %5.3f.' % (energy, R, dR)
    g2 = ROOT.TMultiGraph('g2','resolution'); g2.Add(e3); g2.Add(e4); g2.SetTitle(T2)
    l2 = ROOT.TLegend(0.75, 0.75, 0.98, 0.95, '', 'brNDC');
    l2.AddEntry(e3, '#sigma_{R}', 'lpe')
    l2.AddEntry(e4, '#sigma_{L}', 'lpe')

    cs = ROOT.TCanvas('cs','scale & resolution', 2, 2, 1202, 1002); cs.Divide(1,2)
    cs.cd(1); cs.GetPad(1).SetGrid(); g1.Draw('AP');  l1.Draw()
    cs.cd(2); cs.GetPad(2).SetGrid(); g2.Draw('AP');  l2.Draw()
    for g in [g1,g2]:
      g.GetXaxis().SetTimeDisplay(1); g.GetXaxis().SetTimeFormat(TF); g.GetXaxis().SetLabelOffset(0.03)
      g.GetYaxis().SetDecimals();     g.GetYaxis().SetTitle('keV');   g.GetYaxis().SetTitleOffset(1.2)

    cs.Modified(); cs.Update(); cs.SaveAs(self.outfile + '.pdf')

    print 'Average linear energy    = %8.3f ± %8.3f keV' % (p1, dp1)
    print 'Average corrected energy = %8.3f ± %8.3f keV' % (p2, dp2)
    print 'Average resolution σR    = %8.3f ± %8.3f keV' % (p3, dp3)
    print 'Average resolution σL    = %8.3f ± %8.3f keV' % (p4, dp4)
    raw_input('Have a look, then press <Enter> to exit.')
    self.cv.cd(); self.cv.Clear(); self.cv.Modified(); self.cv.Update() # This is to prevent segmentation fault on exit()

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Check_Pulser(self, energy):
    from numpy import asarray
    d = [] 
    with open(self.outfile,'rb') as fp:
      while True:
        try:   d.append(cPickle.load(fp))
        except EOFError: break
    P_true, n_true, P_fals, n_fals, P_all = ROOT.TGraph(), 0, ROOT.TGraph(), 0, ROOT.TMultiGraph()
    for rec in d:
      t, PB5 = 0.5*(rec['infos'][0] + rec['infos'][1]), rec['infos'][3]
      P0, P1 =      rec['scale'][2],  rec['scale'][3]
#      print '%10d  %8.3f  %8.3f %d' % (t, P0, P1, PB5)
      if PB5:
        P_true.SetPoint(n_true, t, (energy - P0)/P1); n_true  += 1
      else:
        P_fals.SetPoint(n_fals, t, (energy - P0)/P1); n_fals += 1
    self.cv.cd();              self.cv.SetTitle('Pulser Check');   self.cv.SetGrid()
    P_true.SetMarkerStyle(20); P_true.SetMarkerColor( ROOT.kRed);  P_all.Add(P_true)
    P_fals.SetMarkerStyle(24); P_fals.SetMarkerColor(ROOT.kBlue);  P_all.Add(P_fals)
    P_all.Draw('AP'); 
    P_all.GetXaxis().SetLabelOffset(0.03); P_all.GetXaxis().SetTimeDisplay(1); 
    P_all.GetXaxis().SetTimeFormat('#splitline{%b %d}{%H:%M}%F1970-01-01 00:00:00'); 
    P_all.GetYaxis().SetTitle('PB-5 [Volts] for %.1f keV energy' % energy)
    P_all.GetYaxis().SetDecimals()
    self.cv.Modified(); self.cv.Update()
    raw_input('Have a look, then press <Enter> to exit.')
    


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
"""
class LineShape: # monochromatic gamma line response  # class # class # class # class # class # class # class # class # class
  # P0  | P1 | P2      | P3      | P4      | P5         | P6
  # Amp | X0 | Sigma_R | Sigma_L | Compton | Background | kappa
  def __init__(self):
    self.t0 = (0.5 / ROOT.TMath.Pi())**0.5
    self.t1 = (0.5 * ROOT.TMath.Pi())**0.5
  def __call__(self, x, p):
    X = x[0]-p[1]
    if p[4]==100.0: # i.e. for symmetric (pulser) peaks
      return self.t0 * ROOT.TMath.Exp(-0.5*(X/p[2])**2) * p[0] / p[2] + p[5]
    SL = p[2]*p[3]
    if    X >  0.0:       f =                   ROOT.TMath.Exp(-0.5*( X/p[2])**2)
    elif  X < -p[6]*SL:   f = p[4] + (1.0-p[4])*ROOT.TMath.Exp(p[6] * X/SL + 0.5*p[6]**2)
    else:                 f = p[4] + (1.0-p[4])*ROOT.TMath.Exp(-0.5*( X/SL)**2)
    
    nr = ROOT.TMath.Exp(- 0.5*p[6]**2) / p[6]  +  self.t1 * ROOT.TMath.Erf(0.5**0.5 * p[6])
    return p[0] * f / ((1.0-p[4]) * nr * SL + p[2] * self.t1)  +  p[5]
"""

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class LineShape: # monochromatic gamma line response  # class # class # class # class # class # class # class # class # class
  # P0  | P1 | P2      | P3      | P4      | P5         | P6
  # Amp | X0 | Sigma_R | Sigma_L | Compton | Background | kappa
  t0 = (0.5 / ROOT.TMath.Pi())**0.5
  t1 = (0.5 * ROOT.TMath.Pi())**0.5
  pp = 0.5**0.5
  
  def __call__(self, x, p):
    X = x[0]-p[1]
    if p[4]==100.0: # i.e. for symmetric (pulser) peaks
      return self.t0 * ROOT.TMath.Exp(-0.5*(X/p[2])**2) * p[0] / p[2] + p[5]
    SL = p[2]*p[3]
    if    X >  0.0:     f = ROOT.TMath.Exp(-0.5*(X/p[2])**2)
    elif  X < -p[6]*SL: f = ROOT.TMath.Exp(p[6] *X/SL + 0.5*p[6]**2)
    else:               f = ROOT.TMath.Exp(-0.5*(X/SL)**2)
    nr = ROOT.TMath.Exp(- 0.5*p[6]**2) / p[6]  +  self.t1 * ROOT.TMath.Erf( self.pp * p[6])
    try:
      R = p[0] * f / (nr * SL + p[2] * self.t1)  +  0.5 * p[4] * ROOT.TMath.Erfc(self.pp * X/p[2])  +  p[5]
      return R
    except TypeError:
      print self.pp, X, p[2]
      exit()     
    
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class ResolutionModel: # Combine resolution model # class # class # class # class # class # class # class # class # class
  # P0    | P1     | P2     |
  # Noise | Fano_R | Fano_L |
  pce    = 2.96e-3     # electron-hole pair creation energy in Ge, keV

  def __call__(self, x, p):
    X = abs(x[0])
    if x[0]<0.0:
      return 1.e+2 / X * ROOT.TMath.Sqrt(p[0]**2 + p[2]*self.pce*X)
    elif x[0]>0.0:
      return 1.e+2 / X * ROOT.TMath.Sqrt(p[0]**2 + p[1]*self.pce*X)
    else: return 100.0

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
"""
class ResolutionModel: # Combine resolution model # class # class # class # class # class # class # class # class # class
  # P0    | P1   | P2         | P3               |
  # Noise | Fano | Trapping K | Threshold Energy |
  pce    = 2.96e-3     # electron-hole pair creation energy in Ge, keV

  def __call__(self, x, p):
    X = abs(x[0])
    if x[0]<0.0 and p[2]>0.0:
      return 1.e+2 / X * ROOT.TMath.Sqrt(p[0]**2 + p[1]*self.pce*X + p[2]*(X-p[3])**2)
    elif x[0]>0.0:
      return 1.e+2 / X * ROOT.TMath.Sqrt(p[0]**2 + p[1]*self.pce*X                   )
    else: return 100.0
"""

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def fitParameters(fitf):
  pl, np = {}, fitf.GetNumberFreeParameters()
  pl['Chi2'], pl['NDF'], pl['Prob'] = fitf.GetChisquare(),  fitf.GetNDF(),  fitf.GetProb()
  for p in range(np):  pl['p'+str(p)], pl['dp'+str(p)] = fitf.GetParameter(p), fitf.GetParError(p)
  return pl

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def MultiGraphLimits(mg):
  X,Y,eY = [],[],[]
  for g in mg.GetListOfGraphs():
    for n in range(g.GetN()):
      X.append(g.GetX()[n])
      Y.append(g.GetY()[n])
      eY.append(g.GetEY()[n])
  return min(X)-100, max(X)+100, min(Y)-2*max(eY), max(Y)+2*max(eY)



