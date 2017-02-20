# -*- coding: utf-8 -*-
import ROOT, time, sys, ConfigParser
import numpy as np
from atlas   import Atlas
from scipy.interpolate import UnivariateSpline

# http://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.UnivariateSpline.html#scipy.interpolate.UnivariateSpline
class USpline:  # class # class # class # class # class # class # class # class # class # class # class # class # class
  def __init__(self):         self.uspline = UnivariateSpline([0,1,2,3,4],[0,1,2,3,4])
  def Reset(self,X,Y,W):      self.uspline = UnivariateSpline(X, Y, W)
  def __call__(self, x):      return self.uspline(x[0])


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class Isotopes(Atlas): # class # class # class # class # class # class # class # class # class # class # class # class # class
  pce    = 2.96e-3     # electron-hole pair creation energy in Ge, keV
  Colors = [ROOT.kRed+2, ROOT.kBlue+2, ROOT.kGreen+3, ROOT.kGray+2]
  Styles = [20,  21,  22,   23];    Sizes  = [1.2, 1.0, 1.25, 1.25]

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def __init__(self, outfile, cfg_file, action):
    self.InitParameters(cfg_file)    
    self.outfile = outfile
    Atlas.__init__(self)
    if action == 'calibration':
      self.cv   =  ROOT.TCanvas('cv','HPGe calibration', 2, 2, 1002, 1002)
      self.PADS = [ROOT.TPad('PAD5', 'Scale',      0.01, 0.01, 0.49, 0.36, 0, 1),
                   ROOT.TPad('PAD6', 'Resolution', 0.51, 0.01, 0.99, 0.36, 0, 1),
                   ROOT.TPad('PAD0', 'Spectrum',   0.01, 0.66, 0.99, 0.99, 0, 1),
                   ROOT.TPad('PAD1', 'Fit 1',      0.01, 0.37, 0.25, 0.65, 0, 1),
                   ROOT.TPad('PAD2', 'Fit 2',      0.25, 0.37, 0.50, 0.65, 0, 1),
                   ROOT.TPad('PAD3', 'Fit 3',      0.50, 0.37, 0.75, 0.65, 0, 1),
                   ROOT.TPad('PAD4', 'Fit 4',      0.75, 0.37, 0.99, 0.65, 0, 1)]
      for pad in self.PADS:  self.cv.cd(); pad.Draw()

      self.SP = ROOT.TSpectrum()
      # Peaks Fitting Section
      self.asps = ROOT.TF1('asps', LineShape(), 0.0, 1.0, 6);    self.asps.SetLineColor(ROOT.kBlue+2)
      self.asps.SetParNames('Amp','E_{0}, keV', '#sigma_{R}, keV', '#sigma_{L}, keV', 'Compton', 'Background')
      self.p_limits = ((1.0, 1.e+7), (50.0, 1.e+4), (0.5, 10.00),  (0.5, 20.0),  (0.0,0.10), (0.0,1.e+5))
    else:
      self.cv   =  ROOT.TCanvas('cv','HPGe calibration', 2, 2, 1002, 502)
      self.PADS = [ROOT.TPad('PAD5', 'Scale',      0.01, 0.01, 0.49, 0.99, 0, 1),
                   ROOT.TPad('PAD6', 'Resolution', 0.51, 0.01, 0.99, 0.99, 0, 1)]
      for pad in self.PADS:  self.cv.cd(); pad.Draw()
    # Energy Scale Section
    self.NLG      = ROOT.TMultiGraph(); self.NLG.SetTitle('E_{FIT} - E_{REF} [keV]')    
    self.lisc     = ROOT.TF1('lisc', '[0] + [1]*x',                self.emin, self.emax) # linear scale calibration
    self.lipc     = ROOT.TF1('lipc', '[0] + [1]*x',                self.emin, self.emax) # linear pulser calibration
    self.cnst     = ROOT.TF1('cnst', '[0]',                        self.emin, self.emax) 
    self.spline   = USpline()
    self.SplineU  = ROOT.TF1('SplineU', self.spline, self.emin, self.emax, 0)
    self.SplineU.SetLineColor(self.Colors[3]); self.SplineU.SetLineWidth(2); self.SplineU.SetNpx(250); 
    # Energy Resolution Section
    self.ERG      = ROOT.TMultiGraph(); self.ERG.SetTitle('#sigma_{E} / E [%]')
    self.eres_C   = ROOT.TF1('eres_C', ResolutionModel(), -self.emax, self.emax, 4) # combined resolution model [%]
    self.eres_C.SetParameters(1.0,  0.24, 0.000, 100.0); self.eres_C.SetLineColor(self.Colors[0])
    self.eres_C.SetParLimits( 0,    0.40, 10.00);        self.eres_C.SetParLimits(1, 0.10, 0.500)
    self.eres_C.SetParLimits( 2,    0.00, 0.001);        self.eres_C.SetParLimits(3, 0.00, 1000.)
    self.pulser   = ROOT.TF1('pulser', '100.0*([0]/abs(x)+[1])',  -self.emax, self.emax) # pulser resolution [%]
    self.pulser.SetParameters(1.0, 0.0); self.pulser.SetLineColor(self.Colors[2])
    self.InitGraphics()    
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def __del__(self):    Atlas.__del__(self)

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Do(self, utb, ute, pb5, hps, ptype):
    self.UTB, self.UTE, self.hps, self.PB5, self.ptype = utb, ute, hps, pb5, ptype
    self.nbins = self.hps.GetNbinsX()
    self.hps.GetXaxis().SetTitle('E_{#gamma}, keV')
    self.hps.SetBins(self.nbins, self.zero, self.zero + self.gain * self.nbins)
    self.PB5lines()
    self.findPeaks()
    self.identifyPeaks()
    self.ShowScale(0)
    if self.nScalePeaks>2:
      for iteration in range(self.nitr):
        self.fitPeaks(L = self.fitL, R = self.fitR)
        self.Do_Energy_Scale(); 
        self.Show_Energy_Scale()
      self.ShowScale(1)
      self.PeaksTable() 
      self.Show_Energy_Resolution(self.Do_Energy_Resolution())
      self.Save_Calibration()
    return self.zero, self.gain
      #      for i in range(7): self.PADS[i].SaveAs('PAD%1d.eps' % i)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def PB5lines(self):
    self.atlas['PB5'], N = [], len(self.PB5)
    for i in range(N):
      self.atlas['PB5'].append({'Key':str(i), 'W': self.zero_p + self.gain_p*self.PB5[i], 'dW':self.peer_p, 'CC':False, 'RC':False})

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def findPeaks(self):
    rebin = 2
    self.nFoundPeaks, self.FoundPeaks = 0, []
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
    self.KnownPeaks = []
    for p in self.FoundPeaks:
      for nucl,lines in self.atlas.iteritems():
        for line in lines:
          name     = "%5s(%4s)" % (nucl,line['Key'])
          E_t, dE_t  = line['W'], line['dW']
          if abs(p['E']-E_t) < self.erec or ((abs(p['E']-E_t) < 2*self.erec) and ('PB5' in name)):
            self.KnownPeaks.append({'name':name, 'E':E_t, 'dE':dE_t, 'X':p['E'], 'A':p['A'], 'B':p['B'], 'incc':line['CC'], 'inrc':line['RC']})
    self.KnownPeaks.sort(key = lambda peak: peak['E'])
    self.ScalePeaks  = [el for el in self.KnownPeaks if          el['incc']];       self.nScalePeaks = len(self.ScalePeaks)
    self.ShapePeaks  = [el for el in self.KnownPeaks if          el['inrc']];       self.nShapePeaks = len(self.ShapePeaks)
    self.PulsePeaks  = [el for el in self.KnownPeaks if 'PB5' in el['name']];       self.nPulsePeaks = len(self.PulsePeaks)
    self.OtherPeaks  = [el for el in self.KnownPeaks if not el in self.PulsePeaks]
    self.OtherPeaks  = [el for el in self.OtherPeaks if not el in self.ScalePeaks]
    self.nOtherPeaks = len(self.OtherPeaks)
    self.PB5 = self.nPulsePeaks > 5
#    for e in self.KnownPeaks: print e['name'], e['E'], e['X']


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def fitPeaks(self, L = 2.0, R = 2.0):
    for p in self.ScalePeaks + self.PulsePeaks + self.OtherPeaks:
      p1, ALREADY, PB5 = p['X'], p.has_key('shape'), 'PB5' in p['name']
      if ALREADY: # if peak was fitted already:   Amplitude,  Position, Sigma_R, Sigma_L, Compton, Background
        V = p['shape']
        p0,p2,p3,p4,p5 = V['p0'], V['p2'], V['p3'], V['p4'], V['p5']
      else:   # if a peak was not fitted yet:   Amplitude,  Position, Sigma_R, Sigma_L, Compton, Background
        p0,   p2 = p['A'], 0.01 * p1 * (self.pulser.Eval(p1) if PB5 else self.eres_C.Eval(p1))
        p3,p4,p5 = p2,  0.00,  p['B']
      self.asps.SetRange(p1 - L*p3, p1 + R*p2);  self.asps.SetParameters(p0,p1,p2,p3,p4,p5)
      for np in xrange(6): self.asps.SetParLimits(np, self.p_limits[np][0], self.p_limits[np][1])
      if PB5:  
        self.asps.FixParameter(3,p2);  
        self.asps.FixParameter(4,100.0)
#        self.asps.SetRange(p1 - 3*p3, p1 + 4*p2)
      for attempt in xrange(3):
        Result = self.hps.Fit('asps','RSQN')
        if not Result.Status(): break # print Result.Print()
      for np in xrange(6): self.asps.ReleaseParameter(np)
      P = fitParameters(self.asps)
      A_OK = abs(P['dp0']/P['p0']) < 0.50
      E_OK = abs(P['dp1']/P['p1']) < 0.005
      S_OK = abs(P['dp2']/P['p2']) < 0.50
      OK = A_OK and E_OK and S_OK and not Result.Status()
      if OK or not ALREADY: 
        p['shape'] = P
      else:
        if   p in self.ScalePeaks: self.ScalePeaks.remove(p); self.nScalePeaks -= 1
        elif p in self.PulsePeaks: self.PulsePeaks.remove(p); self.nPulsePeaks -= 1
        elif p in self.ShapePeaks: self.ShapePeaks.remove(p); self.nShapePeaks -= 1; self.OtherPeaks.remove(p); self.nOtherPeaks -= 1
        elif p in self.OtherPeaks: self.OtherPeaks.remove(p); self.nOtherPeaks -= 1
      if   self.names[0] in p['name']: self.Show_Peak(p, L, R, 3)
      elif self.names[1] in p['name']: self.Show_Peak(p, L, R, 4)
      elif self.names[2] in p['name']: self.Show_Peak(p, L, R, 5)
      elif self.names[3] in p['name']: self.Show_Peak(p, L, R, 6)
      self.cv.Modified(); self.cv.Update()
    return


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Show_Peak(self, p, VLEVO, VPRAVO, npad):
    C, SR, SL = p['shape']['p1'], p['shape']['p2'], p['shape']['p3']
    L, R, N, M = C-SL*(VLEVO+1), C+SR*(VPRAVO+1), '', ''
    for c in p['name'][0:5]:
      if   c.isdigit(): N+=c
      elif c.isspace(): pass
      else:             M+=c
    name  = '^{%s}%s' % (N,M)    
    self.PADS[npad].cd(); self.PADS[npad].Clear(); self.PADS[npad].SetGrid()
    H = self.hps.DrawCopy();  H.GetXaxis().SetRangeUser(L, R)
    H.SetNameTitle('Spectrum', '%5s (%.3f keV) #chi^{2}/ndf = %.1f/%s' % (name, p['E'], p['shape']['Chi2'], p['shape']['NDF']))
    self.asps.DrawCopy('SAME')
#    self.cv.SaveAs('PEAKs/%.0fkeV.pdf' % peak['E'])
#    raw_input()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Do_Energy_Scale(self):
    K = 'shape'
# 0) LINEAR SCALE CORRECTION WITH CALIBRATION GAMMA LINES
    self.S0.Set(self.nScalePeaks)
    for pk in self.ScalePeaks:
      n = self.ScalePeaks.index(pk)
      self.S0.SetPoint(     n, pk[ 'E'],  pk[K]['p1'] - pk['E'])
      self.S0.SetPointError(n, pk['dE'], (pk['dE']**2 + pk[K]['dp1']**2 + self.tbpa**2)**0.5)
    self.S0.Fit('lisc','QRN')
    k0,  k1 = self.lisc.GetParameter(0), self.lisc.GetParameter(1)
    self.zero = (self.zero-k0)/(1.+k1);     self.gain /= (1.+k1)
    self.hps.SetBins(self.nbins, self.zero, self.zero + self.gain * self.nbins)
    for pk in self.KnownPeaks: pk['X'] = (pk['X']-k0)/(1.+k1)

# 1) SPLINE FOR PULSER LINES
    if self.PB5:
      X, Y, W = [],[],[]
      self.S1.Set(self.nPulsePeaks)
      for pk in self.PulsePeaks:
        n = self.PulsePeaks.index(pk)
        x, y, dx, dy = pk['E'], pk[K]['p1'] - pk['E'], pk['dE'], (pk['dE']**2 + pk[K]['dp1']**2)**0.5
        self.S1.SetPoint(     n,  x,  y);  self.S1.SetPointError(n, dx, dy)
        X.append(x); Y.append(y); W.append(1.0/dy)
      self.spline.Reset(X,Y,W)

# 2) JUST SHOW OTHER GAMMA LINES
    self.S2.Set(self.nOtherPeaks)
    for pk in self.OtherPeaks:
      n = self.OtherPeaks.index(pk)
      self.S2.SetPoint(     n, pk[ 'E'], pk[K][ 'p1'] - pk['E'])
      self.S2.SetPointError(n, pk['dE'], pk[K]['dp1'])

# 3) PULSER AMPLITUDES LINEAR CORRECTION
    if self.PB5:
      R, n = ROOT.TGraphErrors(), 0 
      for pk in self.ScalePeaks:
        R.SetPoint(     n, pk[ 'E'], self.SplineU.Eval(pk['E']) )
        R.SetPointError(n, pk['dE'], (pk['dE']**2 + pk[K]['dp1']**2 + self.tbpa**2)**0.5 ); n+=1
      R.Fit('lipc','QRN'); p0, p1 = self.lipc.GetParameter(0), self.lipc.GetParameter(1)
      for pk in self.PulsePeaks: pk['E'] = (p0-k0+(1.+p1)*pk['E'])/(1.+k1)
      self.zero_p = (p0 - k0 + (1.+p1)*self.zero_p)/(1. + k1);   self.gain_p *= (1. + p1)/(1. + k1)
      del R


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Show_Energy_Scale(self):
    xmin, xmax, ymin, ymax = MultiGraphLimits(self.NLG)
    self.PADS[0].cd();  self.PADS[0].Clear();  self.PADS[0].SetGrid()
    self.NLG.Draw('AP');  
    self.NLG.GetXaxis().SetLimits(   xmin, xmax); self.NLG.GetXaxis().SetTitle('E_{#gamma}, keV')
    self.NLG.GetYaxis().SetRangeUser(ymin, ymax); self.NLG.GetYaxis().SetDecimals()
    self.Lg1.Draw('SAME')
    if self.PB5: self.SplineU.Draw('SAME')
    self.cv.Modified();    self.cv.Update()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Do_Energy_Resolution(self):
    K = 'shape'

# 1. Resolution by pulser (noise level)
    if self.PB5: 
      self.RP.Set(self.nPulsePeaks)
      for p in self.PulsePeaks:
        n = self.PulsePeaks.index(p); 
        self.RP.SetPoint(n, p['E'], 100*p[K]['p2']/p['E']); self.RP.SetPointError(n, p['dE'], 100*p[K]['dp2']/p['E'])
      self.RP.Fit('pulser','QRN');  P = fitParameters(self.pulser); PP = (P['p0'], P['dp0'], 1e+6*P['p1'], 1e+6*P['dp1'], P['Chi2'], P['NDF'])
      print ' ╔ Resolution by pulser: ═══╤══════════════════════════════╤══════════════════════╗'
      print ' ║ Noise: %5.3f ± %5.3f keV │ Slope:  %6.3f ± %5.3f ppm   │ χ2/NDF:   %5.1f/%3d  ║' % PP
      print ' ╚══════════════════════════╧══════════════════════════════╧══════════════════════╝\n'
    
# 2. Resolution by isotopes for calibration
    N = self.nShapePeaks; self.RC.Set(2*N)
    for p in self.ShapePeaks:
      n = self.ShapePeaks.index(p)
      self.RC.SetPoint(N-n-1, -p['E'], 100*p[K]['p3']/p['E']); self.RC.SetPointError(N-n-1, p['dE'], 100*p[K]['dp3']/p['E'])
      self.RC.SetPoint(N+n,    p['E'], 100*p[K]['p2']/p['E']); self.RC.SetPointError(N+n,   p['dE'], 100*p[K]['dp2']/p['E'])

# 3. Resolution for other isotopes
    OP = [p for p in self.OtherPeaks if p not in self.ShapePeaks]; N = len(OP); self.RO.Set(2*N)
    for p in OP:
      n = OP.index(p); 
      self.RO.SetPoint(N-n-1, -p['E'], 100*p[K]['p3']/p['E']); self.RO.SetPointError(N-n-1, p['dE'], 100*p[K]['dp3']/p['E'])
      self.RO.SetPoint(N+n,    p['E'], 100*p[K]['p2']/p['E']); self.RO.SetPointError(N+n,   p['dE'], 100*p[K]['dp2']/p['E'])

    R = self.RC.Fit('eres_C','RSQN')
    P = fitParameters(self.eres_C); quality = '#chi^{2}/NDF = %5.1f/%3d' % (P['Chi2'], P['NDF'])
    PP1 = (P['p0'], P['dp0'], 1e+6*P['p2'], 1e+6*P['dp2'], P['Chi2'], P['NDF'])
    PP2 = (P['p1'], P['dp1'],      P['p3'],      P['dp3'], R.Prob())
    print ' ╔ Resolution by isotopes: ═╤══════════════════════════════╤══════════════════════╗'
    print ' ║ Noise: %5.3f ± %5.3f keV │ Trapping:  %5.3f ± %5.3f ppm │ χ2/NDF:   %5.1f/%3d  ║' % PP1
    print ' ║ Fano fact: %5.3f ± %5.3f │ Threshold: %5.1f ± %5.1f keV │ Probability:  %5.3f  ║' % PP2
    print ' ╚══════════════════════════╧══════════════════════════════╧══════════════════════╝\n'
    return quality

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Show_Energy_Resolution(self, quality):
    xmin, xmax, ymin, ymax = MultiGraphLimits(self.ERG)
    self.Lg2.Clear()
    self.Lg2.AddEntry(self.RC, '#sigma_{L}, #sigma_{R} reference #gamma', 'lpe')
    self.Lg2.AddEntry(self.RO, '#sigma_{L}, #sigma_{R} observed #gamma', 'lpe')
    self.Lg2.AddEntry(self.RP, '#sigma  pulser ', 'lpe')
    self.Lg2.AddEntry(self.eres_C, quality, 'lpe')
    self.PADS[1].cd(); self.PADS[1].Clear(); self.PADS[1].SetGrid()
    self.ERG.Draw('AP'); self.eres_C.DrawCopy('SAME'); self.Lg2.Draw('SAME')
    if self.PB5:         self.pulser.DrawCopy('SAME')
    self.ERG.GetXaxis().SetLimits( -xmax, xmax); self.ERG.GetXaxis().SetTitle('E_{#gamma}, keV');  
    self.ERG.GetYaxis().SetRangeUser(0.0, ymax); self.ERG.GetYaxis().SetDecimals()
    self.cv.Modified(); self.cv.Update()
    return


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def PeaksTable(self):
    K='shape'
    print ' ╔ Pulser: ══╤════════╤══════════╤════════╤════════╤═════════╤════════╤═══════════╗'
    print ' ║ Line Name │ Height │  Eγ,keV  │ σR,keV │ σL,keV │ Compton │ Backgr │   χ2/NDF  ║'
    print ' ╟───────────┼────────┼──────────┼────────┼────────┼─────────┼────────┼───────────╢'
    for p in self.PulsePeaks:
      PP = (p['name'], p[K]['p0'], p[K]['p1'], p[K]['p2'], p[K]['p3'], p[K]['p5'], p[K]['Chi2'], p[K]['NDF'])
      print ' ║%9s│ %6.0f │ %8.2f │ %6.4f │ %6.4f │   ---   │ %6.0f │ %5.1f/%3d ║' % PP
    print ' ╚═══════════╧════════╧══════════╧════════╧════════╧═════════╧════════╧═══════════╝\n'
    print ' ╔ Calibration Isotopes: ════════╤════════╤════════╤═════════╤════════╤═══════════╗'
    print ' ║ Line Name │ Height │  Eγ,keV  │ σR,keV │ σL,keV │ Compton │ Backgr │   χ2/NDF  ║'
    print ' ╟───────────┼────────┼──────────┼────────┼────────┼─────────┼────────┼───────────╢'
    for p in self.ScalePeaks:
      PP = (p['name'], p[K]['p0'], p[K]['p1'], p[K]['p2'], p[K]['p3'], p[K]['p4'], p[K]['p5'], p[K]['Chi2'], p[K]['NDF'])
      print ' ║%9s│ %6.0f │ %8.2f │ %6.4f │ %06.3f │ %6.5f │ %6.0f │ %5.1f/%3d ║' % PP
    print ' ╚═══════════╧════════╧══════════╧════════╧════════╧═════════╧════════╧═══════════╝\n'
    print ' ╔ Other Isotopes: ═══╤══════════╤════════╤════════╤═════════╤════════╤═══════════╗'
    print ' ║ Line Name │ Height │  Eγ,keV  │ σR,keV │ σL,keV │ Compton │ Backgr │   χ2/NDF  ║'
    print ' ╟───────────┼────────┼──────────┼────────┼────────┼─────────┼────────┼───────────╢'
    for p in self.OtherPeaks:
      PP = (p['name'], p[K]['p0'], p[K]['p1'], p[K]['p2'], p[K]['p3'], p[K]['p4'], p[K]['p5'], p[K]['Chi2'], p[K]['NDF'])
      print ' ║%9s│ %6.0f │ %8.2f │ %6.4f │ %06.3f │ %6.5f │ %6.0f │ %5.1f/%3d ║' % PP
    print ' ╚═══════════╧════════╧══════════╧════════╧════════╧═════════╧════════╧═══════════╝\n'
    

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def ShowScale(self,n):
    status = ['═══ Preset',' Corrected']
    print ' ╔%00000009s energy scale: ═══════════════╤═══════════════════════════════════════╗' % status[n]
    print ' ║ γ-lines:                               │ Pulser:                               ║'
    print ' ║ Zero = %5.3f keV, Gain = %6.4f keV/ch │ Zero = %5.3f keV, Gain = %6.1f keV/V ║' % (self.zero, self.gain, self.zero_p, self.gain_p)
    print ' ╚════════════════════════════════════════╧═══════════════════════════════════════╝\n'


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def InitParameters(self, cfg_file):
    self.names = ['Cs137(   0)', 'Co60 (   0)', 'Tl208(  23)', ' O16 (   0)']
    cfg = ConfigParser.ConfigParser(); cfg.read(cfg_file)
    S = 'scale'
    if cfg.has_section(S):
      self.zero = cfg.getfloat(S, 'zero');    self.gain = cfg.getfloat(S, 'gain')
      self.emin = cfg.getfloat(S, 'emin');    self.emax = cfg.getfloat(S, 'emax')
      self.fitL = cfg.getfloat(S, 'fitL');    self.fitR = cfg.getfloat(S, 'fitR')
      self.erec = cfg.getfloat(S, 'erec');    self.tbpa = cfg.getfloat(S, 'tbpa')
      self.nitr = cfg.getint(  S, 'nitr');    self.amin = cfg.getfloat(S, 'amin')
      if cfg.has_option(S, 'name0'): self.names[0] = '%11s' % cfg.get(S, 'name0')
      if cfg.has_option(S, 'name1'): self.names[1] = '%11s' % cfg.get(S, 'name1')
      if cfg.has_option(S, 'name2'): self.names[2] = '%11s' % cfg.get(S, 'name2')
      if cfg.has_option(S, 'name3'): self.names[3] = '%11s' % cfg.get(S, 'name3')
    
    else: print 'Can not read configuration for scale!'; exit()
    S = 'pulser'
    if cfg.has_section(S):
      self.zero_p = cfg.getfloat(S, 'zero_p'); self.gain_p = cfg.getfloat(S, 'gain_p'); self.peer_p = cfg.getfloat(S, 'peer_p')
    else: print 'Can not read configuration for pulser!'


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def InitGraphics(self):
    # Energy Scale Section 
    self.S0   = ROOT.TGraphErrors(); self.S0.SetMarkerColor(self.Colors[0]); self.S0.SetMarkerStyle(self.Styles[0]); self.S0.SetMarkerSize(self.Sizes[0])
    self.S1   = ROOT.TGraphErrors(); self.S1.SetMarkerColor(self.Colors[2]); self.S1.SetMarkerStyle(self.Styles[2]); self.S1.SetMarkerSize(self.Sizes[2])
    self.S2   = ROOT.TGraphErrors(); self.S2.SetMarkerColor(self.Colors[1]); self.S2.SetMarkerStyle(self.Styles[1]); self.S2.SetMarkerSize(self.Sizes[1])
    self.NLG.Add(self.S0); self.NLG.Add(self.S1); self.NLG.Add(self.S2)
    self.Lg1  = ROOT.TLegend(0.62, 0.85, 0.99, 1.00, '', 'brNDC')
    self.Lg1.AddEntry(self.S0, 'reference lines', 'lpe')
    self.Lg1.AddEntry(self.S1, 'pulser lines',    'lpe')
    self.Lg1.AddEntry(self.S2, 'other lines',     'lpe')
    # Energy Resolution Section 
    self.RC   = ROOT.TGraphErrors(); self.RC.SetMarkerColor(self.Colors[0]); self.RC.SetMarkerStyle(self.Styles[0]); self.RC.SetMarkerSize(1.20)
    self.RO   = ROOT.TGraphErrors(); self.RO.SetMarkerColor(self.Colors[1]); self.RO.SetMarkerStyle(self.Styles[1]); self.RO.SetMarkerSize(1.20)
    self.RP   = ROOT.TGraphErrors(); self.RP.SetMarkerColor(self.Colors[2]); self.RP.SetMarkerStyle(self.Styles[2]); self.RP.SetMarkerSize(0.80)
    self.ERG.Add(self.RC); self.ERG.Add(self.RP); self.ERG.Add(self.RO)
    self.Lg2  = ROOT.TLegend(0.62, 0.75, 0.99, 1.00, '', 'brNDC')


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Save_Calibration(self):
    T = '%10d, %10d, "%8s"' % (self.UTB,  self.UTE, self.ptype)
    S = ROOT.TObjString('%.6f, %.6f' % (self.zero, self.gain))
    CARE = ROOT.TList()
    CARE.Add(S)          
    CARE.Add(self.S0) # real peaks TGraphErrors
    CARE.Add(self.S1) # PB-5 peaks TGraphErrors
    CARE.Add(self.RC) # combined: sigma_R if E>0 else sigma_L  TGraphErrors
    CARE.Add(self.RP) # pulser (noise) TGraphErrors
    fp = ROOT.TFile(self.outfile,'UPDATE'); fp.WriteObject(CARE, T); fp.Close()
    del CARE
    print 'Calibration results have been saved to %s' % self.outfile    

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Get_Calibration(self, t = 0, w = 1000.0, pause = False):
    x = np.ndarray(2, 'float64');   x[0], x[1] =  w, -w;   er = np.ndarray(2, 'float64')
    OD = {'T':[], 'dT':[], 'Z':[], 'G':[], 'N':[], 'dW':[], 'C':[], 'dC':[], 'R':[], 'dR':[], 'L':[], 'dL':[]}
    f, CARE  = ROOT.TFile(self.outfile), ROOT.TList()
#    print 'Get calibration from: ', self.outfile  
    for el in [el.GetName() for el in f.GetListOfKeys()]:
      utb, ute, ptype = eval(el)
      if t==0 or utb<t<ute: 
        f.GetObject(el,CARE)
        zero, gain = eval(str(CARE[0]))
        N,X,Y,dX,dY = [],[],[],[],[]
        
        for G in [CARE[1], CARE[2], CARE[3], CARE[4]]:
          Np = G.GetN(); N.append(Np); 
          X.append( [ G.GetX()[n] for n in range(Np)]); 
          dX.append([G.GetEX()[n] for n in range(Np)]); 
          Y.append( [ G.GetY()[n] for n in range(Np)]); 
          dY.append([G.GetEY()[n] for n in range(Np)])
                  
        # 1. Linear Scale Calibration
        self.S0.Set(N[0])
        for p in range(N[0]):  
          self.S0.SetPoint(p,X[0][p],Y[0][p]); self.S0.SetPointError(p,dX[0][p],dY[0][p])
        R = self.S0.Fit('lisc',  'QNRS') 
        R.GetConfidenceIntervals(1, 1, 1, x, er, 0.683, True); dw = er[0]
        
        # 2. Spline Correction & it's uncertainty
        if N[1]>5:
          self.PB5 = True
          for p in range(N[1]):  self.S1.SetPoint(p,X[1][p],Y[1][p]); self.S1.SetPointError(p,dX[1][p],dY[1][p])
          W = [1.0/V for V in dY[1]]; self.spline.Reset(X[1],Y[1],W); c = self.SplineU.Eval(w)
          U = ROOT.TGraphErrors()
          for p in range(N[0]):  
            U.SetPoint(     p,  X[0][p],  Y[0][p] - self.SplineU.Eval(X[0][p]) )
            U.SetPointError(p, dX[0][p], (dY[0][p]**2 + dY[1][p]**2)**0.5)
#            U.SetPointError(p, dX[0][p], dY[0][p])
          R = U.Fit('lipc','QNRS')
          R.GetConfidenceIntervals(1, 1, 1, x, er, 0.683, True); dc = er[0]
          del U

        self.Show_Energy_Scale()

        # 3. Resolution
        self.RC.Set(N[2]) 
        for p in range(N[2]):   self.RC.SetPoint(p,X[2][p],Y[2][p]); self.RC.SetPointError(p,dX[2][p],dY[2][p])

        if self.PB5: 
          self.RP.Set(N[3])
          for p in range(N[3]): self.RP.SetPoint(p,X[3][p],Y[3][p]); self.RP.SetPointError(p,dX[3][p],dY[3][p])
          self.RP.Fit('pulser','QRN')
        
        R = self.RC.Fit('eres_C','QRNS')
        
        quality = '#chi^{2}/NDF = %5.1f/%3d' % (self.eres_C.GetChisquare(),  self.eres_C.GetNDF())
        R.GetConfidenceIntervals(2, 1, 1, x, er, 0.683, True)
        r = self.eres_C.Eval(x[0]); dr = er[0]#;   print r, '±', dr
        l = self.eres_C.Eval(x[1]); dl = er[1]#;   print l, '±', dl
        self.Show_Energy_Resolution(quality)
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
#        time.sleep(0.5)
        if not R.Status(): # i. e. all the fits are O.K.
          OD['T'].append(0.5*(utb+ute));  OD['dT'].append(0.5*(ute-utb))
          OD['Z'].append(zero);           OD['G' ].append(gain)
          OD['N'].append(ptype);          OD['dW'].append(dw)
          OD['C'].append(c);              OD['dC'].append(dc)
          OD['R'].append(0.01*r*w);       OD['dR'].append(0.01*dr*w)
          OD['L'].append(0.01*l*w);       OD['dL'].append(0.01*dl*w)
        if pause: raw_input('Press <Enter> to proceed')
    f.Close()
    return(OD)    
    
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Check_Scale(self, energy, prompt):
    from numpy import asarray
    E = self.Get_Calibration(t=0, w=energy, pause = prompt)
    for k,v in E.iteritems():  E[k] = asarray(v)
    mchan = (energy - E['Z'].mean()) / E['G'].mean()
    NE, W, TF = len(E['T']), [], '#splitline{%b %d}{%H:%M}%F1970-01-01 00:00:00'
    for n in range(NE):  W.append(E['Z'][n] + E['G'][n]*mchan)
    W = asarray(W)
    e1 = ROOT.TGraphErrors(NE, E['T'], W,      E['dT'], E['dW']); e1.SetTitle('CH. %d absolute energy equivalent'           % mchan )
    e2 = ROOT.TGraphErrors(NE, E['T'], E['C'], E['dT'], E['dC']); e2.SetTitle('pulser correction for %5.0f keV #gamma-rays' % energy) 
    e3 = ROOT.TGraphErrors(NE, E['T'], E['R'], E['dT'], E['dR']); e3.SetTitle('#sigma_{RIGHT} for %5.0f keV #gamma-rays'    % energy)
    e4 = ROOT.TGraphErrors(NE, E['T'], E['L'], E['dT'], E['dL']); e4.SetTitle('#sigma_{LEFT} for %5.0f keV #gamma-rays'     % energy)
    lg = ROOT.TLegend(0.75, 0.75, 0.98, 0.95, '', 'brNDC');  lg.AddEntry(e1, self.outfile, 'lpe')
    cs = ROOT.TCanvas('cs','scale & resolution', 2, 2, 1002, 1002); cs.Divide(2,2)
    cs.cd(1); cs.GetPad(1).SetGrid(); e1.Draw('AP'); lg.Draw('SAME')
    cs.cd(2); cs.GetPad(2).SetGrid(); e2.Draw('AP'); lg.Draw('SAME')
    cs.cd(3); cs.GetPad(3).SetGrid(); e3.Draw('AP'); lg.Draw('SAME')
    cs.cd(4); cs.GetPad(4).SetGrid(); e4.Draw('AP'); lg.Draw('SAME')
    for g in [e1,e2,e3,e4]: 
      g.SetMarkerStyle(20);           g.SetMarkerColor(ROOT.kRed);    g.SetLineColor(ROOT.kRed)
      g.GetXaxis().SetTimeDisplay(1); g.GetXaxis().SetTimeFormat(TF); g.GetXaxis().SetLabelOffset(0.03)
      g.GetYaxis().SetDecimals();     g.GetYaxis().SetTitle('keV');   g.GetYaxis().SetTitleOffset(1.2)

    cs.Modified(); cs.Update()
    print 'Average energy      = %8.3f ± %8.3f keV' % (e1.GetMean(2), e1.GetRMS(2))
    print 'Average correction  = %8.3f ± %8.3f keV' % (e2.GetMean(2), e2.GetRMS(2))
    print 'Average right sigma = %8.3f ± %8.3f keV' % (e3.GetMean(2), e3.GetRMS(2))
    print 'Average left  sigma = %8.3f ± %8.3f keV' % (e4.GetMean(2), e4.GetRMS(2))
    raw_input('Have a look, then press <Enter> to exit.')
    self.cv.cd(); self.cv.Clear(); self.cv.Modified(); self.cv.Update() # This is to prevent segmentation fault on exit()

    
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class LineShape: # monochromatic gamma line responce  # class # class # class # class # class # class # class # class # class 
  norm = 1./(2.*ROOT.TMath.Pi())**0.5
  # P0  | P1 | P2      | P3      | P4      | P5         |
  # Amp | X0 | Sigma_R | Sigma_L | Compton | Background |
  def __call__(self, x, p):
    X = x[0]-p[1]
    if p[4]==100.0:
      return self.norm * ROOT.TMath.Exp(-0.5*(X/p[2])**2) * p[0]/p[2] + p[5]
    try:
      if    X >  0.0:            f =                 ROOT.TMath.Exp(-0.5*(X/p[2])**2) 
      else:                      f = p[4] + (1-p[4])*ROOT.TMath.Exp(-0.5*(X/p[3])**2) 
      return 2 * self.norm * f * p[0]/(p[2]+p[3]) + p[5]
    except OverflowError:
      print 'Overflow:', p[0], p[1], p[2], p[3], p[4], p[5]; raw_input()



# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
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


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def fitParameters(fitf):
  pl, np = {}, fitf.GetNumberFreeParameters()
  pl['Chi2'], pl['NDF'] = fitf.GetChisquare(),  fitf.GetNDF()
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



