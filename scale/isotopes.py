# -*- coding: utf-8 -*-
import ROOT, time, os, sys, ConfigParser
import numpy as np
from atlas import Atlas


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class Isotopes(Atlas): # class # class # class # class # class # class # class # class # class # class # class # class # class
  Colors = [ROOT.kRed+2, ROOT.kGreen+3, ROOT.kBlue+2]
  Styles = [20,  22,  24];    Sizes  = [1.25, 1.0, 1.0]
  emin, emax = 100., 10000.  # just for a case when no *.cfg file is present

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def __init__(self, outfile, cfg_file, action):
    self.outfile = outfile
    Atlas.__init__(self)
    if action == 'calibration':
      self.InitParameters(cfg_file)
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
#        if pad > 2:   self.PADS[pad].SetLogy()

      self.SP = ROOT.TSpectrum()
      # Peaks Fitting Section
      self.asps = ROOT.TF1('asps', LineShape(), 0.0, 1.0, 9);    self.asps.SetLineColor(ROOT.kBlue+2)
      self.asps.SetParName(0, 'Amp')
      self.asps.SetParName(1, 'E_{0}, keV')
      self.asps.SetParName(2, '#sigma_{R}, keV')
      self.asps.SetParName(3, '#sigma_{L}, keV')
      self.asps.SetParName(4, '#kappa')
      self.asps.SetParName(5, 'Background')
      self.asps.SetParName(6, 'Amp_{Compton}')
      self.asps.SetParName(7, 'E_{Compton}')
      self.asps.SetParName(8, '#sigma_{Compton}')
      self.p_limits = ((1.0, 1.e+7), (50.0, 1.e+4), (0.1, 10.00),  (0.1, 20.0),  (0.1,10.))
      self.peak = ROOT.TF1('peak', Photopeak(),   0.0, 1.0, 7);    self.peak.SetLineColor(ROOT.kBlack)
      self.asym = ROOT.TF1('asym', Photopeak_A(), 0.0, 1.0, 5);    self.asym.SetLineColor(ROOT.kRed+2)
      self.symm = ROOT.TF1('symm', Photopeak_S(), 0.0, 1.0, 3);    self.symm.SetLineColor(ROOT.kGreen+2)
      self.step = ROOT.TF1('step', Background(),  0.0, 1.0, 4);    self.step.SetLineColor(ROOT.kRed+2)
    else:
      self.cv   =  ROOT.TCanvas('cv','HPGe calibration', 2, 2, 1002, 502)
      self.PADS = [ROOT.TPad('PAD5', 'Scale',      0.01, 0.01, 0.49, 0.99, 0, 1),
                   ROOT.TPad('PAD6', 'Resolution', 0.51, 0.01, 0.99, 0.99, 0, 1)]
      for pad in self.PADS:  self.cv.cd(); pad.Draw()
    # Energy Scale Section
    self.NLG      = ROOT.TMultiGraph(); self.NLG.SetTitle('E_{FIT} - E_{REF} [keV]')
    self.lisc     = ROOT.TF1('lisc', '[0] + [1]*x',                self.emin, self.emax) # linear scale calibration
    self.cnst     = ROOT.TF1('cnst', '[0]',                        self.emin, self.emax)

    # Energy Resolution Section
    self.ERG      = ROOT.TMultiGraph(); self.ERG.SetTitle('#sigma_{E} / E [%]')
    self.eres_C   = ROOT.TF1('eres_C', ResolutionModel(), -self.emax, self.emax, 4) # combined resolution model [%]
    self.eres_C.SetParameters(1.0,  0.24, 0.1e-6, 100.0); self.eres_C.SetLineColor(self.Colors[0])
    self.eres_C.SetParLimits( 0,    0.40, 10.00);         self.eres_C.SetParLimits(1, 0.05, 0.500)
    self.eres_C.SetParLimits( 2,    0.005e-6, 0.15e-6)
    self.InitGraphics()
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def __del__(self):
    Atlas.__del__(self);
    self.cv.cd(); self.cv.Clear()

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
      for iteration in range(self.nitr):
        self.fitPeaks(L = self.fitL, R = self.fitR)
        self.Do_Energy_Scale();
        self.Show_Energy_Scale()
      self.ShowScale(1)
      self.PeaksTable()
      self.Show_Energy_Resolution(self.Do_Energy_Resolution())
    else:
      raw_input('Scale calibration is impossible!')
    return self.zero, self.gain
      #      for i in range(7): self.PADS[i].SaveAs('PAD%1d.eps' % i)

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
    self.ScalePeaks, self.OtherPeaks = [], []
    while len(self.FoundPeaks):
      x = self.FoundPeaks[0]['E']
      name = False
      for nucl,lines in self.atlas.iteritems():
        for line in lines:
          if abs(x - line['W']) < self.erec:
            name = "%5s(%4s)" % (nucl, lines.index(line))
            V = {'name':name, 'E':line['W'], 'dE':line['dW'], 'X':x, 'A':self.FoundPeaks[0]['A'], 'B':self.FoundPeaks[0]['B']}
            if name in self.scale: self.ScalePeaks.append(V)
            else:                  self.OtherPeaks.append(V)
            break
        else:
          continue
        break
      if not name: print 'Unknown line: %.3f keV' % self.FoundPeaks[0]['E']
      self.FoundPeaks.pop(0)

    self.ScalePeaks.sort(key = lambda p: p['E'])
    self.OtherPeaks.sort(key = lambda p: p['E'])


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def fitPeaks(self, L = 2.0, R = 2.0):
    for p in self.ScalePeaks + self.OtherPeaks:
      p1, ALREADY = p['X'], p.has_key('shape')
      if ALREADY: # if peak was fitted already:   Amplitude,  Position, Sigma_R, Sigma_L, kappa, background, Compton
        V = p['shape']
        p0,p1,p2,p3,p4,p5,p6,p7,p8 = V['p0'], V['p1'], V['p2'], V['p3'], V['p4'], V['p5'], V['p6'], V['p7'], V['p8']
      else:   # if a peak was not fitted yet:   Amplitude,  Position, Sigma_R, Sigma_L, kappa, Background, Compton
        re = 0.001 * p['X']
        p0,p1,p2,p3,p4,p5,p6,p7,p8 = 2*re*p['A'], p['X'], re, 1.2*re, 1.0, p['B'], 100.0, p['X']-1.0, re
      self.asps.SetRange(p1 - L*p3, p1 + R*p2);  self.asps.SetParameters(p0,p1,p2,p3,p4,p5,p6,p7,p8)
      for np in xrange(5): self.asps.SetParLimits(np, self.p_limits[np][0], self.p_limits[np][1])

#      self.asps.SetParLimits(8,p2,p3)
#      self.asps.SetParLimits(9,0.0,p0)
#      self.asps.FixParameter(9,3.0)
      for attempt in xrange(3):
        Result = self.hps.Fit('asps','RSQN')
#        print Result.Print()
#        raw_input()
        if not Result.Status(): break
#      for np in range(9,11): self.asps.ReleaseParameter(np)
      P = fitParameters(self.asps)
      A_OK = abs(P['dp0']/P['p0']) < 0.50
      E_OK = abs(P['dp1']/P['p1']) < 0.01
      S_OK = abs(P['dp2']/P['p2']) < 0.50
      OK = A_OK and E_OK and S_OK and not Result.Status()
      if OK or not ALREADY:
        p['shape'] = P
      else:
        if   p in self.ScalePeaks: self.ScalePeaks.remove(p)
        elif p in self.OtherPeaks: self.OtherPeaks.remove(p)
      self.peak.SetRange(p1 - L*p3, p1 + R*p2)
      self.peak.SetParameters(P['p0'], P['p1'], P['p2'], P['p3'], P['p4'], P['p7'], P['p8'])
      self.asym.SetRange(p1 - L*p3, p1 + R*p2)
      self.asym.SetParameters(P['p0'], P['p1'], P['p2'], P['p3'], P['p4'])
      self.symm.SetRange(p1 - L*p3, p1 + R*p2)
      self.symm.SetParameters(P['p0'], P['p7'], P['p8'])
      self.step.SetRange(p1 - L*p3, p1 + R*p2)
      self.step.SetParameters(P['p5'], P['p6'], P['p7'], P['p8'])
      if   self.names[0] in p['name']: self.Show_Peak(p, L, R, 3);  self.Show_Diff(p, L, R, 4)
      elif self.names[1] in p['name']: self.Show_Peak(p, L, R, 5);  self.Show_Diff(p, L, R, 6)
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
    self.asym.DrawCopy('SAME')
    self.symm.DrawCopy('SAME')
    self.cv.Modified(); self.cv.Update()
#    self.PADS[npad].SaveAs('%.0fkeV.png' % p['E'])

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Show_Diff(self, p, VLEVO, VPRAVO, npad):
    C, SR, SL = p['shape']['p1'], p['shape']['p2'], p['shape']['p3']
    L, R, N, M = C-SL*(VLEVO+1), C+SR*(VPRAVO+1), '', ''
    for c in p['name'][0:5]:
      if   c.isdigit(): N+=c
      elif c.isspace(): pass
      else:             M+=c
    name  = '^{%s}%s' % (N,M)
    self.PADS[npad].cd(); self.PADS[npad].Clear(); self.PADS[npad].SetGrid()
    lbin = int((L-self.zero)/self.gain); L = self.zero + self.gain*lbin
    rbin = int((R-self.zero)/self.gain); R = self.zero + self.gain*rbin
    nbin = rbin-lbin
    D    = ROOT.TH1I('D', '', nbin, L, R)
    for ch in range(lbin,rbin+1):
      F = self.peak.Eval(self.zero + (ch-0.5)*self.gain)
      D.SetBinContent(ch-lbin, self.hps.GetBinContent(ch) - F)
      D.SetBinError(ch-lbin, self.hps.GetBinError(ch))
    D.SetNameTitle('Difference', '%5s (%.3f keV) #chi^{2}/ndf = %.1f/%s' % (name, p['E'], p['shape']['Chi2'], p['shape']['NDF']))
    D.DrawCopy()
    self.step.DrawCopy('SAME')
    self.cv.Modified(); self.cv.Update()
#    self.PADS[npad].SaveAs('%.0fkeV.png' % p['E'])

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Save_Peak(self, p, VLEVO, VPRAVO):
    C, SR, SL = p['shape']['p1'], p['shape']['p2'], p['shape']['p3']
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
  def Do_Energy_Scale(self):
    K = 'shape'
# 0) LINEAR SCALE CORRECTION WITH CALIBRATION GAMMA LINES
    self.S0.Set(len(self.ScalePeaks))
    for pk in self.ScalePeaks:
      n = self.ScalePeaks.index(pk)
      self.S0.SetPoint(     n, pk[ 'E'],  pk[K]['p1'] - pk['E'])
      self.S0.SetPointError(n, pk['dE'], (pk['dE']**2 + pk[K]['dp1']**2)**0.5)
    self.S0.Fit('lisc','QRN')
    k0,  k1 = self.lisc.GetParameter(0), self.lisc.GetParameter(1)
    self.zero = (self.zero-k0)/(1.+k1);     self.gain /= (1.+k1)
    self.hps.SetBins(self.nbins, self.zero, self.zero + self.gain * self.nbins)

# 1) JUST SHOW OTHER GAMMA LINES
    self.S2.Set(len(self.OtherPeaks))
    for pk in self.OtherPeaks:
      n = self.OtherPeaks.index(pk)
      self.S2.SetPoint(     n, pk[ 'E'], pk[K][ 'p1'] - pk['E'])
      self.S2.SetPointError(n, pk['dE'], pk[K]['dp1'])

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Show_Energy_Scale(self):
    xmin, xmax, ymin, ymax = MultiGraphLimits(self.NLG)
    self.PADS[0].cd();  self.PADS[0].Clear();  self.PADS[0].SetGrid()
    self.NLG.Draw('AP');
    self.NLG.GetXaxis().SetLimits(   xmin, xmax); self.NLG.GetXaxis().SetTitle(' E_{#gamma} [keV]')
    self.NLG.GetYaxis().SetRangeUser(ymin, ymax); self.NLG.GetYaxis().SetDecimals()
    self.Lg1.Draw('SAME')
#    if self.PB5: self.SplineU.Draw('SAME')
    self.cv.Modified();    self.cv.Update()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Do_Energy_Resolution(self):
    K = 'shape'

# 1. Resolution by isotopes for calibration
    N = len(self.ScalePeaks); self.RC.Set(2*N)
    for p in self.ScalePeaks:
      n = self.ScalePeaks.index(p)
      self.RC.SetPoint(N-n-1, -p['E'], 100*p[K]['p3']/p['E']); self.RC.SetPointError(N-n-1, p['dE'], 100*p[K]['dp3']/p['E'])
      self.RC.SetPoint(N+n,    p['E'], 100*p[K]['p2']/p['E']); self.RC.SetPointError(N+n,   p['dE'], 100*p[K]['dp2']/p['E'])

# 2. Resolution for other isotopes
    OP = [p for p in self.OtherPeaks if p not in self.ScalePeaks]; N = len(OP); self.RO.Set(2*N)
    for p in OP:
      n = OP.index(p);
      self.RO.SetPoint(N-n-1, -p['E'], 100*p[K]['p3']/p['E']); self.RO.SetPointError(N-n-1, p['dE'], 100*p[K]['dp3']/p['E'])
      self.RO.SetPoint(N+n,    p['E'], 100*p[K]['p2']/p['E']); self.RO.SetPointError(N+n,   p['dE'], 100*p[K]['dp2']/p['E'])

    if self.rep3: self.eres_C.FixParameter(3, self.rep3)
    else:         self.eres_C.SetParLimits(3, -500.00, 500.)
#    self.eres_C.FixParameter(2, 0.01e-6)
#    self.eres_C.FixParameter(3, 0.0)
    R = self.RC.Fit('eres_C','RSQN');
#    self.eres_C.ReleaseParameter(2)
#    self.eres_C.ReleaseParameter(3)
    if self.rep3: self.eres_C.ReleaseParameter(3)
    P = fitParameters(self.eres_C); quality = '#chi^{2}/NDF = %5.1f/%3d' % (P['Chi2'], P['NDF'])
    PP1 = (P['p0'], P['dp0'], 1e+6*P['p2'], 1e+6*P['dp2'], P['Chi2'], P['NDF'])
    PP2 = (P['p1'], P['dp1'],      P['p3'],      P['dp3'], R.Prob())
    print ' ╔ Resolution by isotopes: ═╤══════════════════════════════╤══════════════════════╗'
    print ' ║ Noise: %5.3f ± %5.3f keV │ Trapping:  %5.3f ± %5.3f ppm │ χ2/NDF:   %5.1f/%3d  ║' % PP1
    print ' ║ Fano fact: %5.3f ± %5.3f │ Threshold: %5.0f ± %5.0f keV │ Probability:  %5.3f  ║' % PP2
    print ' ╚══════════════════════════╧══════════════════════════════╧══════════════════════╝\n'

#    with open('p3_vs_p2.txt', 'a') as fp:
#      fp.write('%5.3f  %5.3f  %5.0f  %5.0f\n' % (1e+6*P['p2'], 1e+6*P['dp2'], P['p3'], P['dp3']))

    return quality

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Show_Energy_Resolution(self, quality):
    xmin, xmax, ymin, ymax = MultiGraphLimits(self.ERG)
    self.Lg2.Clear()
    self.Lg2.AddEntry(self.RC, '#sigma_{L}, #sigma_{R} reference #gamma', 'lpe')
    self.Lg2.AddEntry(self.RO, '#sigma_{L}, #sigma_{R} observed #gamma', 'lpe')
    self.Lg2.AddEntry(self.eres_C, quality, 'lpe')
    self.PADS[1].cd(); self.PADS[1].Clear(); self.PADS[1].SetGrid()
    self.ERG.Draw('AP'); self.eres_C.DrawCopy('SAME'); self.Lg2.Draw('SAME')
    self.ERG.GetXaxis().SetLimits( -xmax, xmax); self.ERG.GetXaxis().SetTitle(' E_{#gamma} [keV]');
    self.ERG.GetYaxis().SetRangeUser(0.0, ymax); self.ERG.GetYaxis().SetDecimals()
    self.cv.Modified(); self.cv.Update()
    return


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def PeaksTable(self):
    K='shape'
    print ' ╔ Calibration Isotopes: ════════╤════════╤════════╤═════════╤════════╤═══════════╗'
    print ' ║ Line Name │ Height │  Eγ,keV  │ σR,keV │ σL,keV │ Compton │ Backgr │   χ2/NDF  ║'
    print ' ╟───────────┼────────┼──────────┼────────┼────────┼─────────┼────────┼───────────╢'
    for p in self.ScalePeaks:
      PP = (p['name'], p[K]['p0'], p[K]['p1'], p[K]['p2'], p[K]['p3'], p[K]['p6'], p[K]['p5'], p[K]['Chi2'], p[K]['NDF'])
      print ' ║%9s│ %6.0f │ %8.2f │ %6.4f │ %06.3f │ %7.0f │ %6.0f │ %5.1f/%3d ║' % PP
    print ' ╚═══════════╧════════╧══════════╧════════╧════════╧═════════╧════════╧═══════════╝\n'
    print ' ╔ Other Isotopes: ═══╤══════════╤════════╤════════╤═════════╤════════╤═══════════╗'
    print ' ║ Line Name │ Height │  Eγ,keV  │ σR,keV │ σL,keV │ Compton │ Backgr │   χ2/NDF  ║'
    print ' ╟───────────┼────────┼──────────┼────────┼────────┼─────────┼────────┼───────────╢'
    for p in self.OtherPeaks:
      PP = (p['name'], p[K]['p0'], p[K]['p1'], p[K]['p2'], p[K]['p3'], p[K]['p6'], p[K]['p5'], p[K]['Chi2'], p[K]['NDF'])
      print ' ║%9s│ %6.0f │ %8.2f │ %6.4f │ %06.3f │ %7.0f │ %6.0f │ %5.1f/%3d ║' % PP
    print ' ╚═══════════╧════════╧══════════╧════════╧════════╧═════════╧════════╧═══════════╝\n'

    print self.hps.GetTitle().replace('Live','\nLive')
    for p in self.ScalePeaks:
      if 'Co60' in p['name']:
        V = (p[K]['p1'],  p[K]['p2'],  p[K]['p3'],  p[K]['p4'],  p[K]['p5'],  p[K]['p6'],  p[K]['p7'],  p[K]['p8'])
        E = (p[K]['dp1'], p[K]['dp2'], p[K]['dp3'], p[K]['dp4'], p[K]['dp5'], p[K]['dp6'], p[K]['dp7'], p[K]['dp8'])
        print ' ╔═════════╤══════════╤════════╤════════╤═════╤═══════╤═══════╤══════════╤════════╗'
        print ' ║ %000008s│  E₁,keV  │ σR,keV │ σL,keV │  κ  │ Bckgr │ Compt │  E₂,keV  │ σ, keV ║' % p['name'].replace(' ','')
        print ' ╟─────────┼──────────┼────────┼────────┼─────┼───────┼───────┼──────────┼────────╢'
        print ' ║ Value   │ %0008.3f │ %06.3f │ %06.3f │%5.2f│ %005d │ %005d │ %0008.3f │ %06.3f ║' % V
        print ' ║ Error   │ %0008.3f │ %06.3f │ %06.3f │%5.2f│ %005d │ %005d │ %0008.3f │ %06.3f ║' % E
        print ' ╚═════════╧══════════╧════════╧════════╧═════╧═══════╧═══════╧══════════╧════════╝\n'

        CC = 1./(self.gain*self.LiveT)
        ppc, dppc = p[K]['p0']*CC, p[K]['dp0']*CC
        bgc, dbgc = p[K]['p5']*CC, p[K]['dp5']*CC
        scc, dscc = p[K]['p6']*CC, p[K]['dp6']*CC
        R     = scc/ppc
        dR    = R*((dppc/ppc)**2 + (dscc/scc)**2)**0.5
        print 'Integral   = %8f ± %6f [cnts/s]'     % (ppc, dppc)
        print 'Scattering = %8f ± %6f [cnts/keV/s]' % (scc, dscc)
        print 'Ratio      = %8f ± %6f [1/keV]'      % (R,   dR)
        print 'Background = %8f ± %6f [cnts/keV/s]' % (bgc, dbgc)
        outs  = '%8.6f  %7.6f  ' % (ppc, dppc)
        outs += '%8.6f  %7.6f  ' % (scc, dscc)
        outs += '%8.6f  %7.6f  ' % (R,   dR)
        outs += '%8.6f  %7.6f\n' % (bgc, dbgc)
        fname = p['name'].replace(' ','').replace('(','_').replace(')','') + '.dat'
        with open(fname,'a') as f:  f.write(outs)


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def ShowScale(self,n):
    status = ['═══ Preset',' Corrected']
    print ' ╔%00000009s energy scale: ═══════════════╤═══════════════════════════════════════╗' % status[n]
    print ' ║ γ-lines:                               │ Pulser:                               ║'
    print ' ║ Zero = %5.3f keV, Gain = %6.4f keV/ch │ Zero = %5.3f keV, Gain = %6.1f keV/V ║' % (self.zero, self.gain, self.zero, self.gain)
    print ' ╚════════════════════════════════════════╧═══════════════════════════════════════╝\n'


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def InitParameters(self, cfg_file):
    self.scale = ['Cs137(   0)', 'Co60 (   0)', 'Co60 (   1)', 'Tl208(  23)']
    self.names = self.scale
    cfg = ConfigParser.ConfigParser(); cfg.read(cfg_file)

    S = 'scale'
    if cfg.has_section(S):
      self.zero = cfg.getfloat(S, 'zero');    self.gain = cfg.getfloat(S, 'gain')
      self.emin = cfg.getfloat(S, 'emin');    self.emax = cfg.getfloat(S, 'emax')
      self.fitL = cfg.getfloat(S, 'fitL');    self.fitR = cfg.getfloat(S, 'fitR')
      self.erec = cfg.getfloat(S, 'erec');    #self.tbpa = cfg.getfloat(S, 'tbpa')
      self.nitr = cfg.getint(  S, 'nitr');    self.amin = cfg.getfloat(S, 'amin')
      if cfg.has_option(S, 'thre'):           self.rep3 = cfg.getfloat(S, 'thre')
      else:                                   self.rep3 = 0.0
      if cfg.has_option(S, 'escale'): self.scale = cfg.get(S, 'escale').split(', ')
      if cfg.has_option(S, 'toshow'): self.names = cfg.get(S, 'toshow').split(', ')[0:4]
    else:
      print 'Can not read configuration for scale!'; exit()


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def InitGraphics(self):
    # Energy Scale Section
    self.S0   = ROOT.TGraphErrors(); self.S0.SetMarkerColor(self.Colors[0]); self.S0.SetMarkerStyle(self.Styles[0]); self.S0.SetMarkerSize(self.Sizes[0])
    self.S2   = ROOT.TGraphErrors(); self.S2.SetMarkerColor(self.Colors[2]); self.S2.SetMarkerStyle(self.Styles[2]); self.S2.SetMarkerSize(self.Sizes[2])
    self.NLG.Add(self.S0); self.NLG.Add(self.S2)
    self.Lg1  = ROOT.TLegend(0.62, 0.85, 0.99, 1.00, '', 'brNDC')
    self.Lg1.AddEntry(self.S0, 'reference lines', 'lpe')
    self.Lg1.AddEntry(self.S2, 'other lines',     'lpe')
    # Energy Resolution Section
    self.RC   = ROOT.TGraphErrors(); self.RC.SetMarkerColor(self.Colors[0]); self.RC.SetMarkerStyle(self.Styles[0]); self.RC.SetMarkerSize(self.Sizes[0])
    self.RO   = ROOT.TGraphErrors(); self.RO.SetMarkerColor(self.Colors[1]); self.RO.SetMarkerStyle(self.Styles[2]); self.RO.SetMarkerSize(self.Sizes[2])
    self.ERG.Add(self.RC); self.ERG.Add(self.RO)
    self.Lg2  = ROOT.TLegend(0.62, 0.75, 0.99, 1.00, '', 'brNDC')

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class Photopeak_A: # Asymmetric photopeak has 5 parameters:
#  P0  | P1 | P2      | P3      | P4
#  Amp | X0 | Sigma_R | Sigma_L | kappa
  tw, sq  = 0.5**0.5, (0.5*ROOT.TMath.Pi())**0.5
  def __call__(self, x, p):
    X = x[0]-p[1]
    if                  X  >= 0.0:  f = ROOT.TMath.Exp(-0.5*(X/p[2])**2)
    elif   -p[4]*p[3] < X  <  0.0:  f = ROOT.TMath.Exp(-0.5*(X/p[3])**2)
    else:                           f = ROOT.TMath.Exp(p[4]*X/p[3] + 0.5*p[4]**2)
    N1 = p[3]/p[4]*ROOT.TMath.Exp(-0.5*p[4]**2) + self.sq*(p[3]*ROOT.TMath.Erf(p[4]*self.tw) + p[2])
    return 0.75*p[0]*f/N1

class Photopeak_S: # Symmetric photopeak has 3 parameters:
#  P0  | P1 | P2
#  Amp | X0 | Sigma
  isq = (0.5/ROOT.TMath.Pi())**0.5
  def __call__(self, x, p):
    return 0.25*p[0] * ROOT.TMath.Exp(-0.5*((x[0]-p[1])/p[2])**2) * self.isq/p[2]

class Photopeak(Photopeak_A,Photopeak_S): # Complex photopeak has 7 parameters
#  P0  | P1 | P2      | P3      | P4    | P5   | P6
#  Amp | X0 | Sigma_R | Sigma_L | kappa | X0_S | Sigma_S
  pa = ROOT.TF1('pa', Photopeak_A(), 0.0, 1.0, 5)
  ps = ROOT.TF1('ps', Photopeak_S(), 0.0, 1.0, 3)
  def __call__(self, x, p):
    self.pa.SetParameters(p[0], p[1], p[2], p[3], p[4])
    self.ps.SetParameters(p[0], p[5], p[6])
    return self.pa.Eval(x[0]) + self.ps.Eval(x[0])

class Background: # Background has 4 parameters:
#  P0          | P1          | P2         | P3
#  Background  | Amp_Compton | X0_Compton | Sigma_Compton
  isq2   = 0.5**0.5
  def __call__(self, x, p):
    return p[0] + p[1] * 0.5 * ROOT.TMath.Erfc(self.isq2*(x[0]-p[2])/p[3])


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class LineShape(Photopeak,Background):
  """
  Lineshape has 9 parameters:
  P0  | P1 | P2      | P3      | P4    | P5         | P6          | P7         | P8
  Amp | X0 | Sigma_R | Sigma_L | kappa | Background | Amp_Compton | X0_Compton | Sigma_Compton
  """
  pp = ROOT.TF1('pp', Photopeak(),  0.0, 1.0, 7)
  bg = ROOT.TF1('bg', Background(), 0.0, 1.0, 4)

  def __call__(self, x, p):
    self.pp.SetParameters(p[0], p[1], p[2], p[3], p[4], p[7], p[8])
    self.bg.SetParameters(p[5], p[6], p[7], p[8])
    return self.pp.Eval(x[0]) + self.bg.Eval(x[0])



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



