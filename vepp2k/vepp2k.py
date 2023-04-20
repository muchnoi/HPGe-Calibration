# -*- coding: utf-8 -*-
import pickle, os, time, configparser
import numpy as np
import ROOT, sys
from multiprocessing import Process, Queue
from scale.scale import Scale
sql = True
try:
  from .temdbWO import Storage
except ImportError:
  sql = False
  print('database is unavailable')

class Constants:
  c    = 299792458.0      # speed of light          [m/s]
  me   = 0.510998910e+6   # electron rest energy    [eV]
  h    = 4.13566752e-15   # Plank constant          [eV*s]
  hbar = 6.58211928e-16   # Plank constant reduced  [eV*s]
  Bo   = me**2/hbar/c**2  # Bo = me[eV]^2 / hbar[eV*s] / c[m/s]^2 [T]
  RemB = 1.e+6/me/Bo      # i.e. [1 MeV] / Bo [T] / me [eV]

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class EDGE(Constants):

  def __init__(self, scalefile, cfg_file, folder):
    self.CO_Grating = {0  : 5.426463e-6,
                       250: 5.426463e-6,
                       252: 5.426463e-6,
                       322: 5.604325e-6,
                       324: 5.604325e-6,
                       325: 5.617221e-6,
                       327: 5.617221e-6}
    self.exclude = []
    cfg = configparser.ConfigParser(); cfg.read(cfg_file)
    self.MinAmp  = cfg.getfloat(  'edge', 'MinAmpEdge')
    self.Ranger  = cfg.getfloat(  'edge', 'EdgeRanger')
    self.Merger  =   cfg.getint(  'edge', 'BinsMerger')
    self.ETuner  = cfg.getfloat(  'edge', 'EveppTuner')
    self.BTuner  = cfg.getfloat(  'edge', 'BveppTuner')
    self.Radius  = cfg.getfloat(  'edge', 'VeppRadius')
    self.Asymme  = cfg.getboolean('edge', 'Asymmetry')
    self.EdgeP2  = cfg.getboolean('edge', 'EdgePoly2')
    self.SaveDB  = cfg.getboolean('edge', 'SaveForSND') and sql
    self.negate  = cfg.getboolean('edge', 'NegativeCV')

    lwave        = cfg.getfloat(  'edge', 'WaveLength')     # laser wavelength [m]
    Constants.wo = Constants.h*Constants.c/lwave            # laser photon energy [eV]
    Constants.Eo = 0.25e-6 * Constants.me**2 / Constants.wo # (me^2/4wo) [MeV]
#    for i in range(10):
#      b, e = 'ex_from_%1d' % (i), 'ex_upto_%1d' % (i)
#      if cfg.has_option('edge', b) and cfg.has_option('edge', e):
#        self.exclude.append([ cfg.getfloat('edge', b), cfg.getfloat('edge', e) ])
#      else: break

    self.cc      = ROOT.TCanvas('cc','BEMS for VEPP-2000', 800, 600, 800, 600)
    self.const   = ROOT.TF1('const', '[0]')
    self.EdgeSimple = EdgeSimple()
    self.HPGeSpread = HPGeSpread()
    self.simple  = ROOT.TF1('simple', self.EdgeSimple, 0, 1, 7 );     self.simple.SetLineColor(ROOT.kRed)
    self.spreso  = ROOT.TF1('spreso', self.HPGeSpread, 0, 1, 3 );

    self.convol  = ROOT.TF1Convolution('simple', 'spreso', -1, 1); self.convol.SetNofPointsFFT(1000)
    self.comple  = ROOT.TF1('comple', self.convol, 1.0, 2.0, 10);  self.comple.SetNpx(1000)

    self.Legend  = ROOT.TLegend(0.6, 0.6, 0.95, 0.95, '', 'brNDC');
    self.HPGe    = Scale(scalefile, cfg_file, 'application')
    self.plots   = EMSResults(cfg_file, folder)

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def __del__(self):
    self.cc.cd(); self.cc.Clear()

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Go(self, UTB, UTE, hps, filechain, grating, Ie, Ee):

    if 0.15<Constants.wo<0.95:
      if grating in self.CO_Grating:
        print('%3d: %8.6f um' % (grating, 1.e+6*self.CO_Grating[grating]))
        Constants.wo = Constants.h*Constants.c/self.CO_Grating[grating]  # laser photon energy [eV]
        Constants.Eo = 0.25e-6 * Constants.me**2 / Constants.wo # (me^2/4wo) [MeV]
      else:
        print(grating)
        input()
        return False

    
    R = False; #VEPP2K_DB().GetRunInfo(filechain)
    if not R:
      R =  {'E':Ee, 'dE':0.0, 'I':Ie, 'dI':0.0, 'B':0.0, 'dB':0.0}
    if R:
      print(' ╔ VEPP2K conditions: ══════╤═════════════════════════╤═══════════════════════════╗')
      print(' ║  E = %7.2f ± %5.2f MeV │ Ie = %5.1f ± %5.1f mA   │  Bo = %7.5f ± %7.5f T ║' % (R['E'], R['dE'], R['I'], R['dI'], R['B'], R['dB']))
      print(' ╚══════════════════════════╧═════════════════════════╧═══════════════════════════╝\n')
      self.VEPP2K = R
      if abs(self.ETuner)<20.0:     self.Eo = self.ETuner + self.VEPP2K['E']        # [MeV]
      else:                         self.Eo = self.ETuner                           # [MeV]
      if self.Eo<150.: return
      if self.Radius:               self.Bo = 1.e+8*self.Eo/Constants.c/self.Radius # [T]
      elif 'B' not in self.VEPP2K:  self.Bo = 1.e+8*self.Eo/Constants.c/140.0       # [T]
      elif self.VEPP2K['B']<0.1:    self.Bo = 1.e+8*self.Eo/Constants.c/140.0       # [T]
      elif abs(self.BTuner)<0.1:    self.Bo = self.VEPP2K['B'] + self.BTuner        # [T]
      else:                         self.Bo = self.BTuner                           # [T]
    else: self.Eo = 0.0; return
    k    = 4.e+6*self.Eo*Constants.wo/Constants.me**2;                          # [eV*eV / eV**2]
    Wmax = 1.e+3*self.Eo*k/(1.+k)                                               # [keV]
    print("Fit initial beam energy: %7.2f MeV " % self.Eo)
    zero, gain = self.GetCalibrationResults(0.5*(UTB+UTE), Wmax)
    if gain:
      nbins    = hps.GetNbinsX();  hps.SetBins(nbins, zero, zero + gain * nbins)
      self.hps = hps.Clone(); self.hps.Rebin(self.Merger);  self.hps.GetXaxis().SetTitle('E_{#gamma} [keV]')

      # get rid of spikes:
#      self.EP.append(565)
      for spike in self.EP:
        lo = 1 + int((spike - zero - 5 * self.RL)/(gain*float(self.Merger)))
        hi = 1 + int((spike - zero + 7 * self.RR)/(gain*float(self.Merger)))
        for ch in range(lo, hi): self.hps.SetBinContent(ch, 0.0); self.hps.SetBinError(ch, 0.0)

      Wmax = self.fitEdgeSimple(Wmax, self.Ranger)
      self.NicePicture()
      if Wmax:
        Results = self.fitEdgeComple(Wmax, self.Ranger)
        self.NicePicture()
        self.comple.Draw('SAME');      self.Legend.Draw('SAME'); self.cc.Modified(); self.cc.Update()
        #self.hps.Draw('SAME'); self.cc.Modified(); self.cc.Update()
        if Results:
          self.plots.AddPoint(         UTB, UTE, Results)
          if self.SaveDB: Save_for_SND(UTB, UTE, Results)
          Results['T'] = [0.5*(UTB + UTE), 0.5*(UTE - UTB)]
          return Results
      else:
        print('No Edge?')
        self.cc.cd();      self.cc.Clear();    self.cc.SetGrid()
        self.hps.Draw(''); self.cc.Modified(); self.cc.Update()
        return False
    else:
      return False


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def NicePicture(self):
      if self.negate:
        self.cc.SetFillColor(923);            self.cc.SetFrameFillColor(1);          self.cc.SetFrameLineColor(0)
        self.hps.SetMarkerColor(ROOT.kCyan);  self.hps.SetLineColor(ROOT.kCyan)
        self.hps.GetXaxis().SetAxisColor(0);  self.hps.GetXaxis().SetTitleColor(0);  self.hps.GetXaxis().SetLabelColor(0)
        self.hps.GetYaxis().SetAxisColor(0);  self.hps.GetYaxis().SetTitleColor(0);  self.hps.GetYaxis().SetLabelColor(0)
        self.Legend.SetFillColor(923);        self.Legend.SetLineColor(0);           self.Legend.SetTextColor(0)
        self.simple.SetLineColor(ROOT.kOrange)
        self.comple.SetLineColor(ROOT.kPink)

      self.cc.cd(); self.cc.Clear();  self.cc.SetGrid(); self.hps.SetMarkerStyle(24)
      self.hps.Draw(''); self.simple.Draw('SAME');
      self.cc.Modified(); self.cc.Update()

      if self.negate:
        self.cc.FindObject('title').SetFillColor(923)
        self.cc.FindObject('title').SetTextColor(0)
      self.cc.Modified(); self.cc.Update()
#      self.cc.ls()
#      raw_input()

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def fitEdgeSimple(self,W,LBK):

    K = 2.e+3*W/Constants.me; Wmin = W/(1+K); E1, E2 = W - LBK*Wmin, W + LBK*Wmin
    self.hps.GetXaxis().SetRangeUser(E1-10, E2+10)
    self.const.SetRange(1.02*W, E2); self.hps.Fit('const','QRN');  B = self.const.GetParameter(0)
    self.const.SetRange(E1, 0.98*W); self.hps.Fit('const','QRN');  A = self.const.GetParameter(0) - B
    # Eb [MeV] | B [T] | Amplitude | Edge linear | Edge square | Background | Backg. Slope |
    params = np.fromiter([self.Eo,  self.Bo, A, 0.0, 0.0, B, 0.0], np.float)
    self.simple.SetParameters(params); self.simple.SetRange(E1, E2); self.simple.SetNpx(1000)
    self.simple.SetParLimits(2, 0.0, 1.e+6)
    if self.Radius: self.simple.FixParameter(1, self.Bo)
    if self.EdgeP2: self.simple.SetParLimits(4, 0.0, 0.01)
    else:           self.simple.FixParameter(4, 0.0)
    R = self.hps.Fit('simple','RSQN'); OK = False
    if self.Radius: self.simple.ReleaseParameter(1)
    if not self.EdgeP2: self.simple.ReleaseParameter(4)


    if not R.Status():
      E = fitParameters(self.simple)
      tilt = (1e+3*E['p'][3], 1e+3*E['e'][3], 1e+6*E['p'][4], 1e+6*E['e'][4])
      print(' ╔ Simple Edge Fit: ═════════════╤══════════════════════╤═════════════════════════╗')
      print(' ║ Range from %5.0f to %5.0f keV │ E_beam = %7.2f MeV │   W_max = %9.3f keV ║' % (E1, E2, self.Eo, W))
      print(' ╟───────────────────────────────┴────────┬─────────────┴─────────────────────────╢')
      print(' ║ Beam energy: %8.3f ± %5.3f [MeV]    │ Bending field: %6.4f ± %6.4f [T]    ║' % (E['p'][0], E['e'][0], E['p'][1], E['e'][1]))
      print(' ║ Background:  %8.0f ± %5.0f          │ Amplitude:   %8.0f ± %6.0f        ║' % (E['p'][5], E['e'][5], E['p'][2], E['e'][2]))
      print(' ║ Edge tilt pol1: %5.1f ± %5.1f [1/eV]   │ Edge tilt pol2: %5.0f ± %5.0f [1/eV^2]║' % tilt)
      print(' ╟────────────────────────────────────────┼───────────────────────────────────────╢')
      print(' ║         χ²/NDF = %5.1f/%3d             │         Probability: %8.6f         ║' % (R.Chi2(), R.Ndf(), R.Prob()))
      print(' ╚════════════════════════════════════════╧═══════════════════════════════════════╝\n')
      self.Legend.Clear();
      self.Legend.AddEntry(self.simple, '#chi^{2}/NDF = %5.1f/%3d  (Prob: %5.3f)' % (R.Chi2(), R.Ndf(), R.Prob()), 'l')
      if (E['p'][2]>self.MinAmp) and (E['e'][2]/E['p'][2]<0.5):
        k    = 4.e+6*E['p'][0]*Constants.wo/Constants.me**2; # [eV*eV / eV**2]
        return 1.e+3*E['p'][0]*k/(1.+k)                      # Wmax [keV]
    print('Simple fit: bad fit, bad amplitude, spread or χ²');  return 0.0

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def fitEdgeComple(self, W, LBK):
    K = 2.e+3*W/Constants.me; Wmin = W/(1+K);
    E1, E2 = W -     LBK*Wmin, W +     LBK*Wmin
    C1, C2 = W - 0.9*LBK*Wmin, W + 0.9*LBK*Wmin
    def FitP(H, P, q):
      convol = ROOT.TF1Convolution('simple', 'spreso', E1, E2); convol.SetNofPointsFFT(1000)
      comple = ROOT.TF1('comple', convol, C1, C2, 10);  comple.SetNpx(1000)
      comple.SetParameters(np.fromiter(P, np.float))
      if self.Radius: comple.FixParameter(1, self.Bo)
      if self.EdgeP2: comple.SetParLimits(4, 0.0, 0.01)
      else:           comple.FixParameter(4, 0.0)
      comple.FixParameter(7, self.RR)
      comple.FixParameter(8, self.RL)
      comple.SetParLimits(9,0.1,100.0)
      R = H.Fit('comple','RSN')
      if self.Radius: comple.ReleaseParameter(1)
      if not self.EdgeP2: comple.ReleaseParameter(4)
      comple.ReleaseParameter(7);       comple.ReleaseParameter(8)
      q.put([R.Status(), R.Chi2(), R.Ndf(), R.Prob(), fitParameters(comple)])


    self.convol.SetRange(E1, E2);            self.comple.SetRange(C1, C2)
#    self.simple.SetLineColor(ROOT.kBlue);    self.comple.SetLineColor(ROOT.kRed);
    P = fitParameters(self.simple)['p'];     P.extend([self.RR, self.RL, 1.0])
    self.comple.SetParameters(np.fromiter(P, np.float))
    self.comple.Draw('SAME');  self.cc.Modified(); self.cc.Update()

    q = Queue();    p = Process(target=FitP, args=(self.hps, P, q))
    p.start();      R = q.get()
    status, chi2, ndf, prob = R[0:4]
    pV = np.fromiter(R[4]['p'], np.float); self.comple.SetParameters(pV)
    pE = np.fromiter(R[4]['e'], np.float); self.comple.SetParErrors( pE)
    p.join()

    if not status:
      tilt     = (1e+3*pV[3], 1e+3*pE[3], 1e+6*pV[4], 1e+6*pE[4])
      k        =  4e+6*pV[0]*Constants.wo/Constants.me**2; # [eV*eV / eV**2]
      deriv    = (1.+k)**2/k/(2.+k) # dE/dWmax, apply scale correcion to the beam energy:
      BE, dBE  = pV[0]-1e-3*self.SC*deriv,  (pE[0]**2 + (1e-3*deriv*self.dSC)**2)**0.5 # Beam Energy [MeV]
      BF, dBF  = pV[1], pE[1]                                                          # Bending Field [T]
      BS, dBS  = pV[9]*deriv,  pE[9]*deriv                                             # Beam Spread [MeV]
      BR       = 1.e+8*pV[0]/Constants.c/BF;   dBR = BR * ((dBE/BE)**2 + (dBF/BF)**2)**0.5 # Beam Radius [cm]
      print(' ')
      print(' ╔ Convolution Fit: ═════════════╤══════════════════════╤═════════════════════════╗')
      print(' ║ Range from %5.0f to %5.0f keV │ σR = %4.2f ± %4.2f keV │ σL = %5.2f ± %5.2f keV  ║' % (E1, E2, self.RR, self.dRR, self.RL, self.dRL))
      print(' ╟───────────────────────────────┴────────┬─────────────┴─────────────────────────╢')
      print(' ║ Beam energy: %8.3f ± %5.3f [MeV]    │ Bending field: %6.4f ± %6.4f [T]    ║' % (BE, dBE, BF, dBF))
      print(' ║ σ from beam: %8.3f ± %5.3f [keV]    │ Beam spread:   %6.0f ± %6.0f [keV]  ║' % (pV[9], pE[9], BS, dBS))
      print(' ║ Edge tilt pol1: %5.1f ± %5.1f [1/eV]   │ Edge tilt pol2: %5.0f ± %5.0f [1/eV^2]║' % tilt)
      print(' ╟────────────────────────────────────────┼───────────────────────────────────────╢')
      print(' ║         χ²/NDF = %5.1f/%3d             │         Probability: %8.6f         ║' % (chi2, ndf, prob))
      print(' ╚════════════════════════════════════════╧═══════════════════════════════════════╝\n')
      self.Legend.AddEntry(self.comple, '#chi^{2}/NDF = %5.1f/%3d  (Prob: %5.3f)' % (chi2, ndf, prob), 'l')
      self.Legend.AddEntry(0, 'E_{beam} = %8.3f #pm %5.3f [MeV]'   % (BE, dBE), '')
      self.Legend.AddEntry(0, '#sigma_{E} = %6.0f #pm %4.0f [keV]' % (BS, dBS), '')
      self.Legend.AddEntry(0, 'R_{beam} = %6.2f #pm %5.2f [cm]'    % (BR, dBR), '')
      if (pE[0]/pV[0]<0.001) and prob>0.00001:
        return {'BE':[BE,dBE], 'BF':[BF,dBF], 'BS':[BS,dBS], 'BC':[self.VEPP2K['I'], self.VEPP2K['dI']]}
    else:
      print('Status: %d (χ²/NDF = %5.1f/%3d - probability: %8.6f)' % (status, chi2, ndf, prob))
    print('Convolution fit: bad fit (status %d), bad amplitude, spread or χ²' % status)
    return False


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def GetCalibrationResults(self, t, wmax):
    for attempt in range(3):
      Scale = self.HPGe.Get_Calibration(t, wmax)
      L = len(Scale['R'])
      if L: break
      else: print('Waiting 15s for calibration results...');   time.sleep(15)


    if L == 1: c = 0
    elif L >1:
      Q = [Scale['dW'][c] * Scale['dC'][c] * Scale['dR'][c] * Scale['dL'][c] for c in range(L)]
#      for c in range(L): Q.append(Scale['dW'][c] * Scale['dC'][c] * Scale['dR'][c] * Scale['dL'][c])
      if len(Q)==0:    print('No valid calibration was found'); return 0, 0
      c = Q.index(min(Q))
    else:              print('No valid calibration was found'); return 0, 0

    zero, gain         = Scale['Z'][c], Scale['G' ][c]

    self.dW            =                Scale['dW'][c] # linear calibration statistical error, keV
    self.SC, self.dSC  = Scale['C'][c], Scale['dC'][c] # PB-5  scale correction and its error, keV
    self.RR, self.dRR  = Scale['R'][c], Scale['dR'][c] # Right resolution sigma and its error, keV
    if self.Asymme:    self.RL, self.dRL  = Scale['L'][c], Scale['dL'][c] # Left  resolution sigma and its error, keV
    else:              self.RL, self.dRL  = Scale['R'][c], Scale['dR'][c] # Left  resolution sigma and its error, keV
    self.EP            = Scale['X'][c]                 # exclude peaks, keV
    print(' ╔ HPGe calibration: %15s ══════════════════╤══════════════════════════╗' % (self.HPGe.outfile.split('/')[-1]))
    print(' ║  W_max  = %9.3f keV  │ σR = %6.3f ± %5.3f keV  │ σL = %6.3f ± %5.3f keV  ║' % (wmax, self.RR, self.dRR, self.RL, self.dRL))
    print(' ╚══════════════════════════╧══════════════════════════╧══════════════════════════╝\n')
    return zero, gain

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class VEPP2K_DB:
#  cc = ROOT.TCanvas('cc','LOAD [Hz] vs I[mA]', 5, 5, 1000, 1200)
#  NAMEs = { 't':'t', 'Ee':'energy', 'Ie':'e_current', 'Bo':'nmr_field', 'L1':'load1', 'L2':'load2', 'L3':'load3'}
  NAMEs = { 't':'t', 'Ee':'energy', 'Ie':'e_current', 'Bo':'nmr_field'}
  def GetRunInfo(self,filechain):
    R, T = {'E':0.0, 'dE':0.0, 'I':0.0, 'dI':0.0, 'B':0.0, 'dB':0.0}, { k:[] for k in list(self.NAMEs.keys())}
    for el in filechain:
      fn = el.split('.')[0].replace('SPECTRA','LOADS') + '.sta'
      try:
        with open(fn,'rb') as fp: data = pickle.load(fp)
      except:
        print('sta read error'); return R

      OK, LREC = (len(list(data.keys()))==10), len(data['t'])
      if not OK: print('Bad file: %s' % fn); return False
      else:
        for k in list(data.keys()):
          L = len(data[k])
          if L+2 < LREC:
            print('Warning: LREC_t = %d, while LREC_%s = %d' % (LREC,k,L))
            data[k].extend([data[k][-1]]*(LREC-L))
          elif L > LREC:
            print('Error: LREC_t = %d, while LREC_%s = %d' % (LREC,k,L))
            return R
      for k,v in self.NAMEs.items(): T[k].extend(data[v])

    E = np.fromiter(T['Ee'], np.float);  R['E'], R['dE']  = E.mean(), E.std()
    B = np.fromiter(T['Bo'], np.float);  R['B'], R['dB']  = B.mean(), B.std()
#    I = np.fromiter(T['Ie'], np.float);  R['I'], R['dI']  = I.mean(), I.std()
    I = [i for i in T['Ie'] if i > 2.0]
    if len(I):
      I = np.fromiter(I, np.float);      R['I'], R['dI']  = I.mean(), I.std()
    else: R['I'], R['dI']  = 0.0, 0.0

#    n = len(I); dI = np.zeros(n)
#    L = np.fromiter(T['L1'], np.float) + np.fromiter(T['L2'], np.float) + np.fromiter(T['L3'], np.float)
#    L = np.fromiter(T['L2'], np.float)
#    dL = L**0.5
#    self.cc.cd()
#    self.cc.SetGrid()
#    LvsI = ROOT.TGraphErrors(n,I,L,dI,dL)
#    LvsI.Draw('AP')
#    self.cc.Modified(); self.cc.Update()
#    raw_input()
#    LvsI.Delete()

#    print ' ╔ VEPP2K conditions: ══════╤═════════════════════════╤═══════════════════════════╗'
#    print ' ║  E = %7.2f ± %5.2f MeV │ Ie = %5.1f ± %5.1f mA   │  Bo = %7.5f ± %7.5f T ║' % (R['E'], R['dE'], R['I'], R['dI'], R['B'], R['dB'])
#    print ' ╚══════════════════════════╧═════════════════════════╧═══════════════════════════╝\n'

#    if R['I'] < 5.0 : return False
    return R

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class EdgeSimple(Constants): # Ai integral (classical formula) and Klein-Nishina cross section
  ROOT.gROOT.LoadMacro(sys.argv[0].replace(sys.argv[0].split('/')[-1], 'vepp2k/airy_ai_int.C'))
  ROOT.airy_ai_int(np.array([0.0]), np.array([0.0])) # first call => initialization (could be up to 5 sec)
  iAiry = ROOT.TF1("af", ROOT.airy_ai_int, -55, 5, 1)
  # p[0] - Beam Energy    [MeV] | p[1] - Bending field  [T]   | p[2] - Amplitude
  # p[3] - Linear Slope at Wmax | p[4] - Square Slope at Wmax | p[5] and p[6] - Background and its Slope
  def __call__(self, x, p):
    u = x[0]/(1.e+3*p[0]-x[0])       # u     [keV / keV]
    k = p[0]      / Constants.Eo     # kappa [MeV / MeV]
    X = p[0]*p[1] * Constants.RemB   # Xi, free field if comment next line
    D = x[0] - p[0]*k/(1.+k)*1.e+3   # (W - Wmax) [keV]

    I = self.iAiry(ROOT.TMath.Power((u/X),0.6666666666666666) * (1.-k/u))

    return I * p[2]*(1.0 + p[3]*D + p[4]*D**2) + p[5] + p[6]*D

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class HPGeSpread: # Analitical convolution of the HPGe bifurcated Gaussian with normal Gaussian from beam spread
  No = 1./ROOT.TMath.Sqrt(ROOT.TMath.Pi())
# p[0] - HPGe sigma_R [keV] | p[1] - HPGe sigma_L [keV]  | p[2] Beam spread smearing [keV]
  def __call__(self, x, p):
    SR = 0.5/(p[0]*p[0] + p[2]*p[2]); RR = ROOT.TMath.Sqrt(SR)
    SL = 0.5/(p[1]*p[1] + p[2]*p[2]); RL = ROOT.TMath.Sqrt(SL)
    R = ROOT.TMath.Exp(-x[0]*x[0]*SR) * ROOT.TMath.Erfc(-x[0] * RR * p[0]/p[2]) * RR / (1. + p[1]/p[0])
    L = ROOT.TMath.Exp(-x[0]*x[0]*SL) * ROOT.TMath.Erfc( x[0] * RL * p[1]/p[2]) * RL / (1. + p[0]/p[1])
    return self.No * (L + R)
"""
# ЭТО НЕВЕСТЬ КАК ПОЛУЧЕННАЯ НЕПРАВИЛЬНАЯ СВЕРТКА!!! (СОХРАНЕНО ДЛЯ ИСТОРИИ)
class HPGeSpread: # Analitical convolution of the HPGe bifurcated Gaussian with normal Gaussian from beam spread
# p[0] - HPGe sigma_R [keV] | p[1] - HPGe sigma_L [keV]  | p[2] Beam spread smearing [keV]
  C = 1./ROOT.TMath.Sqrt(2*ROOT.TMath.Pi())
  def __call__(self, x, p):
    R, L = p[0]*p[0] + p[1]*p[1], p[0]*p[0] + p[2]*p[2]
    R = p[0] * ROOT.TMath.Exp(-0.5*x[0]**2/R) * ROOT.TMath.Erfc(-x[0]*p[0]/(p[2]*ROOT.TMath.Sqrt(R)))
    L = p[1] * ROOT.TMath.Exp(-0.5*x[0]**2/L) * ROOT.TMath.Erfc( x[0]*p[1]/(p[2]*ROOT.TMath.Sqrt(L)))
    return self.C*(R+L)/(p[0]+p[1])
# ЭТО НЕВЕСТЬ ОТКУДА ПОЛУЧЕННАЯ НЕПРАВИЛЬНАЯ СВЕРТКА!!!
"""
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class EMSResults:
  # Beam Energy, Bending Field, Beam energy Spread, Beam Current
  BE, BF, BS, BC = ROOT.TGraphErrors(), ROOT.TGraphErrors(), ROOT.TGraphErrors(), ROOT.TGraphErrors()
  SI = ROOT.TProfile('SI','spread vs current', 30, 0, 150, 'g')
  headline = ('t, s', 'dt, s', '  E, MeV', 'dE, MeV', 'S, keV',  'dS, keV', '  B, T', 'dB, T', 'I, mA', 'dI, mA')
  headline = '# %8s  %5s  %8s  %5s  %5s  %5s  %6s  %6s  %5s  %6s\n' % headline
  dataline = '%10d  %5d  %8.3f  %7.3f  %6d  %7d  %6.4f  %6.4f  %5.1f  %6.1f\n'

  def __init__(self, cfg_file, folder):

    cfg = configparser.ConfigParser()
    cfg.read(cfg_file)
    self.roofile = folder + cfg.get('scan', 'file')
    self.datfile = folder + cfg.get('scan', 'file') +'.txt'
    self.espread = folder + 'SvsI.pdf'
    self.SI.SetTitle(folder)
    with open(self.datfile,'a') as f:
#      f.write('# R=140cm, assym = false, edgepoly2 = false\n')
      f.write(self.headline)

  def __del__(self):
    if hasattr(self, 'rc'): self.rc.cd(); self.rc.Clear()
    if hasattr(self, 'ic'): self.ic.cd(); self.ic.Clear()

  def SaveGraphs(self):
    fp = ROOT.TFile(self.roofile, 'RECREATE')
    if fp.IsOpen():
      fp.WriteObject(self.BE, 'BE')
      fp.WriteObject(self.BF, 'BF')
      fp.WriteObject(self.BS, 'BS')
      fp.WriteObject(self.BC, 'BC')
      fp.Close()

  def ReadGraphs(self):
    if os.path.isfile(self.roofile):
      fp = ROOT.TFile(self.roofile, 'READ')
      self.BE = fp.BE; n = self.BE.GetN()
      self.BF = fp.BF; m = self.BF.GetN()
      self.BS = fp.BS; k = self.BS.GetN()
      self.BC = fp.BC; l = self.BC.GetN()
      fp.Close()
      return n, m, k, l
    else:
      return 0, 0, 0, 0

  def AddPoint(self, UTB, UTE, R):
    n, m, k, l = self.ReadGraphs()
    t, dt = 0.5*(UTB + UTE), 0.5*(UTE - UTB)
    e, de = R['BE']; self.BE.SetPoint(n, t, e); self.BE.SetPointError(n, dt, de)
    b, db = R['BF']; self.BF.SetPoint(m, t, b); self.BF.SetPointError(m, dt, db)
    s, ds = R['BS']; self.BS.SetPoint(k, t, s); self.BS.SetPointError(k, dt, ds)
    c, dc = R['BC']; self.BC.SetPoint(l, t, c); self.BC.SetPointError(l, dt, dc)
    self.SaveGraphs()
    with open(self.datfile,'a') as f:  f.write(self.dataline % (t,dt, e,de, s,ds, b,db, c,dc))

  def ShowRunInfo(self):
    n, m, k, l = self.ReadGraphs()
    for g in [self.BE, self.BF, self.BS, self.BC]:
      g.GetXaxis().SetTimeDisplay(1); g.GetXaxis().SetTimeFormat('#splitline{%b%d}{%H:%M}%F1970-01-01 00:00:00')
      g.GetXaxis().SetTitle('time');  g.GetXaxis().SetLabelOffset(0.02)
      g.SetMarkerColor(ROOT.kRed); g.SetLineColor(ROOT.kRed)
      g.SetMarkerStyle(20); g.SetMarkerSize(1.25); g.GetYaxis().SetDecimals()
    self.rc = ROOT.TCanvas('rc','BEMS results for VEPP-2000', 0, 0, 1600, 1200)
    self.rc.Divide(1,3)

    self.rc.cd(1); self.rc.GetPad(1).SetGrid();  self.BE.Fit('pol0'); self.BE.Draw('AP');
    self.BE.GetXaxis().SetTitle('time');         self.BE.GetYaxis().SetTitle('Beam energy [MeV]')
    self.rc.cd(2); self.rc.GetPad(2).SetGrid();  self.BS.Fit('pol0'); self.BS.Draw('AP');
    self.BS.GetXaxis().SetTitle('time');         self.BS.GetYaxis().SetTitle('Beam energy spread [keV]')
    self.rc.cd(3); self.rc.GetPad(3).SetGrid();  self.BC.Fit('pol0'); self.BC.Draw('AP');
    self.BC.GetXaxis().SetTitle('time');         self.BC.GetYaxis().SetTitle('Beam current [mA]')
#    self.rc.cd(4); self.rc.GetPad(4).SetGrid();  self.BF.Fit('pol0'); self.BF.Draw('AP');
#    self.BF.GetXaxis().SetTitle('time');         self.BF.GetYaxis().SetTitle('Bending field [T]')

    self.rc.Modified();  self.rc.Update();  self.rc.SaveAs(self.roofile + '.pdf')

  def EnergySpread(self):
    n, m, k, l = self.ReadGraphs()
    xmin, xmax = 100.0, 0.0
    ymin, ymax = 10.0,  0.0
    for p in range(n):
      t = self.BE.GetPointX(p) 
      e = self.BE.GetPointY(p); de = self.BE.GetErrorY(p)
      s = self.BS.GetPointY(p); ds = self.BS.GetErrorY(p)
      i = self.BC.GetPointY(p); di = self.BC.GetErrorY(p)
      if (i>0.0) and (s>0.0) and abs(di/i)<0.5 and abs(ds/s)<0.5:
        S, dS = 10.*s/e, 10.*ds/e
        self.SI.Fill(i,S,1./dS**2)
        xmin = min(xmin, i-di);  xmax = max(xmax, i+di)
        ymin = min(ymin, S-dS);  ymax = max(ymax, S+dS)
    self.SI.GetXaxis().SetTitle('Current [mA]');                  #self.SI.GetXaxis().SetLabelOffset(0.03)
    self.SI.GetXaxis().SetRangeUser(xmin, xmax)
    self.SI.GetYaxis().SetTitle('#sigma_{E}/E #upoint 10^{-4}');  #self.SI.GetYaxis().SetLabelOffset(0.02)
    self.SI.GetYaxis().SetRangeUser(ymin, ymax)
    self.SI.SetMarkerColor(ROOT.kRed); self.SI.SetLineColor(ROOT.kRed)
    self.SI.SetMarkerStyle(20); self.SI.SetMarkerSize(1.25); self.SI.GetYaxis().SetDecimals()
    self.sc = ROOT.TCanvas('sc','BEMS results for VEPP-2000', 0, 0, 1200, 800)
    self.sc.cd(), self.sc.SetGrid()
    self.SI.Draw()
    self.sc.Modified();  self.sc.Update();  self.sc.SaveAs(self.espread)

class Points_Splitter:
  ENERGY_CHANGE = 5.0
  SMALL_CURRENT = 5.0
  E_MIN,  E_MAX = 100., 1000.
  success = []
  failure = []

  def Go(self, flist, todo):
    from hpge import DataFile
    T = DataFile()
    T.ReadHat(flist[0], 100)
    BDATE = str(int(T.HAT['Date']))
    R = VEPP2K_DB().GetRunInfo([flist[0]])
    if R: E0 = R['E']
    else: E0 = 0.0

    while len(flist):
      f = flist.pop(0)
      T.ReadHat(f, 100)
      t, dt = 0.5*(T.utb + T.ute), 0.5*(T.ute - T.utb)
      R = VEPP2K_DB().GetRunInfo([f])
      if T and R:
        E, dE, I, dI = R['E'], R['dE'], R['I'], R['dI']
        SUCCESS = self.E_MIN < E < self.E_MAX
        SUCCESS = dt > 100.  and    dE  < self.ENERGY_CHANGE
        SUCCESS = SUCCESS and abs(E-E0) < self.ENERGY_CHANGE
        SUCCESS = SUCCESS and         I > self.SMALL_CURRENT
        SUCCESS = SUCCESS and len(flist)

        FAILURE  = E < self.E_MIN or  I < self.SMALL_CURRENT
        FAILURE  = FAILURE        or dE > self.ENERGY_CHANGE
        FAILURE  = FAILURE        and len(flist)

        PFINISH   = self.E_MIN < E < self.E_MAX
        PFINISH   = PFINISH and abs(E-E0) > self.ENERGY_CHANGE
        PFINISH   = PFINISH or len(flist)==0

        if   SUCCESS:    self.success.append(f)
        elif FAILURE:    self.failure.append(f)
        elif PFINISH:
          EDATE  = str(int(T.HAT['Date']))
          folder = '%8s-%8s~%3.0fMeV' % (BDATE, EDATE, E0)
          print(folder)
          print('number of good files: ', len(self.success))
          print('number of  bad files: ', len(self.failure))

          if not os.path.exists(folder): os.makedirs(folder)
          with open('%s/success.list' % (folder), 'w') as fp:
            for el in self.success: fp.write(el+'\n')
          with open('%s/failure.list' % (folder), 'w') as fp:
            for el in self.failure: fp.write(el+'\n')

          E0 = R['E']
          BDATE = EDATE
          flist.insert(0,f)
          del self.success[:]
          del self.failure[:]

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def Save_for_SND(UTB, UTE, R):
  t, dt = UTB, UTE - UTB;  E, dE = R['BE'];  B, dB = R['BF'];  S, dS = R['BS']
  PL = {'/EMS/DT':dt, '/EMS/E':E, '/EMS/DE':dE, '/EMS/B':B, '/EMS/DB':dB, '/EMS/S':S, '/EMS/DS':dS}
  print('WRITING TO SND ', end=' ')
  DB = Storage()
  for k,v in PL.items():
    DB.register(k); DB[k]=v; print('.', end=' ')
  try:
    DB.insert(t)
  except:
    try:
      DB.replace(t)
    except:
      print('Warning! Can not write nor replace data!')
      return
  print(' DONE')
  return

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def fitParameters(fitf):
  n, p, e = fitf.GetNumberFreeParameters(), [], []
  for i in range(n):
    p.append(fitf.GetParameter(i))
    e.append(fitf.GetParError(i))
  return {'p':p,'e':e}

