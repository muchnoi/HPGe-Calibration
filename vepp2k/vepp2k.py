# -*- coding: utf-8 -*-
import ROOT, numpy, cPickle, os, time, ConfigParser
from scale.isotopes  import Isotopes
sql = True
try: 
  from temdbWO import Storage
except ImportError:
  sql = False 
  print 'database is unavailable'

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

  def __init__(self, scalefile, cfg_file):
    self.exclude = []
    cfg = ConfigParser.ConfigParser(); cfg.read(cfg_file)
    self.MinAmp  = cfg.getfloat('edge', 'MinAmpEdge')
    self.Ranger  = cfg.getfloat('edge', 'EdgeRanger')
    self.Merger  =   cfg.getint('edge', 'BinsMerger')
    self.ETuner  = cfg.getfloat('edge', 'EveppTuner')
    self.BTuner  = cfg.getfloat('edge', 'BveppTuner')
    self.Radius  = cfg.getfloat('edge', 'VeppRadius')
    self.Asymme  = cfg.getfloat('edge', 'Asymmetry')
    self.SaveDB  = sql and ('True' in cfg.get('edge', 'SaveForSND'))
    lwave        = cfg.getfloat('edge', 'WaveLength')       # laser wavelength [m]
    Constants.wo = Constants.h*Constants.c/lwave            # laser photon energy [eV]
    Constants.Eo = 0.25e-6 * Constants.me**2 / Constants.wo # (me^2/4wo) [MeV]
    for i in range(10):
      b, e = 'ex_from_%1d' % (i), 'ex_upto_%1d' % (i)
      if cfg.has_option('edge', b) and cfg.has_option('edge', e):
        self.exclude.append([ cfg.getfloat('edge', b), cfg.getfloat('edge', e) ])
      else: break
    self.cc      = ROOT.TCanvas('cc','BEMS for VEPP-2000', 800, 600, 800, 600)
    self.const   = ROOT.TF1('const', '[0]')
    self.simple  = ROOT.TF1('simple', EdgeSimple(), 0, 1, 7 ); self.simple.SetLineColor(ROOT.kRed)
    self.spreso  = ROOT.TF1('spreso', HPGeSpread(), 0, 1, 3 ); 
    self.Legend  = ROOT.TLegend(0.6, 0.6, 0.95, 0.95, '', 'brNDC'); 
    self.HPGe    = Isotopes(scalefile, cfg_file, 'application')
    
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def __del__(self):    
    self.cc.cd(); self.cc.Clear() 

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Go(self, UTB, UTE, hps, filechain):
    self.hps = hps.Clone(); nbins = hps.GetNbinsX(); self.hps.GetXaxis().SetTitle('E_{#gamma}, keV')
    self.VEPP2K = VEPP2K_DB().GetRunInfo(filechain)
    if self.VEPP2K:
      if abs(self.ETuner)<10.0: self.Eo = self.ETuner + self.VEPP2K['E']        # [MeV]
      else:                     self.Eo = self.ETuner                           # [MeV]
      if self.Radius:           self.Bo = 1.e+8*self.Eo/Constants.c/self.Radius # [T] 
      else:                     self.Bo = self.VEPP2K['B'] + self.BTuner        # [T]
      k    = 4.e+6*self.Eo*Constants.wo/Constants.me**2;                        # [eV*eV / eV**2]
      Wmax = 1.e+3*self.Eo*k/(1.+k)                                             # [keV]
    else: 
      Wmax = 0.0
    zero, gain = self.GetCalibrationResults(0.5*(UTB+UTE), Wmax)
    if gain:  
      self.hps.SetBins(nbins, zero, zero + gain * nbins); self.hps.Rebin(self.Merger)
      if Wmax:
        for spike in self.exclude:
          lo = 1 + int((spike[0]-zero)/(gain*float(self.Merger)))
          hi = 1 + int((spike[1]-zero)/(gain*float(self.Merger)))
          for ch in range(lo, hi): 
            self.hps.SetBinContent(ch, 0.0); self.hps.SetBinError(ch, 0.0)
        Wmax = self.fitEdgeSimple(Wmax, self.Ranger)
        if Wmax:
          Results = self.fitEdgeComple(Wmax, self.Ranger)
          if Results and self.SaveDB: Save_for_SND(UTB, UTE, Results)
      else:
        print 'No Edge?'
        self.cc.cd();      self.cc.Clear();    self.cc.SetGrid()
        self.hps.Draw(''); self.cc.Modified(); self.cc.Update()
        
           
    
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def fitEdgeSimple(self,W,LBK):
    
    K = 2.e+3*W/Constants.me; Wmin = W/(1+K); E1, E2 = W - LBK*Wmin, W + LBK*Wmin
    self.const.SetRange(1.02*W, E2); self.hps.Fit('const','QRN');  B = self.const.GetParameter(0)
    self.const.SetRange(E1, 0.98*W); self.hps.Fit('const','QRN');  A = self.const.GetParameter(0) - B
    # Eb [MeV] | B [T] | Amplitude | Edge linear | Edge square | Background | Backg. Slope |  
    params = numpy.fromiter([self.Eo,  self.Bo, A, 0.0, 0.0, B, 0.0], numpy.float)
    self.simple.SetParameters(params); self.simple.SetRange(E1, E2); self.simple.SetNpx(1000)

    if self.Radius: self.simple.FixParameter(1,self.Bo)
    R = self.hps.Fit('simple','RSQN'); OK = False
    if self.Radius: self.simple.ReleaseParameter(1)

    if not R.Status():
      E = fitParameters(self.simple)
      tilt = (1e+3*E['p'][3], 1e+3*E['e'][3], 1e+6*E['p'][4], 1e+6*E['e'][4])
      print ' ╔ Simple Edge Fit: ═════════════╤══════════════════════╤═════════════════════════╗' 
      print ' ║ Range from %5.0f to %5.0f keV │ E_beam = %7.2f MeV │   W_max = %9.3f keV ║' % (E1, E2, self.Eo, W)
      print ' ╟───────────────────────────────┴────────┬─────────────┴─────────────────────────╢'
      print ' ║ Beam energy: %8.3f ± %5.3f [MeV]    │ Bending field: %6.4f ± %6.4f [T]    ║' % (E['p'][0], E['e'][0], E['p'][1], E['e'][1])
      print ' ║ Background:  %8.0f ± %5.0f          │ Amplitude:   %8.0f ± %6.0f        ║' % (E['p'][5], E['e'][5], E['p'][2], E['e'][2])
      print ' ║ Edge tilt pol1: %5.3f ± %5.3f [1/eV]   │ Edge tilt pol2: %5.2f ± %5.2f [1/eV^2]║' % tilt
      print ' ╟────────────────────────────────────────┼───────────────────────────────────────╢'
      print ' ║         χ2/NDF = %5.1f/%3d             │         Probability: %8.6f         ║' % (R.Chi2(), R.Ndf(), R.Prob())
      print ' ╚════════════════════════════════════════╧═══════════════════════════════════════╝\n'
      OK = (E['p'][2]>self.MinAmp) and (E['e'][2]/E['p'][2]<0.5) 
    
    self.cc.cd(); self.cc.Clear();  self.cc.SetGrid(); self.hps.SetMarkerStyle(20)
    self.hps.Draw(''); self.simple.Draw('SAME'); self.hps.GetXaxis().SetRangeUser(E1, E2)
    self.cc.Modified(); self.cc.Update()

    if OK: 
      k    = 4.e+6*E['p'][0]*Constants.wo/Constants.me**2; # [eV*eV / eV**2]
      return 1.e+3*E['p'][0]*k/(1.+k)                      # Wmax [keV]
#    if offline: raw_input('Bad Fit?')
    print 'Simple fit: bad fit, bad amplitude, spread or χ2';  return 0.0

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def fitEdgeComple(self,W,LBK):
    
    K = 2.e+3*W/Constants.me; Wmin = W/(1+K); E1, E2 = W - LBK*Wmin, W + LBK*Wmin
    cnvl = ROOT.TF1Convolution('simple', 'spreso'); cnvl.SetNofPointsFFT(1000)
    self.comple = ROOT.TF1('comple',cnvl,E1,E2,10);   self.comple.SetNpx(1000)
    self.simple.SetLineColor(ROOT.kBlue);  self.comple.SetLineColor(ROOT.kRed)
    p = fitParameters(self.simple)['p'];  p[2] *= 0.6;  p[5] *= 0.6;  p.extend([self.RR, self.RL, 1.0])
    self.comple.SetParameters(numpy.fromiter(p, numpy.float))
    self.comple.Draw('SAME');  self.cc.Modified(); self.cc.Update()

    if self.Radius: self.comple.FixParameter(1,self.Bo)
    self.comple.FixParameter(7,self.RR);  self.comple.FixParameter(8,self.RL)
    R = self.hps.Fit('comple','RSQN');    OK = False
    if self.Radius: self.comple.ReleaseParameter(1)
    self.comple.ReleaseParameter(7);      self.comple.ReleaseParameter(8)

    if not R.Status():
      E = fitParameters(self.comple)
      tilt     = (1e+3*E['p'][3], 1e+3*E['e'][3], 1e+6*E['p'][4], 1e+6*E['e'][4])
      k        =  4e+6*E['p'][0]*Constants.wo/Constants.me**2; # [eV*eV / eV**2]
      deriv    = (1.+k)**2/k/(2.+k) # dE/dWmax, apply scale correcion to the beam energy:
      BE, dBE  = E['p'][0]-1e-3*self.SC*deriv,  (E['e'][0]**2 + (1e-3*deriv*self.dSC)**2)**0.5 # Beam Energy [MeV]
      BF, dBF  = E['p'][1], E['e'][1]                                                          # Bending Field [T]  
      BS, dBS  = E['p'][9]*deriv,  E['e'][9]*deriv                                             # Beam Spread [MeV]
      BR       = 1.e+8*E['p'][0]/Constants.c/BF;   dBR = BR * ((dBE/BE)**2 + (dBF/BF)**2)**0.5 # Beam Radius [cm] 
      print ' ╔ Convolution Fit: ═════════════╤══════════════════════╤═════════════════════════╗' 
      print ' ║ Range from %5.0f to %5.0f keV │ σR = %4.2f ± %4.2f keV │ σL = %5.2f ± %5.2f keV  ║' % (E1, E2, self.RR, self.dRR, self.RL, self.dRL)
      print ' ╟───────────────────────────────┴────────┬─────────────┴─────────────────────────╢'
      print ' ║ Beam energy: %8.3f ± %5.3f [MeV]    │ Bending field: %6.4f ± %6.4f [T]    ║' % (BE, dBE, BF, dBF)
      print ' ║ σ from beam: %8.3f ± %5.3f [keV]    │ Beam spread:   %6.0f ± %6.0f [keV]  ║' % (E['p'][9], E['e'][9], BS, dBS)
      print ' ║ Edge tilt pol1: %5.3f ± %5.3f [1/eV]   │ Edge tilt pol2: %5.2f ± %5.2f [1/eV^2]║' % tilt
      print ' ╟────────────────────────────────────────┼───────────────────────────────────────╢'
      print ' ║         χ2/NDF = %5.1f/%3d             │         Probability: %8.6f         ║' % (R.Chi2(), R.Ndf(), R.Prob())
      print ' ╚════════════════════════════════════════╧═══════════════════════════════════════╝\n'
      OK = (E['e'][0]/E['p'][0]<0.001) and (E['e'][1]/E['p'][1]<0.1) 
      self.Legend.Clear(); self.Legend.SetHeader('#chi^{2}/NDF = %5.1f/%3d  (Probability %5.3f)' % (R.Chi2(), R.Ndf(), R.Prob()))
      self.Legend.AddEntry(self.comple, 'E_{beam} = %8.3f #pm %5.3f [MeV]'   % (BE, dBE), 'l')
      self.Legend.AddEntry(self.comple, '#sigma_{E} = %6.0f #pm %4.0f [keV]' % (BS, dBS), 'l')
      self.Legend.AddEntry(self.comple, 'R_{beam} = %6.2f #pm %5.2f [cm]'    % (BR, dBR), 'l')
      self.Legend.Draw('SAME')
    self.comple.Draw('SAME');  self.cc.Modified(); self.cc.Update()
    if OK: return {'BE':[BE,dBE], 'BF':[BF,dBF], 'BS':[BS,dBS]}
    else:  print 'Complex fit: bad fit, bad amplitude, spread or χ2';  return False


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def GetCalibrationResults(self, t, wmax):
    for attempt in xrange(3):
      Scale = self.HPGe.Get_Calibration(t, wmax)
      if len(Scale): break
      else:
        print 'Waiting for scale calibration results...';   time.sleep(20) 
    if not len(Scale):
      print 'No valid calibration was found'; return 0, 0
#    Quality = []
#    for c in range(len(Scale['T'])):
#      Quality.append(Scale['dW'][c] * Scale['dC'][c] * Scale['dR'][c] * Scale['dL'][c])
#    if len(Quality)==0:   
#      print 'No valid calibration was found'; return 0,0
#    c = Quality.index(min(Quality))
    c = 0
    zero, gain         = Scale['Z'][c], Scale['G' ][c]
    self.dW            =                Scale['dW'][c] # linear calibration statistical error, keV
    self.SC, self.dSC  = Scale['C'][c], Scale['dC'][c] # PB-5  scale correction and its error, keV
    self.RR, self.dRR  = Scale['R'][c], Scale['dR'][c] # Right resolution sigma and its error, keV
    self.RL, self.dRL  = Scale['L'][c], Scale['dL'][c] # Left  resolution sigma and its error, keV
#    print Scale['N'][c]
    print ' ╔ HPGe calibration: %15s ══════════════════╤══════════════════════════╗' % (self.HPGe.outfile)
    print ' ║  W_max  = %9.3f keV  │ σR = %6.3f ± %5.3f keV  │ σL = %6.3f ± %5.3f keV  ║' % (wmax, self.RR, self.dRR, self.RL, self.dRL)
    print ' ╚════════════════════════╧════════════════════════════╧══════════════════════════╝\n'
    return zero, gain  



# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class VEPP2K_DB:
  NAMEs = { 't':'t', 'Ee':'energy', 'Ie':'e_current', 'Bo':'nmr_field'}
  def GetRunInfo(self,filechain):
    R, T = {'E':0, 'dE':0.0, 'I':0.0, 'dI':0.0}, { k:[] for k in self.NAMEs.keys()}
    for el in filechain:
      fn = el.split('.')[0].replace('SPECTRA','LOADS') + '.sta'
      try:   
        with open(fn,'rb') as fp: data = cPickle.load(fp)
      except: 
        print 'sta read error'; return R
      
      OK, LREC = (len(data.keys())==10), len(data['t'])
      if not OK: print 'Bad file: %s' % fn; return R
      else: 
        for k in data.keys():
          L = len(data[k])
          if L < LREC:
            print 'Warning: LREC_t = %d, while LREC_%s = %d' % (LREC,k,L)
            data[k].extend([data[k][-1]]*(LREC-L))
          elif L > LREC:
            print 'Error: LREC_t = %d, while LREC_%s = %d' % (LREC,k,L)
            return R
      for k,v in self.NAMEs.iteritems(): T[k].extend(data[v])
        
    A = numpy.fromiter(T['Ee'], numpy.float);  R['E'], R['dE']  = A.mean(), A.std()
    A = numpy.fromiter(T['Ie'], numpy.float);  R['I'], R['dI']  = A.mean(), A.std()
    A = numpy.fromiter(T['Bo'], numpy.float);  R['B'], R['dB']  = A.mean(), A.std()
    if R['E']<100. or R['I'] < 5.0 or R['B'] <1.0: return False 
    print ' ╔ VEPP2K conditions: ══════╤═════════════════════════╤═══════════════════════════╗' 
    print ' ║  E = %6.3f ± %5.3f MeV │ Ie = %5.1f ± %5.1f mA   │  Bo = %7.5f ± %7.5f T ║' % (R['E'], R['dE'], R['I'], R['dI'], R['B'], R['dB'])
    print ' ╚══════════════════════════╧═════════════════════════╧═══════════════════════════╝\n'
    return R

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class EdgeSimple(Constants): # Ai integral (classical formula) and Klein-Nishina cross section
  ROOT.airy_ai_int(numpy.array([0.0]), numpy.array([0.0])) # first call => initialization (could be up to 5 sec) 
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
# p[0] - HPGe sigma_R [keV] | p[1] - HPGe sigma_L [keV]  | p[2] Beam spread smearing [keV] 
  C = 1./ROOT.TMath.Sqrt(2*ROOT.TMath.Pi())
  def __call__(self, x, p):
    R, L = p[0]*p[0] + p[1]*p[1], p[0]*p[0]+ p[2]*p[2]
    R = p[0] * ROOT.TMath.Exp(-0.5*x[0]**2/R) * ROOT.TMath.Erfc(-x[0]*p[0]/(p[2]*ROOT.TMath.Sqrt(R)))
    L = p[1] * ROOT.TMath.Exp(-0.5*x[0]**2/L) * ROOT.TMath.Erfc( x[0]*p[1]/(p[2]*ROOT.TMath.Sqrt(L)))
    return self.C*(R+L)/(p[0]+p[1])

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
"""
class TheResults(cfg_file):
  cfg = ConfigParser.ConfigParser()
  cfg.read(cfg_file)
  roofile = cfg.get('edge', 'file')
  
  def AddPoint(self, UTB, UTE, R):
    pass

  def GetGraph(self, UTB, UTE, R):
    f, MERE  = ROOT.TFile(self.roofile), ROOT.TList()
    f.GetListOfKeys()
    for el in [el.GetName() for el in f.GetListOfKeys()]:
      pass

  def DoReview(self): 
    pass
"""    


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def Save_for_SND(UTB, UTE, R):
  t, dt = UTB, UTE - UTB;  E, dE = R['BE'];  B, dB = R['BF'];  S, dS = R['BS']
  PL = {'/EMS/DT':dt, '/EMS/E':E, '/EMS/DE':dE, '/EMS/B':B, '/EMS/DB':dB, '/EMS/S':S, '/EMS/DS':dS}
  print 'WRITING TO SND ',
  DB = Storage()
  for k,v in PL.iteritems(): 
    DB.register(k); DB[k]=v; print '.', 
  try:
    DB.insert(t)
  except:
    try:
      DB.replace(t)
    except: 
      print 'Warning! Can not write nor replace data!'
      return
  print ' DONE'
  return

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def fitParameters(fitf):
  n, p, e = fitf.GetNumberFreeParameters(), [], []
  for i in range(n):  
    p.append(fitf.GetParameter(i))
    e.append(fitf.GetParError(i))
  return {'p':p,'e':e}

