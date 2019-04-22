# -*- coding: utf-8 -*-

import ROOT, cPickle, os, time, ConfigParser
import numpy as np
from scale.scale import Scale
from hpge import DataFile, Histogram
ssh = True
try:
  import paramiko
except ImportError:
  ssh = False
  print 'paramiko is unavailable'


class EDGE:
  me       = 0.510998910e+3   # electron rest mass, keV
  h        = 4.135667662e-15  # Plank conatsant, eV*s
  c        = 299792458        # speed of light, m/s
  hc       = h*c

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def __init__(self, scalefile, cfg_file):
    cfg = ConfigParser.ConfigParser(); cfg.read(cfg_file)
    self.wo      = 1.e-3*self.hc/cfg.getfloat('edge', 'WaveLength') # laser photon energy [keV]
    self.MinAmp  =               cfg.getfloat('edge', 'MinAmpEdge')
    self.Ranger  =               cfg.getfloat('edge', 'EdgeRanger')
    self.Merger  =                 cfg.getint('edge', 'BinsMerger')
    self.ETuner  =               cfg.getfloat('edge', 'EbepcTuner')
    self.Asymme  =               cfg.getfloat('edge', 'Asymmetry')
    self.SaveDB  =     ssh and cfg.getboolean('edge', 'SaveForDB')

    self.cc      = ROOT.TCanvas('cc','BEMS for BEPC-II', 800, 600, 800, 600); # self.cc.Divide(1,2)
    self.const   = ROOT.TF1('const', '[0]')
    self.simple  = ROOT.TF1('simple', EdgeSimple(), 0, 1, 7); self.simple.SetLineColor(ROOT.kRed)
    self.comple  = ROOT.TF1('comple', EdgeComple(), 0, 1, 9); self.comple.SetLineColor(ROOT.kAzure)
    self.Lg1     = ROOT.TLegend(0.55, 0.71, 0.98, 0.91, '', 'brNDC');
    self.Lg2     = ROOT.TLegend(0.55, 0.49, 0.98, 0.69, '', 'brNDC');
    self.HPGe    = Scale(scalefile, cfg_file, 'application')

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Go(self, UTB, UTE, hps, filechain, grating):
    if   filechain == [e for e in filechain if 'elec' in e]: self.lepton = 'electron'
    elif filechain == [e for e in filechain if 'posi' in e]: self.lepton = 'positron'
    else: print 'Not a valid measurement: mixed spectrum'; return False
    self.BEPC = BEMS_DB().GetRunInfo(filechain)
    if self.BEPC:
      Wmax = self.Wmax_expect()
      zero, gain = self.GetCalibrationResults(0.5*(UTB+UTE), Wmax)
      if gain:
        nbins = hps.GetNbinsX(); hps.SetBins(nbins, zero, zero + gain * nbins)
        self.hps = hps.Clone(); self.hps.Rebin(self.Merger)
        self.hps.GetXaxis().SetTitle('E_{#gamma} [keV]');
        # get rid of spikes:
        for spike in self.EP:
          lo = 1 + int((spike - zero - 3 * self.RR)/(gain*float(self.Merger)))
          hi = 1 + int((spike - zero + 3 * self.RL)/(gain*float(self.Merger)))
          for ch in range(lo, hi): self.hps.SetBinContent(ch, 0.0); self.hps.SetBinError(ch, 0.0)
        Edge = self.fitEdgeSimple(Wmax, self.Ranger)
        if Edge:
         Edge = self.fitEdgeComplex(Edge, self.Ranger)
         if Edge:
           R = self.Beam_Energy(UTB, UTE, Edge)
           if R:
             self.Save_Files((R['t'], R['dt'], R['E'], R['dE'], R['S'], R['dS']))
             if self.SaveDB:
               self.Save_DB((R['t'], R['dt'], R['Eo'], R['dE1'], R['BS'], R['dBS']))
               self.Save_WWW()
           return R
    else: return False

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def fitEdgeSimple(self,W,LBK):
    K = 2.*W/self.me; Wmin = W/(1+K); E1, E2 = W - LBK*Wmin, W + LBK*Wmin; self.simple.SetRange(E1, E2)
    self.const.SetRange(1.02*W, E2); self.hps.Fit('const','QRN');  B = self.const.GetParameter(0)
    self.const.SetRange(E1, 0.98*W); self.hps.Fit('const','QRN');  A = self.const.GetParameter(0) - B
    # Amplitude | Position | Sigma | edge tilt ~x | edge tilt ~x^2 | Background | Background Slope
    self.simple.SetParameters(A, W, self.RR, 0.0, 0.0, B, 0.0)
    self.simple.SetParLimits(0, 1.0, 1.0e+4)
    self.simple.SetParLimits(1, E1,      E2)
    self.simple.SetParLimits(2, 1.0, 1.0e+2)
    R = self.hps.Fit('simple','RSQN'); OK = False
    if not R.Status():
      W = self.simple.GetParameter(1);   E1, E2 = W - LBK*Wmin, W + LBK*Wmin; self.simple.SetRange(E1,E2)
      R = self.hps.Fit('simple','RSQN')
      E = fitParameters(self.simple); Prob = R.Prob()
      print ' ╔ Simple Edge Fit: ═══════════╤══════════════════════╤═══════════════════════════╗'
      print ' ║ Range from %4.0f to %4.0f keV │ E_beam = %7.2f MeV │  W_max = %9.3f keV    ║' % (E1, E2, self.BEPC[self.lepton]['E'], W)
      print ' ╟─────────────────────────────┴───────┬──────────────┴───────────────────────────╢'
      print ' ║ Edge amplitude:  %6.1f ± %6.1f    │ Edge Wmax:       %10.3f ± %10.3f ║' % (E['p0'], E['dp0'], E['p1'], E['dp1'])
      print ' ║ Edge σW, keV:    %6.2f ± %6.2f    │ Background:      %10.3f ± %10.3f ║' % (E['p2'], E['dp2'], E['p5'], E['dp5'])
#      print ' ║ Edge σW, keV:    %10.3f ± %10.3f │ Edge tilt pol1:   %10.7f ± %10.8f ║' % (E['p2'], E['dp2'], E['p5'], E['dp5'])
#      print ' ║ Edge tilt pol2:  %10.8f ± %10.8f │ Background:       %10.3f ± %10.3f ║' % (E['p4'], E['dp4'], E['p5'], E['dp5'])
      print ' ╟─────────────────────────────────────┼──────────────────────────────────────────╢'
      print ' ║         χ²/NDF = %5.1f/%3d          │        Probability: %5.3f                ║' % (E['Chi2'], E['NDF'],Prob)
      print ' ╚═════════════════════════════════════╧══════════════════════════════════════════╝\n'
      OK = (E['p0']>self.MinAmp) and (E['dp0']/E['p0']<1.0) and (E['Chi2']/E['NDF']<10.0) and (E['p2']<100.0)
      self.Lg1.Clear(); self.Lg1.SetHeader('#chi^{2}/NDF = %5.1f/%3d  (Prob: %5.3f)' % (E['Chi2'], E['NDF'], Prob))
      self.Lg1.AddEntry(self.simple, '#omega_{max}  = %7.2f #pm %4.2f keV' % (E['p1'], E['dp1']), 'l')
      self.Lg1.AddEntry(self.simple, '#sigma_{edge} = %7.2f #pm %4.2f keV' % (E['p2'], E['dp2']), 'l')
    self.cc.cd(); self.cc.Clear();  self.cc.SetGrid(); self.hps.SetMarkerStyle(20)
    self.hps.Draw(''); self.hps.GetXaxis().SetRangeUser(E1-10, E2+10); self.simple.Draw('SAME')
    self.Lg1.Draw('SAME'); self.cc.Modified(); self.cc.Update()
    if OK: return E
#    raw_input('Bad Fit?')
    print 'Simple fit: bad fit, bad amplitude, spread or χ²';
    return False


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def fitEdgeComplex(self,E,LBK):
    W = E['p1']
    K = 2.*W/self.me; Wmin = W/(1+K); E1, E2 = W - LBK*Wmin, W + LBK*Wmin; self.comple.SetRange(E1, E2)
    S = E['p2']**2 - self.RR**2; S = 1.0 if S<0.0 else S**0.5 # estimated influence of beam energy spread to the edge width

    Ebeam = self.BEPC[self.lepton]['E'] # MeV
    if self.Asymme>=1.0: self.RL = self.Asymme*self.RR
    # Amplitude | Position | Sigma | edge tilt ~x | edge tilt ~x^2 | Background | Background Slope | Sigma_R | Sigma_L |
    self.comple.SetParameters(E['p0'], E['p1'], S, E['p3'], E['p4'], E['p5'], E['p6'], self.RR, self.RL);
    FIXED = [5,6,7,8];
    for n in FIXED: V = self.comple.GetParameter(n); self.comple.FixParameter(n,V)
    R = self.hps.Fit('comple','RSQN')
    if R.Status(): return False
    for n in FIXED: self.comple.ReleaseParameter(n)

    Ec = fitParameters(self.comple); Prob = R.Prob();

    print ' ╔ Complex Edge Fit: ══════════╤════════════════════════╤═════════════════════════╗'
    print ' ║ Range from %4.0f to %4.0f keV │ Amplitude = %6.1f     │  W_max = %9.3f keV  ║' % (E1, E2, E['p0'], E['p1'])
    print ' ║    σR = %6.3f ± %5.3f keV  │  σL = %5.2f ± %4.2f keV │     σW  = %6.3f keV    ║' % (self.RR, self.dRR, self.RL, self.dRL, S)
    print ' ╟─────────────────────────────┴────────┬───────────────┴─────────────────────────╢'
    print ' ║ Edge amplitude: %8.1f ± %6.1f    │ Edge wmax, keV:   %8.3f ± %7.3f    ║' % (Ec['p0'], Ec['dp0'], Ec['p1'], Ec['dp1'])
    print ' ║ Edge σW, keV:     %6.3f ± %5.3f     │ Background:       %8.1f ± %7.1f    ║' % (Ec['p2'], Ec['dp2'], Ec['p5'], Ec['dp5'])
    print ' ║ HPGe σR, keV:     %6.3f ± %5.3f     │ HPGe σL, keV:  %10.3f ± %10.3f  ║' % (Ec['p7'], Ec['dp7'], Ec['p8'], Ec['dp8'])
    print ' ╟──────────────────────────────────────┼─────────────────────────────────────────╢'
    print ' ║         χ²/NDF = %5.1f/%3d           │          Probability: %5.3f             ║' % (Ec['Chi2'], Ec['NDF'], Prob)
    print ' ╚══════════════════════════════════════╧═════════════════════════════════════════╝\n'

    self.Lg2.Clear(); self.Lg2.SetHeader('#chi^{2}/NDF = %5.1f/%3d  (Prob: %5.3f)' % (Ec['Chi2'], Ec['NDF'], Prob))
    self.Lg2.AddEntry(self.comple, '#omega_{max} =  %7.2f #pm %4.2f keV' % (Ec['p1'], Ec['dp1']), 'l')
    self.Lg2.AddEntry(self.comple, '#sigma_{R} = %5.2f keV;  #sigma_{L} = %5.2f keV' % (self.RR, self.RL), 'l')
    self.Lg2.AddEntry(self.comple, '#sigma_{beam} = %7.2f #pm %4.2f keV' % (Ec['p2'], Ec['dp2']), 'l')
    self.comple.DrawCopy('SAME'); self.Lg2.Draw('SAME'); self.cc.Modified(); self.cc.Update()

#    print ' Wmax: %7.2f ± %4.2f keV (symmetric fit)'      % (E['p1'], E['dp1'])
#    Wmax, dWmax = Ec['p1'],  Ec['dp1']
#    print ' Wmax: %7.2f ± %4.2f keV (asymmetric fit)'     % (Wmax, dWmax)
#    Wmax, dWmax = Wmax-self.SC, (dWmax**2+self.dSC**2)**0.5;
#    print ' Wmax: %7.2f ± %4.2f keV (spline correction )' % (Wmax, dWmax)

#    outs = '%10d  %5d  %4d  %4d   %4d  %4d # point %2d (%s)\n' % (t, dt, int(E[4]), int(E[5]), int(P[4]), int(P[5]), PointN, Points_Names[PointN])
#    outs = '%6.3f  %5.3f  %7.5f  %7.5f  %9.7f  %9.7f\n' % (Ec['p2'], Ec['dp2'], Ec['p3'], Ec['dp3'], Ec['p4'], Ec['dp4'])
#    with open('%8s.dat' % self.lepton,'a') as f: f.write(outs)

#    self.Wmax, self.dWmax = Wmax, dWmax
#    self.BsCt, self.dBsCt = Ec['p2'], Ec['dp2']
    BAD = Ec['p0']<self.MinAmp or Ec['dp0']/Ec['p0']>1.0 or Ec['Chi2']/Ec['NDF']>10.0
    if BAD:  print 'Complex fit: bad amplitude, spread or χ²';  return False
    return Ec

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Beam_Energy(self, utb, ute, edge):
    R = {}
    def ToE(W, dW, SW, dSW):
      E = 0.5 * W * (1.0 + (1.0 + self.me**2/W/self.wo)**0.5)
      x = 4.0 * self.wo * E /self.me**2
      D = (1.0 + x)**2 / ((2.0 + x) * x)
      R['SRC'] = 4.75e-3*(1e-6*E)**4 # Energy Correction due to SR losses [MeV]
      R['E'  ] = 1e-3*E + R['SRC']
      R['dE1'] = 1e-3*D*dW[0] # keV -> MeV
      R['dE2'] = 1e-3*D*dW[1] # keV -> MeV
      R['dE3'] = 1e-3*D*dW[2] # keV -> MeV
      R['S'  ] = D*SW
      R['dS' ] = D*dSW
      return R

    R   = ToE(edge['p1']-self.SC, [edge['dp1'], self.dW, self.dSC], edge['p2'], edge['dp2'])
    st  = 'Measurement time '
    st += time.strftime('from %Y.%m.%d/%H:%M:%S', time.localtime(utb))
    st += time.strftime(' to %Y.%m.%d/%H:%M:%S', time.localtime(ute))
    print '\n ╔ %8s Beam Energy Determination:  ══════════════════════════════════════════╗' % self.lepton
    print ' ║      BEPC beam energy = %8.3f ± %5.3f MeV was taken from database           ║' % (self.BEPC[self.lepton]['E'], self.BEPC[self.lepton]['dE'])
    print ' ║    %66s          ║' % (st)
    print ' ║                                                                                ║'
    print ' ║      BEMS beam energy = %8.3f ± %5.3f ± %5.3f ± %5.3f MeV                   ║'      % (R['E'], R['dE1'], R['dE2'], R['dE3'])
    print ' ║ Correction +%5.3f MeV for half-a-ring synchrotron radiation losses was added   ║'   %  R['SRC']
    print ' ║      BEMS beam energy spread = %4.0f ± %4.0f keV                                 ║' % (R['S'], R['dS'])
    print ' ╚════════════════════════════════════════════════════════════════════════════════╝\n'
    R['t'], R['dt'] = 0.5*(utb + ute), 0.5*(ute - utb)
    R['dE'] = (R['dE1']**2 + R['dE2']**2 + R['dE3']**2)**0.5
    self.Lg2.AddEntry(self.comple, 'E_{beam} = %8.2f #pm %5.2f MeV'   % (R['E'], R['dE']), 'l')
    self.Lg2.AddEntry(self.comple, '#sigma_{E} = %4.0f #pm %4.0f keV' % (R['S'], R['dS']), 'l')
    self.Lg2.Draw('SAME'); self.cc.Modified(); self.cc.Update()
    return R


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def GetCalibrationResults(self, t, wmax):
    for attempt in xrange(3):
      Scale = self.HPGe.Get_Calibration(t, wmax)
      L = len(Scale['R'])
      if L: break
      else: print 'Waiting 20s for calibration results...';   time.sleep(20)
    if L == 1: c = 0
    elif L >1:
      Q = [Scale['dW'][c] * Scale['dC'][c] * Scale['dR'][c] * Scale['dL'][c] for c in range(L)]
#      for c in range(L): Q.append(Scale['dW'][c] * Scale['dC'][c] * Scale['dR'][c] * Scale['dL'][c])
      if len(Q)==0:    print 'No valid calibration was found'; return 0, 0
      c = Q.index(min(Q))
    else:              print 'No valid calibration was found'; return 0, 0
    zero, gain         = Scale['Z'][c], Scale['G' ][c]
    self.dW            =                Scale['dW'][c] # linear calibration statistical error, keV
    self.SC, self.dSC  = Scale['C'][c], Scale['dC'][c] # PB-5  scale correction and its error, keV
    self.RR, self.dRR  = Scale['R'][c], Scale['dR'][c] # Right resolution sigma and its error, keV
    self.RL, self.dRL  = Scale['L'][c], Scale['dL'][c] # Left  resolution sigma and its error, keV
    self.EP            = Scale['X'][c]                 # exclude peaks, keV
#    print Scale['N'][c]
    print ' ╔ HPGe calibration: %15s ══════════════════╤══════════════════════════╗' % (self.HPGe.outfile)
    print ' ║  W_max  = %9.3f keV  │ σR = %6.3f ± %5.3f keV  │ σL = %6.3f ± %5.3f keV  ║' % (wmax, self.RR, self.dRR, self.RL, self.dRL)
    print ' ╚══════════════════════════╧══════════════════════════╧══════════════════════════╝\n'
    return zero, gain

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Wmax_expect(self):
    self.BEPC[self.lepton]['E'] += self.ETuner
    Eo = 1.e+3*self.BEPC[self.lepton]['E']
    k  = 4.*Eo*self.wo/self.me**2;
    return Eo*k/(1.+k)

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Save_Files(self, results):
    print 'Saving (%s) to history file' % (self.lepton)
    S = results + (self.BEPC[self.lepton]['E'], self.BEPC[self.lepton]['dE'], self.BEPC[self.lepton]['I'], self.BEPC[self.lepton]['dI'])
    S = '%10d  %5d  %8.3f  %5.3f  %5.0f  %5.0f  %7.2f  %4.2f  %3.0f  %3.0f\n' % S
    if   'electron' in self.lepton:
      with open('E.results','a') as f:  f.write(S)
    elif 'positron' in self.lepton:
      with open('P.results','a') as f:  f.write(S)

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Save_DB(self, (t, dt, EB, dEB, BS, dBS)):
    print 'Saving (%s) to BEPC database' % (self.lepton)
    from EpicsCA import caput
    if 'lectron' in self.lepton:
      caput('CBSE:E',  EB); caput('CBSE:DE', dEB); caput('CBSE:S', BS); caput('CBSE:DS', dBS); caput('CBSE:BT', t-dt); caput('CBSE:ET', t+dt)
    elif 'ositron' in self.lepton:
      caput('CBSP:E',  EB); caput('CBSP:DE', dEB); caput('CBSP:S', BS); caput('CBSP:DS', dBS); caput('CBSP:BT', t-dt); caput('CBSP:ET', t+dt)
    else:
      pass

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
  def Save_WWW(self):
    print 'Sending to web-server...'
    self.cc.SaveAs('Energy.png')
    os.system('ourweb.py')
    time.sleep(10)
    if   'lectron' in self.lepton: rfile='/home/BEMSystem/public_html/images/E.png'
    elif 'ositron' in self.lepton: rfile='/home/BEMSystem/public_html/images/P.png'
    else:                          rfile='/home/BEMSystem/public_html/images/B.png'
#    if not self.prompt:
    try:
      ssh = paramiko.SSHClient()
      ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
      ssh.connect('docbes3.ihep.ac.cn', username='BEMSystem', password='88230944')
      ftp = ssh.open_sftp()
      print ftp
      print 'OK'
      ftp.put('./Energy.png', rfile)
      ftp.put('./in-time.png', '/home/BEMSystem/public_html/images/in-time.png')
      ftp.put('./index.html',  '/home/BEMSystem/public_html/index.html')
      ftp.close()
      ssh.close()
    except:
      pass
    print 'All Done'
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class EdgeSimple():
  sqrt2pi = (2.*ROOT.TMath.Pi())**0.5
  """
  Edge Simple has 7 parameters:
  Amplitude | Position | Sigma_S | edge tilt ~x | edge tilt ~x^2 | Background | Background Tilt |
  """
  def __call__(self, x, p):
    X  = x[0]-p[1]
    Er = 0.5*ROOT.TMath.Erfc(X/(2.**0.5*p[2]))
    Ex =     ROOT.TMath.Exp(-0.5*(X/p[2])**2)/self.sqrt2pi
    Ec = Er*p[2]**2 + Ex*X*p[2]
    T1 = (1 + p[3]*X + p[4]*X**2) * Er
    T2 = (    p[3] + 2*p[4]*X   ) * Ex * p[2]
    T3 = p[4]*p[2]**2 * (Er + X * Ex)
    return p[0] * (T1 - T2 + T3) + p[5] + X * p[6]
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class EdgeComple():
  """
  Edge Complex has 9 parameters:
  Amplitude | Position | Sigma_S | edge tilt ~x | edge tilt ~x^2 | Background | Background Tilt | Sigma_R | Sigma_L
  F1, F2, F3 have 3 parameters:  Sigma_R | Sigma_L | Sigma_S |
  """
  infty = 20
  D1 = '(2*([0]**2 + [2]**2))'; D2 = '(2*([1]**2 + [2]**2))'
  T  = '[0]*exp(-x**2/%s) * TMath::Erfc(-x*[0]/([2]*sqrt%s)) / sqrt%s + ' % (D1,D1,D1)
  T += '[1]*exp(-x**2/%s) * TMath::Erfc( x*[1]/([2]*sqrt%s)) / sqrt%s'    % (D2,D2,D2)
  T  = '(' + T + ') / (sqrt(pi)*([0]+[1])) '
  F1 = ROOT.TF1('F1',          T,  0.0, 1.0)
  F2 = ROOT.TF1('F2', 'x*'   + T , 0.0, 1.0)
  F3 = ROOT.TF1('F3', 'x*x*' + T , 0.0, 1.0)
  parset = []

  def __call__(self, x, p):
    self.curpar  =        p[7], p[8], p[2]
    self.F1.SetParameters(p[7], p[8], p[2])
    self.F2.SetParameters(p[7], p[8], p[2])
    self.F3.SetParameters(p[7], p[8], p[2])
    X = x[0]-p[1]
    if self.parset == self.curpar:
      self.I1 += self.F1.Integral(X, self.Xo)
      self.I2 += self.F2.Integral(X, self.Xo)
      self.I3 += self.F3.Integral(X, self.Xo)
    else:
      self.parset  = self.curpar
      self.I1  = self.F1.Integral(X, self.infty * p[2])
      self.I2  = self.F2.Integral(X, self.infty * p[2])
      self.I3  = self.F3.Integral(X, self.infty * p[2])
    self.Xo = X
    return p[0] * (self.I1*(1 + p[3]*X + p[4]*X**2) - (p[3] + 2*p[3]*p[4])*self.I2 + p[4]**2*self.I3) + p[5] + X * p[6]
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class BEMS_DB:
  NAMEs = { 't':'t', 'Ee':'e_energy', 'Ep':'p_energy', 'Ie':'e_current', 'Ip':'p_current' }

  def GetRunInfo(self,filechain):
    T = { k:[] for k in self.NAMEs.keys()}
    for el in filechain:
      fn = el.split('.')[0].replace('SPECTRA','LOADS') + '.sta'
      try:
        with open(fn,'rb') as fp: data = cPickle.load(fp)
      except:
        print 'sta read error'; return
      OK, LREC = (len(data.keys())==8), len(data['t'])
      if not OK: print 'Bad file: %s' % fn; return False
      else:
        for k in data.keys():
          L = len(data[k])
          if L < LREC:
            print 'Warning: LREC_t = %d, while LREC_%s = %d' % (LREC,k,L)
            data[k].extend([data[k][-1]]*(LREC-L))
          elif L > LREC:
            print 'Error: LREC_t = %d, while LREC_%s = %d' % (LREC,k,L)
            return False
      for k,v in self.NAMEs.iteritems(): T[k].extend(data[v])

    BAD, LREC = 0, len(T['t'])
    for nr in range(LREC):
      if abs(T['Ee'][nr-BAD]-T['Ep'][nr-BAD]) > 0.01: # i.e. 10 MeV, supresses 'jumps'
        for k in T.keys(): T[k].pop(nr-BAD)
        BAD += 1
    if BAD>0: LREC -= BAD; print '(BAD Energy is in %d records)' % BAD
    if float(BAD)/float(LREC)>0.1: return False
    R = {'electron':{}, 'positron':{}}
    A = np.fromiter(T['Ee'], np.float);  R['electron']['E'], R['electron']['dE']  = 1000.*A.mean(), 1000.*A.std()
    A = np.fromiter(T['Ie'], np.float);  R['electron']['I'], R['electron']['dI']  =       A.mean(),       A.std()
    A = np.fromiter(T['Ep'], np.float);  R['positron']['E'], R['positron']['dE']  = 1000.*A.mean(), 1000.*A.std()
    A = np.fromiter(T['Ip'], np.float);  R['positron']['I'], R['positron']['dI']  =       A.mean(),       A.std()
    return R


def fitParameters(fitf):
  pl, np = {}, fitf.GetNumberFreeParameters()
  pl['Chi2'], pl['NDF'] = fitf.GetChisquare(),  fitf.GetNDF()
  for p in range(np):  pl['p'+str(p)], pl['dp'+str(p)] = fitf.GetParameter(p), fitf.GetParError(p)
  return pl


