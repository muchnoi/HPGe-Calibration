# -*- coding: utf-8 -*-
#import paramiko

import ROOT, ConfigParser
import numpy as np
from bepc import EDGE, BEMS_DB
from hpge import DataFile, Histogram

class Energy_Scan:
  IP = 2.0*ROOT.TMath.Cos(0.011)

  def Go(self, flist, todo):
    T, Emin, Emax = DataFile(), 1830.0, 1850.0
    ngood, TABLE = 0, {'file':[], 't':[], 'dt':[], 'E':[], 'dE':[], 'Se':[], 'Sp':[]}
    cfg = ConfigParser.ConfigParser(); cfg.read(todo.cfg_file)
    if cfg.has_option('scan', 'efrom'):  Emin = cfg.getfloat('scan', 'efrom')
    if cfg.has_option('scan', 'eupto'):  Emax = cfg.getfloat('scan', 'eupto')
    for f in flist:
      T.ReadData(f,16384,100)
      t, dt = 0.5*(T.utb + T.ute), 0.5*(T.ute - T.utb)
      R = BEMS_DB().GetRunInfo([f])
      if T and R:
        E  = 0.5*(R['electron']['E']     + R['positron']['E']    )
        dE = 0.5*(R['electron']['dE']**2 + R['positron']['dE']**2)**0.5
        Se = R['electron']['dE']
        Sp = R['positron']['dE']
        if dt>200 and dE<0.5 and not ('wall' in f) and (Emin < E < Emax):
          ngood += 1
          TABLE['file'].append(f)
          TABLE['t' ].append(  t); TABLE['dt'].append(dt)
          TABLE['E' ].append(  E); TABLE['dE'].append(dE)
          TABLE['Se'].append( Se); TABLE['Sp'].append(Sp)
      else: print "bad records in %s (acq. time is %d s)" % (f, 2*dt)

    self.cfh = ROOT.TCanvas('cfh','Files history', 2, 2, 1002, 1002);  self.cfh.Divide(1,2)
    self.Show_Points(TABLE);   raw_input('All good points...')

    n, PointN = 0, 0
    Points_Names, Energy_Points = [], [ [ TABLE['file'][n] ] ]
    Points_Names.append('%.3f MeV' % TABLE['E'][n])
    while n < ngood-1:
      Condition1 = abs(TABLE['E'][n] - TABLE['E'][n+1]) < (TABLE['dE'][n] + TABLE['dE'][n+1])
      Condition2 = abs(TABLE['t'][n+1] - TABLE['t'][n]) < 4000 # separate points with "large" time gaps
      if Condition1 and Condition2:
        Energy_Points[PointN].append(TABLE['file'][n+1])
      else:
        Energy_Points.append([TABLE['file'][n+1]]); PointN += 1
        Points_Names.append('%.3f MeV' % TABLE['E'][n+1])
      n += 1

    print 'Initially:'
    for PointN in range(len(Energy_Points)):
      print 'Point %2d energy: %s. Number of files: %d' % (PointN, Points_Names[PointN], len(Energy_Points[PointN]))

# Clear TABLE
    p, tnp = 0, len(Energy_Points)
    while p < tnp:
      Ne = len([fn for fn in Energy_Points[p] if 'elec' in fn])
      Np = len([fn for fn in Energy_Points[p] if 'posi' in fn])
      if len(Energy_Points[p]) < 3 or not (Ne and Np):
        for el in Energy_Points[p]:
          rem_rec = TABLE['file'].index(el)
          for key in TABLE.keys(): TABLE[key].pop(rem_rec)
        Energy_Points.pop(p);      Points_Names.pop(p)
        tnp -= 1
      else:
        p += 1

    self.Show_Points(TABLE)

    print 'After selection:'
    for PointN in range(len(Energy_Points)):
      print 'Point %2d energy: %s. Number of files: %d' % (PointN, Points_Names[PointN], len(Energy_Points[PointN]))

    raw_input( '%d energy scan points have been founded' % PointN )

    self.H = Histogram(todo)
    Bs_elec, Bs_posi = [], []
    for PointN in range(len(Energy_Points)):
      print '\nPoint %2d energy: %s. Number of files: %d' % (PointN, Points_Names[PointN], len(Energy_Points[PointN]))
      if todo.edge:
        flist = [el for el in Energy_Points[PointN] if 'elec' in el]; self.H.nfile = len(flist); E = self.H.Go(flist)
        flist = [el for el in Energy_Points[PointN] if 'posi' in el]; self.H.nfile = len(flist); P = self.H.Go(flist)
        if E and P:
          Bs_elec.append(E['S']);   Bs_posi.append(P['S'])
          tmin, tmax = min(E['t'] - E['dt'], P['t'] - P['dt']), max(E['t'] + E['dt'], P['t'] + P['dt'])
          t   , dt   = 0.5*(tmin + tmax), 0.5*(tmax - tmin)
          Ecm   = self.IP *(E['E'] * P['E'])**0.5
          dE1  = (E['dE1']**2 + P['dE1']**2)**0.5
          dE2  = (E['dE2']**2 + P['dE2']**2)**0.5
          dE3  = (E['dE3']**2 + P['dE3']**2)**0.5
          Scm  = (  E['S']**2 +   P['S']**2)**0.5;
          dScm = ((E['S']*E['dS'])**2 + (P['S']*P['dS'])**2)**0.5/Scm
          Scm  = 0.001*Scm
          dScm = 0.001*dScm
          SRC  = E['SRC']+P['SRC']
          outs = (t, dt, Ecm, SRC, dE1, dE2, dE3, Scm, dScm, PointN, Points_Names[PointN])
          outs = '%10d  %5d  %8.3f  %5.3f  %5.3f  %5.3f  %5.3f  %6.3f  %5.3f # point %2d (%s)\n' % outs
          with open('SCAN.results','a') as f: f.write(outs)
#        raw_input('Point %d energy measurement is ready' % PointN)
      else:
        flist = [el for el in Energy_Points[PointN]];  self.H.nfile = len(flist);  self.H.Go(flist)
        print 'Point %d scale calibration is ready (porridge)'  % PointN

    if todo.edge:
      Se = np.fromiter(Bs_elec, np.float)
      Sp = np.fromiter(Bs_posi, np.float)
      outs  = '# The beam energy spreads averaged over the scan points:\n'
      outs += "# e- beam spread = %.0f ± %.0f keV\n" % (Se.mean(), Se.std())
      outs += "# e+ beam spread = %.0f ± %.0f keV\n" % (Sp.mean(), Sp.std())
      with open('SCAN.results','a') as f: f.write(outs)




    raw_input('So what?')


  def Show_Points(self,TABLE):
    TE, TP = {'t':[], 'dt':[], 'E':[], 'dE':[]}, {'t':[], 'dt':[], 'E':[], 'dE':[]}
    for n in range(len(TABLE['file'])):
      if   'elec' in TABLE['file'][n]:
        for k in TE.keys(): TE[k].append(TABLE[k][n])
      elif 'posi' in TABLE['file'][n]:
        for k in TP.keys(): TP[k].append(TABLE[k][n])
    Te, dTe = np.fromiter(TE['t'], np.float), np.fromiter(TE['dt'], np.float)
    Ee, dEe = np.fromiter(TE['E'], np.float), np.fromiter(TE['dE'], np.float)
    Tp, dTp = np.fromiter(TP['t'], np.float), np.fromiter(TP['dt'], np.float)
    Ep, dEp = np.fromiter(TP['E'], np.float), np.fromiter(TP['dE'], np.float)
    GEE = ROOT.TGraphErrors(len(Te),Te,Ee,dTe,dEe); GEE.SetMarkerColor(ROOT.kBlack); GEE.SetLineColor(ROOT.kBlack); GEE.SetMarkerStyle(20)
    GEP = ROOT.TGraphErrors(len(Tp),Tp,Ep,dTp,dEp); GEP.SetMarkerColor(ROOT.kRed);   GEP.SetLineColor(ROOT.kRed);   GEP.SetMarkerStyle(20)
    GEE.GetXaxis().SetTimeDisplay(1);  GEE.GetXaxis().SetTimeFormat('#splitline{%b%d}{%H:%M}%F1970-01-01 00:00:00')
    GEE.GetXaxis().SetLabelSize(0.02); GEE.GetXaxis().SetTitle('China time');        GEE.GetXaxis().SetLabelOffset(0.02)
    GEE.SetTitle('Mean BEPC Energy [MeV]')
    self.GEE = GEE
    self.GEP = GEP

    T,  dT = np.fromiter(TABLE['t'],  np.float), np.fromiter(TABLE['dt'], np.float)
    Se, Sp = np.fromiter(TABLE['Se'], np.float), np.fromiter(TABLE['Sp'], np.float)
    GSE = ROOT.TGraph(len(T),T,Se); GSE.SetMarkerColor(ROOT.kBlack); GSE.SetLineColor(ROOT.kBlack); GSE.SetMarkerStyle(20)
    GSP = ROOT.TGraph(len(T),T,Sp); GSP.SetMarkerColor(ROOT.kRed);   GSP.SetLineColor(ROOT.kRed);   GSP.SetMarkerStyle(20)
    GSE.GetXaxis().SetTimeDisplay(1);  GSE.GetXaxis().SetTimeFormat('#splitline{%b%d}{%H:%M}%F1970-01-01 00:00:00')
    GSE.GetXaxis().SetLabelSize(0.02); GSE.GetXaxis().SetTitle('China time'); GSE.GetXaxis().SetLabelOffset(0.02)
    GSE.SetTitle('R.m.s. BEPC Energy [MeV]')
    self.GSE = GSE
    self.GSP = GSP

    self.cfh.cd(1); self.cfh.SetGrid(); self.GEE.Draw('AP'); self.GEP.Draw('PSAME');
    self.cfh.cd(2); self.cfh.SetGrid(); self.GSE.Draw('APL'); self.GSP.Draw('PLSAME');

    self.cfh.Modified(); self.cfh.Update()

