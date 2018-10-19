#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os, sys, getopt, time, string, gzip, ConfigParser


def Usage():
  print '''
 ╔ Python script to generate HPGe spectra. © 2018 N.Yu.Muchnoi ═══════════════════╗
 ║                                             ╭────────────────────────────────╮ ║
 ║ Usage: %0000000012s [options]               │ Last update: October 01, 2018  │ ║
 ║                                             ╰────────────────────────────────╯ ║
 ║ List of options:                                                               ║
 ║ -h,          --help               : print this help message and exit,          ║
 ║ -c filename, --cfg     = filename : file to read various parameters,           ║
 ║ -n N,        --nf   = N           : put several files into one spectrum.       ║
 ╚════════════════════════════════════════════════════════════════════════════════╝
 ''' % sys.argv[0].split('/')[-1];  sys.exit(0)

if ('-h' in sys.argv) or ('--help' in sys.argv):
  Usage()
else:
  import ROOT
#  ROOT.gROOT.LoadMacro(sys.argv[0].replace(sys.argv[0].split('/')[-1], 'vepp2k/airy_ai_int.C'))
  ROOT.gROOT.SetStyle("Plain")
  ROOT.gROOT.ForceStyle()
  ROOT.gStyle.SetTitleBorderSize(0)
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptFit(0)


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class ToDo:

  def __init__(self):

    cfg = ConfigParser.ConfigParser(); cfg.read('_globals_.cfg')
    if cfg.has_option('globals', 'data_folder'):
      self.root = cfg.get('globals', 'data_folder')
    else:
      exit()

    try:
      Years = os.listdir(self.root); Years.sort()
    except OSError:
      print 'root folder does not exist!'
      sys.exit(1)

    self.nst = 1
    self.cfg = 'MC.sim'
    self.Options(sys.argv)

  def Options(self,argv):
    sopt = "c:n:h"
    lopt = ["cfg=","nst=","help"]
    try:
      opts, args = getopt.getopt(argv[1:], sopt, lopt)
    except getopt.GetoptError:
      print 'Wrong option'; Usage(); sys.exit(2)

    for opt, arg in opts:
      if   opt in ("-n", "--nst")  : self.nst = int(arg)
      elif opt in ("-c", "--cfg")  : self.cfg =     arg



#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class MonteCarlo:
  nbins    = 16384   # number of bins in a histogram
  hps      = ROOT.TH1I('hps','',nbins, 0, nbins)

  def __init__(self,todo):
    from vepp2k.vepp2k import Constants
    self.nst = todo.nst
    cfg   = ConfigParser.ConfigParser(); cfg.read(todo.cfg)
    self.lwave   = cfg.getfloat('edge', 'WaveLength')       # laser wavelength [m]
    self.ebeam   = cfg.getfloat('edge', 'BeamEnergy')       # beam energy [MeV]
    self.sbeam   = cfg.getfloat('edge', 'BeamSpread')       # beam energy spread [MeV]
    self.field   = cfg.getfloat('edge', 'BendDipole')       # bending field induction [T]

    if self.field>0.0:
      try: from vepp2k.vepp2k import EdgeSimple, HPGeSpread
      except ImportError: print 'VEPP2K is not yet implemented'; exit(0)
    else:
      try: from bepc.bepc     import EdgeSimple, EdgeComple
      except ImportError: print 'BEPCII is not yet implemented'; exit(0)

    self.cv = ROOT.TCanvas('cv','HPGe spectrum', 2, 2, 1002, 1002)

    Constants.wo = Constants.h*Constants.c/self.lwave       # laser photon energy [eV]
    Constants.Eo = 0.25e-6 * Constants.me**2 / Constants.wo # (me^2/4wo) [MeV]
    self.k       = self.ebeam/Constants.Eo                  # kappa
    SCALE        = self.k/(1. + self.k)
    self.wmax    = 1000 * self.ebeam * SCALE                # wmax [keV]
    SCALE        = SCALE * (2.0 + self.k)/(1.0 + self.k)
    self.smax    = 1000 * self.sbeam * SCALE                # beam spread blur at wmax [keV]

    print self.smax

    E1, E2 = int(0.95*self.wmax), int(1.01*self.wmax)

    self.simple = ROOT.TF1('simple', EdgeSimple(), 0, 1, 7 )
    self.simple.SetRange(float(E1), float(E2))
    self.spreso = ROOT.TF1('spreso', HPGeSpread(), 0, 1, 3 );
    self.convol = ROOT.TF1Convolution('simple', 'spreso', 0.99*float(E1), 1.01*float(E2));
    self.convol.SetNofPointsFFT(2000)
    self.spread = ROOT.TF1('spread', self.convol, float(E1), float(E2), 10)
    self.comple = ROOT.TF1('comple', self.convol, float(E1), float(E2), 10)

    self.simple.SetNpx(1000); self.simple.SetLineColor(ROOT.kBlue)
    self.spread.SetNpx(1000); self.spread.SetLineColor(ROOT.kBlack)
    self.comple.SetNpx(1000); self.comple.SetLineColor(ROOT.kRed)


  def __del__(self): pass

  def Go(self):

    AMP = 100.0
    # p[0] - Beam Energy [MeV] | p[1] - Bending field [T] | p[2] - Amplitude     | p[3] - Linear Slope at Wmax | p[4] - Square Slope at Wmax
    # p[5] - Background        | p[6] - Background Slope  | p[7] - sigma_R [keV] | p[8] - sigma_L [keV]        | p[9] - Beam spread blur [keV]

    P = [self.ebeam, self.field, AMP, -0.001, 0.00001, 10.0, -0.005, 0.1, 0.1, self.smax]

    """
    self.spreso.SetParameters(0.1, 0.2, 0.2)
    print 'I = ', self.spreso.Integral(-1,1)
    self.spreso.SetRange(-1,1)
    self.cv.cd();       self.cv.SetGrid()
    self.spreso.Draw()
    self.cv.Modified(); self.cv.Update()

#    """
    self.simple.SetParameters(np.fromiter(P[0:7], np.float))

    self.spread.SetParameters(np.fromiter(P,      np.float))

    P[7], P[8] = 1.5, 1.5
    self.comple.SetParameters(np.fromiter(P,      np.float))

    self.cv.cd();       self.cv.SetGrid()
    self.simple.Draw()
    self.spread.Draw('SAME')
    self.comple.Draw('SAME')
    self.cv.Modified(); self.cv.Update()
#    """




def main(argv):
  todo  = ToDo()
  try:
    A = MonteCarlo(todo)
    A.Go()
  except KeyboardInterrupt: print '\nExecution is interrupted'
  raw_input('Quit?')
  exit()


if __name__ == "__main__": main(sys.argv)

