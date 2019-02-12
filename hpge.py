#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys, getopt, time, string, gzip, ConfigParser

def Usage():
  print '''
 ╔ Python script to process HPGe spectra. © 2005-2018 Nickolai Muchnoi ═══════════╗
 ║                                             ╭────────────────────────────────╮ ║
 ║ Usage: %0000000012s [options]               │ Last update: December 4, 2018  │ ║
 ║                                             ╰────────────────────────────────╯ ║
 ║ List of options:                                                               ║
 ║ -h,          --help               : print this help message and exit.          ║
 ║ -i,          --interactive        : if set, script will prompt to proceed.     ║
 ║ -k,          --keV                : use channels not keV.                      ║
 ║ -l,          --list               : just show the list of files.               ║
 ║              --toi                : show the table of known isotopes.          ║
 ║ -f expr,     --file = expr        : file name(s) under specified folder(s).    ║
 ║ -n N,        --nf   = N           : put several files into one spectrum.       ║
 ║ -d YYYYMMDD, --folder  = YYYYMMDD : date to start from (year, month, day).     ║
 ║ -e YYYYMMDD, --efolder = YYYYMMDD : date to  end  with (year, month, day).     ║
 ║ -t HHMMSS,   --time    = HHMMSS   : time from which to start (hour, min, sec). ║
 ║ -c filename, --cfg     = filename : file to read various parameters,           ║
 ║                                     otherwise "online.cfg" is used.            ║
 ║ -s filename, --scale   = filename : file to store/get calibration results,     ║
 ║                                     otherwise "escale.root" is used.           ║
 ║ -v energy,   --verify  = energy   : calibration results for an energy [keV].   ║
 ║              --edge               : try to measure beam energy by Compton edge.║
 ║              --escan              : deal with beam energy scan experiment.     ║
 ║ -g       ,   --generate           : generate subfolders with 'success.list'    ║
 ║                                     and 'failure.list' file containers.        ║
 ║ -p folder,   --point   = folder   : subfolder under the working folder,        ║
 ║                                     where there is the 'success.list' file.    ║
 ╚════════════════════════════════════════════════════════════════════════════════╝
 ''' % sys.argv[0].split('/')[-1];  sys.exit(0)

def List_TOI():
  from scale.atlas import Atlas
  OUT = {}
  for k,v in Atlas().atlas.iteritems():
    for el in v:
      OUT[el['W']] = "%5s(%4s): Eγ = %9.3f ± %5.3f keV" % (k, v.index(el), el['W'], el['dW'])
  for key in sorted(OUT.keys()): print OUT[key]
  sys.exit(0)

if ('-h' in sys.argv) or ('--help' in sys.argv):
  Usage()
elif ('--toi' in sys.argv):
  List_TOI()
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
    if cfg.has_option('globals', 'data_folder'):  self.root = cfg.get('globals', 'data_folder')
    else:                                         exit()
    if cfg.has_option('globals', 'data_origin'):  self.orig = cfg.get('globals', 'data_origin')
    else:                                         exit()

    try:
      Years = os.listdir(self.root);       Years.sort()
    except OSError:
      print 'root folder does not exist!'; exit()

    Days  = os.listdir(self.root+Years[-1])
    if len(Days)==1: Days.extend(os.listdir(self.root+Years[-2]));
    Days.sort()
    self.folder,     self.efolder  = Days[-2], Days[-1]
    self.tfrom,      self.nfile    = 0   , 1
    self.filename  = self.lastfile = ''
    self.tokeV     = self.online   = True
    self.point     = './'
    self.prompt    = self.listonly = self.edge = self.escan = self.verify = self.generate = False
    self.cfg_file  = 'online.cfg'
    self.scalefile = 'escale.root'
    self.Options(sys.argv)

  def Options(self,argv):
    S = D = E = False
    sopt = "c:d:e:f:n:p:s:t:v:ghikl"
    lopt = ["cfg=", "folder=", "efolder=", "file=", "time=", "nfiles=", "scale=", "verify=", "point=",
            "generate", "interactive", "help", "keV", "list", "edge", "escan"]
    try:
      opts, args = getopt.getopt(argv[1:], sopt, lopt)
    except getopt.GetoptError:
      print 'Wrong option'; Usage(); sys.exit(2)

    for opt, arg in opts:
      if   opt in ("-i", "--interactive"): self.prompt    = True;     self.online = False
      elif opt in ("-n", "--nf")         : self.nfile     = int(arg)
      elif opt in ("-d", "--folder")     : self.folder    = arg;      self.online = False; D = True
      elif opt in ("-e", "--efolder")    : self.efolder   = arg;      self.online = False; E = True
      elif opt in ("-f", "--file")       : self.filename  = arg;      self.online = False
      elif opt in ("-t", "--time")       : self.tfrom     = int(arg); self.online = False
      elif opt in ("-k", "--keV")        : self.tokeV     = False
      elif opt in ("-s", "--scale")      : self.scalefile = arg;                           S = True
      elif opt in ("-v", "--verify")     : self.verify    = True;     self.online = False; self.energy = float(arg)
      elif opt in ("-l", "--list")       : self.listonly  = True;
      elif opt in (      "--edge")       : self.edge      = True;                          self.tokeV  = False
      elif opt in (      "--escan")      : self.escan     = True;     self.online = False
      elif opt in ("-g", "--generate")   : self.generate  = True;     self.online = False
      elif opt in ("-c", "--cfg")        : self.cfg_file  = arg;      self.online = False
      elif opt in ("-p", "--point")      : self.point     = arg;      self.online = False
    cfg = ConfigParser.ConfigParser(); cfg.read(self.cfg_file)
    if cfg.has_option('scale', 'file') and not S: self.scalefile = self.point + cfg.get('scale', 'file')
    if not self.online:
      if cfg.has_option('scan', 'bdate') and not D:  self.folder = cfg.get('scan', 'bdate')
      if cfg.has_option('scan', 'edate') and not E: self.efolder = cfg.get('scan', 'edate')

  def GetList(self):
    if self.point != './':
      try:
        with open('%ssuccess.list' % self.point) as fp: file_list = fp.readlines()
      except IOError:
        print 'Folder does not exist!'; exit()
      for i in range(len(file_list)): file_list[i] = file_list[i].strip('\n')
      file_list = [el for el in file_list if self.filename in el]
    else:
      fold_list, file_list = [], []
      if self.online and os.path.exists(self.root + time.strftime("%Y/%Y%m%d/")): self.efolder = time.strftime("%Y%m%d")
      for Y in range(int(self.folder[0:4]), int(self.efolder[0:4])+1): fold_list.extend(os.listdir(self.root+str(Y)+'/'))
      fold_list = [el for el in fold_list if (int(self.folder) <= int(el) <= int(self.efolder))]
      for n in range(len(fold_list)): fold_list[n] = self.root+fold_list[n][0:4]+'/'+fold_list[n]+'/'
      fold_list.sort()
      first_date = True
      for folder in fold_list:
        if os.path.isdir(folder):
          files = [el for el in os.listdir(folder) if self.filename in el]
          if first_date:
            files = [el for el in files if int(string.split(el,'.')[0]) > self.tfrom]
            first_date = False
          for el in files: file_list.append(folder+el)
        else:
          print 'File system warning: %s is not a folder' % folder;   sys.exit(1)
      file_list.sort()

    if self.listonly:
      if self.online:
        for f in file_list[-self.nfile:]: print f
      else:
        for f in file_list: print f;
      sys.exit(0)
    if self.online:
      if self.lastfile != file_list[-1]:
        self.lastfile = file_list[-1]
        return file_list[-self.nfile:]
      else: return []
    else:   return file_list


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class Histogram:
  nbins    = 16384   # number of bins in a histogram
  tmin     = 10.0    # minimum data acquisition time for one file , s
  EME      = False
  hps      = ROOT.TH1I('hps','',nbins, 0, nbins)


  def __init__(self, todo):
    self.nfile, self.prompt, self.tokeV, self.folder = todo.nfile, todo.prompt, todo.tokeV, todo.point
    cfg = ConfigParser.ConfigParser(); cfg.read(todo.cfg_file)
    if cfg.has_option('scale', 'cdif'):  self.cdif = cfg.getfloat('scale', 'cdif')
    else:                                self.cdif = 0.01
    if   todo.edge and 'BEPC'   in todo.orig:
      try: from bepc.bepc     import EDGE;  self.EME = EDGE(todo.scalefile, todo.cfg_file)
      except ImportError: print 'BEPC is not yet implemented'; exit(0)
    elif todo.edge and 'VEPP2K' in todo.orig:
      try: from vepp2k.vepp2k import EDGE;  self.EME = EDGE(todo.scalefile, todo.cfg_file, todo.point)
      except ImportError: print 'V2K is not yet implemented'; exit(0)
    if self.tokeV:
      if 'HPGe' in todo.orig:
        from scale.isotopes import Isotopes
        self.CALIBRATION = Isotopes(todo.scalefile, todo.cfg_file, 'calibration')
      else:
        from scale.scale import Scale
        self.CALIBRATION = Scale(todo.scalefile, todo.cfg_file, 'calibration')
    elif not self.EME:
      self.cv = ROOT.TCanvas('cv','HPGe spectrum', 2, 2, 1002, 1002)

  def __del__(self): pass

  def Go(self,flist):
    SPEC = DataFile()
    while len(flist):
      n, self.LiveT, filechain = 0, 0.0, []
      self.hps.Scale(0)
      while n < self.nfile:
        if SPEC.ReadData(flist[0], self.nbins, self.tmin):
          if n==0: self.UTB = SPEC.utb
          elif (SPEC.utb - self.UTE) > 1800:  break

          self.UTE = SPEC.ute;  self.LiveT += SPEC.tLive;  n += 1;  filechain.append(flist[0])
          for nbin in range(self.nbins): self.hps.SetBinContent(nbin, self.hps.GetBinContent(nbin) + SPEC.DATA[nbin])
          print 'Add:', flist[0]
        else:
          print 'Skip:', flist[0]
        flist.pop(0)
        if len(flist)==0: break

      threshold = 1000
      for nbin in range(threshold):
        self.hps.SetBinContent(nbin, 0); self.hps.SetBinError(nbin, 0)
      for nbin in range(threshold, self.nbins):
        N = self.hps.GetBinContent(nbin); self.hps.SetBinError(nbin, (N + (self.cdif*N)**2)**0.5 )

      H,S = divmod(self.LiveT, 3600); M,S = divmod(S, 60)

      if   len([f for f in filechain if 'elec' in f]) == len(filechain): ptype = 'electron'
      elif len([f for f in filechain if 'posi' in f]) == len(filechain): ptype = 'positron'
      else:                                                              ptype = 'porridge'

      self.sname  = ptype + ': '
      self.sname += time.strftime('%Y.%m.%d [%H:%M:%S -', time.localtime(self.UTB))
      self.sname += time.strftime(' %H:%M:%S] %Y.%m.%d.', time.localtime(self.UTE))
      self.sname += ' Live-time: %d hours %d min %d s (%d files).' % (H, M, S, n)
      self.hps.SetNameTitle('hps',self.sname);  self.hps.SetEntries(self.hps.Integral());  self.hps.SetLineColor(ROOT.kBlack)
      if self.tokeV:
        R = self.CALIBRATION.Do(self.UTB, self.UTE, SPEC.PB5, self.hps, ptype)
      elif self.EME:
        R = self.EME.Go(self.UTB, self.UTE, self.hps, filechain, SPEC.Grate)
      else:
        R = 0
        self.hps.GetXaxis().SetTitle('E_{#gamma}, channels')
        self.cv.cd(); self.cv.SetGrid(); self.hps.Draw('HIST'); self.cv.Modified(); self.cv.Update()
      if self.prompt: raw_input('Press <Enter> to proceed')
    return R


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class DataFile:
  DATA, HAT, PB5         = [], {}, []
  utb, ute, Tlive, pType = 0, 0, 0, 'None'
#  TZ = 8 * 3600 # China TimeZone: GMT +8

  def ReadHat(self,ifile,tmin): # reads hat of spectrum
    def UnixTime(t,d):
      H,S = divmod(t,3600); M,S = divmod(S,60)
      return time.mktime(time.strptime("%8d%02d%02d%02d" % (d, H, M, S), "%Y%m%d%H%M%S")) # - self.TZ
    del self.DATA[:];    self.HAT.clear();    del self.PB5[:]

    with gzip.open(ifile,'rb') as fp: self.DATA = fp.readlines()
    while self.DATA[0][0] == '#':
      T = self.DATA.pop(0)[1:].strip().split()
      if 'PB5' in T[0]:
        for v in T[2:]: self.PB5.append(float(v))
      else:
        try:                self.HAT[T[0]] = float(T[1])
        except ValueError:  self.HAT[T[0]] = T[1]

#    self.PB5 = [0.800, 0.900, 1.100, 1.500, 1.800, 2.000, 2.500, 2.800, 3.200, 3.700, 4.500, 5.200, 5.400, 5.600]

    if self.HAT['Tlive']>tmin:
      self.tLive = self.HAT['Tlive']
      self.utb   = UnixTime(self.HAT['Begin'], self.HAT['Date'])
      self.ute   = UnixTime(self.HAT['End'],   self.HAT['eDate'])
      if self.HAT.has_key('pType'):   self.pType = self.HAT['pType']
      else:                           self.pType = 'None'
      if self.HAT.has_key('GRATING'): self.Grate = int(self.HAT['GRATING'])
      else:                           self.Grate = 0
      return True
    else:
      print 'Data acquisition time is to short, only %d seconds' % int(self.HAT['Tlive'])
      return False


  def ReadData(self,ifile,nbins,tmin): # reads one file spectrum

    if self.ReadHat(ifile,tmin) and (len(self.DATA) == nbins):
      for i in range(len(self.DATA)): self.DATA[i] = int(self.DATA[i])
      return True
    else:  print 'Wrong number of data lines in %s' % ifile; exit(1)


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


def main(argv):
  todo  = ToDo()
  flist = todo.GetList()
  try:
    if todo.escan:
      if 'BEPC' in todo.orig:
        from bepc.energy_scan import Energy_Scan
        A = Energy_Scan()
        A.Go(flist,todo)
      elif 'VEPP2K' in todo.orig:
        from vepp2k.vepp2k import EMSResults
        A = EMSResults(todo.cfg_file, todo.point)
        A.ShowRunInfo()
        A.EnergySpread()
        raw_input()
    elif todo.generate:
      if 'VEPP2K' in todo.orig:
        T = raw_input('Do you understand what you are going to do?')
        if 'Y' in T or 'y' in T:
          from vepp2k.vepp2k import Points_Splitter
          A = Points_Splitter()
          A.Go(flist,todo)
    elif todo.verify:
      from scale.scale import Scale
      A = Scale(todo.scalefile, todo.cfg_file, 'verification')
      A.Check_Scale(todo.energy, todo.prompt)
    else:
      A = Histogram(todo)
      while True:
        A.Go(flist)
        while todo.online:
          time.sleep(15.0)
          flist=todo.GetList()
          if len(flist)>0:
            time.sleep(5.0)
            break
        else:
          break
      raw_input('All done. <Enter> to quit.')
  except KeyboardInterrupt: print '\nExecution is interrupted'
#  del A
  exit()


if __name__ == "__main__": main(sys.argv)

