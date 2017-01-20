#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys, getopt, time, string, gzip, ConfigParser

def Usage():
  print '''
 ╔ Python script to process HPGe spectra. © 2005-2017 N.Yu.Muchnoi ════════════════╗
 ║                                              ╭────────────────────────────────╮ ║
 ║ Usage: %0000000012s [options]                │ Last update: January 18, 2017  │ ║
 ║                                              ╰────────────────────────────────╯ ║
 ║ List of options:                                                                ║
 ║ -h,          --help                : print this help message and exit,          ║
 ║ -i,          --interactive         : if set, script will prompt to proceed,     ║
 ║ -k,          --keV                 : use channels not keV,                      ║
 ║ -l,          --list                : just show the list of files,               ║
 ║ -f expr,     --file = expr         : file name(s) under specified folder(s),    ║
 ║ -n N,        --nf   = N            : put several files into one spectrum,       ║
 ║ -d YYYYMMDD, --folder  = YYYYMMDD  : date to start from (year, month, day),     ║
 ║ -e YYYYMMDD, --efolder = YYYYMMDD  : date to  end  with (year, month, day),     ║
 ║ -t HHMMSS,   --time    = HHMMSS    : time from which to start (hour, min, sec), ║
 ║ -s filename, --scale   = filename  : file to store/restore calibration results, ║
 ║ -c filename, --cfg     = filename  : file to read various parameters,           ║
 ║              --edge                : try to measure beam energy by Compton edge,║
 ║              --escan               : deal with beam energy scan experiment.     ║
 ╚═════════════════════════════════════════════════════════════════════════════════╝
 ''' % sys.argv[0].split('/')[-1];  sys.exit(0)

if ('-h' in sys.argv) or ('--help' in sys.argv): Usage() 
else: import ROOT

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class ToDo:

  def __init__(self):

    cfg = ConfigParser.ConfigParser(); cfg.read('_globals_.cfg')
    if cfg.has_option('globals', 'data_folder'):  self.root = cfg.get('globals', 'data_folder')
    if cfg.has_option('globals', 'data_origin'):  self.orig = cfg.get('globals', 'data_origin')
    else:                                         exit()
    try:
      Years = os.listdir(self.root); Years.sort()
    except OSError:
      print 'root folder does not exist!'
      sys.exit(1)
    Days  = os.listdir(self.root+Years[-1])
    if len(Days)==1: Days.extend(os.listdir(self.root+Years[-2]));  
    Days.sort()
    self.folder,   self.efolder   = Days[-2], Days[-1]
    self.tfrom,    self.nfile     = 0   , 1
    self.filename = self.lastfile = '' 
    self.tokeV    = self.online   = True
    self.prompt   = self.listonly = False
    self.edge     = self.escan    = False
    
    self.scalefile  = 'escale.root'
    self.cfg_file   = 'online.cfg'
    if len(sys.argv)>1: self.Options(sys.argv)

  def Options(self,argv):
    sopt = "d:e:f:t:n:s:c:ihkl"
    lopt = ["folder=","efolder=","file=","time=","nfiles=","scale=","cfg=","interactive","help","keV","list","edge","escan"]
    try:
      opts, args = getopt.getopt(argv[1:], sopt, lopt)
    except getopt.GetoptError:
      print 'Wrong option'; Usage(); sys.exit(2)
    date, time = False, False
    for opt, arg in opts:
      if   opt in ("-h", "--help"):        Usage()
      elif opt in ("-i", "--interactive"): self.prompt    = True;     self.online = False
      elif opt in ("-n", "--nf"):          self.nfile     = int(arg)
      elif opt in ("-d", "--folder"):      self.folder    = arg;      self.online = False; date = True
      elif opt in ("-e", "--efolder"):     self.efolder   = arg;      self.online = False
      elif opt in ("-f", "--file"):        self.filename  = arg;      self.online = False
      elif opt in ("-t", "--time"):        self.tfrom     = int(arg); self.online = False; time = True
      elif opt in ("-k", "--keV"):         self.tokeV     = False
      elif opt in ("-s", "--scale"):       self.scalefile = arg;  
      elif opt in ("-c", "--cfg"):         self.cfg_file  = arg;  
      elif opt in ("-l", "--list"):        self.listonly  = True
      elif opt in (      "--edge"):        self.edge      = True;     self.tokeV = False
      elif opt in (      "--escan"):       self.escan     = True;     self.online = False
    cfg = ConfigParser.ConfigParser(); cfg.read(self.cfg_file)
    if cfg.has_option('scale', 'file'): self.scalefile = cfg.get('scale', 'file'); self.online = False; date = True
    if cfg.has_option('scan', 'bdate'):    self.folder = cfg.get('scan', 'bdate'); self.online = False; date = True
    if cfg.has_option('scan', 'edate'):   self.efolder = cfg.get('scan', 'edate'); self.online = False

  def GetList(self):
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
        print 'File system warning: %s is not a folder' % folder
        sys.exit(1)
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
  tmin     = 100.0   # minimum data acquisition time for one file , s
  ROOT.gROOT.SetStyle("Plain"); ROOT.gROOT.ForceStyle()
  ROOT.gStyle.SetTitleBorderSize(0);  ROOT.gStyle.SetOptStat(0);    ROOT.gStyle.SetOptFit(0)
  hps      = ROOT.TH1I('hps','',nbins, 0, nbins)

  def __init__(self,todo):
    self.nfile, self.prompt = todo.nfile, todo.prompt
    self.tokeV, self.orig, self.EME    = todo.tokeV, todo.orig, todo.edge
    cfg = ConfigParser.ConfigParser(); cfg.read(todo.cfg_file)
    if cfg.has_option('scale', 'cdif'):  self.cdif = cfg.getfloat('scale', 'cdif')
    else:                                self.cdif = 0.01
    if   self.EME and 'BEPC' in todo.orig:
      try:                from bepc.bepc import EDGE;      self.EME = EDGE(todo.scalefile, todo.cfg_file)
      except ImportError: print 'BEPC is not yet implemented'; exit(0)
    elif self.EME and  'VEPP2K' in todo.orig:
      try:                from vepp2k.vepp2k import EDGE;  self.EME = EDGE(todo.scalefile, todo.cfg_file)
      except ImportError: print 'V2K is not yet implemented'; exit(0)
    if self.tokeV: 
      from scale.isotopes import Isotopes
      self.CALIBRATION = Isotopes(todo.scalefile, todo.cfg_file, 'calibration')
    elif not self.EME: 
      self.cv = ROOT.TCanvas('cv','HPGe spectrum', 2, 2, 1002, 1002)

  def __del__(self): pass

  def Go(self,flist):
    DATA     = DataFile()
    while len(flist)>=self.nfile:
      n, self.LiveT, filechain = 0, 0.0, []
      self.hps.Scale(0)
      while n < self.nfile:
        spectrum = DATA.ReadData(flist[0], self.nbins, self.tmin)
        if spectrum:
          if n==0: self.UTB = DATA.utb
          self.UTE = DATA.ute;  self.LiveT += DATA.tLive;  n += 1;  filechain.append(flist[0])
          for nbin in range(self.nbins): self.hps.SetBinContent(nbin, self.hps.GetBinContent(nbin) + spectrum[nbin])
          print 'Add:', flist[0]
        else:
          print 'Skip:', flist[0]
        flist.pop(0)

      for nbin in range(self.nbins):     
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
      if self.tokeV: self.CALIBRATION.Do(self.UTB, self.UTE, DATA.PB5, self.hps, ptype)
      elif self.EME: 
        R = self.EME.Go(self.UTB, self.UTE, self.hps, filechain)
        if self.prompt: raw_input('Press <Enter> to proceed')
        return R
      else:
        self.hps.GetXaxis().SetTitle('E_{#gamma}, channels')
        self.cv.cd(); self.cv.SetGrid(); self.hps.Draw('HIST'); self.cv.Modified(); self.cv.Update()
      if self.prompt: raw_input('Press <Enter> to proceed')
    return 0


#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
class DataFile:

  TZ = 8 * 3600 # China TimeZone: GMT +8
  
  def ReadData(self,ifile,nbins,tmin): # reads one file spectrum

    def UnixTime(t,d): 
      H,S = divmod(t,3600); M,S = divmod(S,60)
      return time.mktime(time.strptime("%8d%02d%02d%02d" % (d, H, M, S), "%Y%m%d%H%M%S")) # - self.TZ

    self.PB5, self.utb, self.ute, self.Tlive, self.pType = [], 0, 0, 0, 'None'
    HAT, counts  = {}, []
    with gzip.open(ifile,'rb') as fp: DATA = fp.readlines()
    for line in DATA:
      if line[0] == '#': 
        T = line[1:].strip().split()
        if 'PB5 VOLTS' in line:
          for v in T[2:]: self.PB5.append(float(v))
        else:    
          try:                HAT[T[0]] = float(T[1])
          except ValueError:  HAT[T[0]] = T[1]
      else: counts.append(int(line))
    

    if HAT['Tlive']>tmin:
      self.tLive = HAT['Tlive']
      self.utb, self.ute = UnixTime(HAT['Begin'], HAT['Date']), UnixTime(HAT['End'], HAT['eDate'])
      if HAT.has_key('pType'): self.pType = HAT['pType']
      else:                    self.pType = 'None' 
      if nbins == len(counts): return counts
      else: print 'Wrong number of data lines in %s' % ifile; exit(1)
    else:
      print 'Data acquisition time is to short, only %d seconds' % int(HAT['Tlive'])
      return False



#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


def main(argv):
  todo  = ToDo()
  flist = todo.GetList()
  if todo.escan:
    from bepc.energy_scan import Energy_Scan
    InTime = Energy_Scan()
    InTime.Go(flist,todo)
  else:
    SP = Histogram(todo)
    try:
      while True:
        SP.Go(flist)
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
  sys.exit(0)
  

if __name__ == "__main__": main(sys.argv)  

