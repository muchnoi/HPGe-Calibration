#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ROOT, time

toff = 0 # -8*3600

def Create_html(le,lp):
  htmlhead = '''<html>
  <head>  
  <meta http-equiv="Content-Type" content="text/html; charset=UTF-8">
  <meta name="author" content="Nickolai Muchnoi">
  <meta http-equiv="refresh" content="30">
  <title>BEPC-II Beam Energy Measurement System Status</title>
  </head>

  <style type="text/css">
  <!--
  BODY      {font-family: Arial, Helvetica, sans-serif; font-size: 24px; font-weight: normal; color: #ffffff}
  TABLE     {font-family: Arial, Helvetica, sans-serif; font-size: 24px; font-weight: normal; color: #ffffff}
  A         {color: #ffffff; background-color: transparent; text-decoration: none;}
  A:visited {color: #ffffff; background-color: transparent; text-decoration: none;}
  A:hover   {color: #ffffff; background-color: transparent; text-decoration: underline;}
  -->
  </style>

  <SCRIPT type=text/javascript>
  function zoom(n,x,y)
  {
  window.open(n, '_blank', 'width='+x+',height='+y+',status=no,menubar=no,toolbar=no,location=no,scrollbars=no,resizable=yes');
  }
  </SCRIPT>

  <body bgcolor="#404040">
  <div align='center'>
  <H4><font color=#ffff40>Beam Energy Measurement System Status </font>
  <font color=#40ffff>(%s)</font></H4>
  ''' % (time.ctime())
  E_bt = time.ctime(le['t']-toff-le['dt']);  E_et = time.ctime(le['t']-toff+le['dt'])
  E_EB = le['E'];          E_dE = le['dE'];  E_BS = le['S'];         E_dS = le['dS']
  P_bt = time.ctime(lp['t']-toff-lp['dt']);  P_et = time.ctime(lp['t']-toff+lp['dt'])
  P_EB = lp['E'];          P_dE = lp['dE'];  P_BS = lp['S'];          P_dS = lp['dS']
 
  htmltable = '''
  <table width=1000 border=1 bordercolor=#ffffff cellpadding=4 cellspacing=0>
  <tr>
  <td>&nbsp;</td>
  <td nosave align=center> <font color=#40ff40><a href='images/E.png' target='_blank'>Electrons</a></font> </td>
  <td nosave align=center> <font color=#ff4040><a href='images/P.png' target='_blank'>Positrons</a></font> </td>
  </tr><tr>
  <td>Energy, MeV:</td>
  <td nosave align=center> <font color=#40ff40>%8.3f ± %8.3f</font> </td>
  <td nosave align=center> <font color=#ff4040>%8.3f ± %8.3f</font> </td>
  </tr><tr>
  <td>Energy spread, keV:</td>
  <td nosave align=center> <font color=#40ff40>%5.0f ± %5.0f</font> </td>
  <td nosave align=center><font color=#ff4040> %5.0f ± %5.0f</font> </td>
  </tr><tr>
  <td>Measured from:</td>
  <td nosave align=center> <font color=#40ff40>%s</font> </td>
  <td nosave align=center> <font color=#ff4040>%s</font> </td>
  </tr><tr>
  <td>Measured until:</td>
  <td nosave align=center> <font color=#40ff40>%s</font> </td>
  <td nosave align=center> <font color=#ff4040>%s</font> </td>
  </tr>

  </table>
  ''' % (E_EB, E_dE, P_EB, P_dE, E_BS, E_dS, P_BS, P_dS, E_bt, P_bt, E_et, P_et)
#  ''' % (E_EB, E_dE, P_EB, P_dE, E_BS, E_dS, P_BS, P_dS, E_bt, P_bt, E_et, P_et)
  htmlfoot = '''
  <img src='images/in-time.png' align='center', border=0, hspace=10, vspace=10>  
  </BODY></html>
  '''

  with open('index.html','w+') as f:
    f.write(htmlhead)
    f.write(htmltable)
    f.write(htmlfoot)
    

def Get_Graph(filename):
  t, dt = int(time.time()), 12*3600
  L, G1, G2, N = [], ROOT.TGraphErrors(), ROOT.TGraphErrors(),  0
  with open(filename,'r')  as f: lines = f.readlines()
#  with open(filename,'w+') as f:
  for line in lines[1:]:
    OK, fields = False, line.strip().split()
    if int(fields[0]) > t-dt or line==lines[-1]:
      L.append({'t':int(fields[0])+toff,'dt':int(fields[1]),
                'E':float(fields[2]),'dE':float(fields[3]),     # BEMS Energy
                'S':float(fields[4]),'dS':float(fields[5]),     # BEMS Spread
                'B':float(fields[6]),'dB':float(fields[7])})    # BEPC Energy
      t_last, e_last = int(fields[0]), float(fields[2])
      
#        f.write(line)
  for R in L: 
    if R['t'] > t-dt:
      G1.SetPoint(N, R['t'],  R['E']); G1.SetPointError(N, R['dt'], R['dE'])
      G2.SetPoint(N, R['t'],  R['B']); G2.SetPointError(N, R['dt'], R['dB']); N+=1
  if N: return (G1, G2, t_last, e_last, L[-1])
  else: return (0,0,0,0,L[-1])

def Energy_Plot():
  EG1, EG2, te, ee, le = Get_Graph('E.results')
  PG1, PG2, tp, ep, lp = Get_Graph('P.results')
  
  MG = ROOT.TMultiGraph()
  if EG2: EG2.SetMarkerStyle(20); EG2.SetMarkerColor(ROOT.kGray+1);  EG2.SetLineColor(ROOT.kGray+1);  EG2.SetLineWidth(1); MG.Add(EG2)
  if PG2: PG2.SetMarkerStyle(20); PG2.SetMarkerColor(ROOT.kGray+1);  PG2.SetLineColor(ROOT.kGray+1);  PG2.SetLineWidth(1); MG.Add(PG2)
  if EG1: EG1.SetMarkerStyle(20); EG1.SetMarkerColor(ROOT.kGreen+2); EG1.SetLineColor(ROOT.kGreen+2); EG1.SetLineWidth(2); MG.Add(EG1)
  if PG1: PG1.SetMarkerStyle(20); PG1.SetMarkerColor(ROOT.kRed+2);   PG1.SetLineColor(ROOT.kRed+2);   PG1.SetLineWidth(2); MG.Add(PG1)
 
  cv = ROOT.TCanvas('cvt','Compton Beam Energy Monitor',0,0,1000,800);  cv.cd();  cv.SetGrid()
  MG.Draw('AP');  MG.GetYaxis().SetDecimals(); MG.GetYaxis().SetTitle('beam energy, MeV'); MG.GetYaxis().SetTitleOffset(1.2)
  MG.GetXaxis().SetTitle('Beijing time');  MG.GetXaxis().SetTimeDisplay(1);  MG.GetXaxis().SetTimeOffset(0);  MG.GetXaxis().SetTimeFormat("%H:%M")
  
  lg = ROOT.TLegend(0.8, 0.8, 0.98, 0.98, '', 'brNDC');  lg.SetFillColor(ROOT.kGreen-9)
  lg.AddEntry(PG2, 'BEPC energy', 'lpe')
  lg.AddEntry(EG1, 'BEMS: e^{-}', 'lpe')
  lg.AddEntry(PG1, 'BEMS: e^{+}', 'lpe')
  lg.Draw('SAME')
  cv.SetFillColor(ROOT.kGray);  cv.GetFrame().SetFillColor(ROOT.kGreen-10)
  cv.Modified();  cv.Update();  cv.SaveAs('in-time.png')
  return le, lp

print 'Creating plot'
le,lp = Energy_Plot()

print 'Creating html'
Create_html(le,lp)


