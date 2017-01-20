#!/usr/bin/env python
# -*- coding: utf-8 -*-
import ROOT
from isotopes import Isotopes
from numpy import asarray

scalefile = 'JPsi_2011.root'
conf_file = 'JPsi_2011.cfg'
testch = int(4150/0.45850)

C = Isotopes(scalefile, conf_file, 'application')
E = C.Get_Calibration(t=0, w=4000.0)

for k,v in E.iteritems():   E[k] = asarray(v)

NE  = len(E['T'])
CE = ROOT.kRed+1
TF = '#splitline{%b %d}{%H:%M}%F1970-01-01 00:00:00'

E1 = ROOT.TGraphErrors(NE, E['T'], E['W'], E['dT'], E['dW']); E1.SetTitle('CH. %d linear energy equivalent' % testch)
E2 = ROOT.TGraphErrors(NE, E['T'], E['C'], E['dT'], E['dC']); E2.SetTitle('CH. %d pulser energy correction' % testch) 
E3 = ROOT.TGraphErrors(NE, E['T'], E['R'], E['dT'], E['dR']); E3.SetTitle('CH. %d resolution - #sigma_{R}'  % testch)
E4 = ROOT.TGraphErrors(NE, E['T'], E['L'], E['dT'], E['dL']); E4.SetTitle('CH. %d resolution - #sigma_{L}'  % testch)
for GE in [E1,E2,E3,E4]: GE.SetMarkerStyle(20); GE.SetMarkerColor(CE); GE.SetLineColor(CE)

LEGEND = ROOT.TLegend(0.75, 0.75, 0.98, 0.95, '', 'brNDC');
LEGEND.AddEntry(E1, scalefile, 'lpe')

cv = ROOT.TCanvas('cv','scale & resolution', 2, 2, 1002, 1002)
cv.Divide(2,2)

cv.cd(1); cv.GetPad(1).SetGrid()
E1.Draw('AP'); E1.GetXaxis().SetTimeDisplay(1); E1.GetXaxis().SetTimeFormat(TF); E1.GetXaxis().SetLabelOffset(0.03)
E1.GetYaxis().SetDecimals(); E1.GetYaxis().SetTitle('keV'); E1.GetYaxis().SetTitleOffset(1.2)
LEGEND.Draw('SAME')

cv.cd(2); cv.GetPad(2).SetGrid()
E2.Draw('AP'); E2.GetXaxis().SetTimeDisplay(1); E2.GetXaxis().SetTimeFormat(TF); E2.GetXaxis().SetLabelOffset(0.03)
E2.GetYaxis().SetDecimals(); E2.GetYaxis().SetTitle('keV'); E2.GetYaxis().SetTitleOffset(1.2)
LEGEND.Draw('SAME')

cv.cd(3); cv.GetPad(3).SetGrid()
E3.Draw('AP'); E3.GetXaxis().SetTimeDisplay(1); E3.GetXaxis().SetTimeFormat(TF); E3.GetXaxis().SetLabelOffset(0.03)
E3.GetYaxis().SetDecimals(); E3.GetYaxis().SetTitle('keV'); E3.GetYaxis().SetTitleOffset(1.2)
LEGEND.Draw('SAME')

cv.cd(4); cv.GetPad(4).SetGrid()
E4.Draw('AP'); E4.GetXaxis().SetTimeDisplay(1); E4.GetXaxis().SetTimeFormat(TF); E4.GetXaxis().SetLabelOffset(0.03 )
E4.GetYaxis().SetDecimals(); E4.GetYaxis().SetTitle('keV'); E4.GetYaxis().SetTitleOffset(1.2)
LEGEND.Draw('SAME')

cv.Modified()
cv.Update()

raw_input()


exit(0)


