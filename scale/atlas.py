class Atlas:
  def __init__(self):
    self.me = 510.998910  # electron rest mass, keV
    self.atlas={}

    self.atlas['Cs137']   =   [{'W': 661.659, 'dW':.003 }] # [{'W': 661.657, 'dW':.003 }]
    self.atlas['Co60 ']   =   [{'W':1173.240, 'dW':.003 }] # [{'W':1173.228, 'dW':.003 }]
    self.atlas['Co60 '].append({'W':1332.508, 'dW':.004 }) # ({'W':1332.492, 'dW':.004 })

    self.atlas['H1   ']   =   [{'W':2223.250, 'dW':.050 }]
    self.atlas['K40  ']   =   [{'W':1460.750, 'dW':.060 }]
    self.atlas['annig']   =   [{'W': 510.998, 'dW':.060 }]

    self.atlas['O16  ']   =   [{'W':6129.12, 'dW':.050 }]
    self.atlas['O16  '].append({'W':self.atlas['O16  '][0]['W']-  self.me, 'dW':.050 })
    self.atlas['O16  '].append({'W':self.atlas['O16  '][0]['W']-2*self.me, 'dW':.050 })

    """
    self.atlas[' O16 ']   =   [{'W':6128.630, 'dW':.040  }]
    self.atlas['Ale81']   =   [{'W':6129.240, 'dW':.050  }]
    self.atlas['She82']   =   [{'W':6129.142, 'dW':.032  }]
    self.atlas['Ken86']   =   [{'W':6129.119, 'dW':.040  }]

    self.atlas['Na24 ']   =   [{'W':1368.625, 'dW':.005 }]
    self.atlas['Na24 '].append({'W':2754.008, 'dW':.011 })
    self.atlas['Na24 '].append({'W':3866.190, 'dW':.011 })
    self.atlas['Al28 ']   =   [{'W':1778.969, 'dW':.012 }]

    self.atlas['Co56 ']   =   [{'W': 846.764, 'dW':.020 }]
    self.atlas['Co56 '].append({'W':1238.274, 'dW':.020 })

    self.atlas['Mn56 ']   =   [{'W': 846.764, 'dW':.020 }]
    self.atlas['Mn56 '].append({'W':1810.726, 'dW':.040 })

    self.atlas['In116'] = []
    self.atlas['In116'].append({'W': 416.860, 'dW':.030, 'I':27.7, 'dI':.120 })
    self.atlas['In116'].append({'W': 818.718, 'dW':.021, 'I':11.5, 'dI':.040 })
    self.atlas['In116'].append({'W':1097.326, 'dW':.022, 'I':56.2, 'dI':.110 })
    self.atlas['In116'].append({'W':1293.558, 'dW':.015, 'I':84.4, 'dI':.170 })
    self.atlas['In116'].append({'W':1507.670, 'dW':.040, 'I':10.0, 'dI':.030 })
    self.atlas['In116'].append({'W':2112.312, 'dW':.022, 'I':15.5, 'dI':.040 })

    self.atlas['FePNG'] = []
    self.atlas['FePNG'].append({'W': 139.750, 'dW':.030, 'I':.005, 'dI':.001 })
    self.atlas['FePNG'].append({'W': 198.380, 'dW':.030, 'I':.005, 'dI':.001 })
    self.atlas['FePNG'].append({'W': 352.332, 'dW':.016, 'I':.273, 'dI':.001 })
    self.atlas['FePNG'].append({'W': 366.737, 'dW':.016, 'I':.050, 'dI':.001 })
    self.atlas['FePNG'].append({'W': 569.852, 'dW':.021, 'I':.014, 'dI':.001 })
    self.atlas['FePNG'].append({'W': 691.914, 'dW':.016, 'I':.137, 'dI':.001 })
    self.atlas['FePNG'].append({'W':1612.770, 'dW':.030, 'I':.153, 'dI':.001 })
    self.atlas['FePNG'].append({'W':1725.255, 'dW':.024, 'I':.181, 'dI':.001 })
    self.atlas['FePNG'].append({'W':5920.250, 'dW':.080, 'I':.225, 'dI':.001 })
    self.atlas['FePNG'].append({'W':6018.290, 'dW':.080, 'I':.227, 'dI':.001 })
    self.atlas['FePNG'].append({'W':7278.830, 'dW':.100, 'I':.137, 'dI':.001 })
    self.atlas['FePNG'].append({'W':7631.050, 'dW':.090, 'I':.653, 'dI':.001 })
    self.atlas['FePNG'].append({'W':7645.480, 'dW':.090, 'I':.549, 'dI':.001 })
    """
# 232-Th decay chain

    self.atlas['Ac228'] = []
    self.atlas['Ac228'].append({'W': 911.209, 'dW':.006, 'I':.265 })
#    self.atlas['Ac228'].append({'W':964.768, 'dW':.007, 'I':.050 })
#    self.atlas['Ac228'].append({'W':968.968, 'dW':.005, 'I':.161 })

    self.atlas['Bi212'] = []
    self.atlas['Bi212'].append({'W': 727.330, 'dW':.010, 'I':6.74 })
    self.atlas['Bi212'].append({'W': 785.370, 'dW':.090, 'I':1.15, 'dI':.00024 })
    self.atlas['Bi212'].append({'W': 893.410, 'dW':.020, 'I':0.39, 'dI':.00024 })
    self.atlas['Bi212'].append({'W':1078.630, 'dW':.100, 'I':0.56, 'dI':.00024 })
    self.atlas['Bi212'].append({'W':1512.700, 'dW':.080, 'I':0.29, 'dI':.00024 })
    self.atlas['Bi212'].append({'W':1620.740, 'dW':.010, 'I':1.52, 'dI':.00040 })
    self.atlas['Bi212'].append({'W':1805.960, 'dW':.100, 'I':0.12, 'dI':.00040 })

    self.atlas['Tl208'] = []
    self.atlas['Tl208'].append({'W': 583.187, 'dW':.005, 'I':.85100, 'dI':.00900 })
    self.atlas['Tl208'].append({'W': 860.560, 'dW':.030, 'I':.12600, 'dI':.00200 })
    self.atlas['Tl208'].append({'W':2614.511, 'dW':.010, 'I':.99900, 'dI':.00400 })
#    self.atlas['Tl208'].append({'W':self.atlas['Tl208'][2]['W']  -self.me, 'dW':.010})
#    self.atlas['Tl208'].append({'W':self.atlas['Tl208'][2]['W']-2*self.me, 'dW':.010})

    self.atlas['Bi214'] = []
    self.atlas['Bi214'].append({'W': 609.312, 'dW':.007, 'I':46.1 })
    self.atlas['Bi214'].append({'W': 768.356, 'dW':.010, 'I':4.94 })
    self.atlas['Bi214'].append({'W': 806.174, 'dW':.018, 'I':1.22 })
    self.atlas['Bi214'].append({'W': 934.061, 'dW':.012, 'I':3.03 })
    self.atlas['Bi214'].append({'W':1120.287, 'dW':.010, 'I':15.1 })
    self.atlas['Bi214'].append({'W':1120.287, 'dW':.010, 'I':15.1 })
    self.atlas['Bi214'].append({'W':1155.190, 'dW':.020, 'I':1.63 })
    self.atlas['Bi214'].append({'W':1238.110, 'dW':.012, 'I':5.79 })
    self.atlas['Bi214'].append({'W':1280.960, 'dW':.020, 'I':1.43 })
    self.atlas['Bi214'].append({'W':1377.669, 'dW':.012, 'I':4.00 })
    self.atlas['Bi214'].append({'W':1401.500, 'dW':.040, 'I':1.27 })
    self.atlas['Bi214'].append({'W':1407.980, 'dW':.040, 'I':2.15 })
#    self.atlas['Bi214'].append({'W':1509.228, 'dW':.015, 'I':2.11 })
#    self.atlas['Bi214'].append({'W':1509.228, 'dW':.015, 'I':2.11 })
#    self.atlas['Bi214'].append({'W':1661.280, 'dW':.020, 'I':1.15 })
#    self.atlas['Bi214'].append({'W':1729.595, 'dW':.015, 'I':2.92 })
#    self.atlas['Bi214'].append({'W':1764.494, 'dW':.014, 'I':15.4 })
    self.atlas['Bi214'].append({'W':1847.420, 'dW':.025, 'I':2.11 })
    self.atlas['Bi214'].append({'W':2118.550, 'dW':.030, 'I':1.14 })
    self.atlas['Bi214'].append({'W':2204.210, 'dW':.040, 'I':5.08 })
    self.atlas['Bi214'].append({'W':2447.860, 'dW':.010, 'I':1.57 })
    """
    self.atlas['Tl208'].append({'W': 211.400, 'dW':.150, 'I':.00018, 'dI':.00001 })
    self.atlas['Tl208'].append({'W': 233.360, 'dW':.150, 'I':.00290, 'dI':.00030 })
    self.atlas['Tl208'].append({'W': 252.610, 'dW':.100, 'I':.00780, 'dI':.00030 })
    self.atlas['Tl208'].append({'W': 277.358, 'dW':.010, 'I':.06370, 'dI':.00120 })
    self.atlas['Tl208'].append({'W': 485.950, 'dW':.150, 'I':.00047, 'dI':.00006 })
    self.atlas['Tl208'].append({'W': 510.770, 'dW':.100, 'I':.22500, 'dI':.00500 })
    self.atlas['Tl208'].append({'W': 587.700, 'dW':.300, 'I':.00042, 'dI':.00017 })
    self.atlas['Tl208'].append({'W': 650.100, 'dW':.300, 'I':.00036, 'dI':.00006 })
    self.atlas['Tl208'].append({'W': 705.200, 'dW':.300, 'I':.00022, 'dI':.00006 })
    self.atlas['Tl208'].append({'W': 722.040, 'dW':.120, 'I':.00240, 'dI':.00020 })
    self.atlas['Tl208'].append({'W': 748.700, 'dW':.200, 'I':.00045, 'dI':.00003 })
    self.atlas['Tl208'].append({'W': 763.130, 'dW':.080, 'I':.01780, 'dI':.00060 })
    self.atlas['Tl208'].append({'W': 821.200, 'dW':.200, 'I':.00042, 'dI':.00003 })
    self.atlas['Tl208'].append({'W': 883.300, 'dW':.200, 'I':.00031, 'dI':.00003 })
    self.atlas['Tl208'].append({'W': 927.600, 'dW':.200, 'I':.00128, 'dI':.00020 })
    self.atlas['Tl208'].append({'W': 982.700, 'dW':.200, 'I':.00190, 'dI':.00020 })
    self.atlas['Tl208'].append({'W':1093.900, 'dW':.200, 'I':.00420, 'dI':.00030 })
    self.atlas['Tl208'].append({'W':1160.800, 'dW':.300, 'I':.00011, 'dI':.00003 })
    self.atlas['Tl208'].append({'W':1185.200, 'dW':.300, 'I':.00017, 'dI':.00006 })
    self.atlas['Tl208'].append({'W':1282.800, 'dW':.300, 'I':.00050, 'dI':.00006 })
    """



  def __del__(self):
    del self.atlas

