class Atlas:
  def __init__(self):
    self.me = 510.998910  # electron rest mass, keV
    self.atlas={}

    self.atlas['  H1 ']   =   [{'W':2223.250, 'dW':.050 }]
    self.atlas[' K40 ']   =   [{'W':1460.750, 'dW':.060 }]
    self.atlas['annig']   =   [{'W': 510.998, 'dW':.060 }]

    self.atlas[' O16 ']   =   [{'W':6129.075, 'dW':.050 }]
#    self.atlas[' O16 ']   =   [{'W':6128.630, 'dW':.040  }]
#    self.atlas['Ale81']   =   [{'W':6129.240, 'dW':.050  }]
#    self.atlas['She82']   =   [{'W':6129.142, 'dW':.032  }]
#    self.atlas['Ken86']   =   [{'W':6129.119, 'dW':.040  }]

    self.atlas['Cs137']   =   [{'W': 661.657, 'dW':.003 }]

    self.atlas['Co60 ']   =   [{'W':1173.228, 'dW':.003 }]
    self.atlas['Co60 '].append({'W':1332.492, 'dW':.004 })

    self.atlas['Na24 ']   =   [{'W':1368.625, 'dW':.005 }]
    self.atlas['Na24 '].append({'W':2754.008, 'dW':.011 })
    self.atlas['Na24 '].append({'W':3866.190, 'dW':.011 })

    self.atlas['Al28 ']   =   [{'W':1778.969, 'dW':.012 }]

    self.atlas['Tl208']   =   [{'W': 238.600, 'dW':.100 }]
    self.atlas['Tl208'].append({'W': 241.000, 'dW':.100 })
    self.atlas['Tl208'].append({'W': 277.400, 'dW':.100 })


    self.atlas['Co56 ']   =   [{'W': 846.764, 'dW':.020 }]
    self.atlas['Co56 '].append({'W':1238.274, 'dW':.020 })

    self.atlas['Mn56 ']   =   [{'W': 846.764, 'dW':.020 }]
    self.atlas['Mn56 '].append({'W':1810.726, 'dW':.040 })


    self.atlas['Pb212'] = []
    self.atlas['Pb212'].append({'W': 238.632, 'dW':.002, 'I':43.3, 'dI':0.40 })
    self.atlas['Pb212'].append({'W': 300.087, 'dW':.010, 'I':3.28, 'dI':0.30 })

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


    self.atlas['Tl208'] = []
    self.atlas['Tl208'].append({'W': 211.400, 'dW':.150, 'I':.00018, 'dI':.00001 })
    self.atlas['Tl208'].append({'W': 233.360, 'dW':.150, 'I':.00290, 'dI':.00030 })
    self.atlas['Tl208'].append({'W': 252.610, 'dW':.100, 'I':.00780, 'dI':.00030 })
    self.atlas['Tl208'].append({'W': 277.358, 'dW':.010, 'I':.06370, 'dI':.00120 })
    self.atlas['Tl208'].append({'W': 485.950, 'dW':.150, 'I':.00047, 'dI':.00006 })
    self.atlas['Tl208'].append({'W': 510.770, 'dW':.100, 'I':.22500, 'dI':.00500 })
    self.atlas['Tl208'].append({'W': 583.187, 'dW':.005, 'I':.85100, 'dI':.00900 })
    self.atlas['Tl208'].append({'W': 587.700, 'dW':.300, 'I':.00042, 'dI':.00017 })
    self.atlas['Tl208'].append({'W': 650.100, 'dW':.300, 'I':.00036, 'dI':.00006 })
    self.atlas['Tl208'].append({'W': 705.200, 'dW':.300, 'I':.00022, 'dI':.00006 })
    self.atlas['Tl208'].append({'W': 722.040, 'dW':.120, 'I':.00240, 'dI':.00020 })
    self.atlas['Tl208'].append({'W': 748.700, 'dW':.200, 'I':.00045, 'dI':.00003 })
    self.atlas['Tl208'].append({'W': 763.130, 'dW':.080, 'I':.01780, 'dI':.00060 })
    self.atlas['Tl208'].append({'W': 821.200, 'dW':.200, 'I':.00042, 'dI':.00003 })
    self.atlas['Tl208'].append({'W': 860.530, 'dW':.020, 'I':.12600, 'dI':.00200 })
    self.atlas['Tl208'].append({'W': 883.300, 'dW':.200, 'I':.00031, 'dI':.00003 })
    self.atlas['Tl208'].append({'W': 927.600, 'dW':.200, 'I':.00128, 'dI':.00020 })
    self.atlas['Tl208'].append({'W': 982.700, 'dW':.200, 'I':.00190, 'dI':.00020 })
    self.atlas['Tl208'].append({'W':1093.900, 'dW':.200, 'I':.00420, 'dI':.00030 })
    self.atlas['Tl208'].append({'W':1160.800, 'dW':.300, 'I':.00011, 'dI':.00003 })
    self.atlas['Tl208'].append({'W':1185.200, 'dW':.300, 'I':.00017, 'dI':.00006 })
    self.atlas['Tl208'].append({'W':1282.800, 'dW':.300, 'I':.00050, 'dI':.00006 })
    self.atlas['Tl208'].append({'W':2614.511, 'dW':.010, 'I':.99900, 'dI':.00400 })

    self.atlas['Bi212'] = []
    self.atlas['Bi212'].append({'W': 39.8580, 'dW':.004, 'I':.01020, 'dI':.00005 })
    self.atlas['Bi212'].append({'W': 288.080, 'dW':.070, 'I':.00314, 'dI':.00011 })
    self.atlas['Bi212'].append({'W': 327.940, 'dW':.060, 'I':.00120, 'dI':.00004 })
    self.atlas['Bi212'].append({'W': 433.700, 'dW':.200, 'I':.00012, 'dI':.00003 })
    self.atlas['Bi212'].append({'W': 452.830, 'dW':.100, 'I':.00333, 'dI':.00018 })
    self.atlas['Bi212'].append({'W': 473.600, 'dW':.200, 'I':.00044, 'dI':.00003 })
    self.atlas['Bi212'].append({'W': 727.200, 'dW':.050, 'I':.00500, 'dI':.00110 })
    self.atlas['Bi212'].append({'W': 785.370, 'dW':.080, 'I':.01094, 'dI':.00024 })
    self.atlas['Bi212'].append({'W': 893.430, 'dW':.090, 'I':.00381, 'dI':.00012 })
    self.atlas['Bi212'].append({'W': 952.100, 'dW':.200, 'I':.00140, 'dI':.00010 })
    self.atlas['Bi212'].append({'W':1078.600, 'dW':.600, 'I':.00630, 'dI':.00040 })
    self.atlas['Bi212'].append({'W':1512.660, 'dW':.070, 'I':.00278, 'dI':.00015 })
    self.atlas['Bi212'].append({'W':1620.735, 'dW':.100, 'I':.01490, 'dI':.00040 })
    self.atlas['Bi212'].append({'W':1679.300, 'dW':.100, 'I':.00060, 'dI':.00010 })
    self.atlas['Bi212'].append({'W':1801.000, 'dW':.500, 'I':.00110, 'dI':.00020 })
    self.atlas['Bi212'].append({'W':1805.800, 'dW':.100, 'I':.00900, 'dI':.00020 })
#    self.atlas['Bi212'].append({'W': 893.408, 'dW':.005, 'I':.00381, 'dI':.00012 })
#    self.atlas['Bi212'].append({'W': 620.400, 'dW':.300, 'I':.00004, 'dI':.00001 })
#    self.atlas['Bi212'].append({'W': 580.500, 'dW':.300, 'I':.00001, 'dI':.00001 })
#    self.atlas['Bi212'].append({'W': 492.700, 'dW':.300, 'I':.00007, 'dI':.00003 })




    """
    for nucl, lines in self.atlas.iteritems():
      n,T = 0,[]
      for line in lines:
        n+=1
        if line['W']>2000.0: # Add single- and double-escape peaks
          T.append({'Key':line['Key']+'se', 'W': line['W']-   self.me, 'dW':line['dW']  })
          T.append({'Key':line['Key']+'de', 'W': line['W']-2.*self.me, 'dW':line['dW']  })
      for elem in T: lines.append(elem)
    """

# The following were measured by using the previous to calibrate HPGe
#    self.atlas['Ac228']   =   [{'W': 911.282, 'dW':.004  , 'Compton':.0099 }]
#    self.atlas['Ac228'].append({'W': 964.799, 'dW':.010   })


#      print nucl
#      for line in lines:
#        print "%8s\t%8.3f" % (line['Key'], line['W'])
#      raw_input()

#        for next in lines[n:]: # Possible summs of two gammas
#          if line['W']+next['W'] > 2615.0:
#            T.append({'Key':line['Key']+'+'+next['Key'], 'W': line['W']+next['W'],
#                     'dW':sqError([line['dW'],next['dW']])  })


  def __del__(self):
    del self.atlas



