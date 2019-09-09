#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os, fnmatch
import sys
from core.curves import make_curve

def main():
  nbCoeff = 6
  #directory = '.' #directory from which the script is executed
  directory = 'raws'
  data_curves = dict()
  for file in fnmatch.filter(os.listdir(directory), '*.raw'):
    processRawFile(directory + "/" + file, nbCoeff, data_curves)
  make_curve(data_curves, directory)

def processRawFile(file, nbCoeff, data_curves):
  data = np.genfromtxt( file , usecols =(1, 2), delimiter=";", skip_header=65, encoding="latin1", dtype=None)

  depla = []
  charge = []
  chargeN = []
  lissage = []
  PPmax = []
  aire = []

  for i, el in enumerate(data):
    depla.append(float(el[0].replace(',','.')))
    charge.append(float(el[1].replace(',','.')))
    chargeN.append(float(el[1].replace(',','.'))*9.80665)

  for i in range(nbCoeff):
    lissage.append(0)
    PPmax.append(0)
    aire.append(0)

  div = 0;
  for i in range(1, nbCoeff+1):
    div += i * 2

  PMax = 0;
  lissage_mean = 0
  for i in range(nbCoeff, len(depla) - nbCoeff):
    tmp = 0
    for j in range(1, nbCoeff+1):
      #print(i, j, nbCoeff - j + 1, chargeN[i - j], chargeN[i + j])
      tmp += (chargeN[i - j] + chargeN[i + j])*(nbCoeff - j + 1)
    tmp = tmp/div
    if tmp > PMax:
      PMax = tmp
    lissage.append(tmp)
    lissage_mean += tmp
  lissage_mean /= len(lissage)


  for i in range(nbCoeff, len(lissage)):
    tmp = lissage[i]*100/PMax
    PPmax.append(tmp)

  for i in range(nbCoeff, len(lissage) - 1):
    tmp = (lissage[i] + lissage[i+1]) / 2 * (depla[i+1] - depla[i])
    aire.append(tmp)


  KBond = 0
  KDbond = 0
  yb = 0
  ya = 0
  xa = 0
  xb = 0
  incr = 0
  for i in range(len(lissage)):
    if PPmax[i] > 30:
      yb = lissage[i]
      xb = depla[i]
      break

  for i in range(len(lissage)):
    if PPmax[i] > 50:
      ya = lissage[i]
      xa = depla[i]
      incr = i
      break

  KBond = (yb - ya) / (xb - xa)

  yb = 0
  xb = 0

  for i in range(incr, len(lissage)):
    if PPmax[i] < 50:
      ya = lissage[i]
      xa = depla[i]
      incr = i
      break

  for i in range(incr, len(lissage)):
    if PPmax[i] < 30:
      yb = lissage[i]
      xb = depla[i]
      break

  KDbond = (yb - ya) / (xb - xa)
  KDbond_p = ya - xa * KDbond
  KDbond_y0 = - KDbond_p / KDbond


  KRes = 0
  yb = 0
  ya = 0
  xa = 0
  xb = 0
  incr = 0

  for i in range(len(lissage)):
    if depla[i] > 5:
      ya = lissage[i]
      xa = depla[i]
      incr = i
      break

  for i in range(incr, len(lissage)):
    if depla[i] > 8:
      yb = lissage[i]
      xb = depla[i]
      incr = i
      break

  KRes = (yb - ya) / (xb - xa)
  KRes_p = ya - xa * KRes
  PRes = KDbond * (KDbond_p - KRes_p) / (KRes - KDbond) + KDbond_p

  d1=0
  d2=0
  d3=0
  incr=0

  for i in range(incr, len(lissage)):
    if lissage[i] > PMax/2:
      d1 = depla[i]
      incr = i
      break

  for i in range(incr, len(lissage)):
    if lissage[i] == PMax:
      d2 = depla[i]
      incr = i
      borne2 = i
      break

  for i in range(incr, len(lissage)):
    if lissage[i] < PMax/2:
      d3 = depla[i]
      break

  s1=2*(d2 - d1)
  s2=2*(d3 - d2)

  borne1=0
  borne3=0

  for i in range(len(lissage)):
    if depla[i] > d2 - s1:
      borne1 = i
      break

  for i in range(borne2, len(lissage)):
    if depla[i] > d2 + s2:
      borne3 = i
      break

  aire1=0
  aire2=0
  for i in range(borne1, borne2):
    aire1 += aire[i]
  for i in range(borne2, borne3):
    aire2 += aire[i]


  depla0=0
  for i in range(50, len(lissage)):
    if lissage[i] < 0.05 * lissage_mean:
      depla0 = depla[i]
      break

  fout = open(os.path.splitext(file)[0] + ".csv", 'w')

  fout.write('PMax;' + str(PMax) + "\n")
  fout.write('Kbond;' + str(KBond) + "\n")
  fout.write('Kdebond;' + str(KDbond) + "\n")
  fout.write('Kdebond_y0;' + str(KDbond_y0) + "\n")
  fout.write('depla0;' + str(depla0) + "\n")
  fout.write('LG;' + str(depla0 - KDbond_y0) + "\n")
  fout.write('KRes;' + str(KRes) + "\n")
  fout.write('PRes;' + str(PRes) + "\n")
  fout.write('aire1;' + str(aire1) + "\n")
  fout.write('aire2;' + str(aire2) + "\n")
  computed_params = dict()
  computed_params["PMax"] = PMax
  computed_params["Kbond"] = KBond
  computed_params["Kdebond"] = KDbond
  computed_params["KRes"] = KRes
  computed_params["PRes"] = PRes
  computed_params["aire1"] = aire1
  computed_params["aire2"] = aire2

  units = dict()
  units["PMax"] = "N"
  units["Kbond"] = "N/mm"
  units["Kdebond"] = "N/mm"
  units["KRes"] = "N/mm"
  units["PRes"] = "N"
  units["aire1"] = "mm*mm"
  units["aire2"] = "mm*mm"


  for i in range(len(aire)):
    fout.write('{};{};{};{};{};{}\n'.format(depla[i], charge[i], chargeN[i], lissage[i], PPmax[i], aire[i]))
  fout.close()
  data_curves[file] = (depla, lissage, computed_params, units)

main()
