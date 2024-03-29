#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os, fnmatch
import sys
import csv
from core.curves import make_curve
from core.compute_params import find_point

def main():
  if len(sys.argv) != 3:
    print("Usage : ")
    print("py", sys.argv[0], "directory parameters_file")
    sys.exit(0)
  nbCoeff = 6
  directory = sys.argv[1] #directory from which the script is executed
  data_curves = dict()
  for file in fnmatch.filter(os.listdir(directory), '*.raw'):
    with open(sys.argv[2]) as csv_param_file:
      csv_reader = csv.reader(csv_param_file, delimiter=';')
      param_file_data = None
      for row in csv_reader:
        if row[0] == file:
          param_file_data = row
      if param_file_data != None:
        processRawFile(directory + "/" + file, nbCoeff, data_curves, param_file_data)
      else:
        print("File", file, "cannot be processed since there is no corresponding entry in", sys.argv[2])
  make_curve(data_curves, directory)

def processRawFile(file, nbCoeff, data_curves, param_file_data):
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
  PMax_index = 0
  for i in range(nbCoeff, len(depla) - nbCoeff):
    tmp = 0
    for j in range(1, nbCoeff+1):
      #print(i, j, nbCoeff - j + 1, chargeN[i - j], chargeN[i + j])
      tmp += (chargeN[i - j] + chargeN[i + j])*(nbCoeff - j + 1)
    tmp = tmp/div
    if tmp > PMax:
      PMax = tmp
      PMax_index = i
    lissage.append(tmp)


  for i in range(nbCoeff, len(lissage)):
    tmp = lissage[i]*100/PMax
    PPmax.append(tmp)

  for i in range(nbCoeff, len(lissage) - 1):
    tmp = (lissage[i] + lissage[i+1]) / 2 * (depla[i+1] - depla[i])
    aire.append(tmp)


  ra = find_point(0, len(lissage), lissage, depla, (lambda i: PPmax[i] > 50))
  rb = find_point(0, len(lissage), lissage, depla, (lambda i: PPmax[i] > 30))
  KBond = (rb["y"] - ra["y"]) / (rb["x"] - ra["x"])

  ra = find_point(ra["index"], len(lissage), lissage, depla, (lambda i: PPmax[i] < 50))
  rb = find_point(ra["index"], len(lissage), lissage, depla, (lambda i: PPmax[i] < 30))
  KDbond = (rb["y"] - ra["y"]) / (rb["x"] - ra["x"])
  KDbond_p = ra["y"] - ra["x"] * KDbond
  KDbond_y0 = - KDbond_p / KDbond

  ra = find_point(0, len(lissage), lissage, depla, (lambda i: depla[i] > 5))
  rb = find_point(ra["index"], len(lissage), lissage, depla, (lambda i: depla[i] > 8))
  KRes = (rb["y"] - ra["y"]) / (rb["x"] - ra["x"])
  KRes_p = ra["y"] - ra["x"] * KRes
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

  depla0 = float(param_file_data[1])

  LG = depla0 - KDbond_y0
  ra = find_point(0, len(lissage), lissage, depla, (lambda i: depla[i] > KDbond_y0 + 0.3 * LG))
  rb = find_point(ra["index"], len(lissage), lissage, depla, (lambda i: depla[i] > KDbond_y0 + 0.5 * LG))
  nKRes = (rb["y"] - ra["y"]) / (rb["x"] - ra["x"])
  nKRes_p = ra["y"] - ra["x"] * nKRes
  nPRes = KDbond * (KDbond_p - nKRes_p) / (nKRes - KDbond) + KDbond_p

  fout = open(os.path.splitext(file)[0] + ".csv", 'w')

  ext_mass = float(param_file_data[2])
  ext_mass_on_depla0 = ext_mass / depla0

  fout.write('PMax;' + str(PMax) + "\n")
  fout.write('Kbond;' + str(KBond) + "\n")
  fout.write('Kdebond;' + str(KDbond) + "\n")
  fout.write('Kdebond_y0;' + str(KDbond_y0) + "\n")
  fout.write('depla0;' + str(depla0) + "\n")
  fout.write('ext_mass;' + str(ext_mass) + "\n")
  fout.write('ext_mass_on_depla0;' + str(ext_mass_on_depla0) + "\n")
  fout.write('LG;' + str(LG) + "\n")
  fout.write('nKRes;' + str(nKRes) + "\n")
  fout.write('nPRes;' + str(nPRes) + "\n")
  fout.write('aire1;' + str(aire1) + "\n")
  fout.write('aire2;' + str(aire2) + "\n")
  fout.write('s1;' + str(s1) + "\n")
  fout.write('s2;' + str(s2) + "\n")
  fout.write('d1;' + str(d1) + "\n")
  fout.write('d2;' + str(d2) + "\n")
  fout.write('d3;' + str(d3) + "\n")
  computed_params = dict()
  computed_params["PMax"] = PMax
  computed_params["Kbond"] = KBond
  computed_params["Kdebond"] = KDbond
  computed_params["nKRes"] = nKRes
  computed_params["nPRes"] = nPRes
  computed_params["aire1"] = aire1
  computed_params["aire2"] = aire2
  computed_params["depla0"] = depla0
  computed_params["ext_mass"] = ext_mass
  computed_params["s1"] = s1
  computed_params["s2"] = s2
  computed_params["d1"] = d1
  computed_params["d2"] = d2
  computed_params["d3"] = d3
  computed_params["ext_mass_on_depla0"] = ext_mass_on_depla0


  units = dict()
  units["PMax"] = "N"
  units["Kbond"] = "N/mm"
  units["Kdebond"] = "N/mm"
  units["nKRes"] = "N/mm"
  units["nPRes"] = "N"
  units["aire1"] = "mm*mm"
  units["aire2"] = "mm*mm"
  units["depla0"] = "mm"
  units["ext_mass"] = "g"
  units["s1"] = "mm"
  units["s2"] = "mm"
  units["d1"] = "mm"
  units["d2"] = "mm"
  units["d3"] = "mm"
  units["ext_mass_on_depla0"] = "g/mm"


  for i in range(len(aire)):
    fout.write('{};{};{};{};{};{}\n'.format(depla[i], charge[i], chargeN[i], lissage[i], PPmax[i], aire[i]))
  fout.close()
  data_curves[file] = (depla, lissage, computed_params, units)

main()
