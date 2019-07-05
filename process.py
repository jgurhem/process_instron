#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import os, fnmatch
import sys



def main():
	#processRawFile("ec0.5ac0.1 1D.is_tens.raw", 6)
	#processRawFile("ec0.5ac0.1 1G.is_tens.raw", 6)
	nbCoeff = 6
	#directory = '.' #directory from which the script is executed
	directory = 'raws'
	for file in fnmatch.filter(os.listdir(directory), '*.raw'):
		processRawFile(directory + "/" + file, nbCoeff)

def processRawFile(file, nbCoeff):
	data = np.genfromtxt( file , usecols =(1, 2), delimiter=";", skip_header=65, encoding=None, dtype=None)
	float(data[1][1].replace(',','.'))

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

	lissageMax = 0;
	for i in range(nbCoeff, len(depla) - nbCoeff):
		tmp = 0
		for j in range(1, nbCoeff+1):
			#print(i, j, nbCoeff - j + 1, chargeN[i - j], chargeN[i + j])
			tmp += (chargeN[i - j] + chargeN[i + j])*(nbCoeff - j + 1)
		tmp = tmp/div
		if tmp > lissageMax:
			lissageMax = tmp
		lissage.append(tmp)


	for i in range(nbCoeff, len(lissage)):
		tmp = lissage[i]*100/lissageMax
		PPmax.append(tmp)

	for i in range(nbCoeff, len(lissage) - 1):
		tmp = (lissage[i] + lissage[i+1]) / 2 * (depla[i+1] - depla[i])
		aire.append(tmp)


	pente = 0
	yb = 0
	ya = 0
	xa = 0
	xb = 0
	for i in range(len(lissage)):
		if PPmax[i] > 30:
			yb = lissage[i]
			xb = depla[i]
			break

	for i in range(len(lissage)):
		if PPmax[i] > 50:
			ya = lissage[i]
			xa = depla[i]
			break

	pente = (yb - ya) / (xb - xa)

	d1=0
	d2=0
	d3=0
	incr=0

	for i in range(incr, len(lissage)):
		if lissage[i] > lissageMax/2:
			d1 = depla[i]
			incr = i
			break

	for i in range(incr, len(lissage)):
		if lissage[i] == lissageMax:
			d2 = depla[i]
			incr = i
			borne2 = i
			break

	for i in range(incr, len(lissage)):
		if lissage[i] < lissageMax/2:
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

	fout = open(os.path.splitext(file)[0] + ".csv", 'w')
	
	fout.write('PMax;' + str(lissageMax) + "\n")
	fout.write('Kbond;' + str(pente) + "\n")
	fout.write('aire1;' + str(aire1) + "\n")
	fout.write('aire2;' + str(aire2) + "\n")


	for i in range(len(aire)):
		fout.write('{};{};{};{};{};{}\n'.format(depla[i], charge[i], chargeN[i], lissage[i], PPmax[i], aire[i]))
	fout.close()

main()
