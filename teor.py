#!/usr/bin/python3
import pandas as pd
import numpy as np
from scipy import interpolate
from array import *
from scipy.interpolate import interp1d
import matplotlib.pyplot as pl
import math
import os, sys
import zipfile
import glob
import shutil
         
fantasy_zip = zipfile.ZipFile('/home/dasha/kmfr-0.95/output/Er.zip')
fantasy_zip.extractall('/home/dasha/kmfr-0.95/output/data/')
fantasy_zip.close()

fantasy_zip2 = zipfile.ZipFile('/home/dasha/seltzer.zip')
fantasy_zip2.extractall('/home/dasha/data2/')
fantasy_zip2.close()

den = 0.958678 #target density
Amass = (2*166.17+3*16)
M = 2*den*6.02*(10**(-4))/Amass

#need_points = np.arange(0.125,55.625+0.25,0.25)
'''A = [196,198,199,200,201,202,204] #HG
ab = [0.0015,0.0997,0.1687,0.231,0.1318,0.2986,0.0687] #abundunce(in parts) of HG in the same order!'''
A = [162,164,166,167,168,170] #Er
ab = [0.0014,0.0161,0.336,0.2295,0.268,0.149] #Er
neutron = [1,2,3,4,5,6,7,0,1,2,3,4]
proton = [0,0,0,0,0,0,0,1,1,1,1,1]

number = []
for i in range(1,100):
	if (i<10):
		i = '0'+str(i)
	number.append(i)
	
for l in range(len(neutron)-1):
	for j in A:
	 
		key = r"/home/dasha/kmfr-0.95/output/data/{id}Er/key".format(id = j)
		df = pd.read_csv(key, delim_whitespace=True, skiprows=0,names=['col1','col2', 'col3', 'col4','col5'])
		k1,k2,k3,k4,k5=df['col1'].values, df['col2'].values, df['col3'].values, df['col4'].values, df['col5'].values
		
		for i in range(0,len(k1)):
			if (float(k4[i]) == proton[l] and float(k5[i]) == neutron[l]): 
				
				crs = r"/home/dasha/kmfr-0.95/output/data/{id}Er/sct{ik}.dat".format(id = j,ik = i+1)
				df = pd.read_csv(crs, delim_whitespace=True, skiprows=0,names=['col1','col2', 'col3', 'col4','col5', 'col6','col7'])
				crs1,crs2,crs3,crs4,crs5,crs6,crs7=df['col1'].values, df['col2'].values, df['col3'].values, df['col4'].values, df['col5'].values, df['col6'].values, df['col7'].values
				
				for x in number: #range of brem. (end point)
				
					bremm = r"/home/dasha/data2/seltzer-{id}.txt".format(id = x)
					df = pd.read_csv(bremm, delim_whitespace=True, skiprows=0,names=['col1','col2','col3','col4','col5','col6','col7','col8'])
					b1,b2,b3,b4,b5,b6,b7,b8=df['col1'].values, df['col2'].values,df['col3'].values, df['col4'].values,df['col5'].values, df['col6'].values,df['col7'].values, df['col8'].values
					
					br = [0 for i in range(3,len(b1)-1)]
					for i in range(3,len(b1)-1):
						br[i-3] = float(b1[i])*float(b2[i])/b4[len(b4)-1]
					norm = sum(br)
				
					points = np.arange(69.375, crs1[0], -0.25) #69.375 is the last point in sct, closest to needed points (last one in sct is 69.4)
					reversed_points = points[::-1]
					
					f = interp1d(crs1, crs2, kind='cubic')
					cs = f(reversed_points)
					for k in range(0,len(cs)):
						if (cs[k]<0):
							cs[k]=0
							
					y = [0 for x in range(0,len(cs))]
					for z in range(3,len(b1)-1): 
						if (float(b1[z]) == reversed_points[0]-0.125):
							Z = z
					for m in range(Z,len(cs)):
						y[m] = cs[m-Z]*float(b2[m])/float(b4[len(b4)-1])
					Ykmfr = sum(y)
					Y = Ykmfr*M*ab[A.index(j)]
					Q = Ykmfr*ab[A.index(j)]*int(x)/norm
					#print(neutron[l],'n',proton[l],'p','from',j,'\t',Y)
					print(neutron[l],'n',proton[l],'p','from',j,'\t', int(x),'\t',Q)
						
shutil.rmtree('/home/dasha/kmfr-0.95/output/data') 
shutil.rmtree('/home/dasha/data2')


