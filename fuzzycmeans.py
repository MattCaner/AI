# -*- coding: utf-8 -*-

"""
IMPLEMENTATION OF FUZZY C-MEANS CLUSTERING ALGORITHM

@author: MattCaner
"""

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import skfuzzy as fuzz
import csv
from skfuzzy import control as ctrl
import pandas
import statistics as stats
import random
from datetime import datetime


random.seed(datetime.now())

fil = pandas.read_csv('data.csv',delimiter = ",",decimal = ".",dtype=np.float64)



class DataBlock:
    def __init__(self,in1,in2,in3,out):
        self.in1 = in1
        self.in2 = in2
        self.in3 = in3
        self.out = out

class Cluster:
    def __init__(self,x,y,z,o):
        self.x = x
        self.y = y
        self.z = z
        self.o = o
        
def euclD(d,c):
    return np.sqrt( (d.in1-c.x)**2 + (d.in2-c.y)**2 + (d.in3-c.z)**2 + (d.out-c.o)**2 )

tabX = []
tabY = []
tabZ = []
tabOut = []


for x in fil['in1']:
    tabX.append(float(x))

for x in fil['in2']:
    tabY.append(float(x))
    
for x in fil['in3']:
    tabZ.append(float(x))

for x in fil['out']:
    tabOut.append(float(x))

DataTab = []

for i in range(0,255):
    DataTab.append(DataBlock(tabX[i],tabY[i],tabZ[i],tabOut[i]))
    
cluster_n = 27
error = 0.005
m = 2
dataSize = 255

ClusterTab = []

for i in range(0,cluster_n):
    ClusterTab.append(Cluster(random.uniform(min(tabX),max(tabX)),random.uniform(min(tabY),max(tabY)),random.uniform(min(tabZ),max(tabZ)),random.uniform(min(tabOut),max(tabOut))))
    print("app")

UMatrix = np.zeros((cluster_n,dataSize)) 
UOld = np.zeros((cluster_n,dataSize))


iternum = 0
iternumMax = 100


while iternum < iternumMax:
    
    for i in range(0,cluster_n):
        for j in range(0,255):
            UOld[i][j] = UMatrix[i][j]
    
    
   #OBLICZANIE NOWEJ MACIERZY U
    for i in range(0,dataSize):
        suma = 0
        for j in range(0,cluster_n):
            suma = 0
            for k in range(0,cluster_n):
                suma+=(euclD(DataTab[i],ClusterTab[j])/euclD(DataTab[i],ClusterTab[k]))**(2/(m-1))
            UMatrix[j][i] = 1/suma

    #KOREKCJA SRODKOW KLASTROW
    for i in range(0,cluster_n):
        suma1 = 0
        suma2 = 0
        for j in range(0,dataSize):
            suma1 +=(UMatrix[i][j]**m) * DataTab[j].in1
            suma2 +=(UMatrix[i][j]**m)
        ClusterTab[i].x = suma1/suma2
        suma1 = 0
        suma2 = 0
        for j in range(0,255):
            suma1 +=(UMatrix[i][j]**m) * DataTab[j].in2
            suma2 +=(UMatrix[i][j]**m)
        ClusterTab[i].y = suma1/suma2
        suma1 = 0
        suma2 = 0
        for j in range(0,255):
            suma1 +=(UMatrix[i][j]**m) * DataTab[j].in3
            suma2 +=(UMatrix[i][j]**m)
        ClusterTab[i].z = suma1/suma2
        suma1 = 0
        suma2 = 0
        for j in range(0,255):
            suma1 +=(UMatrix[i][j]**m) * DataTab[j].out
            suma2 +=(UMatrix[i][j]**m)
        ClusterTab[i].o = suma1/suma2
    

    print("ITERNUM:             ",iternum)
    for i in ClusterTab:
        print("Cluster loc:",round(i.x,4),round(i.y,4),round(i.z,4),round(i.o,4))

    #SZUKAJ MAKSYMALNEJ ROZNICY
    maxim = 0
    for i in range(0,3):
        for j in range(0,255):
            if abs(UOld[i][j]-UMatrix[i][j]) > maxim:
                maxim = abs(UOld[i][j]-UMatrix[i][j])
                #print("MAXIM in ",i,"  ",j," = ",maxim)
    
    #SPRAWDZ, CZY ROZNICA JEST MNIEJSZA OD USTALONEGO BLEDU
    if maxim < error:
        break    

    iternum +=1
    
#SPRAWDZ DO KTOREGO KLASTRA NALEZY PUNKT

MembershipTable = np.zeros(dataSize)

for i in range(0,dataSize):
    t_index = 0
    for j in range(0,cluster_n):
        if UMatrix[j][i] >= UMatrix[t_index][i]:
            t_index = j
    MembershipTable[i] = t_index
    
valx = ctrl.Antecedent(np.arange(min(tabX),max(tabX),0.0001), 'valx')
valy = ctrl.Antecedent(np.arange(min(tabY),max(tabY),0.0001), 'valy')
valz = ctrl.Antecedent(np.arange(min(tabZ),max(tabZ),0.0001), 'valz')
valo = ctrl.Consequent(np.arange(min(tabOut),max(tabOut),0.0001), 'valo')

valx['lo'] = fuzz.zmf(valx.universe,min(tabX),stats.mean(tabX))
valx['med'] = fuzz.gaussmf(valx.universe,stats.mean(tabX),0.05)
valx['hi'] = fuzz.smf(valx.universe,stats.mean(tabX),max(tabX))

valxF = [valx['lo'],valx['med'],valx['hi']]

valy['lo'] = fuzz.zmf(valy.universe,min(tabY),stats.mean(tabY))
valy['med'] = fuzz.gaussmf(valy.universe,stats.mean(tabY),0.05)
valy['hi'] = fuzz.smf(valy.universe,stats.mean(tabY),max(tabY))

valyF = [valy['lo'],valy['med'],valy['hi']]

valz['lo'] = fuzz.zmf(valz.universe,min(tabZ),stats.mean(tabZ))
valz['med'] = fuzz.gaussmf(valz.universe,stats.mean(tabZ),0.05)
valz['hi'] = fuzz.smf(valz.universe,stats.mean(tabZ),max(tabZ))

valzF = [valz['lo'],valz['med'],valz['hi']]

valo['lo'] = fuzz.zmf(valo.universe,min(tabOut),stats.mean(tabOut))
valo['med'] = fuzz.gaussmf(valo.universe,stats.mean(tabOut),0.05)
valo['hi'] = fuzz.smf(valo.universe,stats.mean(tabOut),max(tabOut))

valoF = [valo['lo'],valo['med'],valo['hi']]

CountTable = np.zeros((3,3,3))

CountRuleTable = np.zeros((3,3))


MTab = []

for i in range(0,cluster_n):
    maxO = [0,0,0]
    maxX = [0,0,0]
    maxY = [0,0,0]
    maxZ = [0,0,0]
    for j in range(0,dataSize):
        if MembershipTable[j] == i:
            tmpTabX = []
            tmpTabY = []
            tmpTabZ = []
            tmpTabO = []
            for k in range(0,3):
                tmpTabX.append(fuzz.interp_membership(valx.universe, valxF[k].mf, DataTab[j].in1))
                tmpTabY.append(fuzz.interp_membership(valy.universe, valyF[k].mf, DataTab[j].in2))
                tmpTabZ.append(fuzz.interp_membership(valz.universe, valzF[k].mf, DataTab[j].in3))
                tmpTabO.append(fuzz.interp_membership(valo.universe, valoF[k].mf, DataTab[j].out))
            #print(tmpTabO)
            maxo = 0
            maxx = 0
            maxy = 0
            maxz = 0
            for ii in range(0,3):
                if tmpTabX[ii] > tmpTabX[maxx]:
                    maxx = ii
                if tmpTabY[ii] > tmpTabY[maxy]:
                    maxy = ii
                if tmpTabZ[ii] > tmpTabZ[maxz]:
                    maxz = ii
                if tmpTabO[ii] > tmpTabO[maxo]:
                    maxo = ii
            maxO[maxo] += 1
            maxX[maxx] += 1
            maxY[maxy] += 1
            maxZ[maxz] += 1
    #print(maxO)
    #print(maxX)
    #print(maxY)
    #print(maxZ)
    tabsofc = [maxO,maxX,maxY,maxZ]
    tabofnames = ["outs","xs","ys","zs"]
    tabofgrnames = ["lo","med","hi"]
    print("FOR cluster number ",i,":")
    for ij in range(0,4):
        tmax = max(tabsofc[ij])
        if tmax == 0:
            print("No ",tabofnames[ij])
        else:
            if tmax == tabsofc[ij][0]:
                print(tabofnames[ij]," are ",tabofgrnames[0])
            if tmax == tabsofc[ij][1]:
                print(tabofnames[ij]," are ",tabofgrnames[1])
            if tmax == tabsofc[ij][2]:
                print(tabofnames[ij]," are ",tabofgrnames[2])
    print("     ")
    

rule1 = ctrl.Rule(valx['lo'] & valy['lo'] & valz['lo'],valo['lo'])
rule2 = ctrl.Rule(valx['med'] & valy['med'] & valz['med'],valo['med'])
rule3 = ctrl.Rule(valx['lo'] & valy['med'] & valz['lo'],valo['med'])
rule4 = ctrl.Rule(valx['med'] & valy['med'] & valz['hi'],valo['med'])
rule5 = ctrl.Rule(valx['hi'] & valy['med'] & valz['hi'],valo['hi'])
rule6 = ctrl.Rule(valx['hi'] & valy['hi'] & valz['hi'],valo['hi'])
rule7 = ctrl.Rule(valx['lo'] & valy['med'] & valz['med'],valo['lo'])
rule8 = ctrl.Rule(valx['hi'] & valy['hi'] & valz['med'],valo['hi'])

ctsys = ctrl.ControlSystem([rule1,rule2,rule3,rule4,rule5,rule6,rule7,rule8])

ct = ctrl.ControlSystemSimulation(ctsys)

mse = 0
for x in DataTab:
    ct.input['valx'] = x.in1
    ct.input['valy'] = x.in2
    ct.input['valz'] = x.in3
    ct.compute()
    
    #print(ct.output['valo'],"   ",x.out)
    mse+=(ct.output['valo']-x.out)**2

mse = mse/255
print(mse)
