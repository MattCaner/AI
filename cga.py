"""
Created by M Trzcinski 13:01 26.10.2018

CGA Implementation

Praise Omnissiah.
"""

import numpy as np
import random as rnd
from datetime import datetime
import math

rnd.seed(datetime.now())



class Chrm:
    def __init__(self,size):
        self.size = size
        self.geneTab = np.zeros(size)
    def randomize(self):
        for x in range(self.size):
            self.geneTab[x] = rnd.randint(0,1)
    def validateT(self,vFunc):          #validate (as Table)
        return vFunc(*self.geneTab)
    def toNum(self):
        n = 0
        for i in range(self.size):
            n+=self.geneTab[self.size-i-1]*2**i
        n = n/(4194304/3)
        n = n-1
        return n
    def validateN(self,vFunc):          #validate (as a single Number)
        return vFunc(*(self.toNum()))
    def mutate(self,locus):
        self.geneTab[locus] = not self.geneTab[locus]


class Popul:
    def __init__(self,popsize,chrmsize):
        self.popsize = popsize
        self.popTable = []
        for i in range(popsize):
            self.popTable.append(Chrm(chrmsize))
        for x in self.popTable:
            x.randomize()
    def upgradeToSurvivors(self,SurvTab):
        NewTab = []
        for i in range(self.popsize):
            NewTab.append(self.popTable[SurvTab[i]])
    def randomlyMutate(self,mutc):
        howmany = (self.popsize*self.popTable[0].size)*mutc
        for i in range(math.ceil(howmany)):
            tmp = rnd.randint(0,self.popsize-1)
            self.popTable[tmp].mutate(rnd.randint(0,self.popTable[0].size-1))
            

def RuletteSel(pop,vFunc):
    ChanceTab = []
    TargTab = np.zeros(pop.popsize)
    TabOfSurvivors = []
    suma = 0
    for x in range(pop.popsize):
        evaluated = vFunc((pop.popTable[x].toNum()))
        ChanceTab.append(evaluated)
        suma += evaluated
        if x == 0:
            TargTab[x] = evaluated
        else:
            TargTab[x] = TargTab[x-1] + evaluated
    for x in range(pop.popsize):
        randtarg = rnd.randint(0,int(suma))
        for y in range(pop.popsize):
            if randtarg <= TargTab[y]:
                TabOfSurvivors.append(y)
                break
    #print(TabOfSurvivors)
    return TabOfSurvivors
    
def CrsBreedOneP(par1,par2,locus,siz):
    tmp = Chrm(par1.size)
    for i in range(locus):
        tmp.geneTab[i] = par1.geneTab[i]
    for i in range(locus,par1.size-locus):
        tmp.geneTab[i] = par2.geneTab[i]
    return tmp;

        
def PerformCrossBreedingOneP(pop,bredcoeff):
    #tmppop = copy.deepcopy(pop)
    tmppop = pop
    ntobreed = int(bredcoeff*pop.popsize)
    BreedTableNum = []
    for i in range(ntobreed):
        BreedTableNum.append(rnd.randint(0,pop.popsize-1))
    for j in range(int(ntobreed/2)):
        Par1 = rnd.choice(BreedTableNum)
        P1 = rnd.choice(pop.popTable)
        P2 = rnd.choice(pop.popTable)
        Par2 = rnd.choice(BreedTableNum)
        locusl = rnd.randint(0,tmppop.popTable[1].size)
        Chld1 = CrsBreedOneP(P1,P2,locusl,tmppop.popTable[1].size)
        Chld2 = CrsBreedOneP(P2,P1,locusl,tmppop.popTable[1].size)
        tmppop.popTable[Par1] = Chld1
        tmppop.popTable[Par2] = Chld2
    return tmppop

def Validator(iks):
    return iks*math.sin(10*math.pi*iks) + 1


iterations = 100

mutationcoeff = 0.1
breedingcoeff = 0.5
"""
for x in Pop1.popTable:
    print(x.toNum())
"""

GeneralPop = Popul(100,22)

rememberbest = 0
rememberbesti = 0

print ("HERE IT ALL STARTS")

for iterator in range(iterations):
    
    SurvivorTable = RuletteSel(GeneralPop,Validator)
    GeneralPop.upgradeToSurvivors(SurvivorTable)
    GeneralPop = PerformCrossBreedingOneP(GeneralPop,breedingcoeff)
    GeneralPop.randomlyMutate(mutationcoeff)
    
    for i in GeneralPop.popTable:
        evaluated = Validator(i.toNum())
        if evaluated > rememberbest:
            print (i.toNum())
            print (i.geneTab)
            rememberbest = evaluated
            rememberbesti = i.toNum()
    
    print("Iteration no: ",iterator," Opt: x=",rememberbesti," y=",rememberbest)
        
    
    
