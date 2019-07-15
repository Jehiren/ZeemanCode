# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 04:09:46 2019

@author: Leo
"""

import numpy as np
from matplotlib import pyplot as p
from matplotlib.patches import ConnectionPatch
import os
os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/Main Program")


BmagG = np.linspace(0, .5, 1000) #0 to 3T magnetic field
BmagE = np.linspace(0, .02, 400) #0 to 3T magnetic field

BmagCycling = np.linspace(0, .05, 500) #0 to 3T magnetic field

BmagG = BmagCycling
BmagE = BmagCycling


#Ground State constants for Rb
g_s = 2.00023
g_L = 1
g_I = -0.00029364
mu_b = 9.274e-24
ggJ = 2 #g-lande factor is 2 for the ground state

E_hfs = 3.036e9  #Hyperfine splitting for zero field in MHz

p.close('all')
hbar = 1.0546e-34
hh = 2*np.pi*hbar
mu_b = 9.274e-24


cg = mu_b/(hh) #~14GHz/T
cg2= mu_b/hbar

Ahfs1 = 1.0119e9
gI = 5/2
gm_I = np.array([5/2, 3/2, 1/2, -1/2, -3/2, -5/2])
gS = 1/2
gm_S = [1/2, -1/2]
gL = 0
gm_J = gm_S
#mf states are always conserved.
gmf = []
#define our states using |mI, mJ> basis in an array, this is what the 12x12 matrix is based off of for the linked terms.
gmimj = []
for gss in gm_S:
    for GII in gm_I:
        gmf.append(GII + gss)
        gmimj.append([GII, gss])


PsiG = (np.zeros((12,12)))
Gstates =np.mat(np.identity(12))
#Define the CG coefficients for the ground state of RB85:

#F=3 states

PsiG[0] = Gstates[0]
PsiG[1] = np.sqrt(5/6)*Gstates[2] + np.sqrt(1/6)*Gstates[1]
PsiG[2] = np.sqrt(2/3)*Gstates[4] + np.sqrt(1/3)*Gstates[3]
PsiG[3] = np.sqrt(0.5)*Gstates[6] + np.sqrt(0.5)*Gstates[5]
PsiG[4] = np.sqrt(1/3)*Gstates[8] + np.sqrt(2/3)*Gstates[7]
PsiG[5] = np.sqrt(1/6)*Gstates[10] + np.sqrt(5/6)*Gstates[9]
PsiG[6] = Gstates[11]

#F = 2 states

PsiG[7] = np.sqrt(1/6)*Gstates[2] - np.sqrt(5/6)*Gstates[1]
PsiG[8] = np.sqrt(1/3)*Gstates[4] - np.sqrt(2/3)*Gstates[3]
PsiG[9] = np.sqrt(0.5)*Gstates[6] - np.sqrt(0.5)*Gstates[5]
PsiG[10] = np.sqrt(2/3)*Gstates[8] - np.sqrt(1/3)*Gstates[7]
PsiG[11] = np.sqrt(5/6)*Gstates[10] - np.sqrt(1/6)*Gstates[9]


PsiBG = np.mat(np.zeros((12,12)))


mjMatG = (np.diagflat(6*[gm_J]))
for i in range(12):
    tPsiL = PsiG[i]
    for j in range(12):
        PsiBG[i,j] = PsiBG[i,j] + np.linalg.multi_dot([tPsiL, mjMatG, PsiG[j]])

PsiCG = np.mat(np.zeros((12,12)))


miMatG = np.diagflat((np.kron(gm_I.T, 2*[1])))

for i in range(12):
    tPsiL = PsiG[i]
    for j in range(12):
        PsiCG[i,j] = PsiCG[i,j] + np.linalg.multi_dot([tPsiL, miMatG, PsiG[j]])
    
HtotZG = np.mat(np.zeros((12,12)))
HtotZG = ggJ*PsiBG + g_I*PsiCG

#The Zeeman terms are dependent upon PsiB and PsiC and are actually correct now! Just need to get the rest of the terms now.

JpTG = np.sqrt((0.5 - np.array(gm_J))*(0.5+np.array(gm_J)+1))
JmTG = np.sqrt((0.5 + np.array(gm_J))*(0.5-np.array(gm_J)+1))
IpTG = np.sqrt((2.5 - np.array(gm_I))*(2.5+np.array(gm_I)+1))
ImTG = np.sqrt((2.5 + np.array(gm_I))*(2.5-np.array(gm_I)+1))

JpmatG = np.diagflat(6*[JpTG],1)[1:,1:]
JmmatG = np.diagflat(6*[JmTG], -1)[:-1,:-1]

IpmatG = np.diagflat((np.kron(IpTG.T, 2*[1])),2)[2:,2:]
ImmatG = np.diagflat((np.kron(ImTG.T, 2*[1])), -2)[:-2, :-2]


HPsiPG = np.mat(np.zeros((12,12)))

for i in range(12):
    temppsi = PsiG[i]
    tempIpJm = 0.5*np.linalg.multi_dot([temppsi, IpmatG, JmmatG])
    tIpJm = tempIpJm
    tempImJp = 0.5*np.linalg.multi_dot([temppsi, ImmatG, JpmatG])
    tImJp = tempImJp
#    print(tIpJm)
#    print('')
    TIzJz = np.linalg.multi_dot([temppsi, np.dot(miMatG, mjMatG)])
    tIdotJ = np.mat(tIpJm + tImJp + TIzJz)
    for j in range(12):
        HPsiPG[i,j] = HPsiPG[i,j] + np.dot(tIdotJ, PsiG[j])

        
cd = 1.4e6

HCyclingG = np.zeros_like(PsiG)
HtotCZG = np.zeros_like(PsiG)
HtotCG = np.zeros_like(PsiG)

HtotPG = np.zeros_like(PsiG)



        
        
HCyclingG[0] = Ahfs1*HPsiPG[0] 
HCyclingG[6] = Ahfs1*HPsiPG[6] 

#HCycling[1] = Ahfs2*HPsiP[1] + Hoff2[1]
#HCycling[7] = Ahfs2*HPsiP[7] + Hoff2[7]

HtotCZG[0] = HtotZG[0]
HtotCZG[6] = HtotZG[6]

#HtotCZ[1] = HtotZ[1]
#HtotCZ[7] = HtotZ[7]


evalCG = []
evalPG = []

for bvals2 in BmagG:
    HtotPG =  (mu_b/hh)*bvals2*np.mat(HtotZG) + HPsiPG*Ahfs1    
    HtotCG =  (mu_b/hh)*bvals2*np.mat(HtotCZG) + np.mat(HCyclingG) 
    evalCG.append( (np.linalg.eigvalsh(HtotCG)))  
    evalPG.append( (np.linalg.eigvalsh(HtotPG)))
    
evalTPG = np.transpose(evalPG)
evalTCG = np.transpose(evalCG)

evalCsortG = []
for items in range(len(BmagG)):
    if evalTCG[-2][items] != 0:
        evalCsortG.append((evalTCG[-2][items])/1e6)
    else:
        evalCsortG.append((evalTCG[0][items]/1e6))

#Combine the cycling states into an easier to manage array
CyclingPG = [np.array(evalTCG[-1]/1e6), np.array(evalCsortG)]

p.figure()

for items2 in evalTPG[:5]:
    p.plot(10000*BmagG, items2/1e6, color = 'green', alpha=0.9 , linewidth = .85)
for items2 in evalTPG[5:]:
    p.plot(10000*BmagG, items2/1e6, color = 'blue', alpha=0.9 , linewidth = .85)


#This plots the cycling transitions overlaying the states in the hyperfine structure.
#p.plot(10000*BmagS, CyclingPG[0], color = 'orange', linewidth = 3)
#p.plot(10000*BmagS, CyclingPG[1] , color= 'cyan' , linewidth = 3)

p.show()

#p.title('Zeeman Shift in Rb85 Ground State', fontsize = 18)
p.ylabel('Energy shift (E/h) MHz', fontsize = 16)
p.xlabel('Magnetic Field (G)', fontsize = 16)

import os
os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
#p.savefig('RbHyperfineGround.jpeg', bbox_inches='tight', dpi=800)
p.show() 

# =============================================================================
# Ground state complete! Below is the work for the excited state. It uses the same structure but has 24 possible states.
# =============================================================================
#
#
###Excited state constants and possible states.

Ahfs2 = 25.002e6
#Ahfs2 = 1
#Ahfs2S= 120.53e6
Bhfs2 = 25.88e6
#Bhfs2 = 1
em_I = np.array([5/2, 3/2, 1/2, -1/2, -3/2, -5/2])
#em_I = [-2.5, -1.5, -.5, .5, 1.5, 2.5]
eS = 1/2
em_J = [3/2,0.5,-0.5, -3/2]
#em_J = [-1.5, -.5, .5, 1.5]
eL = 1
gm_L = 0
geJ = 4/3



# =============================================================================
# Using the CG coefficients for the 24 states, we'll define column matricies for the states then I'll cycle through and record the overlaps in a 24x24 matrix
# Then the Zeemanshift will be added, but haven't figured out 100% how yet. (Is it just on the diagonal? Will try that first.
# =============================================================================
# When referencing np.mat, it's Psi[row,column]
Psi = (np.zeros((24,24)))
states =np.mat(np.identity(24))
#Now to define the 24 states using the 35 CG coefficients...this will take a minute.
#These 24 states will be written based on emimj |mI,mJ> such that |1> = |5/2, 3/2>, |2> = |5/2, 1/2>, so we cycle through the mJ's for each mI vaule.
#Lowering as we go.

# F, mF = 4,4
Psi[0] = states[0]

# F, mF = 4, 3
Psi[1] = np.sqrt(5/8)*states[4] + np.sqrt(3/8)*states[1]

#F, mF = 4, 2
Psi[2] = np.sqrt(5/14)*states[8] + np.sqrt(15/28)*states[5] + np.sqrt(3/28)*states[2]

#F, mF = 4, 1
Psi[3] = np.sqrt(5/28)*states[12] + np.sqrt(15/28)*states[9] + np.sqrt(15/56)*states[6] + np.sqrt(1/56)*states[3]

#F, mf = 4, 0
Psi[4] = np.sqrt(1/14)*states[16] + np.sqrt(3/7)*states[13] + np.sqrt(3/7)*states[10] + np.sqrt(1/14)*states[7]

#F, mF = 4, -1
Psi[5] = np.sqrt(5/28)*states[11] + np.sqrt(15/28)*states[14] + np.sqrt(15/56)*states[17] + np.sqrt(1/56)*states[20]

#F, mF = 4, -2
Psi[6] = np.sqrt(5/14)*states[15] + np.sqrt(15/28)*states[18] + np.sqrt(3/28)*states[21]
 
# F, mF = 4, -3
Psi[7] = np.sqrt(3/8)*states[22] + np.sqrt(5/8)*states[19]

# F, mF = 4, -4
Psi[8] = states[23]

##Now for F = 3 states
# F, mF = 3, 3
Psi[9] = np.sqrt(5/8)*states[1] - np.sqrt(3/8)*states[4]

#F, mF = 3, 2
Psi[10] = np.sqrt(5/12)*states[2] + np.sqrt(1/12)*states[5] - np.sqrt(1/2)*states[8]

#F, mF = 3, 1
Psi[11] = np.sqrt(1/8)*states[3] + np.sqrt(49/120)*states[6] - np.sqrt(1/60)*states[9] - np.sqrt(9/20)*states[12]

#F, mf = 3, 0
Psi[12] = np.sqrt(3/10)*states[7] + np.sqrt(1/5)*states[10] - np.sqrt(1/5)*states[13] - np.sqrt(3/10)*states[16]

#F, mF = 3, -1
Psi[13] = -np.sqrt(1/8)*states[20] - np.sqrt(49/120)*states[17] + np.sqrt(1/60)*states[14] + np.sqrt(9/20)*states[11]

#F, mF = 3, -2
Psi[14] = -np.sqrt(5/12)*states[21] - np.sqrt(1/12)*states[18] + np.sqrt(1/2)*states[15]
 
# F, mF = 3, -3
Psi[15] = -np.sqrt(5/8)*states[22] + np.sqrt(3/8)*states[19]

##Now the F=2 states
#F, mF = 2, 2
Psi[16] = np.sqrt(10/21)*states[2] - np.sqrt(8/21)*states[5] + np.sqrt(1/7)*states[8]

#F, mF = 2, 1
Psi[17] = np.sqrt(5/14)*states[3] + np.sqrt(1/42)*states[6] - np.sqrt(25/84)*states[9] + np.sqrt(9/28)*states[12]

#F, mf = 2, 0
Psi[18] = np.sqrt(3/7)*states[7] - np.sqrt(1/14)*states[10] - np.sqrt(1/14)*states[13] + np.sqrt(3/7)*states[16]

#F, mF = 2, -1
Psi[19] = np.sqrt(5/14)*states[20] + np.sqrt(1/42)*states[17] - np.sqrt(25/84)*states[14] + np.sqrt(9/28)*states[11]

 #F, mF = 2, -2
Psi[20] = np.sqrt(10/21)*states[21] - np.sqrt(8/21)*states[18] + np.sqrt(1/7)*states[15]

#Finally wrapping up with the F=1 states
#F, mF = 1, 1
Psi[21] = np.sqrt(1/2)*states[3] - np.sqrt(3/10)*states[6] + np.sqrt(3/20)*states[9] - np.sqrt(1/20)*states[12]

#F, mf = 1, 0
Psi[22] = np.sqrt(1/5)*states[7] - np.sqrt(3/10)*states[10] + np.sqrt(3/10)*states[13] - np.sqrt(1/5)*states[16]

#F, mF = 1, -1
Psi[23] = -np.sqrt(1/2)*states[20] + np.sqrt(3/10)*states[17] - np.sqrt(3/20)*states[14] + np.sqrt(1/20)*states[11]


PsiB = np.mat(np.zeros((24,24)))


mjMat = (np.diagflat(6*[em_J]))
for i in range(24):
    tPsiL = Psi[i]
    for j in range(24):
        PsiB[i,j] = PsiB[i,j] + np.linalg.multi_dot([tPsiL, mjMat, Psi[j]])

PsiC = np.mat(np.zeros((24,24)))


miMat = np.diagflat((np.kron(em_I.T, 4*[1])))

for i in range(24):
    tPsiL = Psi[i]
    for j in range(24):
        PsiC[i,j] = PsiC[i,j] + np.linalg.multi_dot([tPsiL, miMat, Psi[j]])
    
HtotZ=np.mat(np.zeros((24,24)))
HtotZ = geJ*PsiB + g_I*PsiC

#The Zeeman terms are dependent upon PsiB and PsiC and are actually correct now! Just need to get the rest of the terms now.

JpT = np.sqrt((1.5 - np.array(em_J))*(1.5+np.array(em_J)+1))
JmT = np.sqrt((1.5 + np.array(em_J))*(1.5-np.array(em_J)+1))
IpT = np.sqrt((2.5 - np.array(em_I))*(2.5+np.array(em_I)+1))
ImT = np.sqrt((2.5 + np.array(em_I))*(2.5-np.array(em_I)+1))

Jpmat = np.diagflat(6*[JpT],1)[1:,1:]
Jmmat = np.diagflat(6*[JmT], -1)[:-1,:-1]

Ipmat = np.diagflat((np.kron(IpT.T, 4*[1])),4)[4:,4:]
Immat = np.diagflat((np.kron(ImT.T, 4*[1])), -4)[:-4, :-4]


HPsiP = np.mat(np.zeros((24,24)))
for i in range(24):
    temppsi = Psi[i]
    tempIpJm = 0.5*np.linalg.multi_dot([temppsi, Ipmat, Jmmat])
    tIpJm = tempIpJm
    tempImJp = 0.5*np.linalg.multi_dot([temppsi, Immat, Jpmat])
    tImJp = tempImJp
#    print(tIpJm)
#    print('')
    TIzJz = np.linalg.multi_dot([temppsi, np.dot(miMat, mjMat)])
    tIdotJ = np.mat(tIpJm + tImJp + TIzJz)
    for j in range(24):
        HPsiP[i,j] = HPsiP[i,j] + np.dot(tIdotJ, Psi[j])
        

#Defining the offset values as were provided in Dr. Paradis' old matlab code.
Hoff = np.zeros((24,24))
Hoff[0,0] = 18e6
Hoff[1,1] = 18e6
Hoff[2,2] = 18e6
Hoff[3,3] = -2e6
Hoff[4,4] = -2e6
Hoff[5,5]= Hoff[4,4]
Hoff[6,6] = Hoff[4,4]
Hoff[7,7] = Hoff[4,4]
Hoff[8,8] = -14e6
Hoff[9,9] = Hoff[8,8]
Hoff[10,10] = Hoff[8,8]
Hoff[11,11] = Hoff[8,8]
Hoff[12,12] =Hoff[8,8]
Hoff[13,13] = Hoff[8,8]
Hoff[14,14] = Hoff[8,8]
Hoff[15,15] = 6.5e6
Hoff[16,16] = Hoff[15,15]
Hoff[17,17] = Hoff[15,15]
Hoff[18,18] = Hoff[15,15]
Hoff[19,19] = Hoff[15,15]
Hoff[20,20] = Hoff[15,15]
Hoff[21,21] = Hoff[15,15]
Hoff[22,22] = Hoff[15,15]
Hoff[23,23] = Hoff[15,15]

Hoffm = np.mat(Hoff)
        
cd = 1.4e6
HintP = (np.zeros((24,24)))
HintZ = np.zeros((24,24))

HCycling = np.zeros_like(HintP)
HtotCZ = np.zeros_like(HintP)
HtotC = np.zeros_like(HintP)

HtotP = np.zeros_like(HintP)

Hoff2 =  np.zeros_like(HintP)



for ihz in range(24):
    for izh in range(24):
        Hoff2[ihz,ihz] = Hoff[-(ihz+1),-(ihz+1)]
        
        
HCycling[0] = Ahfs2*HPsiP[0] + Hoff2[0]
HCycling[8] = Ahfs2*HPsiP[8] + Hoff2[8]

#HCycling[1] = Ahfs2*HPsiP[1] + Hoff2[1]
#HCycling[7] = Ahfs2*HPsiP[7] + Hoff2[7]

HtotCZ[0] = HtotZ[0]
HtotCZ[8] = HtotZ[8]

#HtotCZ[1] = HtotZ[1]
#HtotCZ[7] = HtotZ[7]


evalC = []
evalP = []

for bvals2 in BmagE:
#    Htot3 = Psi2 + 1*HintP 
    HtotP =  (mu_b/hh)*bvals2*np.mat(HtotZ) + HPsiP*Ahfs2+ Hoff2    
    HtotC =  (mu_b/hh)*bvals2*np.mat(HtotCZ) + np.mat(HCycling) 
    evalC.append( (np.linalg.eigvalsh(HtotC)))  
    evalP.append( (np.linalg.eigvalsh(HtotP)))
    
evalTP= np.transpose(evalP)
evalTC = np.transpose(evalC)

evalCsort = []
for items in range(len(BmagE)):
    if evalTC[-2][items] != 0:
        evalCsort.append((evalTC[-2][items])/1e6)
    else:
        evalCsort.append((evalTC[0][items]/1e6))

#Combine the cycling states into an easier to manage array
CyclingP = [np.array(evalTC[-1]/1e6), np.array(evalCsort)]

p.figure()

#Making the lines thinner so you can see the structure easily
for items2 in evalTP[:6]:
    p.plot(10000*BmagE, items2/1e6, color = 'green', alpha=0.9 , linewidth = .85)
for items2 in evalTP[6:12]:
    p.plot(10000*BmagE, items2/1e6, color = 'blue' , alpha=0.9 , linewidth = .85)
for items2 in evalTP[12:18]:
    p.plot(10000*BmagE, items2/1e6, color = 'red' , alpha=0.9 , linewidth = .85)
for items2 in evalTP[18:]:
    p.plot(10000*BmagE, items2/1e6, color = 'purple' , alpha=0.9 , linewidth = .85)


#This plots the cycling transitions overlaying the states in the hyperfine structure.
#p.plot(10000*Bmag2, CyclingP[0], color = 'orange', linewidth = 3)
#p.plot(10000*Bmag2, CyclingP[1] , color= 'cyan' , linewidth = 3)

p.show()

#p.title('Zeeman Shift in Rb85 Excited State', fontsize = 18)
p.ylabel('Energy shift (E/h) MHz', fontsize = 16)
p.xlabel('Magnetic Field (G)', fontsize = 16)

import os
os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
#p.savefig('RbHyperfineExcited.jpeg', bbox_inches='tight', dpi=800)
p.show() 


Eshift = 384.23e12

# =============================================================================
# Attempting to have a single plot with both the ground and excited states shown.
# =============================================================================
EshiftP = []
EshiftS = []

for items in evalTPG:
    Etemp = [(Pval)/1e6 for Pval in items]
    EshiftS.append(Etemp)
    
    
for items in evalTP:
    Etemp = [(Pval + 0*Eshift)/1e6 for Pval in items]
    EshiftP.append(Etemp)



fig, (ax, ax2) = p.subplots(2, 1, sharex=True)

BmagG = [10000*bx for bx in BmagCycling]

alphaval = 0.45
for items2 in evalTPG[:5]:
    ax2.plot(BmagG, items2/1e6, color = 'green', alpha=0.25 , linewidth = alphaval)
for items2 in evalTPG[5:]:
    ax2.plot(BmagG, items2/1e6, color = 'blue', alpha=0.25 , linewidth = alphaval)



#Making the lines thinner so you can see the structure easily
for items2 in EshiftP[:6]:
    ax.plot(BmagG, items2, color = 'green', alpha=0.25 , linewidth = alphaval)
for items2 in EshiftP[6:12]:
    ax.plot(BmagG, items2, color = 'blue' , alpha=0.25 , linewidth =alphaval)
for items2 in EshiftP[12:18]:
    ax.plot(BmagG, items2, color = 'red' , alpha=0.25 , linewidth = alphaval)
for items2 in EshiftP[18:]:
    ax.plot(BmagG, items2, color = 'purple' , alpha=0.25 , linewidth = alphaval)


#for items in EshiftS:
#    ax2.plot(BmagG, items, color = 'blue')
#
#for items2 in EshiftP:
#    ax.plot(BmagG, items2, color = 'green')

#Ground state plotting
ax2.plot(BmagG, CyclingPG[1], color = 'brown' , label='|F=3, mf=-3>')
ax2.plot(BmagG, CyclingPG[0], color = 'darkorange' , label='|F=3, mf=3>')
#Excited state plotting

#

ax.plot(BmagG, CyclingP[1], color = 'brown', label='|F\'=4, mf\'=-4')
ax.plot(BmagG, CyclingP[0], color = 'darkorange', label='|F\'=4, mf\'=4')




#p.annotate(s='', xy=(200,1550), xytext=(200,6000), arrowprops=dict(arrowstyle='<->'))


# zoom-in / limit the view to different portions of the data
#ax2.set_ylim(-3e3, 2e3)  # most of the data
ax2.set_ylim(0, 2e3)  # most of the data

ax.set_ylim(-1.5e3, 1.5e3)  # most of the data

#ax.set_ylim((0*Eshift - 1.5e9)/1e9, (0*Eshift + 1.5e9)/1e9)  # outliers only

# hide the spines between ax and ax2
ax.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()

p.xlabel('Magnetic Field (G)',  fontsize = 16)
# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1).  Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.

dl = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass to plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((-dl, +dl), (-dl, +dl), **kwargs)        # top-left diagonal
ax.plot((1 - dl, 1 + dl), (-dl, +dl), **kwargs)  # top-right diagonal

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-dl, +dl), (1 - dl, 1 + dl), **kwargs)  # bottom-left diagonal
ax2.plot((1 - dl, 1 + dl), (1 - dl, 1 + dl), **kwargs)  # bottom-right diagonal

# What's cool about this is that now if we vary the distance between
# ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
# the diagonal lines will move accordingly, and stay right at the tips
# of the spines they are 'breaking'
#p.suptitle('Cycling Transition for Rubidium 85' , fontsize = 18)
fig.text(0.03, 0.5, 'Zeeman Shift (E/h) MHz', ha='center', va='center', rotation='vertical' ,  fontsize = 16)
fig.text(0.15, 0.5, r'$\Delta$E= 384.23 THz', ha='left', va='center',  fontsize = 12)

xy = (200, 1550)
xy2 = (200, 660)
coordsA = "data"
coordsB = "data"
con = ConnectionPatch(xyA=xy, xyB=xy2, coordsA=coordsA, coordsB=coordsB,
                      axesA=ax2, axesB=ax,
                      arrowstyle="<->")
ax2.add_artist(con)

xya = (400, 705)
xyb = (400, -1020)
coordsA = "data"
coordsB = "data"
con2 = ConnectionPatch(xyA=xya, xyB=xyb, coordsA=coordsA, coordsB=coordsB,
                      axesA=ax2, axesB=ax,
                      arrowstyle="<->")
ax2.add_artist(con2)

p.savefig('RbCyclingOnly2.jpeg', bbox_inches='tight', dpi=1200)

p.show()