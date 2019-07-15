# -*- coding: utf-8 -*-
"""
Created on Tue Jan  1 18:47:32 2019

@author: Leo
"""
import os,csv
os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/Thesis Version (NO TOUCHING)/Main Program (Clean Edit)/")
import Generate_Loop_Array as GLA
import numpy as np
from matplotlib import pyplot as p
#import BiasFieldHighB as BBias
import Zeeman_Ideal_Import as ZII
#import matplotlib.style
#import matplotlib as mpl
#import pdb
p.close('all')

        
  
#numval = str(106034)
#numval = '02'


#For +/- 1.75
#fld_a = 22501
#fld_b = 1501

#for +/- 2.0
fld_a= 24001   
fld_b = 1001

#for the model to be trimmed based on the fld_a/fld_b values
#for +/- 1.75
#cp_a = 5000
#cp_b-7999

#for +/- 2.0
cp_a = 5000
cp_b = -11199

#This is for the trimmed CP. For L = 0.68, a = 2500, b = -9199
trim_a = 2500
trim_b = -9199
fld = GLA.Loop_Length[-fld_a:-fld_b]

#Need new values for Sr trimmed plot field:
SRtrim_a = 10500
SRtrim_b = -10699
SRpfld = fld[SRtrim_a: SRtrim_b] 

Solenoid_Dep = []
Solenoid_Len = []

Rb_nb_Cur = []
Rb_wb_Cur = []
Sr_nb_Cur = []
Sr_wb_Cur = []



# =============================================================================
# Use this block for importing a new lmfit solution with the numval format
# =============================================================================
#os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/LmFitSolutions/" +numval+"/")
#with open('Zeeman_Params' + numval + '.csv', mode='r', newline='') as readme:
#    reader = csv.reader(readme, delimiter=',', quoting=csv.QUOTE_NONE)
#    itemp = next(reader)
#    wire_diam = float(float(itemp[2]))
#    Slower_rad = float(float(itemp[4]))
#    cp_initial_offset = float(float(itemp[6]))
#    next(reader)
#    for row in reader:
#        Solenoid_Dep.append(int(float(row[1])))
#        Solenoid_Len.append(int(float(row[3])))        
#        Rb_nb_Cur.append(float(float(row[6])))
#        Rb_wb_Cur.append(float((float(row[7]))))
#        Sr_nb_Cur.append(float(float(row[8])))
#        Sr_wb_Cur.append(float(float(row[9])))
#        
# =============================================================================
# This block opens up the current best solution from the main directory.    
# =============================================================================
os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/Thesis Version (NO TOUCHING)/Main Program (Clean Edit)/")

#with open('SigmaPlusSolution2.csv', mode='r', newline='') as readme:
with open('SigmaPlusSolution2.csv', mode='r', newline='') as readme:

    reader = csv.reader(readme, delimiter=',', quoting=csv.QUOTE_NONE)
    itemp = next(reader)
    wire_diam = float(float(itemp[2]))
    Slower_rad = float(float(itemp[4]))
    cp_initial_offset = float(float(itemp[6]))
    next(reader)
    for row in reader:
        Solenoid_Dep.append(int(float(row[1])))
        Solenoid_Len.append(int(float(row[3])))        
        Rb_nb_Cur.append(float(float(row[6])))
        Rb_wb_Cur.append(float((float(row[7]))))
        Sr_nb_Cur.append(float(float(row[8])))
        Sr_wb_Cur.append(float(float(row[9])))


#Now that we've pulled our input data, we'll aim our output location here rather than at each plot I save.
os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
#os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Presentation April19")
wire_diam_steps = int(wire_diam*10000)
CP_Count = len(Solenoid_Len)
#CP_Count = 15
Plot_SLA = GLA.Slower_Loop_Array(wire_diam, Slower_rad)

def offsetgen(cplenin):
    tempoffset = []
    tempy = cp_initial_offset
    for cplen in cplenin:
        tempoffset.append(tempy)
        tempy = tempy + cplen
    return tempoffset
#
def Slower_CoilPacks_Short(scpCoil_Depth, scpCoil_Current, scpCoil_Length, scpCoil_Offset):
    BCP_Totes = np.zeros_like(fld)
    
    Position = np.arange((int(float(scpCoil_Offset))), int(float(scpCoil_Offset+scpCoil_Length)) ,1)
    for pw in range(int(scpCoil_Length)):
        BCP_Totes =  BCP_Totes + scpCoil_Current*(Plot_SLA[int(scpCoil_Depth)][(-fld_a-(Position[pw]*wire_diam_steps)):(-fld_b-(Position[pw]*wire_diam_steps))])
    return BCP_Totes

    
def Slower_CoilPacks_Sparseouter(scpCoil_Depth, scpCoil_Current, scpCoil_Length, scpCoil_Offset):
    BCP_Totes = np.zeros_like(fld)
    
    Position = np.arange((int(float(scpCoil_Offset))), int(float(scpCoil_Offset+scpCoil_Length)) ,1)
    for pw in range(int(scpCoil_Length)):
        if pw%2 ==0:
            BCP_Totes =  BCP_Totes + scpCoil_Current*(Plot_SLA[int(scpCoil_Depth)][(-fld_a-(Position[pw]*wire_diam_steps)):(-fld_b-(Position[pw]*wire_diam_steps))])
        else:
            BCP_Totes =  BCP_Totes + scpCoil_Current*(Plot_SLA[int(scpCoil_Depth-1)][(-fld_a-(Position[pw]*wire_diam_steps)):(-fld_b-(Position[pw]*wire_diam_steps))])
    return BCP_Totes

def CP_Total_B(Coil_Depth, Coil_Current, Coil_Length):
    B_Totals  = np.zeros_like(fld)
    offset = offsetgen(Solenoid_Len)
    offset[1] = 0
    for nn in range(CP_Count):
        B_Totals  = B_Totals  + Slower_CoilPacks_Short(Coil_Depth[nn], Coil_Current[nn], Coil_Length[nn], offset[nn])     
    del offset, Coil_Depth, Coil_Current, Coil_Length, nn
    return B_Totals   


def CP_Total_B_Trimmed(x, Coil_Depth, Coil_Current, Coil_Length, offsetme):
    B_Toteme = np.zeros_like(fld)
    B_Toteme = Slower_CoilPacks_Short(Coil_Depth, Coil_Current, Coil_Length, offsetme)
    del offsetme, Coil_Depth, Coil_Current, Coil_Length
    return B_Toteme[cp_a:cp_b]

def CP_Total_B2(Coil_Depth, Coil_Current, Coil_Length, offsetme):
    B_Toteme = np.zeros_like(fld)
    B_Toteme = Slower_CoilPacks_Short(Coil_Depth, Coil_Current, Coil_Length, offsetme)
    del offsetme, Coil_Depth, Coil_Current, Coil_Length
    return B_Toteme

def CP_Total_B_Full(Coil_Depth, Coil_Current, Coil_Length):
    B_Totals  = np.zeros_like(fld)
    offsetcp = offsetgen(Solenoid_Len)
    offsetcp[1] = 0
    for nn in range(CP_Count):
        B_Totals  = B_Totals  + Slower_CoilPacks_Short(Coil_Depth[nn], Coil_Current[nn], Coil_Length[nn], offsetcp[nn])   
    del offsetcp, Coil_Depth, Coil_Current, Coil_Length, nn
    return B_Totals[trim_a:]   



def CP_Total_B_Plot_Trimmed(Coil_Depth, Coil_Current, Coil_Length):
    B_Totals  = np.zeros_like(fld)
    offsetcp = offsetgen(Solenoid_Len)
    offsetcp[1] = 0
    for nn in range(CP_Count):
        B_Totals  = B_Totals  + Slower_CoilPacks_Short(Coil_Depth[nn], Coil_Current[nn], Coil_Length[nn], offsetcp[nn])   
    del offsetcp, Coil_Depth, Coil_Current, Coil_Length, nn
    return B_Totals[trim_a:trim_b]   

def CP_Total_B_Plot_TrimmedSR(Coil_Depth, Coil_Current, Coil_Length):
    B_Totals  = np.zeros_like(fld)
    offsetcp = offsetgen(Solenoid_Len)
    offsetcp[1] = 0
    for nn in range(CP_Count):
        B_Totals  = B_Totals  + Slower_CoilPacks_Short(Coil_Depth[nn], Coil_Current[nn], Coil_Length[nn], offsetcp[nn])   
    del offsetcp, Coil_Depth, Coil_Current, Coil_Length, nn
    return B_Totals[SRtrim_a:SRtrim_b]   


def CP_Total_B_Plot_TrimmedSparse(Coil_Depth, Coil_Current, Coil_Length):
    B_Totals  = np.zeros_like(fld)
    offsetcp = offsetgen(Solenoid_Len)
    offsetcp[1] = 0
    for nn in range(CP_Count):
        B_Totals  = B_Totals  + Slower_CoilPacks_Sparseouter(Coil_Depth[nn], Coil_Current[nn], Coil_Length[nn], offsetcp[nn])   
    del offsetcp, Coil_Depth, Coil_Current, Coil_Length, nn
    return B_Totals[trim_a:trim_b]   



# =============================================================================
# 
# =============================================================================
#def Slower_CoilPacks_Short2(scpCoil_Depth, scpCoil_Current, scpCoil_Length, scpCoil_Offset):
#    BCP_Totes = np.zeros_like(fld)
#    
#    Position = np.arange((int(float(scpCoil_Offset))), int(float(scpCoil_Offset+scpCoil_Length)) ,1)
#    for pw in range(int(scpCoil_Length)):
#        BCP_Totes =  BCP_Totes + scpCoil_Current*(GLA.SLA2[int(scpCoil_Depth)][(-fld_a-(Position[pw]*wire_diam_steps)):(-fld_b-(Position[pw]*wire_diam_steps))])
#    return BCP_Totes
#
#def CP_Total_B_Plot_Trimmed2(Coil_Depth, Coil_Current, Coil_Length):
#    B_Totals  = np.zeros_like(fld)
#    offsetcp = cp_initial_offset 
#    for nn in range(CP_Count):
#        B_Totals  = B_Totals  + Slower_CoilPacks_Short2(Coil_Depth[nn], Coil_Current[nn], Coil_Length[nn], offsetcp)
#        offsetcp = offsetcp + Coil_Length[nn]       
#    del offsetcp, Coil_Depth, Coil_Current, Coil_Length, nn
#    return B_Totals[2500:-7999]  






# =============================================================================
# 
# =============================================================================
# =============================================================================
# Importing the solutions found from the Optimizer then plotting them rather than re-obtaining them every time you plot them 
# =============================================================================



#CP_Count = len(Solenoid_Len)   
#CP_Count = 15
#
#min_coil_len = 6900 + 0.5*cp_initial_offset*wire_diam_steps
#Solenoid_Len = [int(x) for x in (CP_Count*[1.01*min_coil_len/(CP_Count*wire_diam_steps)])]
#Solenoid_Dep = []
#Rb_nb_Cur = []
#Rb_wb_Cur = []
#for cpn in range(CP_Count):
#    Solenoid_Dep.append(int(5 + 11*np.sqrt(1-cpn/CP_Count)))
#    Rb_nb_Cur.append(-1.475 + 4.9*np.sqrt(1-0.45*cpn/CP_Count))
#    Rb_wb_Cur.append(-4.1 + 7.35*np.sqrt(1-0.45*cpn/CP_Count))
#    
#Rb_wb_Cur[-1] = 0
#print(Rb_wb_Cur)

# ===========================================================


#result_Rb_bias = result_Rb_nobias
#result_Sr_nobias = result_Rb_nobias
#result_Sr_bias = result_Rb_nobias
       

Rb_bias_fit = GLA.Target_Field
Rb_nobias_fit = GLA.Ideal_Field
#Sr_bias_fit = GLA.Target_Sr
#Sr_nobias_fit = GLA.Ideal_Sr

#
#Rb_bias_fit = GLA.Inc_tar
#Rb_nobias_fit = GLA.Ideal_Inc
#

# =============================================================================
# Rubidium Slower without High-B Bias Field
# =============================================================================
#
Rb_nobias = CP_Total_B_Plot_Trimmed(Solenoid_Dep , Rb_nb_Cur, Solenoid_Len)

p.figure()
p.title('Rubidium With No Bias Field' ,fontsize = 18)
#p.title('Ideal Zeeman Slower Filed')
p.ylim(-30.0 , 550)
p.ylabel('Magnetic Field (G)', fontsize = 16)
p.xlabel('Distance along designed slower length (m)', fontsize = 16)
p.plot(GLA.Plot_Range, Rb_nobias, ':', label='Zeeman Coil Solution')
p.plot(GLA.Slower_Length, Rb_nobias_fit, '-' , color='green', label="Ideal Rb Zeeman Field", alpha = 0.6)
#p.plot(GLA.Slower_Length, Rb_bias_fit, '-' , label="Target Field")
p.legend(loc='best')
p.grid(True, which='both')


#p.savefig('Rb_NoBias.pdf', dpi=1200, bbox_inches='tight', pad_inches=0.1)
#p.savefig('Rb_NoBias2.jpeg', dpi=800, bbox_inches='tight', pad_inches=0.1)
p.show()





# =============================================================================
# Rubdium Slower with High-B Bias Field
# =============================================================================


#Rbc = 5*[2.5]
#Rbd = 5*[10]
#Rbl = 5*[75]
#rbo = offsetgen(Rbl)
#p.figure()
#for rn in range(5):
#    Brbt = Slower_CoilPacks_Sparseouter(Rbd[rn], Rbc[rn], Rbl[rn], rbo[rn])
#    Brt2 = Slower_CoilPacks_Short(Rbd[rn], Rbc[rn], Rbl[rn], rbo[rn])
#    p.plot(fld, Brbt)
#    p.plot(fld, Brt2, ':')
#p.show()


#This little loop plots the magnetic field of each component solenoid.
#rboffs = offsetgen(Solenoid_Len)
#p.figure()
#for rn in range(CP_Count):
#    Brbt = Slower_CoilPacks_Short(Solenoid_Dep[rn], Rb_wb_Cur[rn], Solenoid_Len[rn], rboffs[rn])
#    p.plot(fld, Brbt)
#p.show()


Rb_bias = CP_Total_B_Plot_Trimmed(Solenoid_Dep , Rb_wb_Cur, Solenoid_Len)
bfld = fld[4500:13801]

Blonger = GLA.B_Bias_Longer
Rb_bias_long = Rb_bias[2000:] + Blonger

p.figure()
p.title('Rubidium With Bias Field', fontsize = 18)
p.ylim(-30.0 , 550)
p.ylabel('Magnetic Field  (G)', fontsize = 16)
p.xlabel('Distance along designed slower length (m)', fontsize = 16)
p.plot(GLA.Plot_Range, Rb_bias, ':', color = 'magenta',  label='Zeeman Coil Solution')
#p.plot(GLA.Plot_Range, Rb_nobias, ':',  color='magenta' ,label='No bias Field Solution')
#p.xlim(-0.25, 0.85 )
p.plot(bfld, Rb_bias_long, ':', color ='red', label = 'Combined Magnetic Fields')
p.plot(GLA.Slower_Length, GLA.Ideal_Field, color='orange' , label="Ideal Rb Zeeman Field" , alpha = 0.6)
p.plot(GLA.Slower_Length, Rb_bias_fit, '-' , color='green', label="Zeeman Coil Field" , alpha = 0.6)
p.plot(GLA.Slower_Length, GLA.B_Bias, '-.' , color='red', label= "High-B Bias Field")
#p.plot(BBias.zlin, BBias.B_Mag,'-.' , color='red', label= "High-B Bias Field" )
#p.plot(GLA.Slower_Length, Rb_bias_res.best_fit, label="Best Fit")
p.legend(loc='best')
p.grid(True, which='both')
#p.savefig('Rb_With_Bias.pdf', dpi=1200, bbox_inches='tight', pad_inches=0.1)
#p.savefig('Rb_With_Bias.jpeg', dpi=800, bbox_inches='tight', pad_inches=0.1)

p.show()


#os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
#p.savefig('ZeemanTarget.pdf', dpi=1200, bbox_inches='tight', pad_inches=0.1)
#p.savefig('ZeemanNoLegend.png', dpi=1200, bbox_inches='tight', pad_inches=0.1)
#p.savefig('Rb_With_Bias_all' +numval + '.jpeg', dpi=800, bbox_inches='tight', pad_inches=0.1)




# =============================================================================
# This is for a sparse outer coil.
# =============================================================================

#Rb_biassprs = CP_Total_B_Plot_TrimmedSparse(Solenoid_Dep , Rb_wb_Cur, Solenoid_Len)
#bfld = fld[4500:13801]
#
#Blonger = GLA.B_Bias_Longer
#Rb_bias_longsprs = Rb_biassprs[2000:] + Blonger

#p.figure()
#p.title('Rubidium With Bias Field', fontsize = 18)
#p.ylabel('Magnetic Field  (G)', fontsize = 16)
#p.xlabel('Distance along designed slower length (m)', fontsize = 16)
##p.plot(GLA.Plot_Range, Rb_biassprs, ':', label='Zeeman Coil Solution')
##p.plot(GLA.Plot_Range, Rb_nobias, ':',  color='magenta' ,label='No bias Field Solution')
##p.xlim(-0.25, 0.85 )
##p.plot(bfld, Rb_bias_longsprs, ':', color ='red', label = 'Combined Magnetic Fields')
#p.plot(GLA.Slower_Length, GLA.Ideal_Field, color='orange' , label="Ideal Rb Zeeman Field")
#p.plot(GLA.Slower_Length, Rb_bias_fit, '-' , color='green', label="Zeeman Coil Field")
#p.plot(GLA.Slower_Length, GLA.B_Bias, '-.' , color='red', label= "High-B Bias Field")
##p.plot(BBias.zlin, BBias.B_Mag,'-.' , color='red', label= "High-B Bias Field" )
##p.plot(GLA.Slower_Length, Rb_bias_res.best_fit, label="Best Fit")
#p.legend(loc='best')
##p.savefig('Rb_With_Bias.pdf', dpi=1200, bbox_inches='tight', pad_inches=0.1)
##p.savefig('Rb_With_Bias2.png', dpi=800, bbox_inches='tight', pad_inches=0.1)
#
#p.show()


#
#
##=============================================================================
## Strontium Slower with High-B Bias Field
### =============================================================================
#
#Sr_bias = CP_Total_B_Plot_TrimmedSR(Solenoid_Dep , Sr_wb_Cur, Solenoid_Len)
##
#p.figure()
#p.title('Strontium With Bias Field')
#p.ylabel('Magnetic Field Strength (G)')
#p.xlabel('Distance along designed slower length (m)')
#p.plot(SRpfld, Sr_bias, ':', label='Slower Field Solution')
#p.plot(GLA.Srfield, GLA.Ideal_Sr, ':' , label="Ideal Sr Zeeman Field")
#p.plot(GLA.Srfield, Sr_bias_fit, '-' , label="Target Field")
#p.plot(GLA.Srfield, GLA.B_BiasSR, '-.' , label= "High-B Bias Field")
#p.legend(loc='best')
##p.savefig('Sr_With_Bias.png', dpi=800, bbox_inches='tight', pad_inches=0.1)
#p.show()
##
####=============================================================================
#### Strontium Slower without High-B Bias Field
#### =============================================================================
#
#Sr_nobias = CP_Total_B_Plot_TrimmedSR(Solenoid_Dep , Sr_nb_Cur, Solenoid_Len)
#
#p.figure()
#p.title('Strontium With No Bias Field')
#p.ylabel('Magnetic Field Strength (G)')
#p.xlabel('Distance along designed slower length (m)')
#p.plot(SRpfld, Sr_nobias, ':', label='Slower Field Solution')
#p.plot(GLA.Srfield, Sr_nobias_fit, '-' , label="Ideal Sr Zeeman Field")
##p.plot(GLA.Srfield, Rb_nobias_res.best_fit, label='Best Fit')
##p.plot(GLA.Srfield, Rb_nobias_res.init_fit, label='Initial Fit')
#p.legend(loc='best')
##p.savefig('Sr_NoBias.png', dpi=800, bbox_inches='tight', pad_inches=0.1)
##p.savefig('Sr_No_Bias' + numval + '.jpeg', dpi=1200, bbox_inches='tight', pad_inches=0.1)
#p.show()

#rboffs = offsetgen(Solenoid_Len)
#p.figure()
#for rn in range(CP_Count):
#    BrSr = Slower_CoilPacks_Short(Solenoid_Dep[rn], Sr_nb_Cur[rn], Solenoid_Len[rn], rboffs[rn])
#    p.plot(fld, BrSr)
#p.show()


# =============================================================================
# Extending the fields further out to see how the gradient and adiabaticity look, particularly the left oven side of the slower as we can assume the atom velocity is still max there.
# =============================================================================
#
Gradient_int = 1
Bfieldgradient = np.gradient(Blonger, Gradient_int)
B_Bias_longer= GLA.B_Bias_Longer

#rbadbsprs = Rb_biassprs
rbaddb = Rb_bias_long
#rdaddwb = np.gradient(Rb_biassprs, 1)
Rbnb = np.gradient(Rb_nobias[2000:], Gradient_int)
Rbwb = np.gradient(rbaddb, Gradient_int)
dbfld = fld[4500:13801]



p.figure()
p.title('Adiabaticity for Rb Slower', fontsize = 18)
p.ylim(-0.55, 0.55)
p.ylabel('dB/dz', fontsize = 16)
p.xlabel('Distance along designed slower length (m)', fontsize = 16)
p.plot(dbfld, GLA.Rb_Adiabatic2, '-.', color='magenta', label = 'Adiabatic Limit for Rb')
p.plot(dbfld, -GLA.Rb_Adiabatic2, '-.', color='magenta')
#p.plot(dbfld, Bfieldgradient, ':', color='blue', label='High-B Bias')
#p.plot(GLA.Plot_Range, rdaddwb, label = 'Rb_with_bias')
#p.plot(GLA.Slower_Length, GLA.Rb_Adiabatic_limit)
#p.plot(dbfld, Rbnb, label = 'Rb_nobias')
p.plot(bfld, Rbwb, color='orange', label = 'Rb_withbias')
p.legend(loc='best')
p.grid(True, which='both')

xyz = np.linspace(-0.25, 0.25, 50)
xxy = 50*[0]
xxend = 50*[0.68]
p.plot(xxend, xyz, color='green')
p.plot(xxy, xyz, color='red')
#p.savefig('Rb_AdiabaticityDAMOP.jpeg', dpi=800, bbox_inches='tight', pad_inches=0.1)
#p.savefig('Rb_Adiabaticbreak.pdf', dpi=1200, bbox_inches='tight', pad_inches=0.1)
#p.savefig('Rb_Adiabaticity.pdf', dpi=1200, bbox_inches='tight', pad_inches=0.1)
p.show()

# =============================================================================
# Strontium Adiabaticity
# =============================================================================
#B_Bias_longerSR = GLA.B_Bias_LongerSR
#sraddb = Sr_nobias + B_Bias_longerSR
#Srnb = np.gradient(Sr_nobias, 1)
#Srwb = np.gradient(Sr_bias, 1)
#dbfldsr = fld[10750:12051]
#
#p.figure()
#p.title('Adiabaticity for Sr Slower', fontsize = 18)
#p.ylabel('dB/dz')
#p.xlabel('Distance along designed slower length (m)', fontsize = 16)
#p.plot(dbfldsr, GLA.Sr_Adiabatic2, '-.', color='magenta', label = 'Adiabatic Limit for Sr')
#p.plot(dbfldsr, -GLA.Sr_Adiabatic2, '-.', color='magenta')
#p.plot(dbfldsr, Srnb[250:-250], label = 'Sr_nobias')
#p.plot(dbfldsr, Srwb[250:-250], label = 'Sr_withbias')
##p.savefig('Sr_Adiabaticity.jpeg', dpi=800, bbox_inches='tight', pad_inches=0.1)
#p.show()

 



# =============================================================================
# Do a numerical integration of the curve fit field for the velocity and position based on the actual magnetic field.
# We'll do this using fourthe order Runge-Kutta integration
# =============================================================================
#Bverylong = CP_Total_B_Full(Solenoid_Dep , Rb_nb_Cur, Solenoid_Len)[2500:]
Blong2 = Rb_bias_long[500:]
#Blong2 = [(xb+10) for xb in Rb_bias_long[500:]]

# =============================================================================
# This function uses the magnetic field of the solution so that it can be fed based on position into the a_accel2 function for each location.
# This was needed as, prior to running the numerical integration, I don't know what the positions would be.
# =============================================================================
def Bzz(bpos):
    bpv = int(10000*round(bpos, 4)) 
#    print(bpv)
    bzt = Blong2[bpv]
    return bzt

# =============================================================================
# This provides the positional ideal magnetic field for a given eta and detuning.
# =============================================================================
def Bzz2(pos):
    return ZII.Zeeman_B_field(pos, ZII.eta, 150e6)

# =============================================================================
# This should give the velocity based on the magnetic field for a position.
# =============================================================================
def vzz(vpos):
    return (-ZII.delta_o + ZII.cgamma*Bzz(vpos))*ZII.wavelength 

#This one works! Huzzah! But it shouldn't. No really, the units are horribly wrong.
#def a_acel(apos, avel, dh): 
#    zdt = ZII.delta_o - ZII.cgamma*Bzz2(apos) + avel/ZII.wavelength
#    aacel = -ZII.gamma*ZII.s_o*ZII.cgamma**2*ZII.hbar*ZII.k/((8*(zdt**2) +2*ZII.gamma + 2*ZII.s_o*ZII.cgamma**2)*ZII.M_rb)
#    return aacel
    
#This one is directly from the papers quoted everywhere. The only difference that I can tell is that I simplified the expression out front
#therefor removing what apparently was rounding error from the coefficients.
    

def a_acel2(apos, avel, dh): 
    zdt = ZII.delta_o - ZII.cgamma*Bzz(apos) + avel/ZII.wavelength
    aacel2 = - ZII.a_o*ZII.s_o/(1 + ZII.s_o + 4*(zdt/ZII.gamma)**2)
    return aacel2




# =============================================================================
# Define the R4k function for numerical integration
# =============================================================================

def rk4(x, v, a, dt):
    """Returns final (position, velocity) tuple after
    time dt has passed.

    x: initial position (number-like object)
    v: initial velocity (number-like object)
    a: acceleration function a(x,v,dt) (must be callable)
    dt: timestep (number)"""
    x1 = x
    v1 = v
    a1 = a(x1, v1, 0)

    x2 = x + 0.5*v1*dt
    v2 = v + 0.5*a1*dt
    a2 = a(x2, v2, dt/2.0)

    x3 = x + 0.5*v2*dt
    v3 = v + 0.5*a2*dt
    a3 = a(x3, v3, dt/2.0)

    x4 = x + v3*dt
    v4 = v + a3*dt
    a4 = a(x4, v4, dt)

    xf = x + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
    vf = v + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)

    return xf, vf

stepsize = 0.00001


# =============================================================================
# Pull some constants here to avoid having to constantly call them from ZII.
# =============================================================================
k = 2*np.pi/ZII.wavelength
c_d = 7.145e-7 #The constant h/mu_B in units of Gauss/Hz
cgamma = ZII.mu_b/(10000*ZII.hh) #"Zeeman constant" in terms of Hz/Gauss
cgamma2 = ZII.mu_b/(10000*ZII.hbar) #"Zeeman constant" in terms of rad/Gauss*s

# =============================================================================
# Define an array for detuning values that will be looped over.
# =============================================================================

#detunings = [50e6, 75e6, 110e6 , 120e6, 130e6, 140e6, 150e6, 160e6, 170e6]
detunings = [50e6, 110e6 , 130e6, 150e6, 160e6]


# =============================================================================
# Start the plot for the numerically integrated velocity vs position plot for multiple detunings.
# =============================================================================

#mpl.sytle.use('classic')
p.figure()
p.title('Rb Velocity in the Slower' , fontsize = 18)
p.ylabel('Rb Velocity (m/s)' , fontsize = 16)
p.xlabel('Distance along slower' , fontsize = 16)
p.ylim(0, 360)
p.grid(True, which='both')
for detunevals in detunings:

    zzdetune = detunevals
# =============================================================================
#   We define the acceleration function inside the loop so that the proper new detuning value is used in the function. Otherwise, python didn't update the function properly.
# =============================================================================
    def a_acel3(apos, avel, dh): 
        zdt = 2*np.pi*zzdetune - cgamma2*Bzz(apos) + avel*k
    
        aacel = -ZII.gamma*ZII.s_o*ZII.hbar*k/((4*(zdt/ZII.gamma)**2 +1 + ZII.s_o)*2*ZII.M_rb)
        return aacel
    
    va = [350]
    xa = [0.0]
    xf = 0
    vf = 25
    nn = 0
    
    
    ##This loop integrates for position, useful if integration works out well.
#    while xf < 0.695:
#        tempy = rk4(xa[nn], va[nn], a_acel3, stepsize)
#        va.append(tempy[1])
#        xa.append(tempy[0])
#        xf = xa[nn]
#        nn = nn+1
   
    while xf < 0.85:
        tempy = rk4(xa[nn], va[nn], a_acel3, stepsize)
        va.append(tempy[1])
        xa.append(tempy[0])
        vf = va[nn]
        xf = xa[nn]
        nn = nn+1
        if vf <=0:
            break

    #This loop checks the velocity for it's limit and great to check when things aren't working.
    #while vf >=20:
    #    tempy = rk4(xa[nn], va[nn], a_acel3, stepsize)
    #    va.append(tempy[1])
    #    xa.append(tempy[0])
    #    vf = va[nn]
    #    nn = nn + 1
    ##
    

    p.plot(xa, va, label= r'$\Delta\omega$ (MHz) - ' + str(zzdetune/1e6))
#    p.legend(bbox_to_anchor=(1.04, .5), loc='center left')
    p.legend(bbox_to_anchor=(0, 1.06, 1, 0.2), loc='lower center' , ncol = 3)

    p.show()
    
     
p.savefig('RbVelDetuneNum.jpeg', dpi=1200, bbox_inches='tight', pad_inches=0.1)
    
    
    
    
    
    
# =============================================================================
# The original plotting block for a single detuning.  
# =============================================================================
    
#p.figure()
## # p.plot(ifield, v_squared, '.')
#p.title('Rb Velocity in the Slower')
#p.ylabel('Rb Velocity (m/s)')
#p.xlabel('Distance along slower')
#p.plot(ZII.ifield, ZII.v_slower, label='Original Method')
##p.plot(ifield, v3, label='gamma velocity')
#p.plot(xa, va, label='Atom-Frame Numerical Integration')
##p.plot(xa,va2, label='Lab-Frame Numerical Integration')
#p.legend(loc='best')
#p.ylim(0,360)
##p.savefig('Rb_Velocity_Numerical.pdf', dpi=1200, bbox_inches='tight', pad_inches=0.1)
#
#p.show()

#p.close('all')



# =============================================================================
# Calculate the linear portion for the Zeeman shift: This is a rough approximation and doesn't include any information regarding the state of the atom.
# =============================================================================
Zeemanshift = [1.4e6*bb for bb in Blong2]

# =============================================================================
# Calculate the Doppler shift for the velocity curve then find the difference between Zeeman, Doppler and total detuning.
# =============================================================================
Dopplershift = [(1/ZII.wavelength)*vax for vax in va]

Dopplershift_ideal = (1/ZII.wavelength)*ZII.v_slower2[500:]
diff = (Zeemanshift - Dopplershift_ideal - ZII.delta_o)/1e6

bfld2 = bfld[500:]
freqoffs = 8801*[ZII.delta_o*1e-6]  #The offset, delta_o, of the entire process in MHz.

# =============================================================================
# Plot it.
# =============================================================================
p.figure()
p.title('Detuneing along slower')
p.ylabel('Frequency Shift (MHz)')
p.xlabel('Position along slower (m)')
p.ylim(0, 1.1*ZII.delta_o*1e-6)
#p.plot(bfld2, Zeemanshift, label='Z shift')

#p.plot(xa, Dopplershift, label='D shift')
p.plot(bfld2, freqoffs, label= 'Delta_o initial detune')
p.plot(bfld2, diff, label='Zeeman Detuning from resonant Velocity')
#p.plot(bfld2, diffp40, label= 'B+40G')
#p.plot(bfld2, 1e-6*diffm40, label=' B-40G')
p.legend(loc='best')
#p.savefig('Rb_Detuneplot.jpeg', dpi=600, bbox_inches='tight', pad_inches=0.1)

p.show()
#
#
# =============================================================================
# Calculations to determine the amounting of heating this system will be doing
# The Joule heating equation is P=RI^2 for DC currents. R = rho*l/A where rho is resitivity
# l is length of wire, A is cross sectional area for the current. 
# =============================================================================
# =============================================================================
# The first loop here cycles through all of the loops in each solenoid and obtains the total length of wire per solenoid.
# =============================================================================

wire_len = []
for ii in range(CP_Count):
    wiretmp = 0
    for ix in range(Solenoid_Dep[ii]):
        temprad = (Slower_rad + 0.5*(0.5+ix)*wire_diam)
        wiretmp = wiretmp + 2.0*np.pi*temprad  
    wire_len.append(wiretmp*Solenoid_Len[ii])

def JHeat2(current_in):
    B_Heat = []
    JHC = (1.68e-8)*(1/(np.pi*(0.5*wire_diam)**2))
    for ii in range(CP_Count):
        btemp = 0
        for ix in range(Solenoid_Dep[ii]):
            temprad = (Slower_rad + 0.5*(0.5+ix)*wire_diam)
            btemp = btemp + JHC*(current_in[ii]**2)*(2*np.pi*temprad)
        B_Heat.append(btemp*Solenoid_Len[ii])
    return B_Heat
#
Rb_nb_heat2 = JHeat2(Rb_nb_Cur)    
Rb_wb_heat2 = JHeat2(Rb_wb_Cur)
Sr_nb_heat2 = JHeat2(Sr_nb_Cur)
Sr_wb_heat2 = JHeat2(Sr_wb_Cur)
    
    
    
    
def JouleHeat(curin):
    Jouleout = []
    JHC = (1.68e-8)*(1/(np.pi*(0.5*wire_diam)**2))
    for ci in range(CP_Count):
        Jouleout.append(float('%.4f'%((curin[ci]**2)*wire_len[ci]*JHC)))
    return Jouleout



Rb_nb_heat = JouleHeat(Rb_nb_Cur)
Rb_wb_heat = JouleHeat(Rb_wb_Cur)
Sr_nb_heat = JouleHeat(Sr_nb_Cur)
Sr_wb_heat = JouleHeat(Sr_wb_Cur)




# =============================================================================
# Now to determine current density for all four of the slowers. The thought is # of wires per pack*current in that pack / volume of the pack.
# The values are in terms of cm and cm^3 for the results.The power is in terms of mw, so the current density is mW/cm^3
# =============================================================================

def CurrentDensity(Currentin):
    Cur_out = []
    for ci in Currentin:
        Cur_out.append((ci/(wire_diam**2))*(0.01**2))
    return Cur_out

Rb_nb_cd = CurrentDensity(Rb_nb_Cur)
Rb_wb_cd = CurrentDensity(Rb_wb_Cur)
Sr_nb_cd = CurrentDensity(Sr_nb_Cur)
Sr_wb_cd = CurrentDensity(Sr_wb_Cur)



# =============================================================================
# Not sure where I'm going with this, so I'll also grab dP/dV, power per unit volume, dP/dV = J^2*rho
# =============================================================================


def PowerperVol(cd_in):
    ppv = []
    for fdi in range(CP_Count):
        ppv.append((((10000*cd_in[fdi])**2)*1.68E-8)*(.01**3))
    return ppv

Rb_nb_ppv = PowerperVol(Rb_nb_cd)
Rb_wb_ppv = PowerperVol(Rb_wb_cd)
Sr_nb_ppv = PowerperVol(Sr_nb_cd)
Sr_wb_ppv = PowerperVol(Sr_wb_cd)



# =============================================================================
# Expanding on that, we'll take this power density and multiply it by the volume of the solenoid (excluding the central vacuum tube since it doesn't generate power)
# =============================================================================

def PPVtoPower(ppvin):
    TPow = []
    for ppv in range(CP_Count):
        TPow.append((ppvin[ppv]*wire_diam*Solenoid_Len[ppv]*np.pi*(Solenoid_Dep[ppv]*wire_diam)**2)*(100**3))
    return TPow

Rb_nb_TPow = PPVtoPower(Rb_nb_ppv) 
Rb_wb_TPow = PPVtoPower(Rb_wb_ppv)
Sr_nb_TPow = PPVtoPower(Sr_nb_ppv) 
Sr_wb_Tpow = PPVtoPower(Sr_wb_ppv) 


