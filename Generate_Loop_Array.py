##The purpose of this function is to generate an array of current-normalized wrappings.
# The elements of this array are arrays that contain the magnetic field values for each location in range of 1m from each loop of current.
import numpy as np
import os
os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/Thesis Version (NO TOUCHING)/Main Program (Clean Edit)/")

import BiasFieldHighB as BBias
import BiotSavartCalc as BSC
import Zeeman_Ideal_Import as ZII
import Zeeman_Ideal_Import_Sr as ZIISR

from matplotlib import pyplot as p

#Change this to be whatever directory you place all your files in.
# os.chdir(r"C:/coding/Python")

##Create the arrays for the magnetic field to be evalulated over.
# =============================================================================
# Create the main arrays to evaluate the magnetic field. The Loop_Length is the Magnetic Field for a single loop of wire. This wire needs to be
# ~3x longer than the length of the slower so that you don't lose contributions to far edges of the slower.
# The Slower_L and Slower_Length round off the length of the slower so that we avoid the asymptotic nature of the quadratic as z approaches L.
# Plot range takes the Slower_Length and extends the field so that we can see behavior outside the slower length, as well the adiabaticity.
# =============================================================================

Loop_Length = np.linspace(-2.0,1.9, 39001)
Loop_Field = [0,0, Loop_Length]
Slower_L = round(ZII.L, 3) - .01
SL_Steps = round(Slower_L*10000)+1 #stepsize of 0.1 mm
Slower_Length = np.linspace(0,Slower_L, SL_Steps)
field = [(0),(0), (Slower_Length)]
Plotupper = Slower_L + 0.2
Ps = SL_Steps + 4500
Plot_Range= np.linspace(-0.25, Plotupper, Ps)

# =============================================================================
# Strontium Field 
# =============================================================================

Slower_LSR = round(ZIISR.L, 3) 
SL_StepsSR = round(Slower_LSR*10000)+1 #stepsize of 0.1 mm
Slower_LengthSR = np.linspace(0,Slower_LSR, SL_StepsSR)
fieldSR = [(0), (0), (Slower_LengthSR)]

# =============================================================================
# This scaling factor allows the calculations to be done aren't accidentally rounded down as python sometimes hates really small values.
# =============================================================================

Scaling = 10**(-7)
Bconst = Scaling*10000 #Gives Gauss Output

# =============================================================================
# Define the physical parameters of the slower.
# =============================================================================

Slower_radius = 0.0191 #Vacuum tube outer radius (1 inche) 
wire_diameter = 0.0015 #Roughly AWG 15/16

# =============================================================================
# Set the maximum thickness for any particular location along wire. Should never have more than 50.
# =============================================================================

max_wraps = 50

# =============================================================================
# Define a function that takes wire diameter and slower radius input and provides the array of combined magnetic fields for stacking
# n loops of wire at a single position along the slower. This one only does the axial field calculation.
# =============================================================================

def Slower_Loop_Array(wire_d, slow_rad):
    B_Vert=[]
    B_Vtemp = np.zeros_like(Loop_Length)
    for bvn in range(max_wraps):
        B_Vtemp = B_Vtemp + Bconst*(BSC.axialmag(Loop_Field, ([0,0,0]), (slow_rad+(0.5+bvn)*wire_d)))
        B_Vert.append(B_Vtemp)
    return B_Vert

# =============================================================================
# Taking a 3D fieldspace as the input and doing the same total-field calculation for anywhere in space.
# =============================================================================

def Slower_Loop3d(positional_offset, wire_d):
    B_Vert=[]
    B_Vtemp = np.zeros_like(Loop_Length)
    for bvn in range(max_wraps):
        B_Vtemp = B_Vtemp + Bconst*(BSC.BiotSavart_Ring(Loop_Field, (positional_offset), (Slower_radius+(0.5+bvn)*wire_d))[2])
        B_Vert.append(B_Vtemp)
    return B_Vert    

# =============================================================================
# Importing the Bias field arrays as local variables.
# =============================================================================
    
B_Bias = BBias.B_Bias
B_Bias_Longer = BBias.B_Bias_Long



# =============================================================================
# SLA2 allows for calculations that are off-axis as well as different axis cross sectionals. It calculates correctly, but needs work
# before those values could be easily utilized.
# =============================================================================
#SLA2 = Slower_Loop3d(wire_diameter)
b= np.zeros_like(Plot_Range)

Target_Field = (ZII.Zeeman_Ideal(field[2], ZII.eta, ZII.delta_o) -  B_Bias) 
Target_Tesla = [(1/10000)*tf for tf in Target_Field]




Ideal_Field = ZII.Zeeman_Ideal(field[2], ZII.eta, ZII.delta_o)
Ideal_Inc = ZII.Zeeman_Ideal_Increasing(field[2], ZII.eta , -800e6)
Inc_tar = Ideal_Inc - B_Bias

#Srfield = field[2][6000:]
#Srfield2 = [(x-0.6) for x in Srfield]
#Srplotfld = np.linspace(0.55, 0.73, 1801)
#Target_Sr = ZIISR.Zeeman_Ideal(fieldSR[2], ZIISR.eta, ZIISR.delta_o) - B_BiasSR
#Ideal_Sr = ZIISR.Zeeman_Ideal(fieldSR[2], ZIISR.eta, ZIISR.delta_o)

Rb_Adiabatic_limit = ZII.Rb_adiabatic_limit
Rb_Adiabatic2 = ZII.Rb_adiabatic2
#os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")

#Sr_Adiabatic_limit= ZIISR.Sr_adiabatic_limit
SLAPLOT = Slower_Loop_Array(wire_diameter, Slower_radius)

fld_a= 24001   
fld_b = 1001
fld = Loop_Length[-fld_a:-fld_b]
wire_diam_steps = 15
shiftby = 50
Shifted20 = SLAPLOT[19][(-fld_a-(shiftby*wire_diam_steps)):(-fld_b-(shiftby*wire_diam_steps))]

#p.figure()
#p.title('Magnetic Field Contribution for Varying Thickness', fontsize = 18)
#p.xlim(-0.55, 0.55)
#p.ylabel('Magnetic Field (B)', fontsize = 16)
#p.xlabel('Distance from wire loop (m)', fontsize = 16)
##p.plot(dbfld, Bfieldgradient, ':', color='blue', label='High-B Bias')
##p.plot(GLA.Plot_Range, rdaddwb, label = 'Rb_with_bias')
##p.plot(GLA.Slower_Length, GLA.Rb_Adiabatic_limit)
#p.plot(Loop_Length, SLAPLOT[0], label = '1 Layer')
#p.plot(Loop_Length, SLAPLOT[1], label = '2 Layers')
#p.plot(Loop_Length, SLAPLOT[4], label = '5 Layer')
#p.plot(Loop_Length, SLAPLOT[19], label = '20 Layer')
#p.plot(fld,Shifted20, ':' , label='20 Layer - shifted')
#p.plot(Loop_Length, SLAPLOT[-1], label = '50 Layer')
#p.legend(loc='best')
#p.grid(which='both')
#import os
##os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
##p.savefig('Rb_AdiabaticityDAMOP.jpeg', dpi=800, bbox_inches='tight', pad_inches=0.1)
##p.savefig('Rb_Adiabaticbreak.pdf', dpi=1200, bbox_inches='tight', pad_inches=0.1)
##p.savefig('GLAPLOT.jpeg', dpi=1200, bbox_inches='tight', pad_inches=0.1)
#p.show()