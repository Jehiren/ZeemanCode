import numpy as np
from matplotlib import pyplot as p

#p.close('all')

# =============================================================================
# #This is a collection of the constants for testing purposes on how importing works
# These constants are for the Ideal zeeman field calc for Rubidium
# =============================================================================

v_o = 350.0#oven output in m/s
v_f = 18.5#capture velocity v_f = sqrt(2*mu_b*B/M) where B is the trap height
#a_max = 114831 #max accelleration determined by atomic properties
a_max = 113154
T = 356.24 # Temperature for oven output of 350 m/s
kb = 1.3806e-23  #Boltzman thermal constant
#gamma = 2.0*np.pi*6.0666e6 #linewidth of Rubidium
gamma = 2.0*np.pi*5.98e6 #linewidth of Rubidium
M_rb = 1.40999e-25 #Atomic mass of Rb in kg
I_sat = 1.669 #Saturation Intensity mW/cm^2
I_laser = 50*I_sat #laser intensity mW/cm^2, corresponds to 83mw/cm^2
s_o = I_laser/I_sat #saturation parameter
wavelength = 780.24e-9 #wavelength of slower light
gamma_0 = 2*np.pi*5.98e6 #linewidth of Rubidium
L = 0.69 #This value is rounding the above calcuation to a cleaner solution to help maintain adiabaticity.
#eta = (v_o**2 -v_f**2)/(2*L*a_max) #This eta is caluclated based on the Length
eta = 0.8
#eta =0.792
# L = (v_o**2 -v_f**2)/(2*eta*a_max) #Slower Length in m based on eta.
# delta = 0.5*gamma*np.sqrt(1.0+s_o)*np.sqrt((1.0/eta)-1) #total detuning
kb = 1.380649e-23
k = 2*np.pi/(wavelength) #wave number, 2pi/wavelength


a_o = eta*a_max #This is the acceleration we are designing our system to have such that it is less than the maximal allowed.
L_min = eta*L  #minimum length which can achieve adiabacity


#Adjust this detuning to vary the maximum magnetic field. Always Positive.
#More Positive means higher B_Max

delta_o = 150e6 #zero feild detuning between laser and atomic resonance. As delta_o decreases, initial slowed velocity goes UP. 


## Some computational constants to reduce rounding issues with 10^34 factors.
hbar = 1.0546e-34
hh = 2*np.pi*hbar
mu_b = 9.274e-24


c_d = 7.145e-7 #The constant h/mu_B in units of Gauss/Hz
cgamma = mu_b/(10000*hh) #"Zeeman constant" of muB/hh in terms of Hz/Gauss ~1.4 MHz/G
#cg2 = 2*np.pi*cgamma  #Zeeman constant mu_B/hbar


#Ideal for sigma+ (decreasing) field
def Zeeman_Ideal(sfield, zeta, detune):
    Bt = (v_o - v_f)/wavelength
    Bb = v_f/wavelength + detune
    #Need to round and add 0.01 m to prevent issues with discontinutities based on sqrt(1-z/L)
    Lt = 0.69
#    Lt = round((v_o**2 -v_f**2)/(2*zeta*a_max), 2)
    Btemp = (c_d)*(Bb + Bt*np.sqrt(1.-sfield/Lt))
    return Btemp
#Ideal for sigma- (increasing) field    
def Zeeman_Ideal_Increasing(sfield, zeta, detune):
    Bt = (v_o - v_f)/wavelength
    Bb = -(v_f/wavelength + detune)
    #Need to round and add 0.01 m to prevent issues with discontinutities based on sqrt(1-z/L)
    Lt = round((v_o**2 -v_f**2)/(2*zeta*a_max), 2) + .0001
    Btemp = (c_d)*(Bb - Bt*np.sqrt(1.-0.95*sfield/Lt))
    return Btemp

#Function that checks adiabaticity for an input range
def adiabatic_check(field):
    return c_d*a_max/(wavelength*np.sqrt(v_o**2 -2.0*a_o*field))

#adiabatic limit based on velocity    
def adiabatic_limit(va):
    return (1/10000)*c_d*a_max/(wavelength*va)

        
ifield = np.linspace(0,0.68, 6801)
#zz = Zeeman_Ideal(ifield)
zfld = Zeeman_Ideal(ifield, eta, delta_o)
adiabaticity = adiabatic_check(ifield)
#zzdec = Zeeman_Ideal_Decreasing(ifield)
# 
def Zeeman_B_field(x,zeta, detune):
    Bt = (v_o - v_f)/wavelength
    Bb = v_f/wavelength + detune
    #Need to round and add 0.01 m to prevent issues with discontinutities based on sqrt(1-z/L)
    Lt = 0.69
#    print(Lt)
    Btemp = (c_d)*(Bb + Bt*np.sqrt(1.-x/Lt))
    return Btemp    

def Zeeman_Ideal_Perfect(sfield, detune):
    Bt = (v_o - v_f)/wavelength
    Bb = v_f/wavelength + detune
    #Need to round and add 0.01 m to prevent issues with discontinutities based on sqrt(1-z/L)
    Lt = round((v_o**2 -v_f**2)/(2*a_max) + 0.01, 2)+.001
    Btemp = (c_d)*(Bb + Bt*np.sqrt(1.-0.95*sfield/Lt))
    return Btemp

#Zeeman2 = Zeeman_Ideal_Perfect(ifield)

##In case I want to check or plot the velocity
## This function lets me vary the eta parameter and provides the velocity curve for it 
    
def v_slow_max(eta_v):
    return np.sqrt(v_o**2 -2.0*eta_v*a_max*ifield)

def v_slow(eta_vs, vfield, vdetune):
    return wavelength*(Zeeman_Ideal(vfield, eta_vs, vdetune)*cgamma - vdetune)
#    return np.sqrt(v_o**2 -2.0*eta_vs*a_max*vfield)    

def v_therm(velspace, velin):
    vv = velin/np.sqrt(3)
    v_eff = velspace**3/(2*vv**4)*np.exp(-0.5*(velspace/vv)**2)
    return v_eff

v_slower = np.sqrt(v_o**2 - 2.0*a_o*ifield)
v_slower3 = v_slow(eta, ifield, delta_o)
Rb_adiabatic_limit = adiabatic_limit(v_slower)
v_perf = v_slow(1, ifield, delta_o)
v_slower2 = np.concatenate((np.array(500*[v_o]),v_slower , np.array(2000*[20])))
#v_slower2 = 500*[v_o] + v_slower

Rb_adiabatic_limit = adiabatic_limit(v_slower)

Rb_adiabatic2 = adiabatic_limit(v_slower2)

#p.close('all')
#p.figure()

#etalist = np.linspace(0.7, 0.8, 6)
#for etas in etalist:
##    fieldtemp = round((v_o**2 -v_f**2)/(2*etas*a_max) + 0.01, 2)
#    fieldtemp = (v_o**2 -v_f**2)/(2*etas*a_max)
#    fldstp = int(fieldtemp*10000 +1)
#    fldtmp = np.linspace(0,fieldtemp, fldstp )
#    vtemp = v_slow(etas, fldtmp)
#    p.plot(fldtmp, vtemp, label=etas)
#p.legend(loc='best')
#p.show()
#
#vspace1 = np.linspace(0, 600, 6000)
#vspace2 = np.linspace(0, 100, 500)
#fvmax = v_therm(vspace1, v_o)
#fvmin = v_therm(vspace2, v_f)
#fvvmax = max(fvmin)
#p.figure()
#p.title('Effusive Boltzman Velocity Distribution', fontsize = 18)
#p.plot(vspace1, fvmax/fvvmax, label = 'v_ave = 350 m/s')
#p.plot(vspace2, fvmin/fvvmax, label = 'v_ave = 20 m/s')
#p.xlabel('Average Velocity (m/s)', fontsize = 16)
#p.ylim(0, 1)
#p.ylabel('Normalized Boltzman ratio f(v)/fmax', fontsize = 16)
#p.legend(loc='best', fontsize = 12)
#p.show()
#import os
#os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
#p.savefig('ThermalVelBoth.png', bbox_inches='tight', dpi=600)

# =============================================================================
# Plot stuff below.
# =============================================================================
#zperf = Zeeman_Ideal_Perfect(ifield, delta_o)
##
#import os
#os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
#p.figure()
#p.title('Magnetic Field Ideal and broadened Rb Zeeman Slower', fontsize = 18)
#p.plot(ifield, zperf, label = 'Pefectly Ideal')
#p.plot(ifield,zfld, label= r' Broadened Field, $\eta$ = ' + str(eta))
#p.ylabel('Magnetic Field Strength (G)', fontsize = 16)
#p.xlabel('Distance along designed slower length (m)',fontsize = 16)
#p.legend(loc='best')
#p.show()
#p.savefig('IdealRbZeeman.pdf', bbox_inches='tight', dpi=1200)
#

##adbfld2 = np.linspace(-0.05, 0.88, 9301)
##
#p.figure()
## # p.plot(ifield, v_squared, '.')
#p.title('Rb Resonant Velocity in Broadened Zeeman Slower', fontsize = 18)
#p.ylabel('Rb Velocity (m/s)' , fontsize = 16)
#p.xlabel('Distance along slower (m)', fontsize = 16)
#p.plot(ifield, v_perf, label = 'Perfectly Ideal Velocity')
#p.plot(ifield, v_slower, label=r' Broadened Velocity, $\eta$ = ' + str(eta))
##p.plot(ifield, Rb_adiabatic_limit, label= ' Adiabaticity')
##p.plot(adbfld2, Rb_adiabatic2, label='Extended Adiabaticity')
##p.plot(ifield, Rbadbmax)
##p.plot(xa, va, label='Atom-Frame Numerical Integration')
##p.plot(xa,va2, label='Lab-Frame Numeri`cal Integration')
#
##p.legend(loc='best')
#p.show()
#
####
#import os
#os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
#p.savefig('Velocityfin.png', bbox_inches='tight', dpi=800)
## 

# =============================================================================
# Plotting increasing field and decreasing field solution
# =============================================================================

#Zeemaninc = Zeeman_Ideal_Increasing(ifield, eta, -800e6)
#
#p.figure()
#p.plot(ifield, Zeemaninc)
#p.show()

#def magnitude(x):
#    return(int(np.log10(x)))
#ro = 0.01
#db = .25
#da = .47
#Temps = [300, 320, 340, 360 , 380, 400, 420, 440, 460]
#Ninit = 1e15
#Lengths = np.linspace(0.5, 1.25, 200)
#v_o2 = [320, 330, 340, 350, 360, 370 ,380]
#p.figure()
#for temps in Temps:
#    vzmax = np.sqrt(v_f**2 + 2*Lengths*a_o)    
#    vrmax = ro/( (db/vzmax) + (vzmax - v_f)/a_o + da/v_f)
#    
#    alpha = np.sqrt(kb*temps/M_rb)
#    Nout = Ninit*(1-np.exp(-vzmax**2/2*alpha**2))*(1-np.exp(-vrmax**2/alpha**2))
#    p.plot(Lengths, Nout, label=temps)
#    p.show()









