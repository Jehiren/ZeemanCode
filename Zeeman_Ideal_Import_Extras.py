import numpy as np
from matplotlib import pyplot as p

#This is a collection of the constants for testing purposes on how importing works
#These constants are for the Ideal zeeman field calc for Rubidium
#p.close('all')

v_o = 350.0#oven output in m/s
v_f = 20.0#capture velocity v_f = sqrt(2*mu_b*B/M) where B is the trap height
a_max = 114831 #max accelleration determined by atomic properties

T = 356.24 # Temperature for oven output of 350 m/s
kb = 1.3806e-23  #Boltzman thermal constant
gamma = 2.0*np.pi*5.98e6 #linewidth of Rubidium
M_rb = 1.41e-25 #Atomic mass of Rb in kg
I_sat = 1.66 #Saturation Intensity mW/cm^2
I_laser = 50*I_sat #laser intensity mW/cm^2, corresponds to 83mw/cm^2
s_o = I_laser/I_sat #saturation parameter
wavelength = 780.24e-9 #wavelength of slower light

# eta = 0.8#design parameter to ensure adiabacity
L = 0.68 #This value is rounding the above calcuation to a cleaner solution to help maintain adiabaticity.
#eta = (v_o**2 -v_f**2)/(2*L*a_max) 
eta =0.782
# L = (v_o**2 -v_f**2)/(2*eta*a_max) #Slower Length in m
# delta = 0.5*gamma*np.sqrt(1.0+s_o)*np.sqrt((1.0/eta)-1) #total detuning
#delta = 0


a_o = eta*a_max #This is the acceleration we are designing our system to have such that it is less than the maximal allowed.
## Change the sign here for increasing or decreasing field slower.
## Plus is for increasing field, minus is decreasing. Comes from mu_B for sigma+/- transition.


L_min = eta*L  #minimum length which can achieve adiabacity

#Adjust this detuning to vary the maximum magnetic field. ???Always negative.
#More negative means higher B_Max

delta_o = 150e6 #zero feild detuning between laser and atomic resonance
##
hbar = 1.0546e-34
hh = 2*np.pi*hbar
mu_b = 9.274e-24
k = 2*np.pi/wavelength

c_d = 7.145e-7 #The constant h/mu_B in units of Gauss/Hz
cgamma = mu_b/(10000*hh) #"Zeeman constant" in terms of Hz/Gauss
cgamma2 = mu_b/(10000*hbar) #"Zeeman constant" in terms of rad/Gauss*s

def Zeeman_Ideal(sfield, zeta, detune):
    Bt = (v_o - v_f)/wavelength
    Bb = v_f/wavelength + detune
    #Need to round and add 0.01 m to prevent issues with discontinutities based on sqrt(1-z/L)
    Lt = round((v_o**2 -v_f**2)/(2*zeta*a_max), 2) + .0001
    Btemp = (c_d)*(Bb + Bt*np.sqrt(1.-0.95*sfield/Lt))
    return Btemp
    
def Zeeman_Ideal_Increasing(sfield, zeta, detune):
    Bt = (v_o - v_f)/wavelength
    Bb = -(v_f/wavelength + detune)
    #Need to round and add 0.01 m to prevent issues with discontinutities based on sqrt(1-z/L)
    Lt = round((v_o**2 -v_f**2)/(2*zeta*a_max), 2) + .0001
    Btemp = (c_d)*(Bb - Bt*np.sqrt(1.-0.95*sfield/Lt))
    return Btemp


def adiabatic_check(field):
    return c_d*a_max/(wavelength*np.sqrt(v_o**2 -2.0*a_o*field))
    
def adiabatic_limit(va):
    return (1/10000)*c_d*a_max/(wavelength*va)

        
ifield = np.linspace(0,0.68, 6801)
#zz = Zeeman_Ideal(ifield)
zfld = Zeeman_Ideal(ifield, eta, delta_o)
adiabaticity = adiabatic_check(ifield)
#zzdec = Zeeman_Ideal_Decreasing(ifield)
# 


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


# =============================================================================
# I need to numerical integrate the accel to obtain the velocity profile and ensure it looks right. ALL OF THIS HAS BEEN MOVED TO PLOT_SOLUTIONS NOW
# =============================================================================

# =============================================================================
# Numerically integrating the acceleration to get the position using 2nd order Runge-Kutta method for a velocity dependent force
# =============================================================================
#
#
#ztune = 150e6
#
#def F_atom():
#    zf = Zeeman_Ideal(ifield, eta, ztune)
#    zv = v_slower
#    zd = (1/c_d)*zf + zv/wavelength
#    F = (0.5*hbar*k*gamma)*(s_o/(1 + s_o +4*(zd/gamma)**2))
#    return F
#
#
#def func_rk2(func,q,h):
#      k1 = h*func(q)
#      k2 = h*func(q+k1/2)
#      return q+k2
#
#def Zeeman_B_field(x,zeta, detune):
#    Bt = (v_o - v_f)/wavelength
#    Bb = v_f/wavelength + detune
#    #Need to round and add 0.01 m to prevent issues with discontinutities based on sqrt(1-z/L)
#    Lt = round((v_o**2 -v_f**2)/(2*zeta*a_max) + 0.01, 2)
##    print(Lt)
#    Btemp = (c_d)*(Bb + Bt*np.sqrt(1.-x/Lt))
#    return Btemp    
#
#zzeta = eta
#
#
#
##def Bzz(pos):
##    return Zeeman_B_field(pos, zzeta, zzdetune)
#
##def a_acel(apos, avel, dh): 
##    zdt = (1/c_d)*Bzz(apos) - avel/wavelength
##    aacel = -(0.5*hbar*k*gamma/M_rb)*(s_o/(1 + s_o +4*(zdt/gamma)**2))
###    print(zdt, aacel)
##    return aacel
#
##def a_acel(apos, avel, dh): 
###    zdt = (1/c_d)*Bzz(apos) - avel/wavelength
##    zdt = 80e6
##    aacel = -(0.5*hbar*k*gamma/M_rb)*(s_o/(1 + s_o +4*(zdt/gamma)**2))
###    print(zdt, aacel)
##    return aacel
#
#
##This one somehow worked but the units don't match. Crap.
##def a_acel2(apos, avel, dh): 
##    zdt = zzdetune - cgamma*Bzz(apos) + avel/wavelength
##
##    aacel = -gamma*s_o*cgamma**2*hbar*k/((8*(zdt**2) +2*gamma + 2*s_o*cgamma**2)*M_rb)
##    return aacel
#
##This works once zdt is angular. There is still a strange factor of 2pi though...
#
#
##def a_acelBarret(apos, avel, dh): 
##    zdt = zzdetune - cgamma*Bzz(apos) + avel/wavelength
##
##    aacel = -gamma*s_o*hh/((4*(zdt/gamma)**2 +1 + s_o)*2*M_rb*wavelength)
##    return aacel
#
#def rk4(x, v, a, dt):
#    """Returns final (position, velocity) tuple after
#    time dt has passed.
#
#    x: initial position (number-like object)
#    v: initial velocity (number-like object)
#    a: acceleration function a(x,v,dt) (must be callable)
#    dt: timestep (number)"""
#    x1 = x
#    v1 = v
#    a1 = a(x1, v1, 0)
#
#    x2 = x + 0.5*v1*dt
#    v2 = v + 0.5*a1*dt
#    a2 = a(x2, v2, dt/2.0)
#
#    x3 = x + 0.5*v2*dt
#    v3 = v + 0.5*a2*dt
#    a3 = a(x3, v3, dt/2.0)
#
#    x4 = x + v3*dt
#    v4 = v + a3*dt
#    a4 = a(x4, v4, dt)
#
#    xf = x + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
#    vf = v + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)
#
#    return xf, vf
#
#stepsize = 0.00001
#
#
#
#zzdetune = 120e6
#
#def Bzz(pos):
#    return Zeeman_B_field(pos, zzeta, ztune)
#
#def a_acel3(apos, avel, dh): 
#    zdt = 2*np.pi*zzdetune - cgamma2*Bzz(apos) + avel*k
#
#    aacel = -gamma*s_o*hbar*k/((4*(zdt/gamma)**2 +1 + s_o)*2*M_rb)
#    return aacel
#
##This one somehow worked but the units don't match. Crap.
##def a_acel2(apos, avel, dh): 
##    zdt = zzdetune - cgamma*Bzz(apos) + avel/wavelength
##
##    aacel = -gamma*s_o*cgamma**2*hbar*k/((8*(zdt**2) +2*gamma + 2*s_o*cgamma**2)*M_rb)
##    return aacel
#
#
#va = [v_o]
#vf = 330
#xa = [0.0]
#xf = 0
#n = 0
#
#while vf > 15:
#    tempy = rk4(xa[n], va[n], a_acel3, stepsize)
#    va.append(tempy[1])
#    xa.append(tempy[0])
#
#    vf = va[n]
#    n = n +1
#
#
#va2 = []
#for xn in range(len(xa)):
#    va2.append(v_o*np.sqrt(1-xa[xn]/L))
#
##p.figure()
### # p.plot(ifield, v_squared, '.')
##p.title('Rb Velocity in the Slower')
##p.ylabel('Rb Velocity (m/s)')
##p.xlabel('Distance along slower')
##p.plot(ifield, v_slower, label='Original Method')
##p.plot(ifield, v3, label='gamma velocity')
##p.plot(xa, va, label='Atom-Frame Numerical Integration')
#p.plot(xa, va, label='Detuning from Resonance (MHz)  - ' + str(zzdetune/1e6))
##p.plot(xa,va2, label='Lab-Frame Numerical Integration')
##p.plot(ifield, vwhynot, label='Surewtf')
#p.legend(loc='best')
#p.show()
#
#
#
#
##This loop integrates for position, useful if integration works out well.
##while xf < 0.6875:
##    tempy = rk4(xa[n], va[n], a_acel3, stepsize)
##    va.append(tempy[1])
##    xa.append(tempy[0])
##
###    print(xa[n], va[n])
###    vf = va[n]
##    xf = xa[n]
##    n = n +1
#
##This loop integrates for a velocity, useful if things aren't working right. (never reaches end of slower...)
##while vf > 20:
##    tempy = rk4(xa[n], va[n], a_acel3, stepsize)
##    va.append(tempy[1])
##    xa.append(tempy[0])
##
##    vf = va[n]
##    n = n +1
#
##This iteration plots the function while it iterates over the etas in etalist
##p.figure()
##for etavals in etalist:
##    zzdetune = etavals
##    print(zzdetune)
##    va = [v_o]
##    vf = 330
##    xa = [0.0]
##    xf = 0
##    n = 0
##    while vf > 20:
##        tempy = rk4(xa[n], va[n], a_acel3, stepsize)
##        va.append(tempy[1])
##        xa.append(tempy[0])
##        vf = va[n]
##        n = n +1
##    p.plot(xa, va, label = etavals)
#
#
#
#
#
## =============================================================================
## This loop takes the numerically integrated solution for position and plugs it into V_R = v_o*sqrt(1 - z_R/L), where z_R is 
## the atomic position in this decelerating frame.
## =============================================================================
##va2 = []
##for xn in range(len(xa)):
##    va2.append(v_o*np.sqrt(1-xa[xn]/L))
##
###p.figure()
#### # p.plot(ifield, v_squared, '.')
###p.title('Rb Velocity in the Slower')
###p.ylabel('Rb Velocity (m/s)')
###p.xlabel('Distance along slower')
###p.plot(ifield, v_slower, label='Original Method')
###p.plot(ifield, v3, label='gamma velocity')
###p.plot(xa, va, label='Atom-Frame Numerical Integration')
##p.plot(xa, va, label='Detuning from Resonance (MHz)  - ' + str(zzdetune/1e6))
###p.plot(xa,va2, label='Lab-Frame Numerical Integration')
###p.plot(ifield, vwhynot, label='Surewtf')
##p.legend(loc='best')
##p.show()
###
###import os
###os.chdir(r"F:/Zeeman Plots and Diagrams/Nov13")
###p.savefig('Velocity.png', bbox_inches='tight', dpi=1200)
## 
##
### =============================================================================
### Plot the ideal broadened slower and perfect ideal slower
### =============================================================================
##
##p.figure()
##p.title('Ideal Zeeman Magnetic Field Profile for Rubidium')
###p.plot(ifield, Zeeman2, ':' , label = r'$\eta$ = 1.0' )
##p.plot(ifield, zzdec, '-', label = r'$\eta$ = 0.77')    
##p.xlabel('Distance along slower (m)')
##p.ylabel('Magnetic Field B (G)')
##p.legend(loc='best')
##p.show()
###import os
###os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
###p.savefig('IdealRbZeeman.pdf', bbox_inches='tight', dpi=1200)
#
