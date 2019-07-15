import numpy as np


# =============================================================================
# Set the minimum distance to calculate for the magnetic field to prevent asypmtotes.
# =============================================================================

mindist = 1e-7

# =============================================================================
# ###Magnetic field calculation using the Biot Savart Law. The mu factor is scaled
# ### to prevent issues with small values.
# =============================================================================

mu_scaled=4.0*np.pi

# =============================================================================
# Define a function to calculate the distance between the source and field points, using mindist to prevent asymptotes
# =============================================================================

def Bdist(r_field, r_loop_offset, phi, a):
    Tempd = (((a*a)+((r_field[0]-r_loop_offset[0])**2.0)+
            ((r_field[1]-r_loop_offset[1])**2.0)+
            ((r_field[2]-r_loop_offset[2])**2.0)-
            2.0*a*((r_loop_offset[1]-r_field[1])*np.sin(phi)+
            (r_loop_offset[0]-r_field[0])*np.cos(phi)))**(3./2.)) 
            
    np.place(Tempd,Tempd<mindist,mindist)
    # if tempd<mindist:   
    #     tempd==mindist   
    return Tempd
    
# =============================================================================
# The three components of the magnetic field. Calculated aligning the loop to the xy-plane, but allowing for any offset.
# We are calculating the differential, dB form of the Biot-Savart Law. This means we sum over discrete number of current elements (n)
# =============================================================================
    
def B_x(r_field, r_loop_offset,phi,a):
    return ((r_field[2]-r_loop_offset[2])*a*np.cos(phi))
    
def B_y(r_field, r_loop_offset,phi,a):
    return (a*(r_field[2]-r_loop_offset[2])*np.sin(phi))

def B_z(r_field, r_loop_offset,phi,a):
    return ((a**2.0)-a*((r_field[1]-r_loop_offset[1])*np.sin(phi)+
            (r_field[0]-r_loop_offset[0])*np.cos(phi)))

n_loop = 20 ## This is the number of pieces the loop is broken into to evaluate

# =============================================================================
# Combine the components of the dB field into a single function that ignores the non-axial components for computational simplicity.    
# =============================================================================

def axialmag(r_field,r_loop_offset,a):
    
    dp = 0.0
    Bz_temp = 0.0
    dphi = 2.0*np.pi/n_loop
    while n_loop>dp:
        fielddist=Bdist(r_field,r_loop_offset,dp*(dphi),a)
        Bz_temp = Bz_temp + dphi*B_z(r_field,r_loop_offset,dp*dphi,a)/fielddist
        dp=dp+1
    return Bz_temp    
  
# =============================================================================
# Compute all three components of the magnetic field by summing over n discrete "pieces" of current.
# =============================================================================
  
def BiotSavart_Ring(r_field,r_loop_offset,a):
    
    dp = 0.0
    Bz_temp = 0.0
    By_temp = 0.0
    Bx_temp = 0.0
    dphi = 2.0*np.pi/n_loop
    while n_loop>dp:
        fielddist=Bdist(r_field,r_loop_offset,dp*(dphi),a)
        Bz_temp = Bz_temp + dphi*B_z(r_field,r_loop_offset,dp*dphi,a)/fielddist
        Bx_temp = Bx_temp + dphi*B_x(r_field,r_loop_offset,dp*dphi,a)/fielddist
        By_temp = By_temp + dphi*B_y(r_field,r_loop_offset,dp*dphi,a)/fielddist
        dp=dp+1
    return [Bx_temp, By_temp, Bz_temp]       
    