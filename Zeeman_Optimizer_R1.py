import os,csv
import numpy as np
from matplotlib import pyplot as p
#import pdb
from lmfit import Parameters, Model
#from lmfit.model import save_modelresult

#Change this to be whatever directory you place all your files in.
os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/Thesis Version (NO TOUCHING)/Main Program (Clean Edit)/")
import Generate_Loop_Array as GLA



p.close('all')



        
 
# =============================================================================
# This section should contain the only parameters that need to be changed when operating this code. 
# =============================================================================
CP_Count = 20

wire_diam_steps = int(GLA.wire_diameter*10000)

cp_initial_offset = -35 #This is how many coil windings will be placed to the left of the "ideal" slower length to get the B-field where it needs to be.     

# =============================================================================
#  Set limits for total slower length and min/max limits for individual component solenoid lengths.
# =============================================================================
min_coil_len = 6900 + 0.5*abs(cp_initial_offset)*wire_diam_steps

cp_lower_len = int(min_coil_len/(wire_diam_steps*(1.1*CP_Count)))
cp_upper_len = int(min_coil_len/(wire_diam_steps*(0.6*CP_Count)))

cp_current_lim = 6 #Maximum current in amps
Rb_cur_min = 0 #Minimum current in amps
Sr_cur_min = -5.0
cp_depth_max = 30  #Maximium number of loops deep.
cp_depth_min = 1   #Want at least one wire in each location.
BL_min = min_coil_len/wire_diam_steps - 0.75*cp_lower_len
BL_max =BL_min + 1.505*cp_lower_len

# =============================================================================
# End of parameter input
# =============================================================================
# =============================================================================
# Creating the Slower Loop Array (SLA) for the wire diam and Slower_radius listed above.
# =============================================================================

SLA = GLA.Slower_Loop_Array(GLA.wire_diameter, GLA.Slower_radius)

# =============================================================================
# Make a loop that will create coil packs and a current that has the A + B*sqrt(1-n/s) where n is the nth solenoid out of s total.
# =============================================================================   
CP_Length = [int(x) for x in (CP_Count*[2.0*min_coil_len/((2*CP_Count-1)*wire_diam_steps)])]
CP_Depth = []
CP_Current = []
for cpn in range(CP_Count):
    CP_Depth.append(int(3 + 14.5*(1-cpn/CP_Count)**0.875))
#    CP_Current.append(-4.1 + 7.35*np.sqrt(1-0.45*cpn/CP_Count))
    CP_Current.append(3.25)
    

# =============================================================================
# If you are doing Current Only Optimization, you need to uncomment this following section and comment out the above section.
# This section pulls a previous solution for the shape parameters.
# =============================================================================
#    
#Solenoid_Dep = []
#Solenoid_Len = []
#
#with open('SigmaPlusSolution.csv', mode='r', newline='') as readme:
#    reader = csv.reader(readme, delimiter=',', quoting=csv.QUOTE_NONE)
#    itemp = next(reader)
#    wire_diam = float(float(itemp[2]))
#    Slower_rad = float(float(itemp[4]))
#    cp_initial_offset = float(float(itemp[6]))
#    next(reader)
#    for row in reader:
#        Solenoid_Dep.append(int(float(row[1])))
#        Solenoid_Len.append(int(float(row[3])))        
#
#CP_Length = Solenoid_Len
#CP_Depth = Solenoid_Dep      
#CP_Current = []
#for cpn in range(CP_Count):
#    CP_Depth.append(int(5 + 11*np.sqrt(1-cpn/CP_Count)))
#    CP_Current.append(-4.1 + 7.35*np.sqrt(1-0.45*cpn/CP_Count))
#  
# =============================================================================
# Use the wire diamter input to determine how many steps to move per coil loop. The scale is in 0.1 mm per step, so 3.8 mm is 38 steps
# =============================================================================
    
Prefixes = []
for ip in range(1,CP_Count):
    Prefixes.append('B%d_' %(ip))

Prefixes.append('BLast_')
paramnames = ['Coil_Depth', 'Coil_Current', 'Coil_Length', 'offsetme']

print(cp_initial_offset)

############################
# =============================================================================
# 
# Now we will combine the Slower_Loop_Array with the parameters imported to make a B-field plot.
# This function takes the parameter value for # of wrappings and pulls the corresponding field, shifting it to where it needs to be.
# The 17500 and +10000 are values that are based off the difference between the slower field range and the Biot-Savart Calc loop. The 38 corresponds to how many step sizes
#  each wire diameter corresponds to. Our step sizes are consistently 0.1 mm throughout, so for 3.8 mm wire, it's 38 steps.
#  It is important to ensure that the original field for SLA is long enough for the function shifting.
# 
# 
# The index multiple times 38 is based on the wire diamter and is the number of step sizes corresponding to that diameter, with each step size being 0.1 mm.
# 
# ##This function takes input values of [How many coils windings at this location], [how much current through the coil at each location], [how many coils long is the slower]
# =============================================================================

def offsetgen(cplenin):
    tempoffset = []
    tempy = cplenin[0] 
    for cplen in cplenin:
        tempoffset.append(tempy)
        tempy = tempy + cplen
    return tempoffset


# =============================================================================
# Define slicing values for different total lengths of Loop lengths and field ranges.
# =============================================================================

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

fld = GLA.Loop_Length[-fld_a:-fld_b]

srcp_a = 11000
srcp_b = -11199


# =============================================================================
# Slower_CoilPacks_Short generates the magnetic field for a single solenoid based on length, depth, current, and offset. Summed up in other functions.
# =============================================================================

def Slower_CoilPacks_Short(scpCoil_Depth, scpCoil_Current, scpCoil_Length, scpCoil_Offset):
    BCP_Totes = np.zeros_like(fld)    
    Position = np.arange((int(float(scpCoil_Offset))), int(float(scpCoil_Offset+scpCoil_Length)) ,1)
    for pw in range(int(scpCoil_Length)):
        BCP_Totes =  BCP_Totes + scpCoil_Current*(SLA[int(scpCoil_Depth)][(-fld_a-(Position[pw]*wire_diam_steps)):(-fld_b-(Position[pw]*wire_diam_steps))])
    return BCP_Totes

# =============================================================================
# CP_Total_B sums over all contributing solenoids and provides the full Loop_Length (-2m to +2m range) field. 
# =============================================================================

def CP_Total_B(Coil_Depth, Coil_Current, Coil_Length):
    B_Totals  = np.zeros_like(fld)
    offset = cp_initial_offset 
    for nn in range(CP_Count):
        B_Totals  = B_Totals  + Slower_CoilPacks_Short(Coil_Depth[nn], Coil_Current[nn], Coil_Length[nn], offset)
        offset = offset + Coil_Length[nn]       
    del offset, Coil_Depth, Coil_Current, Coil_Length, nn
    return B_Totals   

# =============================================================================
# CP_Total_B_Trimmed does the same as CP_Total_B except that it's trimmed down to the exact length of the slower and no more.
# For Rb, this corresponds to being 0.69m long, so 6900 steps.
# =============================================================================

def CP_Total_B_Trimmed(x, Coil_Depth, Coil_Current, Coil_Length, offsetme):
    B_Toteme = np.zeros_like(fld)
    B_Toteme = Slower_CoilPacks_Short(Coil_Depth, Coil_Current, Coil_Length, offsetme)
    del offsetme, Coil_Depth, Coil_Current, Coil_Length
    return B_Toteme[cp_a: cp_b]

def CP_Total_B_TrimmedSR(x, Coil_Depth, Coil_Current, Coil_Length, offsetme):
    B_Toteme = np.zeros_like(fld)
    B_Toteme = Slower_CoilPacks_Short(Coil_Depth, Coil_Current, Coil_Length, offsetme)
    del offsetme, Coil_Depth, Coil_Current, Coil_Length
    return B_Toteme[srcp_a: srcp_b]

# =============================================================================
# These functions below create the CombinedModel that will add several Coil Packs together.
#  They are done in such a way that it allows me to automate the addition of how many coil packs I have.
# =============================================================================
    
def Make_Model(Modeltype, prefix):
    TempModel = Model(Modeltype, prefix = prefix) 
    return TempModel  

# =============================================================================
# Now we create the Combomodel for a single solenoid component and then loop over the number we want based on CP_Count. 
# =============================================================================
    
Combomod = Make_Model(CP_Total_B_Trimmed, Prefixes[0])
for ii in range(1, CP_Count):
    Combomod = Combomod + Make_Model(CP_Total_B_Trimmed, Prefixes[ii])    

params = Parameters()

# =============================================================================
# This creates the combomodel used for Strontium
# =============================================================================

#Combomod = Make_Model(CP_Total_B_TrimmedSR, Prefixes[0])
#for ii in range(1, CP_Count):
#    Combomod = Combomod + Make_Model(CP_Total_B_TrimmedSR, Prefixes[ii])    
#  
#params = Parameters()

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Set the Parameters that will be varied as well as setting the initial vaues for each of them. If you don't set initial values
# then you will have inf as the min, max, and value by default.
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# =============================================================================
# Set the parameters for the first coil, which is referred to as the "input coil". B1_offsetme having the expression definition ensures that
# we have no issues with the offset not matching the length of the first coil (which happend a few times...). This coil gets us up to the max field quickly.
# =============================================================================

params.add('B1_Coil_Length', value = cp_initial_offset, vary = True, min = 20, max = abs(cp_initial_offset) +10)
params.add('B1_Coil_Current', value = CP_Current[0], vary=True, min = Rb_cur_min, max = 7.5)
params.add('B1_Coil_Depth', value = CP_Depth[0], vary = True , min = 1, max = 40)
params.add('B1_offsetme', value = -20, vary=False , min = -cp_upper_len, max = -10)
params['B1_offsetme'].expr = '-B1_Coil_Length'
params.add('TotalCurrent', value = 2.9, vary = True, min = 0, max = 10)
params.add('BLast_Coil_Current', value = -0.5, vary = True, min = -5, max = 10)


# =============================================================================
# This loop generates the Length, Current, Depth, and offsetme parameters for the main coils of the slower. The first and last coil
# are separate for they have additional requirements to their parameters that were easier to limit this way. 
# =============================================================================

for cm in range(1, CP_Count-1):
    params.add('B%s_Coil_Length' %str(cm+1) , value = CP_Length[cm], vary = True, min = cp_lower_len, max = cp_upper_len)
#    params.add('B%s_Coil_Current'  % str(cm+1) , value = CP_Current[cm], vary=False, min = Rb_cur_min, max = cp_current_lim)
    params.add('B%d_Coil_Current'  %(cm+1) ,value = 2.9, min = Rb_cur_min, max = 10, expr = 'TotalCurrent')
    params.add('B%s_Coil_Depth' % str(cm+1), value = CP_Depth[cm], vary = True , min = 5, max = cp_depth_max)
    params.add('B%s_offsetme' % str(cm+1), value = 0, vary=False , min = 0)

# =============================================================================
# These two commands define the offset expression which requires the m+1 coil to have an offset value defined by the length_m + offset_m
# =============================================================================
    
for i in range(1, CP_Count-1):
    params['B%s_offsetme' % str(i+1)].expr = '(B%s_offsetme + B%s_Coil_Length)' %(str(i), str(i))

#This does that for the last coil since BLast doesn't follow the same naming sheme as the other coils.
exprstr = 'B%d_offsetme + B%d_Coil_Length' % ((CP_Count-1) , (CP_Count-1) ) 

# =============================================================================
# This last section defines the last coil, which has the additional limitations involving the total length of the slower.
# =============================================================================

Total_min = (min_coil_len + abs(0.85*cp_initial_offset*wire_diam_steps))/wire_diam_steps
Total_max = (min_coil_len + abs(2.25*cp_initial_offset*wire_diam_steps))/wire_diam_steps

params.add('BLast_Coil_Length', value = CP_Length[-1], vary = True, min = 0.25*cp_lower_len, max =  cp_upper_len)
#params.add('BLast_Coil_Current', value = CP_Current[-1], vary=False, min = Rb_cur_min, max = cp_current_lim)

params.add('BLast_Coil_Depth', value = CP_Depth[-1], vary = True , min = cp_depth_min, max = 35)
params.add('BLast_offsetme', value = 0, vary=False , min= 0, max = Total_min + cp_initial_offset)

# =============================================================================
# Limit the total length of the slower, should work for any general wire diameter/total slower length. 
# =============================================================================

params['BLast_offsetme'].expr = exprstr 

params.add('Total_Length', value = sum(CP_Length), vary = True, min = Total_min, max =Total_max )
params['BLast_Coil_Length'].expr = 'Total_Length - BLast_offsetme'

# =============================================================================
# Define a function that will change whether a parameter name will vary. This assumes a fairly uniform setup for parameter names
# including prefixes as defined in the lmfit documentation. Also create a function that modifies all values for a param name within
# the combomodel.
# =============================================================================

def Bulk_Param_Vary(paramsin, paramname, isvary):
    for pref in Prefixes:
        paramsin[(pref + str(paramname))].set(vary=isvary)
    return None   
             
def Bulk_Param_Change(paramsin, paramname, paramval):
    for fx in range(CP_Count):
        paramsin[(Prefixes[fx] + str(paramname))].set(paramval[fx])  
    return None

# =============================================================================
# The lmfit package provides answers in the form of a set of dictionary values, this function pulls them into standard lists for iteration and printing into csv
# =============================================================================

def result_grabber(resultsout):
    dictresval = resultsout.best_values    
    rescur = []
    reslen = []
    resdep = []
    resoffset = []
    for ix in range(CP_Count):
        rescur.append(float('%.6f'%(float(dictresval[(Prefixes[ix] + paramnames[1])]))))
        reslen.append(int(round(dictresval[(Prefixes[ix] + paramnames[2])]))    )
        resdep.append(int(round(dictresval[(Prefixes[ix] + paramnames[0])])))
        resoffset.append(int(round(dictresval[(Prefixes[ix] + paramnames[3])])))    
    return [resdep, rescur, reslen, resoffset]

# =============================================================================
# Set params with param['key'].set(value) for individual changes, or the bulk change functino above for modifying all the values at onces.
# =============================================================================

Rb_bias_fit = GLA.Target_Field
Rb_nobias_fit = GLA.Ideal_Field

# =============================================================================
# Starting up the optimizer routine. It's current divided into two sections:
# 1.) The program optimizes to the ideal Zeeman field (without bias). This was done as I found it to provide the best ability
# for the coil pack setup to match for the other slower designs.
# 2.) Those values for length, depth, and current are then fed into a second parameter which will serve for the other 3 models, with the difference of only allowing the current to change value.
# The solutions are then saved into timestamped folders for ease of organization.
# =============================================================================

print('Starting up the Full Optimizer now...')

## =============================================================================
## Making a copy of the optimizer with the order reversed so that it makes the "best" fit for the case of the high-B field active.
## Since that is going to be the primary use, that would be the preferred approach. If the "normal" zeeman slower 
## doesn't operate at the best efficienty, it is less important.
## =============================================================================

Rb_bias_interim_res0 = Combomod.fit(Rb_bias_fit, params,x=0 ) #, method = 'ampgo')
print('Rough Fit done')
paramsT = Rb_bias_interim_res0.params
Rb_bias_interim_res = Combomod.fit(Rb_bias_fit, paramsT, x=0 )#, method='ampgo')
params2 = Rb_bias_interim_res.params
interim_Rb_res = result_grabber(Rb_bias_interim_res)
print("First pass done")

#for ii in range(2):
#    paramsTemp = Rb_bias_interim_res.params
#    Rbbtemp = Combomod.fit(Rb_bias_fit, paramsTemp, x=0 , method='ampgo')
#    Rb_bias_interim_res = Rbbtemp
#    print(ii)
#print('Multi-loop done')
## =============================================================================
## This section takes the solution originally obtained by allowing all three components to vary and does a second pass of finer tuning of just current
## If you want to just do a current optimization for a given coil pack (say for a different background field, uncomment params2=params and comment the initial fit in the section above.)
## =============================================================================

Bulk_Param_Change(params2, 'Coil_Depth', interim_Rb_res[0])
Bulk_Param_Change(params2, 'Coil_Length', interim_Rb_res[2])
Bulk_Param_Vary(params2, 'Coil_Depth', False )
Bulk_Param_Vary(params2, 'Coil_Length', False )
params2['Total_Length'].set(vary=False)

Rb_bias_res =Combomod.fit(Rb_bias_fit, params2,x=0 )#, method = 'ampgo') 

result_Rb_bias = result_grabber(Rb_bias_res)
print('Finished with Rb_bias')
print(result_Rb_bias)   

Rb_nobias_res = Combomod.fit(Rb_nobias_fit, params2, x=0 )#, method='ampgo')
result_Rb_nobias = result_grabber(Rb_nobias_res)
print('Finished with Rb_nobias')


#  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# =============================================================================
# This section is CURRENT ONLY OPTIMIZATION FOR BOTH RUBIDIUM AND STRONTIUM
# =============================================================================
#
#print('Starting up the Current Only Optimizer now...')
#
## =============================================================================
## Turn off the ability to vary Depth, Length, and Total_Length
## =============================================================================
#Bulk_Param_Vary(params, 'Coil_Depth', False )
#Bulk_Param_Vary(params, 'Coil_Length', False )
#params['Total_Length'].set(vary=False)
#
##Rb_bias_res =Combomod.fit(Rb_bias_fit, params,x=0, method = 'ampgo') 
##
##result_Rb_bias = result_grabber(Rb_bias_res)
##print('Finished with Rb_bias')
## 
##Rb_nobias_res = Combomod.fit(Rb_nobias_fit, params, x=0 , method='ampgo')
##result_Rb_nobias = result_grabber(Rb_nobias_res)
##print('Finished with Rb_nobias')
##
## =============================================================================
## There is no reason to have any current in the coils indexed 0-11 for the Strontium solution as they won't impact the actual slower field.
## So we need a little loop to set those param values to 0 and vary=False
## =============================================================================
#
#for i in range(1, 11):
#    params['B%s_Coil_Current' % str(i)].set(value=0, vary=False)
#
#Sr_bias_res = Combomod.fit(Sr_bias_fit, params ,x=0, method='ampgo')
#result_Sr_bias = result_grabber(Sr_bias_res)
#print('Finished with Sr_bias')
#
#Sr_nobias_res = Combomod.fit(Sr_nobias_fit, params, x=0, method='ampgo')
#result_Sr_nobias = result_grabber(Sr_nobias_res)
#print('Finished with Sr_Noias!')
#result_Rb_bias = np.zeros((4,20))
#result_Rb_nobias = np.zeros((4,20))

## $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# =============================================================================
# If the Strontium Solution is not being run, include this section so that the output file doesn't throw an error.
# =============================================================================

result_Sr_bias = np.zeros((4,20))
result_Sr_nobias = np.zeros((4,20))

# =============================================================================
# Place the results into a csv file
# =============================================================================

wire_len = []
for ii in range(CP_Count):
    wiretmp = 0
    for ix in range(result_Rb_bias[0][ii]):
        temprad = (GLA.Slower_radius + 0.5*(0.5+ix)*GLA.wire_diameter)
        wiretmp = wiretmp + 2.0*np.pi*temprad  
    wire_len.append(wiretmp*result_Rb_bias[2][ii])

# =============================================================================
# Importing time so that every solution generated will have a uniquely identifying value based on a ~45 day internal windows system clock.
# =============================================================================

import time

timeval = str(int(time.monotonic()))
print(timeval)

# =============================================================================
# Create a new directory at this location for the lmfit solution to be saved into, keeping the solutions separate.
# =============================================================================
#os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/Thesis Version (NO TOUCHING)/Main Program (Clean Edit)/LmFitSolutions")
#os.makedirs('./'+timeval+'/')
#os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/Thesis Version (NO TOUCHING)/Main Program (Clean Edit)/LmFitSolutions/" +timeval+"/")

# =============================================================================
# Save the values to a csv file with the wire_diameter, slower_radius, and initial offset, as well as the timeval unique ID in the first row before the header row for the data.
# =============================================================================

#with open('Zeeman_Params' + timeval + '.csv', 'w', newline='') as f:
#    writer = csv.writer(f)
#    writer.writerow([timeval,'wire diameter (m)' , GLA.wire_diameter, 'Slower_radius' , GLA.Slower_radius , 'initial offset = ', cp_initial_offset]) #result_Rb_bias[3][0]])
#    writer.writerow(['Solenoid Number', 'Solenoid Thickness (# wires)','Solenoid Outer Radius (m)' , 'Solenoid Length (# wires)','Solenoid Length (m)', 'Wire length per Solenoid (m)' ,'Rb_nobias_current' , 'Rb_withbias_current'])
#    for rows in range(CP_Count-1):
#        writer.writerow([rows, result_Rb_bias[0][rows], GLA.Slower_radius+result_Rb_bias[0][rows]*GLA.wire_diameter ,result_Rb_bias[2][rows] ,result_Rb_bias[2][rows]*GLA.wire_diameter, wire_len[rows]   ,  result_Rb_nobias[1][rows] , result_Rb_bias[1][rows]])


p.figure()
#p.title('Rubidium With Bias Field')
#p.plot(GLA.Plot_Range, Rb_bias_Int, ':', label='Slower Field Solution')
p.plot(GLA.Slower_Length, GLA.Ideal_Field, ':' , label="Ideal Rb Zeeman Field")
p.plot(GLA.Slower_Length, Rb_bias_fit, '-' , label="Target Field")
#p.plot(GLA.Slower_Length, GLA.B_Bias_Trimmed, '-.' , label= "High-B Bias Field")
p.plot(GLA.Slower_Length, Rb_bias_interim_res.init_fit, label="Init Fit")
p.plot(GLA.Slower_Length, Rb_bias_interim_res.best_fit, label="Best Fit")

#Rb_interim.plot_fit()
#Rb_bias_interim_res.plot_fit()
#Rb_bias_res.plot_fit()
#p.plot(GLA.Slower_Length, Rb_interim.best_fit, label="Interim Best Fit")
p.plot(GLA.Slower_Length, Rb_nobias_res.best_fit, label="Nobias Fit")
p.plot(GLA.Slower_Length, Rb_bias_res.best_fit, label="Best Fit")
p.legend(loc='best')
p.show()