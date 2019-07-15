import numpy as np
import matplotlib.pyplot as p
import os, csv
os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/Thesis Version (NO TOUCHING)/github version/")


# =============================================================================
# Import the data from the file BiasBIn.txt. This file should be in the same directory as the rest of the files.
# =============================================================================
#data = np.genfromtxt('BiasBin.txt',usecols=2)
#B_Mag = [10000*x for x in data]
#B_Mag_Tesla = [x for x in data]
#lin = len(B_Mag)/2
#zlin = np.linspace(-1.5,1.5, len(B_Mag))


# =============================================================================
# This section imports the data from a csv file instead
# =============================================================================
B_Mag = []
zlin = []
with open('BiasBIn.csv', mode='r', newline='') as readme:

    reader = csv.reader(readme, delimiter=',', quoting=csv.QUOTE_NONE)
    next(reader)
    for row in reader:
        zlin.append(float(float(row[1])))
        B_Mag.append(float(float(row[1])))

# =============================================================================
# #Slice off the area of the field we are interested in for a 0.68m slower.
# =============================================================================


BBL = 2000
BBR = 8801
B_Bias = B_Mag[BBL:BBR]
B_Bias_Long = B_Mag[(BBL-500):(BBR+2000)]
#zlinbias = zlin[BBL:BBR]
#zz_bias = zlin[920:10920]

#B_Bias2 = B_Mag_Tesla[920:10920]

# =============================================================================
# Plot the Bias Field
# =============================================================================

#p.close('all')
#p.figure()
#p.title('High-B Field for 2.6T Trap', fontsize = 18)
#p.plot(zz_bias[1000:9000],B_Bias2[1000:9000], color='red' , linewidth = 4.0)
#p.plot(zlin, B_Mag_Tesla , ':' )
#p.ylabel('High-B Magnetic Field (T)', fontsize = 16)
#p.xlabel('Distance from Trap Center (m)' , fontsize = 16)
#
#os.chdir(r"C:/LN_Thesis_EMU/Zeeman Plots and Diagrams/Final Collection")
#p.savefig('HighBField2.png', dpi=800, bbox_inches='tight', pad_inches=0.1)
#p.show()

#os.chdir(r"C:/LN_Thesis_EMU/Coding/Python/Main Program")
 
# =============================================================================
# Create an excel spreadsheet with the bias field values
# =============================================================================

# os.chdir(r"F:/Zeeman Plots and Diagrams/Nov13/Solution_V2")
# p.savefig('Trimmed_Bias_Field.png', bbox_inches='tight', dpi=1200)
# with open('Trimmed_Bias_Field.csv', 'w', newline='') as f:
#     writer = csv.writer(f)
#     writer.writerow(['Position (m)', 'Bz field (G)'])
#     for rows in range(len(B_Bias)):
#         writer.writerow([zz_bias[rows] , B_Mag[rows]])