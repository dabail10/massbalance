#--------------------------------------------------------------------
# Example script to demonstrate how to read/write mass budget ASCII files
# The script will:
#  1) Read in saved budget data from example file
#  2) Plot the data to demonstrate that the budget balances
#  3) Write back out to ascii as a deomnstration
#
# Ann Keen & Ed Blockley, Feb 2019
#--------------------------------------------------------------------
 
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Read in the saved data and put into arrrays
#--------------------------------------------

# define input ascii data file
#indatafile='NCAR_CESM2f09g17_HIST_001_ice.txt'
indatafile='NCAR_CESM2f09g17_BL99_001_ice.txt'
#indatafile='NCAR_CESM2f09g17_MUSHY_001_ice.txt'

# how many lines to skip in the header?
nskip=4  # default=4
data = np.genfromtxt(indatafile, skip_header=nskip)

# uncomment for diagnostics info on data file array if required
#print type(data)
#print data.shape

# obtain length of dataset (no. of monthly values)
ndates=len(data[:,0])
# define arrays
years = np.empty([ndates],dtype='int')
months = np.empty([ndates],dtype='int')
icearea = np.empty([ndates],dtype='float')
icemass = np.empty([ndates],dtype='float')
massbudget = np.empty([ndates,9],dtype='float')

# populate arrays (dates, area & mnass)
years = data[:,0]
months = data[:,1]
icearea = data[:,2]
icemass = data[:,3]
# populate buudget array
massbudget[:,0] = data[:,4] # growthbot
massbudget[:,1] = data[:,5] # growthwat
massbudget[:,2] = data[:,6] # melttop
massbudget[:,3] = data[:,7] # meltbot
massbudget[:,4] = data[:,8] # lat
massbudget[:,5] = data[:,9] # si
massbudget[:,6] = data[:,10] # evapsubl
massbudget[:,7] = data[:,11] # dyn
massbudget[:,8] = data[:,12] # total

# Calculate scaled budget in units of Kg/month:
seconds_in_month = 60*60*24*30
massbudgetmonthly = massbudget * seconds_in_month

# ---- Plotting colours etc ----------

lineColours=['red',        # congel/sidmassgrowthbot
             'magenta',    # frazil/sidmassgrowthwat
             'blue',       # meltt/sidmassmelttop
             'lime',       # meltb/sidmassmeltbot
             'cyan',       # meltl/sidmasslat
             'orange',     # snowice/sidmasssi
             'olive',      # evap_ai/sidmassevalsubl
             'darkgrey',   # dvidtd/sidmassdyn
             'black']      # total
Var_labels = ['sidmassgrowthbot',    
                 'sidmassgrowthwat', 
                 'sidmassmelttop',    
                 'sidmassmeltbot',    
                 'sidmasslat',    
                 'sidmasssi',    
                 'sidmassevapsubl',   
                 'sidmassdyn',   
                 'total'   ]
lineMarkers=["^",    # triangle-up
             "<",    # triangle-left
             "v",    # triangle-down
             "8",    # octagon
             "o",    # circle
             "D",    # diamond
             ">",    # triangle-right
             "h",    # hexagon
             "s"]    # square

# use font size=15 always
matplotlib.rcParams.update({'font.size':15})

# ----------- Plot seasonal cycle of budget terms for the first year of data ------
# Using the monthly scaled budgets

# define figure and axis
fig1 = plt.figure(figsize = (8,5)) # Landscape
ax1 = fig1.add_subplot(111)
ax1.set_title('Arctic sea ice mass budget')
ax1.set_ylabel('Mass flux (kg month$^{-1}$)')
ax1.set_xlabel('Month')

# plot budget terms
ax1.plot([0,13],[0.0,0.0],color="black")
for nvar,name in enumerate(Var_labels):
    thisLine = ax1.plot(months[588:599],
                        massbudgetmonthly[588:599,nvar], 
                        color=lineColours[nvar], 
                        linestyle='-',
                        marker=lineMarkers[nvar],
                        label=name,
                        linewidth=1)
ax1.set_xlim([0, 13])
ax1.legend(loc="best", fontsize = 'x-small')    
    
fig1.savefig('Budget_scycle_ice.png', dpi=300)
#plt.show()

# ------ Plots demonstrating the balance for the first year of data -------------------------

# Set up the dmass arrays
# we do this twice for the monthly values and the mid-point between each month
# this helps to account for the discrepancy between monthly-mean and instantaneous mass
mass_midmonths = np.empty([11],dtype=float)
dmass = np.empty([11],dtype=float)
dmass_midmonths = np.empty([10],dtype=float)

for month in range(0,11):
    mass_midmonths[month] = 0.5*(icemass[month+1] + icemass[month])
    dmass[month] = icemass[month+1] - icemass[month]

for month in range(0,10):
    dmass_midmonths[month] = mass_midmonths[month+1] - mass_midmonths[month]

midmonths=months[588:599]+0.5  
    
fig2 = plt.figure(figsize = (15,5)) # Landscape

# -------- Seasonal cycle of Actual Ice Mass  -------------------------
ax2 = fig2.add_subplot(121) # (rows,columns,relative position)
ax2.set_title('Arctic ice mass (kg)')
ax2.set_ylabel('Mass (Kg)')
ax2.set_xlabel('Month')
ax2.set_xlim([0, 13])

ax2.plot(months[588:599],icemass[588:599],color='purple',marker=lineMarkers[6],label='mass') 
ax2.plot(midmonths,mass_midmonths,color='olive',marker=lineMarkers[5],label='mass_midmonths') 

ax2.legend(loc="best", fontsize = 'x-small')

#--------- Mass change from budget terms ------------------
ax3 = fig2.add_subplot(122) # (rows,columns,relative position)
ax3.set_xlim([0,13])
ax3.set_title('Arctic ice mass budget')
ax3.set_ylabel('Mass flux (kg month$^{-1}$)')
ax3.set_xlabel('Month')

ax3.plot(months[588:599],massbudgetmonthly[588:599,-1],color='black',marker=lineMarkers[7],label='total') 
ax3.plot(midmonths,dmass,color='purple',marker=lineMarkers[6],label='dmass') 
ax3.plot(months[1:11],dmass_midmonths,color='olive',marker=lineMarkers[6],label='dmass_midmonths') 
ax3.plot([0,13],[0.0,0.0],color="black")

ax3.legend(loc="best", fontsize = 'x-small')

fig2.savefig('Budget_balance_ice.png', dpi=300)
#plt.show()

# --- Re-generate the data file --------------------------------
# Using original budget terms in units of Kg s-1:

# define demo ascii data file
outdatafile='demo_ice.txt'

# define ascii data strcutures/formats
title_format = "%1s %4s %6s %14s %12s %17s %17s %15s %15s %12s %12s %15s %12s %12s"
data_format = "%6i %6i %14.5e %12.5e %17.5e %17.5e %15.5e %15.5e %12.5e %12.5e %15.5e %12.5e %12.5e"
headers = ('#','Year','Month','Area (Km**2)', 'Mass (Kg)', 
                 'sidmassgrowthbot',    
                 'sidmassgrowthwat', 
                 'sidmassmelttop',    
                 'sidmassmeltbot',    
                 'sidmasslat',    
                 'sidmasssi',    
                 'sidmassevapsubl',   
                 'sidmassdyn',   
                 'total')

# creat file and populate header
data_fileh = open(outdatafile,'w')
data_fileh.write('# Contact: Ann Keen ann.keen@metoffice.gov.uk')
data_fileh.write("\n")
data_fileh.write('# Corresponding HIST file: n/a')
data_fileh.write("\n")
data_fileh.write('# Components of the Arctic sea ice mass budget (Kg s-1):')
data_fileh.write("\n")
data_fileh.write(title_format % headers)
data_fileh.write("\n")

# write out mass budget for each month
ndates=len(years)
for index in range(ndates):
    data_fileh.write(data_format % (years[index],
                    months[index],
                    icearea[index],
                    icemass[index],
                    massbudget[index,0],
                    massbudget[index,1],
                    massbudget[index,2],
                    massbudget[index,3],
                    massbudget[index,4],
                    massbudget[index,5],                                   
                    massbudget[index,6],                                   
                    massbudget[index,7],
                    massbudget[index,8]))                                                                      
    data_fileh.write("\n")

# close ascii file
data_fileh.close()



