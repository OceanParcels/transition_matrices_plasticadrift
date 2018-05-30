import random
import sys
import os, os.path
os.environ["TF_CPP_MIN_LOG_LEVEL"]="3"
import csv
import numpy as np
from scipy import sparse
import calendar

import createrawmatrix

class structtype():
    pass

def ravel(lat, lon) : #function for conversion between latitudes and longitudes (user friendly) to index of model
    ind = np.ravel_multi_index((lon-M.lons[-1][0],lat-M.lats[-1][0]),(M.nx[-1],M.ny[-1]),order='F')
    return ind

region = 'global'
inputFolder = 'TransitMatrices'
outputFolder = 'CSV_files'

beaching = True #whether to include beaching into results or not
beachingmodel = 'frac' #other option: 'adj'. Set beaching model to either the Fractional Coast Model or Adjacent Coast Model, see paper.
scaling = 0.3 #scaling factor for beaching model, in procent

M = createrawmatrix.M #import createrawmatrix for parameters of model (such as range longitudes, latitudes, dx, dt)
landpoints = np.genfromtxt('Landpoints/landpoints.csv', delimiter=',') #array of 1's for land and 0's for sea for every grid cell in model

maxyears = 10 #duration of simulation
minplotval = 1e-12 #minimal value included in results. If some location has less tracer in it than this amount, the tracer is removed from the ocean.
offset = np.arange(M.nt) #array of all starting timesteps that are evaluated
numints = maxyears*M.nt + 1 #number of computational iterations

#release origins to be evaluated. Input is in indices of model, but can be in latitude and longitude with ravel()
indices = ravel(np.array([51,52,53,53,53]) , np.array([3, 4, 5, 6, 7])) #indices for Dutch coast

#rename target area(s) from createrawmatrix for brevity
lon = np.array(M.lons[-1])
lat = np.array(M.lats[-1])

transitMatrix = [[] for x in range(M.nt)]
for x in range(M.nt):
    transitMatrix[x] = sparse.load_npz(os.path.join(inputFolder, 'hybridtransitmatrix_'+M.dir+'_month_'+str(x)+'.npz'))

if beaching: #load beaching data corresponding to set model and determine the chance of particle tracer ending on shore
    beachdata = np.genfromtxt(os.path.join('beachingCoastlines', beachingmodel + '_coast_dx_' + str((M.dd[-1] * 100).astype(int)) + '_region_'+region+'.csv'),delimiter=',',skip_header=1)
    beachdata = np.transpose(beachdata)
    beach = 1 - beachdata[3] * scaling #chance of tracer ending on shore

    if beach.size != M.nc[-1] :
        print("Resolution of beaching grid and transition matrix don't match")
        sys.exit()

for i in indices:
    if landpoints[i] == 0 :
        folder = outputFolder #reset the output folder
    
        if beaching: #create subfolder within index folder to differentiate beaching scaling sizes 
            try :
                os.mkdir(os.path.join(folder, "beach_" + str(scaling * 100) + "_perc"))
                folder = folder + "//beach_" + str(scaling * 100) + "_perc"
            except:
                folder = folder + "//beach_" + str(scaling * 100) + "_perc"
        else:
            try : #create subfolder per index, unless it already exists
                os.mkdir(os.path.join(folder, "no_beaching"))
                folder = folder + "//no_beaching"
            except:
                folder = folder + "//no_beaching"

        for t in offset :
            m = int(np.floor(t * 12/M.nt)) #calendar month on which the iteration starts (first day of month m)

            if beaching : #save tracer per location per timestep in "savefile" and save the amount of tracer ending on shore per location per timestep in "totalbeachedfile"
                savefile = os.path.join(folder, str(i) + "_startsin" + calendar.month_name[m+1] + "_beaching_" + str(scaling * 100) + ".csv") #m+1 because month_index 0 is month 1 (January)
                totalbeachedfile = os.path.join(folder, str(i) + "_startsin" + calendar.month_name[m+1] + "_beaching_" + str(scaling * 100) + "totalbeached.csv")
            else: #save tracer per location per timestep in "savefile"
                savefile = os.path.join(folder, str(i) + "_startsin" + calendar.month_name[m+1] + ".csv") #m+1 because month_index 0 is month 1 (January)
            f1 = open(savefile,'w')
            f1.write('year,month,lat,lon,prob\n') #save index header
            if beaching:
                f2 = open(totalbeachedfile,'w')
                f2.write('year,month,lat,lon,beachedtracer\n') #save index header

            v = np.zeros(M.nc[-1])
            v[i] = 1 #set release location to 1, all other locations 0 for first iteration.
            count = 0 #counter denoting current iteration

            for n in np.arange(numints) + t :
                bm = n % M.nt #timestep of current iteration, modulated by number of timesteps to account for offset t
                mon = np.floor((count % M.nt) * 12/M.nt).astype(int) #calendar month of current timestep
                year = np.floor(count/M.nt).astype(int) #year of current timestep

                #The vector-matrix multiplication
                nums, =np.where(v > minplotval)
                pX, pY = np.unravel_index(nums,(M.nx[-1],M.ny[-1]),order='F') #convert index to geospatial coordinates

                vals = v[nums] #save all tracer amounts > minplotval to file f1
                for j in np.arange(vals.size) :
                    val = vals[j]
                    f1.write(str(year) + ',' + str(mon) + ',' + str(lat[pY[j]]) + ',' + str(lon[pX[j]]) + ',' + str(val) + '\n')
                count += 1
                v = v * transitMatrix[bm] #multiply v with transit matrix for next iteration


                if beaching:
                    initialv = v
                    v = v * beach #multiply v with chance of particles beaching at that location
                    diffv = initialv - v 

                    nums, =np.where(diffv > minplotval)
                    pX, pY = np.unravel_index(nums,(M.nx[-1],M.ny[-1]),order='F')

                    vals = diffv[nums]
                    for j in np.arange(pX.size) : #save the total amount of beached tracer per location to file f2
                        val = vals[j]
                        f2.write(str(year) + ',' + str(mon) + ',' + str(lat[pY[j]]) + ',' + str(lon[pX[j]]) + ',' + str(val) + '\n')
            f1.close()
            if beaching:
                f2.close()
    else:
        print('Index ' + str(i) + ' is a point on land.')
