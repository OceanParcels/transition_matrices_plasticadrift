from __future__ import division
import random
import sys
import os, os.path
os.environ["TF_CPP_MIN_LOG_LEVEL"]="3"
import csv
import numpy as np
from scipy import sparse
import time
start_time = time.time()
import scipy.io as sio
import datetime
import warnings 
with warnings.catch_warnings():
  warnings.filterwarnings("ignore",category=FutureWarning) 
  import h5py

class structtype():
    pass

def timeToMonths(timeList) :
  months = []
  for x in timeList:
    time = datetime.datetime.fromordinal(int(x)-366).strftime('%Y-%m-%d %H:%M:%S')
    months.append(int(time[5:7]))
  return(np.array(months))

#form M: Lons/lats consists of lists for target areas. The last list in M.lons, M.lats is the full global scale, 
#and all lists preceding are smaller target areas with a finer resolution, e.g. a region like Europe.
#dd is the resolution of the grid, i.e. the step size for lons/lats for the different regions. M.nd is used for the amount of loops over the code corresponding to the amount of areas considered. 

M = structtype()

M.dd = np.array([1])
M.nd = M.dd.size
M.lons, M.lats = [], []
for x in np.arange(M.nd):
    M.lons.insert(0,[])
    M.lats.insert(0,[])

#define areas for matrix. Make sure the global area is the last entry.
M.xStarts = np.array([-180])
M.xEnds   = np.array([180])

M.yStarts = np.array([-90]) 
M.yEnds   = np.array([90]) 

for i in np.arange(len(M.lons)) :
  M.lons[i] = np.linspace(M.xStarts[i], M.xEnds[i],  ((M.xEnds[i] - M.xStarts[i])/ M.dd[i] + 1), dtype='i4')
  M.lats[i] = np.linspace(M.yStarts[i], M.yEnds[i],  ((M.yEnds[i] - M.yStarts[i])/ M.dd[i] + 1), dtype='i4')

# nx, ny number of elements in lons, lats per region. nc is the number of total grid cells used in later construction of raw matrix.
# Loop for all target areas
M.nc = np.empty(np.shape(M.dd))
M.nx, M.ny, M.nc = [[] for x in M.dd],[[] for x in M.dd],[[] for x in M.dd]

for d in range(M.nd):
    M.nx[d] = np.size(M.lons[d])
    M.ny[d] = np.size(M.lats[d])
    M.nc[d] = int(M.nx[d] * M.ny[d] + np.sum(M.nc[0:d]))
M.tau=60 #days. This is the timestep used for integration.
M.nt = int(365/M.tau) #total amount of timesteps

M.dir = 'fwd' #for dispersion forward in time, other option is 'bwd' for backward dispersion
types= ['buoy','ofes'] #['test']

# import all trajectory files corresponding to type, hence * in path.
# Different types are combined in createtransitmatrix.py

if __name__ == "__main__" :         #in createtransitmatrix.py this file is imported as module to copy the initial values set above. This conditional prevents running the full code again when imported as module.
  for tp in range(np.size(types)):
    if types[tp] == 'buoy' :
      filepath=[os.path.join('Buoys','buoydata.mat')]
      nFiles = 1
    elif types[tp] == 'ofes' :
      folder='expt_ofes01pd_reinit'
      nFiles = 19
      filepath = []
      for i in range(1,nFiles+1) :
        filename = os.path.join(folder,'traj_' + str(i).zfill(2) + '.mat')
        filepath.append(filename)
    elif types[tp] == 'test' :
      filepath=[os.path.join('TestfilesTraj','testfile.mat')]
      nFiles = 1
    else :
      print('Type unknown')
    
    T = structtype()
    #buoy data and OFES data are in different Matlab version formats, 
    #hence they have to be imported with seperate methods. The resulting arrays in this code are the same format.
    for f in range(nFiles):
      if types[tp] == 'buoy' :
        mat_contents = sio.loadmat(filepath[f])
        T.date   = np.swapaxes(mat_contents['date'],1,0)[0]
        T.lon    = np.swapaxes(mat_contents['lon'],1,0)[0]
        T.lat    = np.swapaxes(mat_contents['lat'],1,0)[0]
        T.starts = np.swapaxes(mat_contents['starts'],1,0)[0]
        T.ends   = np.swapaxes(mat_contents['ends'],1,0)[0]
        
      if types[tp] == 'ofes' :
        datafile = h5py.File(filepath[f],'r')
      
        T.date   = np.array(datafile.get('date'))[0]
        T.lon    = np.array(datafile.get('lon'))[0]
        T.lat    = np.array(datafile.get('lat'))[0]
        T.starts = np.array(datafile.get('starts'))[0]
        T.ends   = np.array(datafile.get('ends'))[0]  

      if types[tp] == 'test' :
        mat_contents = sio.loadmat(filepath[f])
        T.date   = np.swapaxes(mat_contents['date'],1,0)[0]
        T.lon    = np.swapaxes(mat_contents['lon'],1,0)[0]
        T.lat    = np.swapaxes(mat_contents['lat'],1,0)[0]
        T.starts = np.swapaxes(mat_contents['starts'],1,0)[0]
        T.ends   = np.swapaxes(mat_contents['ends'],1,0)[0]

      #adjust 0 < lon < 360 to -180 < lon < 180
      if M.lons[0][0] < 0 :
        T.lon[np.where(T.lon > 180)] = T.lon[np.where(T.lon > 180)]-360
      if M.dir == 'fwd' :
        dm = int(M.tau/(T.date[2]-T.date[1]))
      elif M.dir == 'bwd':
        dm = int(M.tau/(T.date[1]-T.date[2]))

      #ni functions as index through nc for those grid values that are in the target area considered. 
      #The next few lines set the index to negative values so that they are not selected by later conditionals, unless they are changed to positive values because they are in the target area considered.
      T.ni = np.empty(np.shape(T.lon))
      T.ni[:] = -9999

      for d in np.arange(M.nd):
        D = np.where((T.lon>=M.lons[d][0]) & (T.lon <= M.lons[d][-1]) & (T.lat >= M.lats[d][0]) & (T.lat <= M.lats[d][-1]))[0] #take data lons and lats between defined bounds of Model lons and lats.
        xi= np.floor( (T.lon[D]-M.lons[d][0])/M.dd[d] ).astype(int) #iterative start values for lon/lat within area are data values lon/lat minus first lon/lat values divided by stepsize
        yi= np.floor( (T.lat[D]-M.lats[d][0])/M.dd[d] ).astype(int)  #see above. If lon of data and lon of model target match, xi[0] = 0, otherwise e.g. lon.data = 4, lon.target = 1, 5 samples per 1 lat -> xi starts at ((4-1)/0.2) = 15. Note: xi/yi is all eligible data values, not just start value

        D_empty_T = D
        #shifts indices for new regions, see definition M.nc
        if d > 0 :
          D += M.nc[d-1]

        #creates matrix with number of lat/lon samples as dimensions, and takes elements xi, yi from it, and stores those indices in T.ni.
        for i in np.arange(np.size(xi)):
          T.ni[D[i]] = np.ravel_multi_index((xi[i],yi[i]),(M.nx[d],M.ny[d]),order='F') #PROBLEM: so far no solution for crossing paths between regions as defined at starts. Buoys crossing the boundary will be rejected.

        T.lon[D_empty_T]=-9999
        T.lat[D_empty_T]=-9999

      #convert Matlab date format for each datapoint to the corresponding month, and then bin those months in the timestep index.
      #Only works when timestep >= 1 month.
      T.monthi = timeToMonths(T.date)                         #e.g. T.monthi = 3 for datapoint recorded in April
      T.monthi = np.floor((T.monthi-1)/(12/M.nt)).astype(int) #e.g. T.monthi = 1 for same datapoint with M.tau = 60 days (0 for Jan/Feb, 1 for Mar/Apr , ... , 5 for Nov/Dec)

      #initialize couplets origin/destination indi and indj for every timestep in one year. Its function is index in result sparse matrix
      J = structtype()
      J.indi , J.indj = np.array([np.zeros(int(np.round(np.size(T.lon)))) for _ in np.arange(M.nt) ],dtype='int'), np.array([np.zeros(int(np.round(np.size(T.lon)))) for _ in np.arange(M.nt) ],dtype='int')

      tel=np.zeros(M.nt,dtype='int')

      #This loops for every seperate buoy. For every datapoint that has another datapoint at exactly one timestep later,
      #store the original location in indi and the destination after the timestep in indj, in the list for the correct month.  
      for d in range(np.size(T.starts)):
        II = np.arange(T.starts[d]-1,T.ends[d]-1,1,dtype='int')
        ni = T.ni[II]
        ti = T.date[II]
        mi = T.monthi[II]

        for i in np.arange(np.size(ni)-1) :
          j = i+dm
          if j < np.size(ni) and j>0 :
            if np.absolute(ti[j]-ti[i]).astype(int) == M.tau :
              if ni[i] > -9900 and ni[j] > -9900: #PROBLEM: so far no solution for crossing paths between regions as defined at starts. Buoys crossing the boundary will be rejected.
                J.indi[mi[i]][tel[mi[i]]] = ni[i]
                J.indj[mi[i]][tel[mi[i]]] = ni[j]
                tel[mi[i]] += 1
      M.nCross = [[] for _ in np.arange(M.nt) ]
      J.indicrop = [[] for _ in np.arange(M.nt) ]
      J.indjcrop = [[] for _ in np.arange(M.nt) ]

      #delete empty entries in indi and indj
      for m in np.arange(M.nt):
        J.indicrop[m] = J.indi[m][0:tel[m]]
        J.indjcrop[m] = J.indj[m][0:tel[m]]

        #sparse matrix per timestep with every origin and destination as coordinates in matrix, with all possible locations M.nc on i and j axis. 
        sparseV = np.ones(np.shape(J.indicrop[m]))
        sparseI = J.indicrop[m]
        sparseJ = J.indjcrop[m]
        M.nCross[m] = sparse.coo_matrix((sparseV,(sparseI,sparseJ)),shape=(M.nc[-1],M.nc[-1]))
        sparse.save_npz(os.path.join('rawMatrices/',str(types[tp]) + 'transitmatrix_' + M.dir + '_raw_' + str(f) + '_monthindex' + str(m) + '.npz'), M.nCross[m])
  print("--- %s seconds ---" % (time.time() - start_time))
