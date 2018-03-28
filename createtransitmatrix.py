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
import createrawmatrix as raw

class structtype():
    pass

region = 'global'
folder = 'rawMatrices'

Mbuoy = structtype()
Mmodel = structtype()
M = raw.M     #import initial values (such as M.nt, M.nc) from createrawmatrix.py

#target area(s):
lon = M.lons[-1]
lat = M.lats[-1]

M.buoyenhancefactor = 10

Mbuoy.P = [[] for x in range(M.nt)]
for x in range(M.nt):
    Mbuoy.P[x] = sparse.load_npz(os.path.join(folder, 'buoytransitmatrix_'+M.dir+'_raw_0_monthindex'+str(x)+'.npz'))
 
numberOFES = 19

Mmodel.P  = [sparse.coo_matrix((M.nc[-1],M.nc[-1])) for x in range(M.nt)]
M.P  = [sparse.coo_matrix((M.nc[-1],M.nc[-1])) for x in range(M.nt)]
M.nCross  = [sparse.coo_matrix((M.nc[-1],M.nc[-1])) for x in range(M.nt)]

for i in range(numberOFES) :
    for x in range(M.nt):
        b = sparse.load_npz(os.path.join(folder, 'ofestransitmatrix_'+M.dir+'_raw_'+str(i)+'_monthindex'+str(x)+'.npz'))
        Mmodel.P[x] = Mmodel.P[x] + b

M.buoyenhancefactor = 10
#weigh buoy observational data with OFES model data
for m in range(M.nt):
    if np.isfinite(M.buoyenhancefactor) and M.buoyenhancefactor > 0:
        totbuoy = 0
        totofes = 0
        for x in range(M.nt):
            totbuoy = totbuoy + Mbuoy.P[x].sum()
            totofes = totofes + Mmodel.P[x].sum()
        factor = M.buoyenhancefactor*totofes/totbuoy
        for x in range(M.nt):
            M.nCross[x] = Mbuoy.P[x] + Mmodel.P[x]/factor

A = structtype()
#normalizing matrix
for m in range(M.nt):
  A.raw = M.nCross[m].tocoo()
  A.prob = np.zeros(A.raw.row.size) #initialize probability result list
  A.rows = A.raw.row #sparse matrix coordinate i, location particle moves from
  A.cols = A.raw.col #sparse matrix coordinate j, location particle moves to
  nextlocation = (np.array(np.where(np.ediff1d(A.rows)!=0))+1).flatten() #make list of where the row coordinates change, i.e. the next origin location
  A.starts = np.insert(nextlocation,0,0) 
  A.ends = np.append(A.starts[1:] - 1, np.size(A.rows)-1) #all above values minus one, so making iterative range like 1:105, 106:220, 221:332 etc with starts:ends. Last value appended, for taking last value in A.rows

  for x in range(np.size(A.starts)):
    I = np.arange(A.starts[x],A.ends[x]+1,1) #determinant for one origin coordinate, similar to createrawmatrix
    A.prob[I] = A.raw.data[I]/A.raw.data[I].sum() #normalizing by dividing by sum of all crossing counters per origin location value
  M.P[m] = sparse.coo_matrix((A.prob,(A.rows,A.cols)),shape=(M.nc[-1],M.nc[-1]))

#fix points that change into land points
good = np.zeros((M.nc[-1],M.nt))
np.savetxt(os.path.join('Good', 'goodoriginal.csv'),good.sum(axis=1),delimiter=",")

indi = [[] for x in range(M.nt)]
for m in range(M.nt):
    good[:,m] = M.P[m].sum(axis=1).transpose() #creates vector with value-sum of probability per origin location, one vector per month. Total in all of those months should be 1 due to normalization, total of all months should be 6.
    indi[m] = np.zeros(M.nc[-1]) #coordinate list per month
bad = np.where((good.sum(axis=1)<M.nt-0.5) & (good.sum(axis=1)>0.5)) #sum of 'good' finds sum of all prob-sums for all months, per location. Returns number of values equal to elements in sparse matrix
tel = np.zeros(M.nt,dtype='int')                                     #bad' gives indices for which locations have particles in them for some period in the year, but not for the full year. good.sum = 0 never sees a buoy, i.e. land, good.sum = M.nt has buoys going in and out all year -> this conditional takes everything in between.

#for every location that doesn't have crossings all through the year, find the indices for which month the probability sum is 0.
mon = np.where(good[bad[0]]==0) #mon's shape: first array is index in "bad" array that has 0-entry, second array is the months in which it is 0. E.g. first "bad" entry has 0-entry in month 4 and 5, mon[0][0:2] = [0,0] and mon[1][0:2] = [4,5]
unique, index, counts = np.unique(mon[0], return_index=True, return_counts=True) #unique is every unique element, index where it starts, counts how many month entries are 0 for that index. 
for b in unique:
    for m in mon[1][index[b]:index[b]+counts[b]]:
        indi[m][tel[m]] = bad[0][b] #store bad location coordinate in J.indi for corresponding month
        tel[m] += 1

indicrop = [[] for _ in np.arange(M.nt) ]
for m in np.arange(M.nt):
    indicrop[m] = indi[m][0:tel[m]]
    M.P[m]=  M.P[m] + sparse.coo_matrix((np.ones(indicrop[m].size),(indicrop[m],indicrop[m])),shape=(M.nc[-1],M.nc[-1]))
    M.P[m].tocoo()

#fix points where you can get in but not get out
for m in range(M.nt):
    m2 = m%(M.nt-1) + 1 #for next monthstep - basically m+1 with correction for last timestep in year
    indi = np.zeros(M.nc[-1])
    indj = np.zeros(M.nc[-1])
    vals = np.zeros(M.nc[-1])
    tel = 0

    condition1 = np.squeeze(np.asarray(M.P[m].sum(axis=0)))             #condition1 takes sum over all origin probability values per destination coordinate for m, second sum takes destination values per origin for m+1
    condition2 = np.squeeze(np.asarray(M.P[m2].sum(axis=1)).transpose())#i.e.: take coordinates in a timestep that have buoys/particles in it at that timestep, but do not leave that location at all in the following timestep. Origin probability sum = 0 means that the particle will not leave, since no crossing to another cell has been determined starting from that point.
    rows, = np.where((condition1>0) & (condition2==0))
    
    for r in range(rows.size):
        cols = np.nonzero(M.P[m].tocsr()[:,rows[r]])[0] #pick all origins in this month that move to the chosen coordinate (rows[r]) in the next month.
        for c in range(len(cols)):
            indi[tel] = rows[r]
            indj[tel] = cols[c]
            vals[tel] = M.P[m].tocsr()[int(indj[tel]),int(indi[tel])] #Reflective boundary condition, i.e.: if particles can get in the cell, going out happens in the same way. This line takes the probility value for month m and reverses its direction.
            tel +=1

    indicrop = indi[0:tel]
    indjcrop = indj[0:tel]
    valscrop = vals[0:tel]
    M.P[m2]  = M.P[m2] + sparse.coo_matrix((valscrop,(indicrop,indjcrop)),shape=(M.nc[-1],M.nc[-1])) #add the reflection to next month

#renormalize
for m in range(M.nt):
    A.raw = M.P[m].tocoo()
    A.prob = np.zeros(A.raw.row.size)
    A.rows = A.raw.row
    A.cols = A.raw.col
    nextlocation = (np.array(np.where(np.ediff1d(A.rows)!=0))+1).flatten()
    A.starts = np.insert(nextlocation,0,0) #make list of where the row coordinates change, i.e. the next origin location
    A.ends = np.append(A.starts[1:] - 1, np.size(A.rows) - 1) #all above values minus one, so making iterative range like 1:105, 106:220, 221:332 etc with starts:ends. Last value appended, for taking last value in A.rows

    for x in range(np.size(A.starts)):
        I = np.arange(A.starts[x],A.ends[x]+1,1)
        A.prob[I] = A.raw.data[I]/A.raw.data[I].sum()
    M.P[m] = sparse.coo_matrix((A.prob,(A.rows,A.cols)),shape=(M.nc[-1],M.nc[-1]))


#update 'good' with new normalized values
for m in range(M.nt):
    good[:,m] = M.P[m].sum(axis=1).transpose()
np.savetxt(os.path.join('Good', 'goodupdate.csv'),good.sum(axis=1),delimiter=",")

landpoints = good.sum(axis=1)
landpoints[good.sum(axis=1)>0] = 0
landpoints[good.sum(axis=1)==0] = 1
np.savetxt(os.path.join('Landpoints', 'landpoints.csv'),landpoints,delimiter=",")

for m in range(M.nt):
    sparse.save_npz(os.path.join('TransitMatrices/', 'hybridtransitmatrix_' + M.dir + '_month_'+str(m) +'.npz'), M.P[m])
