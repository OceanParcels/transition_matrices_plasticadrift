import random
import sys
import os
import csv
import numpy as np
from numpy import genfromtxt

import createrawmatrix

class structtype():
    pass

M = createrawmatrix.M

region = 'global'

def northCoast(row, col) :
    try:
        earthGrid[row + 1][col]
    except: #index is out of bounds for northernmost latitudes
        return False
    else:
        if earthGrid[row + 1][col] < 9998:
            return True
        else :
            return False

def southCoast(row, col) :
    try:
        if earthGrid[row - 1][col] < 9998:
            return True
        else :
            return False
    except: #index is out of bounds for southernmost latitudes
        return False

def eastCoast(row, col) :
    try:
        if earthGrid[row][col + 1] < 9998:
            return True
        else :
            return False
    except: #index is out of bounds for prime meridian, so check other end of array
        if earthGrid[row][0] < 9998:
            return True
        else :
            return False

def westCoast(row, col) :
    try:
        if earthGrid[row][col - 1] < 9998:
            return True
        else :
            return False
    except: #index is out of bounds for prime meridian, so check other end of array
        if earthGrid[row][-1] < 9998:
            return True
        else :
            return False

earthGrid = genfromtxt(os.path.join(os.getcwd(),'elevationdata.csv'), delimiter=',',dtype=float)
earthGrid = np.array(earthGrid)
#earthGrid is a 3601x1801 matrix, where the first column is latitude (89.95 up until -89.95, from North to South with increasing index)
#and the first row being latitude from -179.95 up to 179.95

lats = earthGrid[1::,0]
lons = earthGrid[0,1::]

#strip grid of index column and row and rearrange the matrix to fit index order of createrawmatrix (and other files)
earthGrid = np.delete(earthGrid,0,1) #delete latitude index column
shift = np.where(earthGrid[0] > 0)[0][0]
earthGrid = np.roll(earthGrid,shift,axis=1) #shift matrix from -180 < lon < 180 to 0 < lon < 360
earthGrid = np.delete(earthGrid,0,0) #delete longitude index row
earthGrid = np.flipud(earthGrid) #reverse matrix from 90 > lon > -90 to -90 < lon < 90

#rearrange index order
lons = np.roll(lons,shift)
lats = np.fliplr([lats])[0]

lons[np.where(lons < 0)] += 360

coastalGrid = np.empty(np.shape(earthGrid),dtype=int)

for row in np.arange(lats.size):
    for col in np.arange(lons.size):
        if earthGrid[row][col] > 9998: #sea
            coastalGrid[row][col] = 0
            if northCoast(row,col): #coastal check
                coastalGrid[row][col] += 1
            if southCoast(row,col):
                coastalGrid[row][col] += 1
            if eastCoast(row,col):
                coastalGrid[row][col] += 1
            if westCoast(row,col):
                coastalGrid[row][col] += 1              
        else: #land
            coastalGrid[row][col] = 10

#convert input grid to flat indices
inputMesh = np.meshgrid(lons,lats)

nx = inputMesh[0].flatten()
ny = inputMesh[1].flatten()
coastalGrid = coastalGrid.flatten()

ni = np.empty(coastalGrid.shape,dtype='int')
ni[:] = -9999

#compare coastalGrid to createrawmatrix resolution 
for d in np.arange(M.nd):
    D = np.where((nx >= M.lons[d][0]) & (nx <= M.lons[d][-1]) & (ny >= M.lats[d][0]) & (ny <= M.lats[d][-1]))[0] #take data lons and lats between defined bounds of Model lons and lats.
    xi = np.floor( (nx[D]-M.lons[d][0])/M.dd[d] ).astype(int) #iterative start values for lon/lat within area are data values lon/lat minus first lon/lat values divided by stepsize
    yi = np.floor( (ny[D]-M.lats[d][0])/M.dd[d] ).astype(int)  #see above. If lon of data and lon of model target match, xi[0] = 0, otherwise e.g. lon.data = 4, lon.target = 1, 5 samples per 1 lat -> xi starts at ((4-1)/0.2) = 15. Note: xi/yi is all eligible data values, not just start value

    for i in np.arange(xi.size):
        ni[D[i]] =  np.ravel_multi_index((xi[i],yi[i]),(M.nx[d],M.ny[d]),order='F')

    resultMesh = np.meshgrid(M.lons,M.lats) #meshgrid of coordinates for model resolution (createrawmatrix)
    resultIndex = np.arange(0,M.nc[d]) #indexation identical to model (createrawmatrix)
    resultCoast = np.zeros(M.nc[d]) #array for summing land/sea/coast data per index
    resultNorm = np.zeros(M.nc[d]) #amount of cells from input grid per cell for output grid
    
    for i in np.arange(coastalGrid.size) :
        if ni[i] > -9900:
            if coastalGrid[i] < 9 : #sea
                resultCoast[ni[i]] += coastalGrid[i]
                resultNorm[ni[i]] += 4
            else : #land
                resultNorm[ni[i]] += 0

    resultNorm[np.where(resultNorm == 0)] = 100 #edge indices (prime meridian, poles) don't have coastalGrid cells in them (so will remain 0) but also do not alter results. The same goes for open ocean. This statements prevents division by zero
    resultCoast = resultCoast / resultNorm

    f = open("adj_coast_dx_" + str((M.dd[d] * 100).astype(int)) + "_region_" + region + ".csv",'w') #also clear it if it exists
    f.write('index,lat,lon,numberCoasts\n')

    for i in np.arange(resultIndex.size) :
        x, y = np.unravel_index(resultIndex[i],(M.nx[d],M.ny[d]),order='F')
        f.write(str(resultIndex[i]) + ',' + str(M.lats[d][y]) + ',' + str(M.lons[d][x]) + ',' + str(resultCoast[i]) + '\n')
    f.close()