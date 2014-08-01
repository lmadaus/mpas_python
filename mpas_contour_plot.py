#!/usr/bin/env python

"""
Script to do simple line contours of MPAS model output using netCDF4, Basemap
and matplotlib.pyplot
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def main():
    mpas_file = './output.nc'
    ncfile = Dataset(mpas_file,'r')

    contour_hplane(ncfile,var='sst')


def contour_hplane(ncfile, var, plevel=None):
    """ Function to contour plot variable "var" on the optional level "plevel"
    """
    from scipy.interpolate import griddata


    # Check to see that this variable exists
    if var not in ncfile.variables.keys():
        print "ERROR: Cannot find variable:", var
        exit(1)

    # For a simple contour plot, easiest thing currently is to interpolate to a
    # regular grid at the highest "equivalent" resolution of the MPAS grid. Sort
    # that out here
    lon_interp, lat_interp = regular_interpolate(ncfile)


    # Set up the plot
    m, xgrid, ygrid = make_map(lon_interp, lat_interp, ncfile)


    # Load the field and get its dimensions
    field = ncfile.variables[var]
    
    field_dims = field.dimensions
    field_shape = field.shape


    # Main deal with MPAS data is that u,v are defined on edges of cells while
    # other variables are based on the cells.  Need to sort out what our
    # defining unit is
    if 'nCells' in field_dims:
        lats = ncfile.variables['latCell'][:]
        lons = ncfile.variables['lonCell'][:]
    elif 'nEdges' in field_dims:
        lats = ncfile.variables['latEdge'][:]
        lons = ncfile.variables['lonEdge'][:]
    else:
        print "Unable to find lat/lon data for var:", var
        exit(1)

    cell_coords = np.vstack((lons,lats))

    # Coordinates are in radians; convert to degrees
    cell_coords *= (180./np.pi)


    # If we have nVertLevels, figure out which vertical level we want
    vert_lev = 0

    # If we have times, loop through all times
    if 'Time' not in field_dims:
        timeloop = [0]
    else:
        timeloop = range(field_shape[field_dims.index('Time')])
    for time in timeloop:
        # Get the slice of "field" that we want
        if 'Time' in field_dims:
            try:
                curfield = field[time, :, vert_lev]
            except:
                curfield = field[time,:]
        else:
            try:
                curfield = field[:, vert_lev]
            except:
                curfield = field[:]

        # We now have lats, lons and values
        # For a simple contour plot, we just interpolate this to a regular grid
        # at the highest resolution
        gridded = griddata(cell_coords.T, curfield, (lon_interp,lat_interp),
                           method='nearest')


        # Now plot
        fig = plt.figure()
        plt.contour(xgrid, ygrid, gridded)
        m.drawcoastlines()
        plt.title('%s   Time: %03d' % (var, time))
        plt.show()
        plt.close()






def regular_interpolate(mpasfile, dx=0.25):
    """ Create a global, rectangular grid based on the specified resolution in
    degrees """
    nlons = 360.0 / dx + 1.0
    nlats = 180 / dx + 1.0
    lat_list = np.linspace(-90,90,nlats)
    lon_list = np.linspace(0,360,nlons)
    lon_grid, lat_grid = np.meshgrid(lon_list,lat_list)
    return lon_grid, lat_grid


def make_map(lons,lats,mpasfile):
    """ Create a basemap object and projected coordinates for a global map """
    m = Basemap(projection='cyl',llcrnrlat=-90, urcrnrlat=90, llcrnrlon=0,
                urcrnrlon=360, resolution='c')
    xgrid, ygrid = m(lons,lats)
    return m, xgrid, ygrid
    





if __name__ == '__main__':
    main()

