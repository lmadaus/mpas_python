#!/usr/bin/env python

"""
Script to do  pcolor-style MPAS model output using netCDF4, Basemap
and matplotlib.pyplot
"""
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

def main():
    mpas_file = '../output.nc'
    ncfile = Dataset(mpas_file,'r')

    pcolor_fill(ncfile,var='temperature_500hPa')


def pcolor_fill(ncfile, var, plevel=None):
    """ Function to do a pcolor-mesh-type plot of "var" on the optional level "plevel"
    """
    import cPickle
    # Check to see that this variable exists
    if var not in ncfile.variables.keys():
        print "ERROR: Cannot find variable:", var
        exit(1)


    # Set up the plot
    m = make_map(ncfile, lons=None, lats=None)


    # Load the field and get its dimensions
    field = ncfile.variables[var]
    
    field_dims = field.dimensions
    field_shape = field.shape

    # Going to create a collection of patches in matplotlib to enable speedy
    # plotting.  Call this function to do that if there is no archived patch
    # collection pickle file already.  It takes several minutes to generate the initial
    # patch collection (could probably be parallelized using multiprocessing or
    # optimized some other way; will leave that to future), so be warned.  I
    # like to generate the collection once and call the saved file every
    # time...makes it go a  lot faster.
    try:
        patchfile = open('mpas_paths.pckl','rb')
        p = cPickle.load(patchfile)
        patchfile.close()
    except:
        p = mpas_grid_to_patches(ncfile, m)



    # Main deal with MPAS data is that u,v are defined on edges of cells while
    # other variables are based on the cells.  Need to sort out what our
    # defining unit is
    if 'nEdges' in field_dims:
        # May just convert this back to a gridded contour plot
        # Not implemented for now
        print "Edge data not implemented for pcolor-type plot"
        exit(1)



    # If we have nVertLevels, figure out which vertical level we want
    vert_lev = plevel 

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


        # Now plot
        fig = plt.figure(figsize=(12,9))
        ax = plt.gca()
        # Set the plotting array of our patch collection to be curfield
        p.set_array(curfield)
        # Set this to 'none' to not plot any of the edges
        p.set_edgecolors('none')
        # Set this to False to not delineate the cells at all (even setting
        # edgecolors to 'none' will show faint cell delineations). The False
        # setting here gives a more "blended" look.
        p.set_antialiaseds(False)
        # Set the colormap here
        p.set_cmap(matplotlib.cm.spectral)
        # Can normalize the colormap here if desired
        #p.set_norm()
        # Add the collection with its colorized properties to the plot
        ax.add_collection(p)
        m.drawcoastlines()
        plt.colorbar(p)
        plt.title('%s   Time: %03d' % (var, time))
        plt.show()
        plt.close()



def mpas_grid_to_patches(mpasfile, bmap):
    """ Function to create a collection of patches in Basemap plotting
    coordinates that define the cells of the MPAS domain """
    print "Defining Path Collection on MPAS Grid"


    # Get the number of cells
    nCells = len(mpasfile.dimensions['nCells'])
    nEdgesOnCell = mpasfile.variables['nEdgesOnCell']
    verticesOnCell = mpasfile.variables['verticesOnCell']
    latVertex = mpasfile.variables['latVertex']
    lonvertex = mpasfile.variables['lonVertex']

    patches = [None] * nCells
    print "    Total num cells:", nCells
    # Need a collection of lats and lons for each vertex of each cell
    for c in xrange(nCells):
        if c % 5000 == 0:
            print "        On:", c
        # Each collection of vertices has a length of maxEdges.  Need to figure
        # out how many vertices are ACTUALLY on the cell, as the rest is just
        # padded with junk data
        # These are the vertex indices
        cell_verts = verticesOnCell[c,:nEdgesOnCell[c]]
        # Add on the final point to close
        cell_verts = np.append(cell_verts,cell_verts[0:1])
        # Subtract one
        cell_verts -= 1
        #cell_verts = mpasfile.variables['indexToVertexID'][cell_vert_index]
        # Get the latitudes and longitudes of these and convert to degrees
        vert_lats = mpasfile.variables['latVertex'][cell_verts] * 180./np.pi
        vert_lons = mpasfile.variables['lonVertex'][cell_verts] * 180./np.pi
        # Check for overlap of date line
        diff_lon = np.subtract(vert_lons, vert_lons[0])
        vert_lons[diff_lon > 180.0] = vert_lons[diff_lon > 180.0] - 360.0
        vert_lons[diff_lon < -180.0] = vert_lons[diff_lon < -180.0] + 360.0
        # Convert to projected coordinates
        vert_x, vert_y = bmap(vert_lons, vert_lats)
        coords = np.vstack((vert_x, vert_y))
        # Now create a path for this
        # Codes follow same format
        cell_codes = np.ones(nEdgesOnCell[c]+1) * mpath.Path.LINETO
        cell_codes[0] = mpath.Path.MOVETO
        cell_codes[-1] = mpath.Path.CLOSEPOLY
        cell_path = mpath.Path(coords.T, codes=cell_codes, closed=True, readonly=True)
        # Convert to a Patch and add to the list
        patches[c] = mpatches.PathPatch(cell_path)

    # Now crate a patch collection
    p = matplotlib.collections.PatchCollection(patches)
    # Archive for future use
    print "    Archiving paths..."
    import cPickle
    outfile = open('mpas_paths.pckl','wb')
    cPickle.dump(p, outfile)
    outfile.close()

    return p




def make_map(mpasfile, lons, lats):
    """ Create a basemap object and projected coordinates for a global map """
    # MPAS cell lats and lons are projected using a cylindrical projection;
    # other projections are not recommended...
    m = Basemap(projection='cyl',llcrnrlat=-90, urcrnrlat=90, llcrnrlon=0,
                urcrnrlon=360, resolution='c')
    #m = Basemap(projection='cyl',llcrnrlat=20, urcrnrlat=70, llcrnrlon=180,
    #           urcrnrlon=340, resolution='c')


    if lons == None or lats == None:
        # Just return the map
        return m
    else:
        xgrid, ygrid = m(lons,lats)
        return m, xgrid, ygrid
    





if __name__ == '__main__':
    main()

