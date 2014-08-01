mpas_python
===========

A few Python scripts for plotting output from the MPAS weather model
All scripts are in the scripts subdirectory.


--> mpas_contour_plot.py --> This is a "cheap" way to do a contour plot and could be adapted to do a contourf plot as well.  It simply sets up a rectangular lat-lon grid at a user-specified resolution and interpolates the MPAS model values at cells (or edges for winds) to that rectangular grid. It's a simple matter to plot after that.

--> mpas_pcolor_plot.py --> This makes more expensive, but rather pleasant looking pcolor-type plots.  A large part of the script invoves looping through all of the cells and defining a collection of patches in Basemap-projected coordinates (takes several minutes for fine meshes with large numbers of cells).  This collection is archived to a Pickled file so that it can be read in for future plots much more quickly.  The patch collection is then filled with colors scaled to the data to produce a "blended" pcolor-type map.
