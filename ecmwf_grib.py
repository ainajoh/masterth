from scipy.stats.kde import gaussian_kde #test
import scipy.ndimage as ndimage
import argparse

import pygrib                           #To read grib files. 
import matplotlib
import matplotlib.pyplot as plt     
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap #To make map
import numpy as np
import math
import spharm  # Conda has it. Converts from spectral grid to gaussian grid, http://code.google.com/p/pyspharm
import warnings



#cdo

def sh( values ):
    """
    ------------------------------------------------------------------
    Handling spectral grid = spherical harmonics.
    NOTE: This is not needed for my new data file since it only contains gg grid.
    But this function is nice to have incase I get a sh grid in somewhere. 
     
    Instead of downloading in totaly raw format(which will contain both sh and gg grid), 
            it can be downloaded as already converted gg grid from ECMWF. Else you also have
            the option in using cdo on the server to convert it one by one to gg.. So this 
            function is now a third backup kind of.. But it is absolutely proven to work well,
            its just time consuming. 
    
    
    uses:       spharm.spectogrd() to transform grid to gausian grid(gg)
    
    returns:    lats, lons, data
    
    Source:     pygrib test files on github
    
    Issues:  
                nlats = 180 in original, but this dosent work for me. 
                nlats = 1280 = size.values, works.
    
    -------------------------------------------------------------------
    """
    # ECMWF normalizes the spherical harmonic coeffs differently than NCEP.
    # (m=0,n=0 is global mean, instead of sqrt(2)/2 times global mean)
    fld = 2.*values/np.sqrt(2.)
    
    #------SPLITTING IT UP IN AN IMAGARY AND REAL PART--------
    fldr = fld[ 0::2 ]                  #annenhver verdi fra 0
    fldi = fld[ 1::2 ]                  #annenhver verdi fra 1
    fldn = np.zeros( fldr.shape, 'F' )  #blir halvparten så stor som orginale fld
    fldn.real = fldr                    #legges da til i fldn vectoren
    fldn.imag = fldi
    #----------------------------------------------------------
     
    nlons = 360                         #Have a feeling it probably is number of values like grid val
    nlats = 1280                        #web sais it shourld be 180.. wellwell, seems to work
    s = spharm.Spharmt( nlons, nlats )  
    
    data = s.spectogrd( fldn )           #Hvis nlats = 180, så feiler denne delen pga hvordan formelen fungerer..
    
    lons = ( 360./nlons ) * np.arange( nlons )
    lats = 90.-( 180./( nlats - 1 ) ) * np.arange( nlats )
    lons, lats = np.meshgrid( lons, lats )
    
    #stack grids side-by-side (in longitiudinal direction), so
    # any range of longitudes (between -360 and 360) may be plotted on a world map.
    lons = np.concatenate(( lons - 360, lons ), 1 )
    lats = np.concatenate(( lats, lats ), 1 )
    data = np.concatenate(( data, data ), 1 )
    
    return lats, lons, data
   
def get_data( obj, prm, lev, date, timelevel=0 ):
    """
    ------------------------------------------------------------------
    fatching the data and finds it latitude and longditude. 
    parameters: - obj       = object of the open file. 
                - prm       = parameter name that we want to find
                - lev       = level. States which level we want to take or data from, feks 1000 hpa level.
                - date      = which date to get
                - timelevel = which time we want. 0: time 0000UTC, 1: time 1200UTC, 2: time 1800UTC
    
    uses:       Pygrib for reading grib file. 
    
    returns:    - lat, lon  = lat,lon
                - data      = values of the parameter. 
    issues: Dont understand why I need to do lon-180 for it to use the whole map. 
            Maybe this can be fixed when defining the drawing of the map or m object..?
            For now lon-180 works fine, but think it can be avoided. 
        
    -------------------------------------------------------------------
    """
    
    parameter = obj( name = prm, level = lev, dataDate = date )[ timelevel ]
    print( parameter.dataDate )
    
    #-----Checking grit type----------------------------------------------
    if parameter.gridType == "sh":
        lat, lon, data = sh( parameter.values )
    elif parameter.gridType == "reduced_gg":
        lat, lon = parameter.latlons() #very easy implementastion with a gg
        lon = lon - 180.        #else it only draws on half the map
        data = parameter.values
    elif parameter.gridType == "regular_gg":
        lat, lon = parameter.latlons() #very easy implementastion with a gg
        lon = lon - 180.        #else it only draws on half the map
        data = parameter.values
    else: 
        print ( parameter.gridType )
           
    return lat, lon, data

def plot_contourf( m, lat, lon, data, C, contur_val=None ):
    """
    ------------------------------------------------------------------
    Plot the data in a filed contour with color C
    Plots the colorbar
    
    parameters: - m         = mab "object" from basemap.
                - lat       = latitude
                - lon       = longditude
                - data      = the data values of our parameter
                - C         = which Color to use.
                -contur_val = int:How many contours to plot.
    Alternatice colors: -#CS = m.contourf(x,y,data,15,cmap=plt.cm.jet)
                        - See version 2.5

     
    -------------------------------------------------------------------
    """
   
    x, y = m( lon,lat )
    min = data.min()
    if contur_val.any():
        CS = m.contourf( x, y, data, contur_val, colors = C, vmin = 264 , vmax=384 )
    else:
        CS = m.contourf( x, y, data, colors = C, vmin = 264 , vmax=384 )
    
    plt.colorbar( drawedges = True )            # draw colorbar       
     
def plot_contour( m, lat, lon, data, contour, clr = 'k' ):
    """
    ------------------------------------------------------------------
    Plot the data in a contour with color xlr
    
    parameters: - m         = mab "object" from basemap.
                - lat       = latitude
                - lon       = longditude
                - data      = the data values of our parameter
                - contour
                - clr       = Color to use.
    Issues:     Had problems finding which contour values I should pot to make it nice. 
                I did try different contour values for different levels, but didnt get any better
                Ended up with taking the "average" over the vertical levels and then filter it horizontally(see def DT),
                but its very time consuming, but looks better.
                Still not great though.. 
    
    -------------------------------------------------------------------
    """
    
    x, y = m( lon,lat )
    
    CS = m.contour( x, y, data, contour , linewidths = 0.8, colors = clr )
    
def map_area( m ):
    """
    ------------------------------------------------------------------
    Draws the background map based on m
    
    parameters: - m    = mab "object" from basemap.
                       It describes witch projection and latitudes to use for everything 
    
    -------------------------------------------------------------------
    """
    
    
    m.drawcoastlines( linewidth = 1.5, linestyle = 'solid', color = [ 75./255., 75/255., 75/255. ] )	
    # ------draw parallels----------------
    circles = np.arange( -90., 90. + 30, 30. ) #delat = 30.
    m.drawparallels( circles, labels = [ 1, 0, 0, 0 ] )
    
    # -------draw meridians---------------
    meridians = np.arange( 0., 360, 60. ) #delon = 60.
    m.drawmeridians( meridians, labels = [ 0, 0, 0, 1 ] )    

def plot_wind_bar( m, lat, lon, u, v ):
    """
    ------------------------------------------------------------------
    Plot wind bars for u and v
    
    parameters: - m    = mab "object" from basemap.
                       It describes witch projection and latitudes to use for everything 
                - lat
                - lon
                - u
                - v
    source:     #http://basemaptutorial.readthedocs.io/en/latest/plotting_data.html
    -------------------------------------------------------------------
    """

    x, y = m( lon, lat )
    
    #---Have to limit number of windbarbs plotted-----
    yy = np.arange( 0, y.shape[ 0 ], 30 ) #skips over every 100th value
    xx = np.arange( 0, x.shape[ 1 ], 30 ) #skips over every 100th value
    points = np.meshgrid( yy, xx )
    m.barbs( x[points], y[points], u[points], v[points], length = 5.5, pivot='middle', linewidth = 1., barbcolor = '#333333' )
   
def DT(time_lvl = 0, date = 160924 ):
    """
    ------------------------------------------------------------------
    Generates the Dynamical Tropopause(DT) map
    
    Parameter:  -C: Found from measuring the color bar at the source.
    
    Uses: Functions: - map_area(m)
                     - get_data( gfile, prm, lev, timelevel=0 )
                     - plot_contourf(m, lat, lon, data, C)
                     - plot_contour(m, lat, lon, data, C)
    
    Issues:   Tried using another basemap, one that is more curved, which is what I want, but the vorticity gets   
                really weird in the North. I dont understand why, so I cant use it yet..
    
    Source:  http://www.atmos.albany.edu/student/abentley/realtime.html
    
    -------------------------------------------------------------------
    """
    
    #-------Customised color in RGB ------------
    C = [[232,232,230],#grey
        [203,203,203], #grey
        [161,161,161], #grey
        [130,130,130], #grey
        [149,53,229],  #lillac, 39	64	197	149,53,229
        [39,64,197],   #blue dark,7,67,194
        [15,110,229],  #blue
        [80,149,240],  #blue
        [74,192,243],  #blue
        [152,219,248], #blue
        [183,237,247], #blue
        [251,217,198], #redish
        [255,197,166], #redish
        [255,172,164], #redish
        [253,139,142], #redish
        [253,101,105], #redish
        [255,66,74],   #redish
        [238,13,28],   #red
        [214,78,166],  #pink
        [214,102,201], 
        [217,155,210],
        [216,181,211]]
    C = np.array( C )
    C = np.divide( C, 255. )  # RGB has to be between 0 and 1 in python
    #-----------------------------------------------------------
    
    fig = plt.figure()
    
    
    #-----Setting our map area and projection of interest-------
    m = Basemap( llcrnrlon = -90., llcrnrlat = 0., urcrnrlon = 50., urcrnrlat=70.,\
               resolution = 'l', area_thresh = 10000., projection = 'merc' )
    #m = Basemap(width=11500000,height=8500000,resolution='l',projection='eqdc',\
    #            lat_1=07.,lat_2=40,lat_0=44,lon_0=-30.)
    #m = Basemap(width=190000,height=2200000,resolution='l', projection='tmerc',lon_0=-30,lat_0=44)
    
    map_area( m )            # ploting background
    path = "gribs/"
    file = path +"DT_var.grib"
    obj = pygrib.open( file )
    
    #-FETCHING ALL THE VALUES----------------------------------------
    #-----Potential temperature---------------------------------------
    lat, lon, data = get_data( obj,'Potential temperature', 2000, date,  timelevel = time_lvl )
    contour_val = np.linspace( 264, 384, 22 ) #contours for potential tempeature
    plot_contourf( m, lat, lon, data, C, contour_val )
    
    #-----Relative vorticity, diff level------------------------------
    contour=[ 2.8E-4, 3.5E-4, 4.5E-4, 6.5E-4, 7.E-4, 7.5E-4, 8.E-4 ] #1.5E-4,2.5E-4]#
    lat, lon, data925 = get_data( obj, 'Vorticity (relative)', 925, date, timelevel = time_lvl )
    lat, lon, data900 = get_data( obj, 'Vorticity (relative)', 900, date, timelevel = time_lvl )
    lat, lon, data850 = get_data( obj, 'Vorticity (relative)', 850, date, timelevel = time_lvl )
    
    #->--->---->--mean value over height and filtering----------------
    data = np.sqrt( data900**2 + 2*data850**2 + data925**2 ) #Vertical "average", weightet values at 850hpa double.
    footprint = np.array([[0,0,0,1,1,1,1,0,0,0],             #footprint=np.ones((3,10))
                          [0,0,1,1,1,2,1,1,0,0],
                          [1,1,1,2,2,1,2,1,1,1],
                          [0,1,1,1,1,2,1,1,1,0],
                          [0,0,1,1,1,1,1,1,0,0]])
    
    data = ndimage.generic_filter( data, np.mean, footprint = footprint, mode='wrap' )
    plot_contour( m, lat,lon, data,contour, clr = 'k' )
    
    #-----Wind barbs----------------------------------------------------
    lat, lon, data_u = get_data( obj , 'U component of wind', 2000, date, timelevel = time_lvl )
    lat, lon, data_v = get_data( obj , 'V component of wind', 2000, date, timelevel = time_lvl  )
    plot_wind_bar( m, lat, lon, data_u, data_v )
    #-----------------------------------------------
    #-----------------------------------------------
    
    
    
    #-SAVE AND CLOSE----------------------------------------------------
    #------------------------------------------------------------------
    obj.close()
    if time_lvl == 0:
        t = "0000"
    elif time_lvl == 1:
        t = "1200"
    elif time_lvl == 2:
        t = "1800"
    else: 
        t = "t_not_set"
    
    fig_name = "DT/DT_" + str( date ) + "_" + str( t )+ ".TIFF"   
    
    ax = plt.gca( )
    plt.rc( 'font', size = 6 )
    fig.set_size_inches( 12.80, 7.15 )
    
    fig.savefig( fig_name, dpi = 600 )
    plt.close( )
    #plt.show()
    #--------------------------
    #----------------------------
 
def user_interface():
    parser = argparse.ArgumentParser(description='Process some integers.')

    t = []          #timelevel
    
    parser.add_argument( "--date", type = int,
        choices= [ 20160920,20160921,20160922,20160923,20160924,20160925,20160926,20160927, 20160928, 20160929, 20160930, 20161001],
        help = "the date you want. " )
    parser.add_argument( "time", type = str,
        help = "the times you want. ex 160924, 160926, 160922, 160930", 
        default = "all" )
        
    args = parser.parse_args( )
    
    date = [ 20160920, 20160921, 20160922, 20160923, 20160924, 20160925, 20160926, 20160927, 20160928, 20160929, 20160930, 20161001 ]
    if args.date:
        date = [args.date]
    
    if args.time == "all":
        t = [ 0, 1, 2 ]
    elif args.time == "0000" or args.time == "00" or args.time == "0":
        t = [ 0 ]
    elif args.time == "1200" or args.time == "12" or args.time == "1":
        t = [ 1 ]
    elif args.time == "1800" or args.time == "18" or args.time == "2":
        t = [ 2 ]
    
    print ("\n-------------------------------------------------------")
    print ("you chosed, time: "+args.time+"corresponding to timelevel: "+str(t))
    print ("you chosed, date: "+str(args.date) ) 
    print ("---------------------------------------------------------")
    print("Creating a dynamic tropopause....\n ...")
    
    #------Setting time and date for plots----------------------
    for d in date:
        print ("one day made")
        for n in t:
            DT(n, d)
            print("one time made")
            
    #------------------------------------------------------------    

#matplotlib warns me that I am using a function(hold = on) that is not supported in updatet version of matplotlib. 
#Though I am not using it, so probably a package that uses it. Might be a problem if I dont update all the packages correctly. 
warnings.filterwarnings("ignore",category=matplotlib.mplDeprecation) 
user_interface()  
    
    
            