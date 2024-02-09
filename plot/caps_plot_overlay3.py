from __future__ import print_function
import traceback
import sys
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import numpy as np
import cartopy.crs as crs
import cartopy.feature as cfeature
import matplotlib.colors as mc
from matplotlib.colors import ListedColormap
from matplotlib.patches import Rectangle
from scipy.ndimage import gaussian_filter
from scipy.ndimage import maximum_filter, minimum_filter

import netCDF4
import os
from eccodes import *
from datetime import datetime, timedelta
import argparse 
import time
import multiprocessing

VERBOSE = 1  # verbose error reporting

def example():
    def parse_case(case_str):
        return datetime.strptime(case_str, "%Y%m%d%H")

    ap = argparse.ArgumentParser(description="Plotting CAPS Ensemble Results")
    ap.add_argument('--case', dest='case', type=parse_case, required=True, help="Date/time to plot in YYYYMMDDHH format")
    ap.add_argument('--fhr',  dest='fhr',  type=int, required=True, help="forecast lead time (hours)")
    ap.add_argument('--opt',  dest='opt',  type=int, required=True, help="Alignment option")
    ap.add_argument('--mask', dest='mask', type=int, required=False, default=0, help="obs mask")

    args = ap.parse_args()
    npass = 2

    #finput_name = '/archive-temp/2023_HMT_summer/{case:%Y%m%d%H}/caps_ens_{case:%Y%m%d%H}_f{fhr:03d}.nc'.format(case=args.case, fhr=args.fhr)
    if (args.opt == 1):
      finput_name = '/non/clee/ens/data/{case:%Y%m%d}/caps_ens_{case:%Y%m%d%H}_f{fhr:03d}.nc'.format(case=args.case, fhr=args.fhr)
    elif (args.opt == 2):
      finput_name = '/non/clee/ens/data2/{case:%Y%m%d}/caps_ens_{case:%Y%m%d%H}_f{fhr:03d}.nc'.format(case=args.case, fhr=args.fhr)
 
    f  = netCDF4.Dataset(finput_name) 

    value1 = f.variables['ensmean'][:]
    value2 = f.variables['enspm'][:]
    value3 = f.variables['enslpm'][:]
    value4 = f.variables['ensshfmean'][:]
    value5 = f.variables['ensshfpm'][:]
    value6 = f.variables['ensshflpm'][:]
    lats = f.variables['lat'][:]
    lons = f.variables['lon'][:]
    f.close()

    value1_ext = maximum_filter(value1, 500, mode='nearest')
    value2_ext = maximum_filter(value2, 500, mode='nearest')
    value3_ext = maximum_filter(value3, 500, mode='nearest')
    value4_ext = []
    value5_ext = []
    value6_ext = []

    value1_ext2 = maximum_filter(value1, 200, mode='nearest')
    value2_ext2 = maximum_filter(value2, 200, mode='nearest')
    value3_ext2 = maximum_filter(value3, 200, mode='nearest')
    value4_ext2 = []
    value5_ext2 = []
    value6_ext2 = []

    mxy1, mxx1 = np.where((value1_ext == value1) & (value1_ext >= 0.5))
    mxy2, mxx2 = np.where((value2_ext == value2) & (value2_ext >= 0.5))
    mxy3, mxx3 = np.where((value3_ext == value3) & (value3_ext >= 0.5))

    dmxy1, dmxx1 = np.where((value1_ext2 == value1) & (value1_ext2 >= 0.5))
    dmxy2, dmxx2 = np.where((value2_ext2 == value2) & (value2_ext2 >= 0.5))
    dmxy3, dmxx3 = np.where((value3_ext2 == value3) & (value3_ext2 >= 0.5))

    mxy4 = [] 
    mxx4 = []
    mxy5 = [] 
    mxx5 = []
    mxy6 = [] 
    mxx6 = []

    dmxy4 = [] 
    dmxx4 = []
    dmxy5 = [] 
    dmxx5 = []
    dmxy6 = [] 
    dmxx6 = []

    for x in range(npass):
        value4_ext.append([])
        value5_ext.append([])
        value6_ext.append([])

        value4_ext2.append([])
        value5_ext2.append([])
        value6_ext2.append([])

        mxy4.append([])
        mxx4.append([])
        mxy5.append([])
        mxx5.append([])
        mxy6.append([])
        mxx6.append([])

        dmxy4.append([])
        dmxx4.append([])
        dmxy5.append([])
        dmxx5.append([])
        dmxy6.append([])
        dmxx6.append([])

        value4_ext[x] = maximum_filter(value4[x,:,:], 500, mode='nearest')
        value5_ext[x] = maximum_filter(value5[x,:,:], 500, mode='nearest')
        value6_ext[x] = maximum_filter(value6[x,:,:], 500, mode='nearest')

        value4_ext2[x] = maximum_filter(value4[x,:,:], 200, mode='nearest')
        value5_ext2[x] = maximum_filter(value5[x,:,:], 200, mode='nearest')
        value6_ext2[x] = maximum_filter(value6[x,:,:], 200, mode='nearest')

        mxy4[x], mxx4[x] = np.where((value4_ext[x] == value4[x,:,:]) & (value4_ext[x] >= 0.5))
        mxy5[x], mxx5[x] = np.where((value5_ext[x] == value5[x,:,:]) & (value5_ext[x] >= 0.5))
        mxy6[x], mxx6[x] = np.where((value6_ext[x] == value6[x,:,:]) & (value6_ext[x] >= 0.5))

        dmxy4[x], dmxx4[x] = np.where((value4_ext2[x] == value4[x,:,:]) & (value4_ext2[x] >= 0.5))
        dmxy5[x], dmxx5[x] = np.where((value5_ext2[x] == value5[x,:,:]) & (value5_ext2[x] >= 0.5))
        dmxy6[x], dmxx6[x] = np.where((value6_ext2[x] == value6[x,:,:]) & (value6_ext2[x] >= 0.5))


    fobs_name = '/home/clee/METplus/Tutorial/data/stage4_conus/st4_conus.{date:%Y%m%d%H}.06h.grb2'.format(date=args.case + timedelta(hours=args.fhr))

    index_keys = ["shortName"]
    iid = codes_index_new_from_file(fobs_name, index_keys)

    codes_index_select(iid,'shortName','tp')
    
    igrib = codes_new_from_index(iid)
    nx_grib = codes_get(igrib,'Nx')
    ny_grib = codes_get(igrib,'Ny')
    missing = codes_get(igrib,'missingValue')
    value_temp = codes_get_array(igrib,'values',float)
    lats_grib = np.reshape(codes_get_array(igrib,'latitudes',float),(ny_grib,nx_grib))
    lons_grib = np.reshape(codes_get_array(igrib,'longitudes',float),(ny_grib,nx_grib))
    value_grib = np.reshape(value_temp,(ny_grib,nx_grib))
    mask_grib = np.ma.getmaskarray(value_grib)
    codes_release(igrib)

    fmask = netCDF4.Dataset("./mask_conus_b50km_stage4.nc")
    mask = fmask.variables['conus'][:] 
    fmask.close()

    if (args.mask > 0):
        for y in range(ny_grib):
            for x in range(nx_grib):
                if ((value_grib[y,x] == missing) or (mask[y,x] == 0)):
                    value_grib[y,x] = 0.
                    mask_grib[y,x] = True 
    else:
        for y in range(ny_grib):
            for x in range(nx_grib):
                if (value_grib[y,x] == missing):
                    value_grib[y,x] = 0.
                    mask_grib[y,x] = True 

    value_grib = np.ma.masked_array(value_grib, mask = mask_grib)

    value_grib_ext = maximum_filter(value_grib, 250, mode='nearest')
    mxy_grib, mxx_grib = np.where((value_grib_ext == value_grib) & (value_grib_ext >= 0.5))

    value_grib_ext2 = maximum_filter(value_grib, 175, mode='nearest')
    dmxy_grib, dmxx_grib = np.where((value_grib_ext2 == value_grib) & (value_grib_ext2 >= 0.5))



    levels = [0, 0.25, 0.5, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
    levels_str = ['0', '0.25', '0.5', '1', '2', '4', '6', '8', '10', '15', '20', '25', '30', '35', '40', '45', '50', '60', '70', '80', '90', '100']
    cmap = ListedColormap([(1,1,1), (0,236/255,236/255), (0,200/255,240/255), (0,160/255,255/255), (0,60/255,255/255), (0,255/255,0), (0,220/255,0),
                           (0,190/255,0), (0,141/255,0), (255/255,255/255,0), (240/255,210/255,0), (231/255,180/255,0), (200/255,120/255,0),
                           (255/255,160/255,160/255), (255/255,60/255,60/255), (230/255,0,0), (180/255,0,0), (255/255,0,255/255),
                           (217/255,0,217/255), (164/255,0,164/255), (120/255,0,120/255)])
    norm = mc.BoundaryNorm(levels, cmap.N)
    data_crs = crs.PlateCarree()

    sectors = ['conus', 'mw', 'ne', 'nw', 'sw', 'se', 'cp', 'ma', 'np', 'sp']
    extents = [[-120,-70,25,50],[-105,-80,32.5,45],[-95,-70,37.5,47.5],[-125,-100,39,49],[-122.5,-93.5,28,42.5],[-105,-77,24,38],[-110,-87,31.5,44],
               [-95,-73,34,44],[-110,-86,38,50],[-110,-86,24,39]]

    processes = []
    pool = multiprocessing.Pool(processes=len(sectors))

    for k in range(len(sectors)):
        p = pool.apply_async(plot_chart, (k, args, npass, sectors, extents, levels, levels_str, lons, lats, lons_grib, lats_grib, 
               value1, value2, value3, value4, value5, value6, value_grib,
               mxx1, mxy1, mxx2, mxy2, mxx3, mxy3, mxx4, mxy4, mxx5, mxy5, mxx6, mxy6, mxx_grib, mxy_grib,
               dmxx1, dmxy1, dmxx2, dmxy2, dmxx3, dmxy3, dmxx4, dmxy4, dmxx5, dmxy5, dmxx6, dmxy6, dmxx_grib, dmxy_grib, cmap, norm, data_crs))
        processes.append(p)
 
    for process in processes:
        process.get()

    pool.close()
    pool.join()
 
def main():
    try:
        example()
    except:

    #except CodesInternalError as err:
        if VERBOSE:
            traceback.print_exc(file=sys.stderr)
        else:
            sys.stderr.write(err.msg + '\n')
 
        return 1

def plot_maxmin_points(ax, lon, lat, data, extrema, nsize, symbol, color='k',
                       plotValue=True, transform=None):
    """
    This function will find and plot relative maximum and minimum for a 2D grid. The function
    can be used to plot an H for maximum values (e.g., High pressure) and an L for minimum
    values (e.g., low pressue). It is best to used filetered data to obtain  a synoptic scale
    max/min value. The symbol text can be set to a string value and optionally the color of the
    symbol and any plotted value can be set with the parameter color
    lon = plotting longitude values (2D)
    lat = plotting latitude values (2D)
    data = 2D data that you wish to plot the max/min symbol placement
    extrema = Either a value of max for Maximum Values or min for Minimum Values
    nsize = Size of the grid box to filter the max and min values to plot a reasonable number
    symbol = String to be placed at location of max/min value
    color = String matplotlib colorname to plot the symbol (and numerica value, if plotted)
    plot_value = Boolean (True/False) of whether to plot the numeric value of max/min point
    The max/min symbol will be plotted on the current axes within the bounding frame
    (e.g., clip_on=True)
    """
    from scipy.ndimage import maximum_filter, minimum_filter

    if (extrema == 'max'):
        data_ext = maximum_filter(data, nsize, mode='nearest')
    elif (extrema == 'min'):
        data_ext = minimum_filter(data, nsize, mode='nearest')
    else:
        raise ValueError('Value for hilo must be either max or min')

    mxy, mxx = np.where(data_ext == data)

    for i in range(len(mxy)):
#        ax.text(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]], symbol, color=color, size=24,
#                clip_on=True, horizontalalignment='center', verticalalignment='center',
#                transform=transform)
        value = np.int64(data[mxy[i], mxx[i]])
        if (value > 0):
            text = ax.text(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]],
                '\n' + str(value),
                color=color, size=10, clip_on=True,
                horizontalalignment='center', verticalalignment='center', transform=transform)
            text.set_clip_box(ax.bbox)

def plot_maxmin_points2(ax, lon, lat, data, mxx, mxy, symbol, color='k',
                       plotValue=True, transform=None):

    for i in range(len(mxy)):
        text = ax.text(lon[mxy[i], mxx[i]], lat[mxy[i], mxx[i]],
                '\n' + str(np.int64(data[mxy[i], mxx[i]])),
                color=color, size=10, clip_on=True,
                horizontalalignment='center', verticalalignment='center', transform=transform)

        text.set_clip_box(ax.bbox)


def plot_chart(k, args, npass, sectors, extents, levels, levels_str, lons, lats, lons_grib, lats_grib, 
               value1, value2, value3, value4, value5, value6, value_grib,
               mxx1, mxy1, mxx2, mxy2, mxx3, mxy3, mxx4, mxy4, mxx5, mxy5, mxx6, mxy6, mxx_grib, mxy_grib,
               dmxx1, dmxy1, dmxx2, dmxy2, dmxx3, dmxy3, dmxx4, dmxy4, dmxx5, dmxy5, dmxx6, dmxy6, dmxx_grib, dmxy_grib, cmap, norm, data_crs):

    fig, axs = plt.subplots(3,4,figsize=(20,12.5), subplot_kw={"projection": crs.LambertConformal(central_longitude=262.5, central_latitude=38.5)})
    fig.subplots_adjust(hspace=0.15, wspace=0.1, right=0.95, top=0.9, left=0.05)
    if (args.opt == 1):
        fig.suptitle("VALID: {valid:%Y.%m.%d.%H}UTC (+{fhr:03d}H) / RUN: {case:%Y.%m.%d.%H}UTC / 6hr Precipitation / Align Method A".format(valid=args.case + timedelta(hours=args.fhr), 
          fhr=args.fhr, case=args.case), fontsize=16.)
    elif (args.opt == 2):
        fig.suptitle("VALID: {valid:%Y.%m.%d.%H}UTC (+{fhr:03d}H) / RUN: {case:%Y.%m.%d.%H}UTC / 6hr Precipitation / Align Method B".format(valid=args.case + timedelta(hours=args.fhr), 
          fhr=args.fhr, case=args.case), fontsize=16.)

    contours = axs[0,0].pcolormesh(lons_grib, lats_grib, value_grib[:,:], cmap=cmap, norm=norm, transform=data_crs)
    axs[0,0].set_title('Obs - Stage IV', fontsize=12.)
    axs[0,0].tick_params(left = False, right = False , labelleft = False ,
            labelbottom = False, bottom = False)
    axs[0,0].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[0,0].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    if (sectors[k] == 'conus'):
        axs[0,0].text(0, -0.065, "max: {:.2f}mm".format(np.max(value_grib)), transform=axs[0,0].transAxes, fontsize=12.)
        plot_maxmin_points2(axs[0,0], lons_grib, lats_grib, value_grib, mxx_grib, mxy_grib, symbol='', color='black',  transform=data_crs)
    else:
        plot_maxmin_points2(axs[0,0], lons_grib, lats_grib, value_grib, dmxx_grib, dmxy_grib, symbol='', color='black',  transform=data_crs)
    axs[0,0].set_extent(extents[k], data_crs)
    axs[0,0].set_facecolor((0.6,0.6,0.6))

    contours = axs[0,1].pcolormesh(lons, lats, value1, cmap=cmap, norm=norm, transform=data_crs)
    axs[0,1].set_title('Ens mean', fontsize=12.)
    axs[0,1].tick_params(left = False, right = False , labelleft = False ,
            labelbottom = False, bottom = False)
    axs[0,1].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[0,1].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    if (sectors[k] == 'conus'):
        axs[0,1].text(0, -0.065, "max: {:.2f}mm".format(np.max(value1[:,:])), transform=axs[0,1].transAxes, fontsize=12.)
        plot_maxmin_points2(axs[0,1], lons, lats, value1[:,:], mxx1, mxy1, symbol='', color='black',  transform=data_crs)
    else:
        plot_maxmin_points2(axs[0,1], lons, lats, value1[:,:], dmxx1, dmxy1, symbol='', color='black',  transform=data_crs)
    axs[0,1].set_extent(extents[k], data_crs)

    contours = axs[1,1].pcolormesh(lons, lats, value2, cmap=cmap, norm=norm, transform=data_crs)
    axs[1,1].set_title('Ens PM', fontsize=12.)
    axs[1,1].tick_params(left = False, right = False , labelleft = False ,
            labelbottom = False, bottom = False)
    axs[1,1].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[1,1].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    if (sectors[k] == 'conus'):
        axs[1,1].text(0, -0.065, "max: {:.2f}mm".format(np.max(value2[:,:])), transform=axs[1,1].transAxes, fontsize=12.)
        plot_maxmin_points2(axs[1,1], lons, lats, value2[:,:], mxx2, mxy2, symbol='', color='black',  transform=data_crs)
    else:
        plot_maxmin_points2(axs[1,1], lons, lats, value2[:,:], dmxx2, dmxy2, symbol='', color='black',  transform=data_crs)
    axs[1,1].set_extent(extents[k], data_crs)

    contours = axs[2,1].pcolormesh(lons, lats, value3, cmap=cmap, norm=norm, transform=data_crs)
    axs[2,1].set_title('Ens LPM', fontsize=12.)
    axs[2,1].tick_params(left = False, right = False , labelleft = False ,
            labelbottom = False, bottom = False)
    axs[2,1].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[2,1].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    if (sectors[k] == 'conus'):
        axs[2,1].text(0, -0.065, "max: {:.2f}mm".format(np.max(value3[:,:])), transform=axs[2,1].transAxes, fontsize=12.)
        plot_maxmin_points2(axs[2,1], lons, lats, value3[:,:], mxx3, mxy3, symbol='', color='black',  transform=data_crs)
    else:
        plot_maxmin_points2(axs[2,1], lons, lats, value3[:,:], dmxx3, dmxy3, symbol='', color='black',  transform=data_crs)
    axs[2,1].set_extent(extents[k], data_crs)

    for x in range(npass):
        contours = axs[0,x+2].pcolormesh(lons, lats, value4[x,:,:], cmap=cmap, norm=norm, transform=data_crs)
        axs[0,x+2].set_title("Ens SAM ({} pass)".format(x+1), fontsize=12.)
        axs[0,x+2].tick_params(left = False, right = False , labelleft = False ,
               labelbottom = False, bottom = False)
        axs[0,x+2].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
        axs[0,x+2].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
        if (sectors[k] == 'conus'):
            axs[0,x+2].text(0, -0.065, "max: {:.2f}mm".format(np.max(value4[x,:,:])), transform=axs[0,x+2].transAxes, fontsize=12.)
            plot_maxmin_points2(axs[0,x+2], lons, lats, value4[x,:,:], mxx4[x], mxy4[x], symbol='', color='black',  transform=data_crs)
        else:
            plot_maxmin_points2(axs[0,x+2], lons, lats, value4[x,:,:], dmxx4[x], dmxy4[x], symbol='', color='black',  transform=data_crs)
        axs[0,x+2].set_extent(extents[k], data_crs)

        contours = axs[1,x+2].pcolormesh(lons, lats, value5[x,:,:], cmap=cmap, norm=norm, transform=data_crs)
        axs[1,x+2].set_title("Ens SAM-PM ({} pass)".format(x+1), fontsize=12.)
        axs[1,x+2].tick_params(left = False, right = False , labelleft = False ,
               labelbottom = False, bottom = False)
        axs[1,x+2].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
        axs[1,x+2].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
        if (sectors[k] == 'conus'):
            axs[1,x+2].text(0, -0.065, "max: {:.2f}mm".format(np.max(value5[x,:,:])), transform=axs[1,x+2].transAxes, fontsize=12.)
            plot_maxmin_points2(axs[1,x+2], lons, lats, value5[x,:,:], mxx5[x], mxy5[x], symbol='', color='black',  transform=data_crs)
        else:
            plot_maxmin_points2(axs[1,x+2], lons, lats, value5[x,:,:], dmxx5[x], dmxy5[x], symbol='', color='black',  transform=data_crs)
        axs[1,x+2].set_extent(extents[k], data_crs)

        contours = axs[2,x+2].pcolormesh(lons, lats, value6[x,:,:], cmap=cmap, norm=norm, transform=data_crs)
        axs[2,x+2].set_title("Ens SAM-LPM ({} pass)".format(x+1), fontsize=12.)
        axs[2,x+2].tick_params(left = False, right = False , labelleft = False ,
               labelbottom = False, bottom = False)
        axs[2,x+2].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
        axs[2,x+2].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
        if (sectors[k] == 'conus'):
            axs[2,x+2].text(0, -0.065, "max: {:.2f}mm".format(np.max(value6[x,:,:])), transform=axs[2,x+2].transAxes, fontsize=12.)
            plot_maxmin_points2(axs[2,x+2], lons, lats, value6[x,:,:], mxx6[x], mxy6[x], symbol='', color='black',  transform=data_crs)
        else:
            plot_maxmin_points2(axs[2,x+2], lons, lats, value6[x,:,:], dmxx6[x], dmxy6[x], symbol='', color='black',  transform=data_crs)
        axs[2,x+2].set_extent(extents[k], data_crs)

    fig.delaxes(axs[1,0])
    fig.delaxes(axs[2,0])

    cax = fig.add_axes([0.05,0.05,0.9,0.02])
    cbar = fig.colorbar(contours, cax=cax, ticks=levels, orientation='horizontal')
    cbar.ax.set_xticklabels(levels_str)
    cbar.ax.set_xlabel('[mm]')

    fout_name = '/non/clee/ens/img/{case:%Y%m%d}/caps_ens_opt{opt:01d}_{case:%Y%m%d%H}_f{fhr:03d}_{sector}_overall.png'.format(opt=args.opt, case=args.case, fhr=args.fhr, sector=sectors[k])
    if os.path.exists(fout_name):
        os.unlink(fout_name)

    plt.savefig(fout_name)
    print(fout_name)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
