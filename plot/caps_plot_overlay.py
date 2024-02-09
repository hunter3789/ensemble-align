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
from scipy.ndimage import gaussian_filter
from scipy.ndimage import maximum_filter, minimum_filter

import netCDF4
import os
from eccodes import *
from datetime import datetime, timedelta
import argparse
import math
import time
import multiprocessing

VERBOSE = 1  # verbose error reporting

def example():
    def parse_case(case_str):
        return datetime.strptime(case_str, "%Y%m%d%H")

    ap = argparse.ArgumentParser(description="Plotting HREF Ensemble Results")
    ap.add_argument('--case', dest='case', type=parse_case, required=True, help="Date/time to plot in YYYYMMDDHH format")
    ap.add_argument('--fhr',  dest='fhr',  type=int, required=True, help="forecast lead time (hours)")
    ap.add_argument('--mask', dest='mask', type=int, required=False, default=0, help="obs mask")
    ap.add_argument('--season', dest='season', type=str, required=False, default="summer", help="season")
    ap.add_argument('--obs', dest='obs', type=int, required=False, default=1, help="obs plot")

    args = ap.parse_args()

    #finput_name = '/archive-temp/2023_HMT_summer/{case:%Y%m%d%H}/caps_ens_{case:%Y%m%d%H}_f{fhr:03d}.nc'.format(case=args.case, fhr=args.fhr)
    finput_name = '/non/clee/ens/data2/{case:%Y%m%d}/caps_ens_{case:%Y%m%d%H}_f{fhr:03d}.nc'.format(case=args.case, fhr=args.fhr)
    f  = netCDF4.Dataset(finput_name)   

    lats = f.variables['lat'][:]
    lons = f.variables['lon'][:]
    value = f.variables['precipitation'][:]
    value_u = f.variables['shiftvec_u'][:]
    value_v = f.variables['shiftvec_v'][:]
    value_shf = f.variables['shift_prcp'][:]
    nx = value_u.shape[3]
    ny = value_u.shape[2]
    nmembers = value_u.shape[1]
    npass = value_u.shape[0]
    value1 = f.variables['ensmean'][:]
    value3 = f.variables['enspm'][:]
    value4 = f.variables['enslpm'][:]
    memname = f.variables['members'][:]
    f.close()

    if (args.season == "winter"):
        args.obs = 0
        value1[:,:] *= 100.
        value3[:,:] *= 100.
        value4[:,:] *= 100.
        value[:,:,:] *= 100.

        #for y in range(ny):
        #    for x in range(nx):
        #        value1[y,x] *= 100.
        #        value3[y,x] *= 100.
        #        value4[y,x] *= 100.
        #        for k in range(nmembers):
        #            value[k,y,x] *= 100.

    value1_ext = maximum_filter(value1, 500, mode='nearest')
    value3_ext = maximum_filter(value3, 500, mode='nearest')
    value4_ext = maximum_filter(value4, 500, mode='nearest')
    mxy1, mxx1 = np.where((value1_ext == value1) & (value1_ext >= 0.5))
    mxy3, mxx3 = np.where((value3_ext == value3) & (value3_ext >= 0.5))
    mxy4, mxx4 = np.where((value4_ext == value4) & (value4_ext >= 0.5))

    value1_ext2 = maximum_filter(value1, 200, mode='nearest')
    value3_ext2 = maximum_filter(value3, 200, mode='nearest')
    value4_ext2 = maximum_filter(value4, 200, mode='nearest')
    dmxy1, dmxx1 = np.where((value1_ext2 == value1) & (value1_ext2 >= 0.5))
    dmxy3, dmxx3 = np.where((value3_ext2 == value3) & (value3_ext2 >= 0.5))
    dmxy4, dmxx4 = np.where((value4_ext2 == value4) & (value4_ext2 >= 0.5))

    no_value = np.zeros(lats.shape)   
    nmembers = value.shape[0]

    labels = []
    for i in range(nmembers):
        labels.append("")
        for x in range(memname.shape[1]):
          char = memname[i][x].decode('utf-8')
          if (char != " "):
               labels[i] = labels[i] + char


    if (args.obs == 1):
        fobs_name = '/archive-temp/2023_HMT_summer/stage4/st4_conus.{date:%Y%m%d%H}.06h.grb2'.format(case=args.case, date=args.case + timedelta(hours=args.fhr))

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

        #codes_close_file(fobs_name)

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
    else:
        lons_grib = 0
        lats_grib = 0
        value_grib = 0
        mxx_grib = 0
        mxy_grib = 0
        dmxx_grib = 0
        dmxy_grib = 0


    labels_ref = ['fv3lam-m0b0l0_p','fv3lam-m0b0l0_pi','fv3lam-m0b0l2_p','fv3lam-m0b0l2_pi','fv3lam-m0b1l0_pi','fv3lam-m0b2l1_p','fv3lam-m0b2l1_pi',
                  'fv3lam-m0b2l2_pi','fv3lam-m1b0l0_p','fv3lam-m1b0l0_pi','fv3lam-m1b0l2_pi','fv3lam-m1b1l0_pi','fv3lam-m1b2l1_pi',
                  'fv3lam-m1b2l2_p','fv3lam-m1b2l2_pi']

    if (args.season == "winter"):
        levels = [0, 0.25, 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 20, 25, 30, 40, 50, 60, 70, 80, 100]
        levels_str = ['0', '0.25', '1', '2', '3', '4', '5', '6', '8', '10', '12', '15', '18', '20', '25', '30', '40', '50', '60', '70', '80', '100']
        cmap = ListedColormap([(1,1,1), (179/255,180/255,221/255), (255/255,234/255,110/255), (255/255,219/255,32/255), (224/255,185/255,0/255), (204/255,170/255,1/255), (202/255,251/255,177/255),
                           (103/255,239/255,77/255), (53/255,204/255,28/255), (135/255,217/255,255/255), (7/255,171/255,255/255), (1/255,141/255,221/255), (0/255,119/255,179/255),
                           (218/255,135/255,255/255), (194/255,62/255,254/255), (173/255,8/255,254/255), (127/255,0,191/255), (250/255,133/255,133/255),
                           (246/255,62/255,62/255), (218/255,0,0/255), (188/255,0,0/255)])
    else:
        levels = [0, 0.25, 0.5, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
        levels_str = ['0', '0.25', '0.5', '1', '2', '4', '6', '8', '10', '15', '20', '25', '30', '35', '40', '45', '50', '60', '70', '80', '90', '100']
        cmap = ListedColormap([(1,1,1), (0,236/255,236/255), (0,200/255,240/255), (0,160/255,255/255), (0,60/255,255/255), (0,255/255,0), (0,220/255,0),
                           (0,190/255,0), (0,141/255,0), (255/255,255/255,0), (240/255,210/255,0), (231/255,180/255,0), (200/255,120/255,0),
                           (255/255,160/255,160/255), (255/255,60/255,60/255), (230/255,0,0), (180/255,0,0), (255/255,0,255/255),
                           (217/255,0,217/255), (164/255,0,164/255), (120/255,0,120/255)])

    sectors = ['conus', 'mw', 'ne', 'nw', 'sw', 'se', 'cp', 'ma', 'np', 'sp']
    extents = [[-120,-70,25,50],[-105,-80,32.5,45],[-95,-70,37.5,47.5],[-125,-100,39,49],[-122.5,-93.5,28,42.5],[-105,-77,24,38],[-110,-87,31.5,44],
               [-95,-73,34,44],[-110,-86,38,50],[-110,-86,24,39]]

    data_crs = crs.PlateCarree()
    norm = mc.BoundaryNorm(levels, cmap.N)

    processes = []
    pool = multiprocessing.Pool(processes=len(sectors))

    for k in range(len(sectors)):
        p = pool.apply_async(plot_chart, (k, args, labels_ref, labels, sectors, extents, levels, levels_str, lons, lats, lons_grib, lats_grib, value, value1, value3, value4, value_grib, no_value, 
               mxx1, mxy1, mxx3, mxy3, mxx4, mxy4, mxx_grib, mxy_grib, dmxx1, dmxy1, dmxx3, dmxy3, dmxx4, dmxy4, dmxx_grib, dmxy_grib, cmap, norm, data_crs))
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
 

def plot_chart(k, args, labels_ref, labels, sectors, extents, levels, levels_str, lons, lats, lons_grib, lats_grib, value, value1, value3, value4, value_grib, no_value, 
               mxx1, mxy1, mxx3, mxy3, mxx4, mxy4, mxx_grib, mxy_grib, dmxx1, dmxy1, dmxx3, dmxy3, dmxx4, dmxy4, dmxx_grib, dmxy_grib, cmap, norm, data_crs):

    fig, axs = plt.subplots(4,5,figsize=(20,12.5), subplot_kw={"projection": crs.LambertConformal(central_longitude=262.5, central_latitude=38.5)})
    fig.subplots_adjust(hspace=0.1, wspace=0.1, right=0.95, top=0.9, left=0.05)
    if (args.season == "winter"):
        fig.suptitle("VALID: {valid:%Y.%m.%d.%H}UTC (+{fhr:03d}H) / RUN: {case:%Y.%m.%d.%H}UTC / 6hr Snow".format(valid=args.case + timedelta(hours=args.fhr), 
          fhr=args.fhr, case=args.case), fontsize=16.)
    else:
        fig.suptitle("VALID: {valid:%Y.%m.%d.%H}UTC (+{fhr:03d}H) / RUN: {case:%Y.%m.%d.%H}UTC / 6hr Precipitation".format(valid=args.case + timedelta(hours=args.fhr), 
          fhr=args.fhr, case=args.case), fontsize=16.)

    for x in range(len(labels_ref)):
        if (labels_ref[x] in labels):
            n = labels.index(labels_ref[x])
            contours = axs[int(x/5),x%5].pcolormesh(lons, lats, value[n,:,:], cmap=cmap, norm=norm, transform=data_crs)
            axs[int(x/5),x%5].tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
            axs[int(x/5),x%5].set_title(labels_ref[x], fontsize=12.)

            if (sectors[k] == 'conus'):
                axs[int(x/5),x%5].text(0, -0.1, "max: {:.2f}mm".format(np.max(value[n,:,:])), fontsize=12., transform=axs[int(x/5),x%5].transAxes)

            axs[int(x/5),x%5].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
            axs[int(x/5),x%5].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
            axs[int(x/5),x%5].set_extent(extents[k], data_crs)
        else:
            contours = axs[int(x/5),x%5].pcolormesh(lons, lats, no_value, cmap=cmap, norm=norm, transform=data_crs)
            axs[int(x/5),x%5].tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
            axs[int(x/5),x%5].set_title(labels_ref[x], fontsize=12.)
            axs[int(x/5),x%5].text(0.4, 0.5, "No Data", fontsize=12., transform=axs[int(x/5),x%5].transAxes) 
            axs[int(x/5),x%5].set_extent(extents[k], data_crs)

    contours = axs[3,1].pcolormesh(lons, lats, value1, cmap=cmap, norm=norm, transform=data_crs)
    contours = axs[3,2].pcolormesh(lons, lats, value3, cmap=cmap, norm=norm, transform=data_crs)
    contours = axs[3,3].pcolormesh(lons, lats, value4, cmap=cmap, norm=norm, transform=data_crs)
    if (args.obs == 1):
        contours = axs[3,4].pcolormesh(lons_grib, lats_grib, value_grib[:,:], cmap=cmap, norm=norm, transform=data_crs)

    axs[3,1].tick_params(left = False, right = False , labelleft = False ,
            labelbottom = False, bottom = False)
    axs[3,2].tick_params(left = False, right = False , labelleft = False ,
            labelbottom = False, bottom = False)
    axs[3,3].tick_params(left = False, right = False , labelleft = False ,
            labelbottom = False, bottom = False)
    if (args.obs == 1):
        axs[3,4].tick_params(left = False, right = False , labelleft = False ,
            labelbottom = False, bottom = False)

    axs[3,1].set_title('Ens mean', fontsize=12.)
    axs[3,2].set_title('Ens PM', fontsize=12.)
    axs[3,3].set_title('Ens LPM', fontsize=12.)
    if (args.obs == 1):
        axs[3,4].set_title('Obs - Stage IV', fontsize=12.)

    if (sectors[k] == 'conus'):
        axs[3,1].text(0., -0.1, "max: {:.2f}mm".format(np.max(value1)), fontsize=12., transform=axs[3,1].transAxes)
        axs[3,2].text(0., -0.1, "max: {:.2f}mm".format(np.max(value3)), fontsize=12., transform=axs[3,2].transAxes)
        axs[3,3].text(0., -0.1, "max: {:.2f}mm".format(np.max(value4)), fontsize=12., transform=axs[3,3].transAxes)
        if (args.obs == 1):
            axs[3,4].text(0., -0.1, "max: {:.2f}mm".format(np.max(value_grib)), fontsize=12., transform=axs[3,4].transAxes)

    axs[3,1].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[3,1].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    axs[3,2].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[3,2].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    axs[3,3].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[3,3].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    if (args.obs == 1):
        axs[3,4].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
        axs[3,4].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")

    if (sectors[k] == 'conus'):
        plot_maxmin_points2(axs[3,1], lons, lats, value1, mxx1, mxy1, symbol='', color='black',  transform=data_crs)
        plot_maxmin_points2(axs[3,2], lons, lats, value3, mxx3, mxy3, symbol='', color='black',  transform=data_crs)
        plot_maxmin_points2(axs[3,3], lons, lats, value4, mxx4, mxy4, symbol='', color='black',  transform=data_crs)
        if (args.obs == 1):
            plot_maxmin_points2(axs[3,4], lons_grib, lats_grib, value_grib, mxx_grib, mxy_grib, symbol='', color='black',  transform=data_crs)
    else:
        plot_maxmin_points2(axs[3,1], lons, lats, value1, dmxx1, dmxy1, symbol='', color='black',  transform=data_crs)
        plot_maxmin_points2(axs[3,2], lons, lats, value3, dmxx3, dmxy3, symbol='', color='black',  transform=data_crs)
        plot_maxmin_points2(axs[3,3], lons, lats, value4, dmxx4, dmxy4, symbol='', color='black',  transform=data_crs)
        if (args.obs == 1):
            plot_maxmin_points2(axs[3,4], lons_grib, lats_grib, value_grib, dmxx_grib, dmxy_grib, symbol='', color='black',  transform=data_crs)

    axs[3,1].set_extent(extents[k], data_crs)
    axs[3,2].set_extent(extents[k], data_crs)
    axs[3,3].set_extent(extents[k], data_crs)
    if (args.obs == 1):
        axs[3,4].set_extent(extents[k], data_crs)
        axs[3,4].set_facecolor((0.6,0.6,0.6))
    else:
        fig.delaxes(axs[3,4])

    fig.delaxes(axs[3,0])

    cax = fig.add_axes([0.05,0.05,0.9,0.02])
    cbar = fig.colorbar(contours, cax=cax, ticks=levels, orientation='horizontal')
    cbar.ax.set_xticklabels(levels_str)
    if (args.season == "winter"):
        cbar.ax.set_xlabel('[cm]')
    else:
        cbar.ax.set_xlabel('[mm]')

    fout_name = '/non/clee/ens/img/{case:%Y%m%d}/caps_ens_{case:%Y%m%d%H}_f{fhr:03d}_{sector}_member.png'.format(case=args.case, fhr=args.fhr, sector=sectors[k])

    if os.path.exists(fout_name):
        os.unlink(fout_name)

    plt.savefig(fout_name)
    print(fout_name)
    plt.close(fig)

if __name__ == "__main__":
    sys.exit(main())
