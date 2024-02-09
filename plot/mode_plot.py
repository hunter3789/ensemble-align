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
from scipy.spatial import ConvexHull
import matplotlib.patheffects as path_effects

import netCDF4
import os
from eccodes import *
from datetime import datetime, timedelta
import argparse
import math
import multiprocessing
from scanf import *

#INPUT = 'fcst.nc'
VERBOSE = 1  # verbose error reporting

def example():
    def parse_case(case_str):
        return datetime.strptime(case_str, "%Y%m%d%H")

    ap = argparse.ArgumentParser(description="Plotting MODE verification  Results")
    ap.add_argument('--case', dest='case', type=parse_case, required=True, help="Date/time to plot in YYYYMMDDHH format")
    ap.add_argument('--fhr',  dest='fhr',  type=int, required=True, help="forecast lead time (hours)")
    ap.add_argument('--opt',  dest='opt',  type=int, required=True, help="Alignment option")

    args = ap.parse_args()

    vdate = args.case + timedelta(hours=args.fhr)

    prefix_arr = ['ensmean', 'ensshfmean1', 'ensshfmean2', 'enslpm', 'ensshflpm1', 'ensshflpm2']
    thresh_arr = ['5','10','15','20','25']

    sectors = ['conus', 'mw', 'ne', 'nw', 'sw', 'se', 'cp', 'ma', 'np', 'sp']
    extents = [[-120,-70,25,50],[-105,-80,32.5,45],[-95,-70,37.5,47.5],[-125,-100,39,49],[-122.5,-93.5,28,42.5],[-105,-77,24,38],[-110,-87,31.5,44],
               [-95,-73,34,44],[-110,-86,38,50],[-110,-86,24,39]]

    levels_precip = [0, 0.25, 0.5, 1, 2, 4, 6, 8, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
    levels_str = ['0', '0.25', '0.5', '1', '2', '4', '6', '8', '10', '15', '20', '25', '30', '35', '40', '45', '50', '60', '70', '80', '90', '100']
    levels = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
    cmap = ListedColormap([(1,1,1), (0,236/255,236/255), (0,200/255,240/255), (0,160/255,255/255), (0,60/255,255/255), (0,255/255,0), (0,220/255,0),
                           (0,190/255,0), (0,141/255,0), (255/255,255/255,0), (240/255,210/255,0), (231/255,180/255,0), (200/255,120/255,0),
                           (255/255,160/255,160/255), (255/255,60/255,60/255), (230/255,0,0), (180/255,0,0), (255/255,0,255/255),
                           (217/255,0,217/255), (164/255,0,164/255), (120/255,0,120/255)])

    data_crs = crs.PlateCarree()

    norm_precip = mc.BoundaryNorm(levels_precip, cmap.N)
    norm = mc.BoundaryNorm(levels, cmap.N)

    for prefix in prefix_arr:
        for thresh in thresh_arr:  
            finput_name = '/non/clee/ens/vrfy2/{case:%Y%m%d}/mode_caps_{prefix}_{thresh}_{fhr:02d}0000L_{vdate:%Y%m%d_%H}0000V_000000A_obj.nc'.format(case=args.case, 
                prefix=prefix, thresh=thresh, fhr=args.fhr, vdate=vdate)

            print(finput_name)
 
            f  = netCDF4.Dataset(finput_name)   

            lats = f.variables['lat'][:]
            lons = f.variables['lon'][:]
            value_f = f.variables['fcst_raw'][:] 
            value_o = f.variables['obs_raw'][:]
            value_fo = f.variables['fcst_obj_id'][:]
            value_fc = f.variables['fcst_clus_id'][:]
            value_oo = f.variables['obs_obj_id'][:]
            value_oc = f.variables['obs_clus_id'][:]
            nx = f.variables['fcst_obj_id'].shape[1]
            ny = f.variables['fcst_obj_id'].shape[0]
            n_fcst_simp = n_obs_simp = n_fcst_clus = n_obs_clus = 0


            if 'fcst_simp_hull_start' in f.variables:
                n_fcst_simp = f.variables['fcst_simp_hull_start'].shape[0]
                fcst_simp_hull_start = f.variables['fcst_simp_hull_start']
                fcst_simp_hull_npts = f.variables['fcst_simp_hull_npts']
                fcst_simp_hull_x = f.variables['fcst_simp_hull_x']
                fcst_simp_hull_y = f.variables['fcst_simp_hull_y']

                n_fcst_bdy = f.variables['fcst_simp_bdy_start'].shape[0]
                fcst_simp_bdy_start = f.variables['fcst_simp_bdy_start']
                fcst_simp_bdy_npts = f.variables['fcst_simp_bdy_npts']
                fcst_simp_bdy_lat = f.variables['fcst_simp_bdy_lat']
                fcst_simp_bdy_lon = f.variables['fcst_simp_bdy_lon']


            if 'obs_simp_hull_start' in f.variables:
                n_obs_simp = f.variables['obs_simp_hull_start'].shape[0]
                obs_simp_hull_start = f.variables['obs_simp_hull_start']
                obs_simp_hull_npts = f.variables['obs_simp_hull_npts']
                obs_simp_hull_x = f.variables['obs_simp_hull_x']
                obs_simp_hull_y = f.variables['obs_simp_hull_y']

                n_obs_bdy = f.variables['obs_simp_bdy_start'].shape[0]
                obs_simp_bdy_start = f.variables['obs_simp_bdy_start']
                obs_simp_bdy_npts = f.variables['obs_simp_bdy_npts']
                obs_simp_bdy_lat = f.variables['obs_simp_bdy_lat']
                obs_simp_bdy_lon = f.variables['obs_simp_bdy_lon']


            if 'fcst_clus_hull_start' in f.variables:
                n_fcst_clus = f.variables['fcst_clus_hull_start'].shape[0]
                fcst_clus_hull_start = f.variables['fcst_clus_hull_start']
                fcst_clus_hull_npts = f.variables['fcst_clus_hull_npts']
                fcst_clus_hull_x = f.variables['fcst_clus_hull_x']
                fcst_clus_hull_y = f.variables['fcst_clus_hull_y']
                fcst_clus_hull_lat = f.variables['fcst_clus_hull_lat']
                fcst_clus_hull_lon = f.variables['fcst_clus_hull_lon']


            if 'obs_clus_hull_start' in f.variables:
                n_obs_clus = f.variables['obs_clus_hull_start'].shape[0]
                obs_clus_hull_start = f.variables['obs_clus_hull_start']
                obs_clus_hull_npts = f.variables['obs_clus_hull_npts']
                obs_clus_hull_x = f.variables['obs_clus_hull_x']
                obs_clus_hull_y = f.variables['obs_clus_hull_y']
                obs_clus_hull_lat = f.variables['obs_clus_hull_lat']
                obs_clus_hull_lon = f.variables['obs_clus_hull_lon']

################ Area size Check
            ftext_name = '/non/clee/ens/vrfy2/{case:%Y%m%d}/mode_caps_{prefix}_{thresh}_{fhr:02d}0000L_{vdate:%Y%m%d_%H}0000V_000000A_obj.txt'.format(case=args.case, 
                prefix=prefix, thresh=thresh, fhr=args.fhr, vdate=vdate)

            f_text  = open(ftext_name, "r") 
            areas = []
            area_thresh = 200
            for strs in f_text:
                lists = scanf("%s %s %d %d %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s", strs)
                if (lists is not None):
                    area = lists[31]
                    if (area != "NA"):
                        areas.append(int(area))

            f_text.close()
################

            paths_f_array = []
            cx_f_array = []
            cy_f_array = []
            label_f_array = []

            for k in range(n_fcst_simp):
                if (areas[k] < area_thresh):
                    continue

                npts = fcst_simp_hull_npts[k]
                start = fcst_simp_hull_start[k]
                points = np.zeros([npts, 2])
                for idx in range(npts):
                    points[idx,0] = fcst_simp_hull_x[start+idx]
                    points[idx,1] = fcst_simp_hull_y[start+idx]

                npts2 = fcst_simp_bdy_npts[k]
                start2 = fcst_simp_bdy_start[k]
                paths = np.zeros([npts2+1, 2])

                for idx in range(npts2):
                    paths[idx,0] = fcst_simp_bdy_lon[start2+idx]
                    paths[idx,1] = fcst_simp_bdy_lat[start2+idx]

                paths[npts2,0] = fcst_simp_bdy_lon[start2]
                paths[npts2,1] = fcst_simp_bdy_lat[start2]
                
                hull = ConvexHull(points)
                cx = int(round(np.mean(hull.points[hull.vertices,0]),0))
                cy = int(round(np.mean(hull.points[hull.vertices,1]),0))

                paths_f_array.append(paths)
                cx_f_array.append(cx)
                cy_f_array.append(cy)
                label_f_array.append(k)

                del points, paths            
                del hull


            paths_o_array = []
            cx_o_array = []
            cy_o_array = []
            label_o_array = []

            for k in range(n_obs_simp):
                if (areas[n_fcst_simp+k] < area_thresh):
                    continue

                npts = obs_simp_hull_npts[k]
                start = obs_simp_hull_start[k]
                points = np.zeros([npts, 2])
                for idx in range(npts):
                    points[idx,0] = obs_simp_hull_x[start+idx]
                    points[idx,1] = obs_simp_hull_y[start+idx]

                npts2 = obs_simp_bdy_npts[k]
                start2 = obs_simp_bdy_start[k]
                paths = np.zeros([npts2+1, 2])

                for idx in range(npts2):
                    paths[idx,0] = obs_simp_bdy_lon[start2+idx]
                    paths[idx,1] = obs_simp_bdy_lat[start2+idx]

                paths[npts2,0] = obs_simp_bdy_lon[start2]
                paths[npts2,1] = obs_simp_bdy_lat[start2]

                hull = ConvexHull(points)
                cx = int(round(np.mean(hull.points[hull.vertices,0]),0))
                cy = int(round(np.mean(hull.points[hull.vertices,1]),0))

                paths_o_array.append(paths)
                cx_o_array.append(cx)
                cy_o_array.append(cy)
                label_o_array.append(k)

                del points, paths
                del hull


            f.close()
################

            processes = []
            pool = multiprocessing.Pool(processes=len(sectors))
 
            for isector in range(len(sectors)):
                p = pool.apply_async(plot_chart, (isector, args, thresh, prefix, sectors, extents, levels, levels_str, levels_precip, lons, lats, value_f, value_o, 
                    paths_f_array, cx_f_array, cy_f_array, label_f_array, paths_o_array, cx_o_array, cy_o_array, label_o_array, cmap, norm, norm_precip, data_crs))
                processes.append(p)

            for process in processes:
                process.get()

            pool.close()
            pool.join()


################
def plot_chart(isector, args, thresh, prefix, sectors, extents, levels, levels_str, levels_precip, lons, lats, value_f, value_o, 
               paths_f_array, cx_f_array, cy_f_array, label_f_array, paths_o_array, cx_o_array, cy_o_array, label_o_array, cmap, norm, norm_precip, data_crs):

    match prefix:
        case 'ensmean':
            label = 'Ens mean'
        case 'ensshfmean1':
            label = 'SAM (1pass)'
        case 'ensshfmean2':
            label = 'SAM (2pass)'
        case 'enslpm':
            label = 'Ens LPM'
        case 'ensshflpm1':
            label = 'SAM-LPM (1pass)'
        case 'ensshflpm2':
            label = 'SAM-LPM (2pass)'
        case default:
            label = prefix

    fig, axs = plt.subplots(2,2,figsize=(10,7.5), subplot_kw={"projection": crs.LambertConformal(central_longitude=262.5, central_latitude=38.5)})
    fig.subplots_adjust(hspace=0, wspace=0, right=0.95, top=0.9, left=0.05)
    #fig.suptitle("{case:%Y.%m.%d.%H}UTC (+{fhr:03d}H) - MODE: {prefix} (>={thresh}mm/6hr)".format(case=args.case, 
    #            fhr=args.fhr, prefix=prefix, thresh=thresh), fontsize=16.)
    if (args.opt == 1):
        fig.suptitle("VALID: {valid:%Y.%m.%d.%H}UTC (+{fhr:03d}H) / RUN: {case:%Y.%m.%d.%H}UTC / MODE: {label} (>={thresh}mm/6hr) / Align Method A".format(valid=args.case + timedelta(hours=args.fhr), 
          fhr=args.fhr, case=args.case, label=label, thresh=thresh), fontsize=11.)
    elif (args.opt == 2):
        fig.suptitle("VALID: {valid:%Y.%m.%d.%H}UTC (+{fhr:03d}H) / RUN: {case:%Y.%m.%d.%H}UTC / MODE: {label} (>={thresh}mm/6hr) / Align Method B".format(valid=args.case + timedelta(hours=args.fhr), 
          fhr=args.fhr, case=args.case, label=label, thresh=thresh), fontsize=11.)

    contours = axs[0,0].pcolormesh(lons, lats, value_f, cmap=cmap, norm=norm_precip, transform=data_crs)
    contours = axs[0,1].pcolormesh(lons, lats, value_o, cmap=cmap, norm=norm_precip, transform=data_crs)

    axs[0,0].tick_params(left = False, right = False, labelleft = False, labelbottom = False, bottom = False)
    axs[0,1].tick_params(left = False, right = False, labelleft = False, labelbottom = False, bottom = False)
    axs[1,0].tick_params(left = False, right = False, labelleft = False, labelbottom = False, bottom = False)
    axs[1,1].tick_params(left = False, right = False, labelleft = False, labelbottom = False, bottom = False)

    axs[0,0].set_title('Forecast', fontsize=12.)
    axs[0,1].set_title('Observation', fontsize=12.)

    axs[0,0].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[0,0].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    axs[0,1].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[0,1].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    axs[1,0].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[1,0].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")
    axs[1,1].add_feature(cfeature.COASTLINE.with_scale('50m'), linewidth=0.8, edgecolor="grey")
    axs[1,1].add_feature(cfeature.STATES.with_scale('50m'), linewidth=0.5, edgecolor="grey")

    ncount = 0
    for k in range(len(paths_f_array)):
        ncount = ncount + 1
        axs[1,0].plot(paths_f_array[k][:,0], paths_f_array[k][:,1], linewidth=0.5, color='black', transform=data_crs)
        axs[1,0].fill(paths_f_array[k][:,0], paths_f_array[k][:,1], linewidth=0.5, color=cmap(ncount), transform=data_crs)
        axs[0,0].plot(paths_f_array[k][:,0], paths_f_array[k][:,1], linewidth=1, color='red', transform=data_crs)

        text = axs[1,0].text(lons[cy_f_array[k], cx_f_array[k]], lats[cy_f_array[k],cx_f_array[k]], str(label_f_array[k]+1), color='black', size=10, clip_on=True,
                horizontalalignment='center', verticalalignment='center', transform=data_crs)
        text.set_path_effects([path_effects.withStroke(linewidth=1, foreground='w')])
        text.set_clip_box(axs[1,0].bbox) 

    ncount = 0
    for k in range(len(paths_o_array)):
        ncount = ncount + 1
        axs[1,1].plot(paths_o_array[k][:,0], paths_o_array[k][:,1], linewidth=0.5, color='black', transform=data_crs)
        axs[1,1].fill(paths_o_array[k][:,0], paths_o_array[k][:,1], linewidth=0.5, color=cmap(ncount), transform=data_crs)
        axs[0,1].plot(paths_o_array[k][:,0], paths_o_array[k][:,1], linewidth=1, color='red', transform=data_crs)

        text = axs[1,1].text(lons[cy_o_array[k], cx_o_array[k]], lats[cy_o_array[k],cx_o_array[k]], str(label_o_array[k]+1), color='black', size=10, clip_on=True,
            horizontalalignment='center', verticalalignment='center', transform=data_crs)
        text.set_path_effects([path_effects.withStroke(linewidth=1, foreground='w')])
        text.set_clip_box(axs[1,1].bbox)


    axs[0,0].set_extent(extents[isector], data_crs)
    axs[0,1].set_extent(extents[isector], data_crs)
    axs[1,0].set_extent(extents[isector], data_crs)
    axs[1,1].set_extent(extents[isector], data_crs)

    axs[0,0].set_facecolor((0.6,0.6,0.6))
    axs[0,1].set_facecolor((0.6,0.6,0.6))
    cax = fig.add_axes([0.05,0.07,0.9,0.02])
    cbar = fig.colorbar(contours, cax=cax, ticks=levels_precip, orientation='horizontal')
    cbar.ax.set_xticklabels(levels_str)
    cbar.ax.set_xlabel('[mm]')

    fout_name = '/non/clee/ens/img/{case:%Y%m%d}/caps_ens_opt{opt:01d}_{case:%Y%m%d%H}_f{fhr:03d}_MODE_{prefix}_{thresh}_{sector}.png'.format(opt=args.opt, case=args.case, 
        fhr=args.fhr, prefix=prefix, thresh=thresh, sector=sectors[isector])
    if os.path.exists(fout_name):
        os.unlink(fout_name)
    plt.savefig(fout_name)
    plt.close(fig)
    print("--- Finish {fname} ---".format(fname=fout_name))

 
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

 
if __name__ == "__main__":
    sys.exit(main())
