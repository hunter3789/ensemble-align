fdate=$1
fhr=$2

YY=`echo $1 | cut -c1-4`
MM=`echo $1 | cut -c5-6`
DD=`echo $1 | cut -c7-8`
HH=`echo $1 | cut -c9-10`

finput=$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}")
echo $finput

fobs_date=`date "--date=${YY}-${MM}-${DD} ${HH}:00 ${fhr}hour" "+%Y%m%d%H"`
fobs_name=$(printf '/archive-temp/2023_HMT_summer/stage4/st4_conus.%s.06h.grb2' "${fobs_date}")
echo $fobs_name

regrid_data_plane \
$finput \
$fobs_name \
$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_mean_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
-field 'name="ensmean"; level="0";' \
-field 'name="enslpm"; level="0";' \

regrid_data_plane \
$finput \
$fobs_name \
$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf1_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
-field 'name="ensshfmean"; level="(0,*,*)";' \
-field 'name="ensshflpm"; level="(0,*,*)";' \

regrid_data_plane \
$finput \
$fobs_name \
$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf2_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
-field 'name="ensshfmean"; level="(1,*,*)";' \
-field 'name="ensshflpm"; level="(1,*,*)";' \


mkdir -p /non/clee/ens/vrfy2/$YY$MM$DD

grid_stat \
$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_mean_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
$fobs_name \
GRIDConfig_mean \
-outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

grid_stat \
$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_mean_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
$fobs_name \
GRIDConfig_lpm \
-outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

grid_stat \
$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf1_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
$fobs_name \
GRIDConfig_shfmean1 \
-outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

grid_stat \
$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf1_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
$fobs_name \
GRIDConfig_shflpm1 \
-outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

grid_stat \
$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf2_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
$fobs_name \
GRIDConfig_shfmean2 \
-outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

grid_stat \
$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf2_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
$fobs_name \
GRIDConfig_shflpm2 \
-outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \
