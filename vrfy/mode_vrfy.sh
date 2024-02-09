fdate=$1
fhr=$2

YY=`echo $1 | cut -c1-4`
MM=`echo $1 | cut -c5-6`
DD=`echo $1 | cut -c7-8`
HH=`echo $1 | cut -c9-10`

finput=$(printf '/non/clee/ens/data2/%s%s%s/caps_ens_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}")
echo $finput

fobs_date=`date "--date=${YY}-${MM}-${DD} ${HH}:00 ${fhr}hour" "+%Y%m%d%H"`
fobs_name=$(printf '/home/clee/METplus/Tutorial/data/stage4_conus/st4_conus.%s.06h.grb2' "${fobs_date}")
echo $fobs_name

mkdir -p /non/clee/ens/vrfy2/$YY$MM$DD

thresh=('5' '10' '15' '20' '25')
for i in "${thresh[@]}"
do
  mode \
  $(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_mean_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
  $fobs_name \
  MODEConfig_APCP_mean_$i \
  -outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

  mode \
  $(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_mean_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
  $fobs_name \
  MODEConfig_APCP_lpm_$i \
  -outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

  mode \
  $(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf1_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
  $fobs_name \
  MODEConfig_APCP_shfmean1_$i \
  -outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

  mode \
  $(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf1_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
  $fobs_name \
  MODEConfig_APCP_shflpm1_$i \
  -outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

  mode \
  $(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf2_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
  $fobs_name \
  MODEConfig_APCP_shfmean2_$i \
  -outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

  mode \
  $(printf '/non/clee/ens/data2/%s%s%s/caps_ens_regridded_shf2_%s_f%03d.nc' "${YY}" "${MM}" "${DD}" "${fdate}" "${fhr}") \
  $fobs_name \
  MODEConfig_APCP_shflpm2_$i \
  -outdir $(printf '/non/clee/ens/vrfy2/%s%s%s' "${YY}" "${MM}" "${DD}") \

done
