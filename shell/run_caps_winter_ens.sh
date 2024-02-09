__conda_setup="$('/usr/local/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/usr/local/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/usr/local/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/usr/local/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

# activate your environment setting
conda activate myenv

date_arr=('2024010600' '2024010800' '2024010900')
#dirname=$(printf '/non/clee/ens/data2/')
fnamelist=$(printf '/home/clee/projects/2023_HMT_summer/input/ensalign_winter_snowfall_opt2.input')
fnamelist2=$(printf '/home/clee/projects/2023_HMT_summer/input/ensalign_winter_precipitation_opt2.input')

for i in "${date_arr[@]}"
do
  YY=`echo $i | cut -c1-4`
  MM=`echo $i | cut -c5-6`
  DD=`echo $i | cut -c7-8`
  HH=`echo $i | cut -c9-10`
  #mkdir -p $dirname$YY$MM$DD

  fhr=6 
  while [ $fhr -le 84 ]
  do
    echo $YY$MM$DD$HH $fhr
    mpirun -np 40 /home/clee/projects/2023_HMT_summer/ensalign $YY$MM$DD$HH $fhr $fnamelist
    mpirun -np 40 /home/clee/projects/2023_HMT_summer/ensalign $YY$MM$DD$HH $fhr $fnamelist2
    cp -f $(printf '/archive-temp/2024_HMT_winter/%s%s%s%s/*f%03d.nc' "${YY}" "${MM}" "${DD}" "${HH}" "${fhr}") /archive-temp/2024_HMT_winter/ForTransfer/$YY$MM$DD$HH

    fhr=$(($fhr+6))
  done
done
