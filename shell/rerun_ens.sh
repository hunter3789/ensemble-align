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

#--real time setting
if [ $# -eq 0 ]
then
  d=`date "+%Y%m%d%H"` 

  YY=`echo $d | cut -c1-4`
  MM=`echo $d | cut -c5-6`
  DD=`echo $d | cut -c7-8`
  HH=$(printf '00')

#--history time setting
else
  YY=`echo $1 | cut -c1-4`
  MM=`echo $1 | cut -c5-6`
  DD=`echo $1 | cut -c7-8`
  HH=`echo $1 | cut -c9-10`
fi

dirname=$(printf '/archive-temp/2023_HMT_summer/')
fnamelist=$(printf '/home/clee/projects/2023_HMT_summer/input/ensalign_opt2.input')
logfile=$(printf '/non/clee/ens/logs/ensalign.%s%s%s%s' "${YY}" "${MM}" "${DD}" "${HH}")
filehead=$(printf 'caps_ens')

#mkdir -p $dirname$YY$MM$DD$HH

# activate your environment setting
conda activate myenv

if [ $# -eq 2 ]
then
  fhr=`echo $2`
  echo $YY$MM$DD$HH $fhr
  mpirun -np 40 /home/clee/projects/2023_HMT_summer/ensalign $YY$MM$DD$HH $fhr $fnamelist
  /home/clee/projects/2023_HMT_summer/rewrite $YY$MM$DD$HH $fhr $fnamelist
  sh /home/clee/projects/2023_HMT_summer/run_plot.sh $YY$MM$DD$HH $fhr
else
  memname=('fv3lam-m0b0l0_p' 'fv3lam-m0b0l0_pi' 'fv3lam-m0b0l2_p' 'fv3lam-m0b0l2_pi' 'fv3lam-m0b1l0_pi' 'fv3lam-m0b2l1_p' 'fv3lam-m0b2l1_pi' 'fv3lam-m0b2l2_pi' 'fv3lam-m1b0l0_p' 'fv3lam-m1b0l0_pi' 'fv3lam-m1b0l2_pi'
  'fv3lam-m1b1l0_pi' 'fv3lam-m1b2l1_pi' 'fv3lam-m1b2l2_p' 'fv3lam-m1b2l2_pi')

  plot=0
  fhr=6 
  while [ $fhr -le 84 ]
  do
    str=`echo $fhr | awk '{printf "%03i", $1}'`

    filename=$(printf '%s/%s%s%s%s/%s_%s%s%s%s_f%s.nc' "${dirname}" "${YY}" "${MM}" "${DD}" "${HH}" "${filehead}" "${YY}" "${MM}" "${DD}" "${HH}" "${str}")
    if [ -f $filename ]
    then
      dump=`ncdump -h ${filename}`
      nmembers=`awk -F'nmembers =|;' '{print $2}' <<< "$dump"`
      nmembers=`echo $nmembers | awk '{printf "%i", $1}'`
    else
      nmembers=0
    fi

    kmembers=0
    for i in "${memname[@]}"
    do
      filename=$(printf '%s/%s%s%s%s/%s_%s%s%s%sf%s.grib2' "${dirname}" "${YY}" "${MM}" "${DD}" "${HH}" "${i}" "${YY}" "${MM}" "${DD}" "${HH}" "${str}")
      if [ -f $filename ]
      then
        kmembers=$(($kmembers+1))
      fi
    done

    echo used members: $nmembers, available members: $kmembers
    echo used members: $nmembers, available members: $kmembers >> $logfile

    if [ $nmembers -ne $kmembers ]
    then
      plot=$(($plot+1))
      echo $YY$MM$DD$HH $fhr
      echo $YY$MM$DD$HH $fhr >> $logfile
      mpirun -np 40 /home/clee/projects/2023_HMT_summer/ensalign $YY$MM$DD$HH $fhr $fnamelist >> $logfile
      #wait

      #sh /home/clee/projects/2023_HMT_summer/rewrite_plot.sh $YY$MM$DD$HH $fhr $plot &
    fi
    fhr=$(($fhr+6))
  done
  wait
fi