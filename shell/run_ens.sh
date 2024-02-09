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
fnamelist=$(printf '/home/clee/projects/2023_HMT_summer/input/ensalign.input')

logfile=$(printf '/non/clee/ens/logs/ensalign.%s%s%s%s' "${YY}" "${MM}" "${DD}" "${HH}")

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
  fhr=6
  while [ $fhr -le 84 ]
  do
    echo $YY$MM$DD$HH $fhr
    echo $YY$MM$DD$HH $fhr >> $logfile
    mpirun -np 40 /home/clee/projects/2023_HMT_summer/ensalign $YY$MM$DD$HH $fhr $fnamelist >> $logfile
    wait

    sh /home/clee/projects/2023_HMT_summer/rewrite_plot.sh $YY$MM$DD$HH $fhr &
    fhr=$(($fhr+6))
  done
  wait
fi