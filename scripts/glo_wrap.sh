#!/bin/bash

# Source useful scripts
. $HOME/lib/glo/lock_or_die.sh
# Check to see if we are being run under PBS/Torque
if [[ $PBS_O_WORKDIR ]]
then
   cd $PBS_O_WORKDIR
   m=$(( PBS_VNODENUM + 1 ))
else
   m=$1
fi

. glo_wrap.conf

# Construct the list of all Rc that we are going to compute
allRcStar=`seq -$RcStarMax $RcStarStep $RcStarMax`

# Loop over the list of constructed tau's
for RcStar in $allRcStar
do
   Rc=`echo "scale=3;$RcStar/$Le/1000.000"|bc|sed -e 's/^\(-?\)\./\10./' -e 's/^0$/0.000/'`e3
   base="lo4-tau=$tau-Le=${Le}-Rc=$Rc-m=$m"
   in="$base.conf"
   out="$base.dat"

   if [[ $force ]]; then unlock_or_die $base; fi
   # If the output file or a lock file already exist do nothing!
   lock_or_die $base || continue
   if [[ ! $force &&  -f $out ]]
   then
      continue
   fi
   rm -vf $out
   cat - > $in << EOF
Calculation = -1
Symmetry = 2
eta = 0.35
VariablePar = Rt
tau = $tau
Rt = $RtMin
Pt = $Pt
Le = $Le
Rc = $Rc
m0 = $m
Truncation = 8
AbsParameterError = 1.00E-08
RelativeGRError = 1.00E-05
MaxIterations = 1000
StepSize = -100
UpperLimit = $RtMax
EOF
   #  
   # Run the program
   $HOME/bin/glo $in $out
   rm $in

   unlock_or_die $base || exit 6
done #Rc
