#!/bin/bash
#
if [[ $# -ne 10 ]]
then
   echo "Usage:"
   echo "$0 -l <Le> -r <Rc> -t <Ta> -R <Ra> -m <m>"
   echo ""
   echo " <Le> - fixed Lewis number for the calculation."
   echo " <Rc> - fixed compositional Rayleigh number for the calculation."
   echo " <Ta> - starting value of the Taylor number (will vary up to 1.0e5)"
   echo " <Ra> - initial estimate of the thermal Rayleigh number."
   echo " <Le> - initial estimate of the most unstable wave number."
   exit 1
fi
while getopts "l:r:t:R:m:" opt
do
   case $opt in
      l)
         Le=$OPTARG;
         ;;
      r)
         Rc=$OPTARG;
         ;;
      t)
         Ta=$OPTARG;
         ;;
      R)
         Ra=$OPTARG;
         ;;
      m)
         m=$OPTARG;
         ;;
   esac
done
#make lo
#========================================================== 
#   input for lo
#* ------- Physical Parameters: ------
#* RA=RAYLEIGH NUMBER,  TAU=Coriolis NUMBER, PR=PRANTEL NUMBER,
#* ETA=RATIO OF RADII,  NT= TRUNCATION,  MODE=WAVE NUMBER,
#* pL = Lewis number,  Rconc=Rayleigh number due to concentration
#* 
#* NE= SYMMETRY PARAMETER, NE=0 : UNDEFINED SYMMETRY,
#* NE=2 : EQUATORIAL SYMMETRY, NE=1 : EQUATORIAL ANTISYMMETRY.
#* DRIFT C IS DEFINED as (phi-c*t).
#*
#* LCALC=-1: Most basic case: find the most unstable growth rate at all parameters fixed
#* LCALC=0 : Critical Rayleigh number at fixed other parameters
#* LCALC=1 : All eigenvalues at fixed parameters including the Rayleigh number
#* LCALC=2 : Onset determined for constant wavenumber M 
#*           (by searching root of grothrate in R, using pegasus.f).
#* LCALC=3 : Onset determined by varying Rayleigh number R and wavenumber M.
#* LCALC=4 : Eigenvector determined for one set of parameters at onset.
#*           (Solution can be then plotted)
#* LCALC=5 : Increment along M
#* LCALC=6 : Vary Le and calculate critical R at fixed P, tau, eta, M
#*
#* LO.F calculates R (crit. Rayleighn.) and Omega (and M) in the
#* range TAU=TTA to TAU=TTF.
#========================================================== 
#
# Initial guesses for each Lewis number
   count=0
   base="lo1-Le=${Le}-Rc=$Rc"
   in="$base.conf"
   out="$base.dat"
   rm -vf $out $base-*.dat
cat > $in << EOT
Symmetry = 2
Calculation = 3
VariablePar = Rt
Rt = ${Ra}
tau = ${ta}
Pt = 1.0d0
eta = 0.35d0 
Le = ${Le} 
Rc = ${Rc}
Truncation = 10 
m0 = ${m}
AbsParameterError = 1.00E-08 
RelativeGRError = 1.00E-03 
MaxIterations = 1000
StepSize = -100 
UpperLimit = 1.00d2
EOT
   #  
   res=10
   # Command
   # Keep running while the output status is 10
   # or we have gone through 40 iterations
   while [[ $res -eq 10 && $count -le 30 ]]
   do
     # Backup config file as it is potentially going to be 
     # rewritten by ./lo
     cp $in $base-$count.conf
     # Run the program
     ./lo $in $out
     res=$?
     # Backup the data file so it won't get rewritten
     mv $out $base-$count.dat
     count=$((count+1))
   done
   count=$((count-1))
   if [[ ${count} -gt 0 ]]
   then
      mv $base-0.dat $out
      for c in `seq 1 $count`
      do
         cat $base-$c.dat| grep -v '#' >> $out
#         rm $base-$c.dat
      done
   fi
   ./listAll.sh
