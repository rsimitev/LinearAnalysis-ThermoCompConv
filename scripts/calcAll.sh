#!/bin/bash
#
i=$1
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
Le=(  0.1     0.3   0.5   1.0   2.0   5.0    10.0     20.0     )
Ra=(  1000   1000  1000  1000  1000  1000    1000     1000      )
m=(   2         2     2     2     2     2       2       2       )
for Rc in 1.0e2 5.0e2 1.0e3 5.0e3 1.0e4 5.0e4 1.0e5
do
   count=0
   base="lo1-Le=${Le[$i]}-Rc=$Rc"
   in="$base.conf"
   out="$base.dat"
   rm -vf $out $base-*.dat
cat > $in << EOT
  NE (0/1/2) | LCALC (-1/0/1/2/3/4/5/6) |
       2               3
 |  RAYLEIGH  |  TAU     |  PRANTEL  |  ETA  |  Lewis |   Rconc   |
     ${Ra[$i]}   1.0e2        1.0       0.35  ${Le[$i]}   $Rc
 |   NTRUNC (>=1) | MODE |
        10         ${m[$i]}
 |   DRA   |  ABSERR  |  RELERR  | NMAX |
  100.00     1.00E-08   0.00E-03   1000
 |   incr_STEP | Tau_END
       100       1.0e5
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
     cp $out $base-$count.dat
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
done # Rc
