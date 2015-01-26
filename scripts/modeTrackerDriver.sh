#!/bin/bash
# ***********************************
# This script drives the mode tracker.
# The mode tracker takes a set of files 
# each representing a set of lines at a value of x
# and reorders those lines so that continuity is imposed.
# ***********************************
# Adimensional parameters
tau=1.2e3
Rc=4.0e3
Le=10.0
prefix=lo6

# ***********************************
# function to extract the value of a parameter from a filename
# constructed as ${label}-${var1}=${value1}-${var2}=$}value2}-...
function extractPar() 
{
   local label=$1; # the prefix of the filename.
   local var=$2; # The parameter whose value we want.
   while read data
   do 
      echo $data|grep -q ${var} && \
      echo $data| sed -e 's/'${label}'//' \
                      -e 's/.*-'${var}'=\(-[0-9.eE]*\).*/\1/g' \
                      -e 's/.*-'${var}'=\([0-9.eE]*\).*/\1/g' \
                      -e 's/\.$//' \
                      -e 's/-$//'
   done 
}

# ***********************************
# Function to sort filenames constructed as ${label}-${var1}=${value1}-${var2}=$}value2}-...
# This function does not deal well with values in exponential notation that
# have different exponents :(
function sortBy()
{
   local par=$1
   local template=$2
   local loc=`echo ${template}|sed -e 's/=/\n/g'|nl|awk '/'${par}'/ {print $1+1}'`
   sort -t '=' -k${loc},${loc} -g
}

function debug()
{
   local n=$1
   for i in `seq 1 $n`
   do
      echo -n '-'
   done
   echo -n ' '
   echo -en "`cat -`\r"
}
# ***********************************
# Set a template for the fileames
template="${prefix}-tau=1-Le=1-Rc=1-Rt=1-m=1.dat"
# We have sets of lines for several m's.
# Loop over the m's
tmp=`mktemp`
tmp0=`mktemp`
for m in `seq 1 5`
do
   debug 0 <<< "Processing m=$m"
   # Remove previously existing sets of lines
   rm -vf ${prefix}-tau=${tau}-Le=${Le}-Rc=${Rc}-Rt=*-m=$m.dat-sorted
   # Prepare the input for the mode tracker
   # 1. Sort fiels by Rt.
   ls -U ${prefix}-tau=${tau}-Le=${Le}-Rc=${Rc}-Rt=*-m=$m.dat|sortBy Rt $template > $tmp0
   # 2. Prepend value of Rt to the file list
   cat /dev/null > $tmp
   for fileName in `cat $tmp0`
   do
      Rt=`echo $fileName | extractPar ${prefix} Rt`
      echo "$Rt   $fileName" >> $tmp
   done 
   # 3. Count the total number of modes
   totalNModes=`head -n 1 $tmp0| xargs grep -c -v '#'`
   # 4. Prepend the total number of modes to consider in these files.
   echo $totalNModes |cat - $tmp > myList
   ./modeTracker2 myList
   # cleanup previously existing positive eigen mode files
   rm -vf ${prefix}-tau=$tau-Le=$Le-Rc=$Rc-m=$m-pem=*.dat
   positiveMode=0
   # For each eigen mode create a file with the whole line.
   # The growth rate is also given the correct sign here.
   for em in `seq 1 ${totalNModes}`
   do 
      debug 1 <<< "Preparing eigenmode $em.\n"
      # Files sorted by Rt
      lo6files=`ls -U *-m=$m.dat-sorted|sortBy Rt $template`
      lo6out="${prefix}-tau=$tau-Le=$Le-Rc=$Rc-m=$m-em=$em.dat"
      cat /dev/null > $lo6out
      for fileName in $lo6files
      do
         Rt=`echo $fileName | extractPar ${prefix} Rt`
         debug 3 <<< "Parsing $fileName."
         # awk does not cope well with fortran style doble precision format
         head -n ${totalNModes} $fileName | sed -e 's/D/e/g' > $tmp
         awk '$1=='$em' {printf "%15.2f %18.5e \n", '$Rt', -$3}' $tmp >> $lo6out
      done
      # Select the modes that become positive out of all the modes an save them to their own file.
      if [[ `awk 'BEGIN {positive=0}; $2>0.0 {positive=positive+1}; END {print positive}' $lo6out` -gt 0 ]]
      then
         debug 2 <<< "Saving eignemode $em which has positive eigen values.\n"
         positiveMode=$((positiveMode+1))
         cp ${prefix}-tau=$tau-Le=$Le-Rc=$Rc-m=$m-em=$em.dat ${prefix}-tau=$tau-Le=$Le-Rc=$Rc-m=$m-pem=$positiveMode.dat
      fi
   done
done

