# Takes one positional parameter that is going to be the 
# base name for the lock file. Returns "True" if the file 
# could be created and "False" otherwise.
function lock_or_die()
{
   local base=$1
   local lock="$base.lock"
   if [[ -f $lock ]] 
   then
      echo "Already computing or previous computation failed."
      echo "ABORTING!"
      false
   else
      touch $lock
      true
   fi
}

# Takes one positional parameter that is going to be the 
# base name for the lock file. Returns "True" if the file 
# could be removed and "False" otherwise.
function unlock_or_die()
{
   local base=$1
   local lock="$base.lock"
   if [[ -f $lock ]] 
   then
      rm $lock && true || false
   else
      true
   fi

}
