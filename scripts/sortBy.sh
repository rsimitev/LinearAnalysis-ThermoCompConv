# Function that sorts strings obeying a certain template according to 
# a certain substring.
function sortBy() 
{    
   local par=$1;    
   local template=$2;    
   local loc=`echo ${template}|sed -e 's/=/\n/g'|nl|awk '/'${par}'/ {print $1+1}'`;    
   sort -t '=' -k${loc},${loc} -g ;
}
