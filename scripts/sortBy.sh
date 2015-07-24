# Function that sorts strings obeying a certain template according to 
# a certain substring. The substring must be of the form <par>=<val>.
# It takes two positional parameters:
# 1. the parameter string, <par> to sort by;
# 2. the template string <par> apears in.
function sortBy() 
{    
   local par=$1;    
   local template=$2;    
   local loc=`echo ${template}|sed -e 's/=/\n/g'|nl|awk '/'${par}'/ {print $1+1}'`;    
   sort -t '=' -k${loc},${loc} -g ;
}
