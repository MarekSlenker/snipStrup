#!/bin/bash

# Author: Marek Šlenker
# License: GNU General Public License 3.0, https://www.gnu.org/licenses/gpl-3.0.html



SAMPLEPLOIDYFILE="$1"


MAXPLOIDY="0"

while read SAMPLE; do 
	read PLOIDY;
	
	if [ "$PLOIDY" -gt "$MAXPLOIDY" ]; then
	MAXPLOIDY="$PLOIDY"
	fi

done < "$SAMPLEPLOIDYFILE"

# return MAXPLOIDY value
echo "$MAXPLOIDY" 

exit

