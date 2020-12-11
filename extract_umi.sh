#!/bin/bash
#### Part of TRAD-Seq from Devin King devking@stanford.edu ####

[ $# -ge 1 -a -f "$1" ] && input="$1" || input="-"

awk '{
      if(substr($1,1,1)=="@")
      {
        print $0
      } 
      else
      {
        split($2,a,":");$2=""; print $0,"RX:Z:"a[4]
      } 
      
    }' OFS="\t" $input

