#!/bin/bash

#for I in {1..$1}
#do
#	cd $I/
#	bsub < job$I
#	cd ..
#done

 i=1  
 while [ $i -lt $1+1 ]  
 do  
     echo $i
	 cd $i/
	 bsub < job$i
	 cd ../
     let i++  
 done  

