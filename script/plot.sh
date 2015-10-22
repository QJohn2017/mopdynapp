#!/bin/bash

#for I in {1..$1}
#do
#	cd $I/
#	bsub < job$I
#	cd ..
#done

 i=1  
 while [ $i -lt $1 ]  
 do  
     echo $i
	 cd $i/
	 mopdynapp 2 $i
	 mdplot d
	 cd result/
	 gnuplot *.plot
	 echo '**************'$i'*****************'
	 cat collision_$i.dob
	 echo '************************************'
	 cd ../
	 mopdyn2xyz collision_$i.coord ./result/collision_$i.xyz
	 cd ../
     let i++  
 done  

