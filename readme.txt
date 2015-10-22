
1  generate gnuplot scripts of K.E., P.E., T.E. and T as a function of time

2  generate gnuplot scripts of K.E., P.E., and T.E. as a function of time

3  *.q3d -> *.xyz

4* which bonds to plot the bond order

5* get accelerations, velocities and coordinates of center of mass of each monomer

6* generate collision.dynput from the results of relaxation: coordinates, velocity and so on.

7  the script to run bsub mopdyn13 thermal.dynput, followed by mopdyn13 collision.dynput

8  generate dynput with a range of parameters in different folder

9  y_0

================================================
0 thermal the monomer for 20ps, and then relax the monomer for 20ps

1 put relax.dynput/coord/veloc... in a folder, and we name that folder ./relax/
    !!!the file name must be relax.*, the folder name should be any one.!!!

2 put *.sample in ./sample/,  autoini *.sample working directory: the same directory as *.sample
  [generate ColDynGen_xx.ini files]
  and then cp ColDynGen_* ../relax/

3 if you want to debug, mopdynapp 1 xx with LColDynGenDebug=.True. to generate debug files in ./relax/Debug/

4 dynputscript in directory ./
  [mopdynapp 1 xx to generate collision_xx.dynput in folder ./relax/]
    !!!set the pathmopdynapp='' before install the program!!!

5 mvdynput in directory ./
  [mkdir ./xx/  mv collision_xx.dynput ColDynGen_xx.ini jobxx ./xx/]

6 job.sh or ...  ./
  [mopdyn13 collision_xx.dynput]

7 after MD simulation

8 plot.sh
  [mopdynapp 2 xx to generate bond order, CM information and mdplot d to generate gnuplot scripts in ./xx/result/]
  [gnuplot *.plot, and mopdyn2xyz *.coord *.plot]
    !!!the monomer should be an integer multiply of 3 in bondorder calculation!!!
    working directory ./xx/
