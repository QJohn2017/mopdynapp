compiler = gfortran
flag = -O2

objects = mopdynapp.o collisiondynput.o rotation.o velocity.o masscendyn.o bondorder.o
objectscom = gaussiancom.o collisiondynput.o rotation.o velocity.o

all:
	make mopdynapp
	make gaussiancom
	make mvdynput
	make dynputscript
	make autoini
	make mdplot

mopdynapp: $(objects)
	$(compiler) $(flag) $(objects) -o $@ 

gaussiancom: $(objectscom)
	$(compiler) $(flag) $(objectscom) -o $@
	
mvdynput: mvdynput.o
	$(compiler) $(flag) mvdynput.o -o $@
	
dynputscript: dynputscript.o
	$(compiler) $(flag) dynputscript.o -o $@
	
autoini: autoini.o
	$(compiler) $(flag) autoini.o -o $@
	
mdplot: plot.o
	$(compiler) $(flag) plot.o -o $@
	    
%.o: %.f90
	$(compiler) $(flag) -c $*.f90

#help:
	
clean: 
	rm -f *.o core *.mod
