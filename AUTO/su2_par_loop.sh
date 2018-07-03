#!/bin/bash


cd ./meshes

for i in $( ls ); do
	echo item: $i
	if [ -f ../airfoil.su2 ];
	then
		echo Removing old mesh...		
		rm ../airfoil.su2
	else
		echo No existing airfoil mesh...
	fi
	
	cp $i ../airfoil.su2
	cd ..

	# Wait for system to copy files
	sleep 1

	# Run SU2
	python $SU2_RUN/parallel_computation.py -f inv_airfoil.cfg -n 30
	#$SU2_RUN/SU2_CFD inv_airfoil_new.cfg

	# Copy Results	
	cp forces_breakdown.dat ./forces/$i.dat
	cp flow.vtk ./flows/$i.vtk
	cd ./meshes
done
