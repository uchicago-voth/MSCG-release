#!/bin/bash
MAX_ITERATION=$1
REF_TRAJECTORY=$2
CG_TRAJECTORY=$3
LAMMPS=$4

ITERATION=1
PREV_ITERATION=0

echo "MAX_ITERATION $MAX_ITERATION"
echo "REF_TRAJECTORY $REF_TRAJECTORY"
echo "CG_TRAJECTORY $CG_TRAJECTORY"
echo "LAMMPS $LAMMPS"


until [ $ITERATION -gt $MAX_ITERATION ]; do
    echo $ITERATION
    mkdir iter_$ITERATION
    let PREV_ITERATION=$ITERATION-1
    cp iter_$PREV_ITERATION/b-spline.out iter_$ITERATION/b-spline-previous.out
    cp iter_$PREV_ITERATION/control.in iter_$ITERATION/control.in
    cp iter_$PREV_ITERATION/rmin.in iter_$ITERATION/rmin.in
    cp iter_$PREV_ITERATION/rmin_b.in iter_$ITERATION/rmin_b.in
    cp iter_$PREV_ITERATION/top.in iter_$ITERATION/top.in
    cp iter_$PREV_ITERATION/newrem.x iter_$ITERATION/newrem.x
    cp iter_$PREV_ITERATION/CG.in iter_$ITERATION/CG.in
    cp iter_$PREV_ITERATION/CG.data iter_$ITERATION/CG.data
    cp iter_$PREV_ITERATION/$LAMMPS iter_$ITERATION/$LAMMPS
    cd iter_$ITERATION/
    ./newrem.x -l_ref ../$REF_TRAJECTORY -l_cg ../iter_$PREV_ITERATION/$CG_TRAJECTORY
    cat fmec* > fmec.table
    mpirun ./$LAMMPS < CG.in
    cd ..
    let ITERATION=$ITERATION+1
done


#//pseudo code

#//mkdir for the current iteration
#//copy basis set coefficients from n-1 iteration
#//copy other inputs for newrem.x code
#//run the newrem.x code
#//format the outputs from newrem so can be used by LAMMPS
#//run the new cg model

#//requirements

#//a zereoth iteration directory with a cg trajectory and basis set coefficients used to create it. 
#//reference trajectory
#//control.in, top.in, rmin.in rmin_b.in
#//code for formatting the newrem.x outputs
#//


