#!/bin/bash
#
# Xiang 'Anthony' Chen
#
# Evaluating voxelization algorithms

make
numTest=6
fileName="result-$(date +"%Y%m%d-%T").csv"

meshes=("../data/sphere/sphere.obj" "../data/teapot/teapot.obj" "../data/teddy/teddy.obj") # "../data/bunny/bunny.obj")
dims=(16 32 64)
mulrays=("-m" "-")
prescreen=("-p" "-")
repetitino=10

touch $fileName
echo number of triangles,voxel grid dimension,multiple rays,triangle prescreening,time, >> $fileName

cnt=1
total=$((${#meshes[@]} * ${#dims[@]} * ${#mulrays[@]} * ${#prescreen[@]}))

for mesh in ${meshes[@]}; do
	for dim in ${dims[@]}; do
		for m in ${mulrays[@]}; do
			for p in ${prescreen[@]}; do
				printf "[$cnt/$total] running ./voxelizer -d $dim $m $p -l $mesh null\n"
				./voxelizer -d $dim $m $p -l $mesh null >> $fileName
				cnt=$((cnt+1))
				# printf "finished $cnt / $total evaluation(s) ... \n"
			done
		done
	done
done

printf "\ndone! \n"

cat $fileName

mv $fileName ../performance-results/