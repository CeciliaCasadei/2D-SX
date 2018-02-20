#!/bin/bash

cd ./Output_r0195/UnassembledImageProcessing 
for i in *.txt
	
	do 
	echo $i
	diff -q $i /mnt/das-gpfs/home/casadei_c/work/casadei/Backup_tilted_Nov2016/Output_r0195/UnassembledImageProcessing/$i
done


