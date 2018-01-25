for OUTPUT in $(ls /mnt/das-gpfs/home/casadei_c/work/casadei/UNIX_@_LCLS/r0195-images/data1)
do
  echo $OUTPUT
  h5diff -v /mnt/das-gpfs/home/casadei_c/work/casadei/UNIX_@_LCLS/r0195-images/data1/$OUTPUT /afs/psi.ch/group/0620/casadei/2D-MX/UNIX_@_LCLS/r0195-images/data1/$OUTPUT
done
#h5diff -v /mnt/das-gpfs/home/casadei_c/work/casadei/UNIX_@_LCLS/r0199-images/data1/LCLS_2013_Nov04_r0199_120223_1f72b.h5 /afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5
#h5diff -v /mnt/das-gpfs/home/casadei_c/work/casadei/Geometry/geometry.h5 /afs/psi.ch/group/0620/casadei/2D-MX/Geometry/geometry.h5
