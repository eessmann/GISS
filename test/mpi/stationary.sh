#! /bin/sh
#
# Test MPI run with 4 process
# 

modules=$HOME/adds/lib/gerris

input=`mktemp /tmp/input.XXXXXX`

output=`mktemp /tmp/output.XXXXXX`
cat > $input <<EOF
1 2 FlDomain FlBox FlGEdge {} {
  FlTime { end = 1.0 }
  FlRefine 6
  GModule $modules/libperiodic_flow.so
  FlInitStationary {} 2
  FlOutputStationaryError { step = 0.05 } stdout 2
}
FlBox {}
1 1 0
1 1 2
EOF
../../src/gerris $input > $output

output_mpi=`mktemp /tmp/output.XXXXXX`
cat > $input <<EOF
4 8 FlDomain FlBox FlGEdge {} {
  FlTime { end = 1.0 }
  FlRefine 5
  GModule $modules/libperiodic_flow.so
  FlInitStationary {} 1
  FlOutputStationaryError { step = 0.05 } stdout 1
}
FlBox {pid=0}
FlBox {pid=1}
FlBox {pid=2}
FlBox {pid=3}
1 2 0
2 1 0
4 3 0
3 4 0
1 4 2
4 1 2
2 3 2
3 2 2
EOF
mpirun -np 4 ../../src/gerris $input > $output_mpi

rm -f $input

if ./compare.sh $output $output_mpi 1.; then
    rm -f $output $output_mpi
    exit 1
fi
rm -f $output $output_mpi
