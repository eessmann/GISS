#! /bin/sh
#
# Tests the fl_domain_split() function
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

input_split=`mktemp /tmp/input.XXXXXX`
../../src/gerris -s 1 $input > $input_split
output_split=`mktemp /tmp/output.XXXXXX`
../../src/gerris $input_split > $output_split
rm -f $input $input_slit

if ./compare.sh $output $output_split 0.02; then
    rm -f $output $output_split
    exit 1
fi
rm -f $output $output_split
