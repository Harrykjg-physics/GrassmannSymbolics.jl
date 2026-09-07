#!/bin/sh
starttime=`date +'%Y-%m-%d %H:%M:%S'`
echo "--------- Results --------"
julia --project=/home/jgkong/software/MyProject/Grassmann_new "$1" > compare.out 2>&1
rc=$?
echo "JULIA_EXIT:$rc" > compare.exit
endtime=`date +'%Y-%m-%d %H:%M:%S'`
start_seconds=$(date --date="$starttime" +%s)
end_seconds=$(date --date="$endtime" +%s)
echo "--------- Running Time -------"
echo "Running time of the scripts: "$((end_seconds-start_seconds))"s"
exit $rc