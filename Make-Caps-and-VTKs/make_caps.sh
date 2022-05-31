prefix=gld422

pid=$(ls pid*.cfg | awk '{ print $NF}')
echo $pid

time_file=gld422.timese
n=$(wc -l ${time_file} | awk '{print $1}')

echo $n

list=$(awk '{ print $1 }' ${time_file} | awk -f transpose.awk)

echo $list $pid
./autocombine.py localhost $pid $list

