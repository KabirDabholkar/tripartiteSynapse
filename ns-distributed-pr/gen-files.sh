i=15
echo $i
iter=$((i+20))

while read -a line
do
    vd="${line[3]}"
    #echo $vd
    sed "s/tic = 40/tic = $vd/; s/isi = \"20\"/isi = \"$i\"/; s/50e3/${iter}e3/" NSI20V40.mdl > NSI${i}V${vd}dist.mdl

done < pr-dist.dat
