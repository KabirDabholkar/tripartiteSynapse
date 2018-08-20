i=400
iter=$((i+30))
echo $iter
while read -a line
do
    vd="${line[3]}"
    echo $vd
    sed "s/tic = 40/tic = $vd/; s/isi = \"20\"/isi = \"$i\"/; s/50e3/${iter}e3/" RSI20V40.mdl > RSI${i}V${vd}dist.mdl

done < pr-dist.dat
