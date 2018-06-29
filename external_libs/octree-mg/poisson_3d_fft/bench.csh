#!/bin/tcsh

if ( $# < 4 ) then
    echo usage: `basename $0` start step end  geocode datacode
    exit 1
endif

set start = $1
set step = $2
set end = $3
set scf = 0
set per = 128
set geocode = $4
set datacode = $5
set h = $start

#echo $geocode

make PSolver

if ($geocode == 'P') then


    rm acc_P.$start-{$end}.txt
while ($h <= $end)

    rm fort.99
    ompirun -np 2 PSolver $h $h $h $scf P $datacode

    head -1 fort.99 >> acc_P.$start-{$end}.txt

    @ h = $h + $step

end

else if ($geocode == 'S') then

    rm acc_S-$per.$start-{$end}.txt
while ($h <= $end)

    rm fort.99
    ompirun -np 2 PSolver $per $h $per $scf S $datacode

    head -1 fort.99 >> acc_S-$per.$start-{$end}.txt

    @ h = $h + $step

end

else if ($geocode == 'F') then

    rm acc_F.$start-{$end}_{$scf}.txt
while ($h <= $end)

    rm fort.99
    ompirun -np 2 PSolver $h $h $h $scf F $datacode

    head -1 fort.99 >> acc_F.$start-{$end}.txt

    @ h = $h + $step

end

endif
