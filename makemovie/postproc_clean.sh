#!/bin/bash
 
count=0
step=1/10
start=0
end=10
#count=6    # count*step=starting time
#step=2     # time inteval between .tga frames
#start=0    # starting .tga frame
#end=1000   # end .tga frame

echo $count

directory=`pwd`

for ((i=$start;i<=$end;i++))
do
  if (( $i < 10 )); then
    file="000"$i.tga
  elif (( $i < 100 )); then
    file="00"$i.tga
  elif (( $i < 1000 )); then
    file="0"$i.tga
  else
    file=$i.tga
  fi
  filename=${file%}      #  Strip ".mac" suffix off filename
                            #+ ('.*c' matches everything
			    #+ between '.' and 'c', inclusive).

  num1=$(( $count * $step ))
  cero=0
  num2=$num1$cero
        num3=$(( $count - $num2))
  num4=$num1.$num3
  if (( "$step" >= 1 )); then
    num4=$num1
  fi
  echo $filename
  convert -fill white -draw "rectangle 10,5 100,25" $file $filename 
  convert -draw "text 20,20 '$num4 (fs)' " $file $filename

  count=$(($count + 1))


done



