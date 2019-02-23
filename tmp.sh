for i in 53 57 41 59; do
   echo $i
   root -b -q "check_UE.c($i)"
done
