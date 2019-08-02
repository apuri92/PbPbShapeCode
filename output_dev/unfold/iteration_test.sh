# systematics=(0 1 2 4 5 6 7 8 9 10 11 14 15 16 17 18 19 20 21 22 23 42 44 45) 
n_iter=(0 1 2 3 4 5 6 7 8 9) 


function unfold_draw_DpT {
	echo "dataset_type: $1" > iter_test.cfg
	echo "isMC: $2" >> iter_test.cfg
	echo "sys_mode: $3" >> iter_test.cfg
	echo "verbose: 1" >> iter_test.cfg
	echo "draw_mode: 0" >> iter_test.cfg
	echo "n_unfold: $4" >> iter_test.cfg

	if [[ $4 == 0 ]]; then
		echo "getting UE and posCorr factors"
		# root -b -q "get_posCorr.c(\"iter_test.cfg\")"
		# root -b -q "UE_factors.c(\"iter_test.cfg\")"
	fi

	# ./analysis iter_test.cfg
	root -b -q "draw_ChPS.c(\"iter_test.cfg\")"
}

for i in ${n_iter[@]}; do
	echo "************************************"
	echo "****** Running Iteration $i *******"
	echo "************************************"
	unfold_draw_DpT PbPb 0 0 $i
	echo "                                    "

done