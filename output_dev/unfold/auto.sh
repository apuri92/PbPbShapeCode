function auto {
	echo "dataset_type: $1" > auto.cfg
	echo "isMC: $2" >> auto.cfg	

	echo "diagnostic_mode: true" >> auto.cfg
	echo "apply_UE_uncert: 1" >> auto.cfg
	echo "sys_mode: 102" >> auto.cfg
	echo "verbose: 0" >> auto.cfg

	if [[ $mode == "unfold" ]]; then
		echo "getting UE and posCorr factors"
		root -b -q "get_posCorr.c(\"auto.cfg\")" 
		echo $dataset_type 
		if [[ $dataset_type=="PbPb" ]]; then
			root -b -q "UE_factors.c(\"auto.cfg\")" 
		fi

		./analysis auto.cfg

	fi
	root -b -q "draw_ChPS.c(\"auto.cfg\")"
	root -b -q "draw_spectra.c(\"auto.cfg\")"
	rm auto.cfg
}



mode=$1 #if running full unfolding as well

auto PbPb 1
auto pp 1
root -b -q "comp_ChPS.c(1)"

auto PbPb 0
auto pp 0
root -b -q "comp_ChPS.c(0)"
