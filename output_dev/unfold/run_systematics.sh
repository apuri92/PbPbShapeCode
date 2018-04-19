systematics=(105) # 101 102 106 107 108 200)


function unfold_draw_DpT {
	echo "dataset_type: $1" > unfold_draw_DpT_auto.cfg
	echo "isMC: $2" >> unfold_draw_DpT_auto.cfg
	echo "sys_mode: $3" >> unfold_draw_DpT_auto.cfg
	echo "verbose: 0" >> unfold_draw_DpT_auto.cfg
	echo "draw_mode: 0" >> unfold_draw_DpT_auto.cfg
	
	if [[ $mode == "unfold" ]]; then
		echo "getting UE and posCorr factors"
		root -b -q "get_posCorr.c(\"unfold_draw_DpT_auto.cfg\")" 
		root -b -q "UE_factors.c(\"unfold_draw_DpT_auto.cfg\")" 
		./analysis unfold_draw_DpT_auto.cfg
	fi

	root -b -q "draw_ChPS.c(\"unfold_draw_DpT_auto.cfg\")"
	root -b -q "draw_spectra.c(\"unfold_draw_DpT_auto.cfg\")"

	rm unfold_draw_DpT_auto.cfg
}

function get_RDpT {
	echo "isMC: $1" > RDpT_config_auto.cfg
	echo "sys_mode: $2" >> RDpT_config_auto.cfg
	echo "verbose: 0" >> RDpT_config_auto.cfg
	echo "draw_mode: 0" >> RDpT_config_auto.cfg
	root -b -q "comp_ChPS.c(\"RDpT_config_auto.cfg\")"

	rm RDpT_config_auto.cfg

}

mode=$1

for i in ${systematics[@]}; do
	echo "************************************"
	echo "****** Running Systematic $i *******"
	echo "************************************" 
	
	unfold_draw_DpT PbPb 0 $i
	unfold_draw_DpT pp 0 $i
	get_RDpT 0 $i

	unfold_draw_DpT PbPb 1 $i
	unfold_draw_DpT pp 1 $i
	get_RDpT 1 $i
done


