	
# systematics=(0 1 2 4 5 6 7 8 9 10 11 14 15 16 17 18 19 20 21 22 23 42 44 45) 
systematics=(0) 


function unfold_draw_DpT {
	echo "dataset_type: $1" > unfold_draw_DpT_auto.cfg
	echo "isMC: $2" >> unfold_draw_DpT_auto.cfg
	echo "sys_mode: $3" >> unfold_draw_DpT_auto.cfg
	echo "verbose: 1" >> unfold_draw_DpT_auto.cfg
	echo "draw_mode: 1" >> unfold_draw_DpT_auto.cfg

	if [[ $mode == "unfold" ]]; then
		echo "getting UE and posCorr factors"
		root -b -q "get_posCorr.c(\"unfold_draw_DpT_auto.cfg\")"
		root -b -q "UE_factors.c(\"unfold_draw_DpT_auto.cfg\")"
		./analysis unfold_draw_DpT_auto.cfg
	fi

	root -b -q "draw_ChPS.c(\"unfold_draw_DpT_auto.cfg\")"
	root -b -q "draw_spectra.c(\"unfold_draw_DpT_auto.cfg\")"

	# rm unfold_draw_DpT_auto.cfg
}

function get_RDpT {
	echo "isMC: $1" > RDpT_config_auto.cfg
	echo "sys_mode: $2" >> RDpT_config_auto.cfg
	echo "verbose: 0" >> RDpT_config_auto.cfg
	echo "draw_mode: 1" >> RDpT_config_auto.cfg
	root -b -q "comp_ChPS.c(\"RDpT_config_auto.cfg\")"

	rm RDpT_config_auto.cfg

}

mode=$1

for i in ${systematics[@]}; do
	echo "************************************"
	echo "****** Running Version $i *******"
	echo "************************************"

	sys_path="c"
	if [[ $i -lt 51 && $i -gt -1 ]]; then
		sys_path="sys"
	fi

	mkdir -p output_pdf_$sys_path$i/root
	mkdir -p output_pdf_$sys_path$i/PbPb
	mkdir -p output_pdf_$sys_path$i/pp


	unfold_draw_DpT PbPb 0 $i
	# unfold_draw_DpT pp 0 $i
	get_RDpT 0 $i

	unfold_draw_DpT PbPb 1 $i
	# unfold_draw_DpT pp 1 $i
	get_RDpT 1 $i


done




