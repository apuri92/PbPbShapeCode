function run {
	echo "mode: $1" > tmp.cfg
	echo "dataset_type: $2" >> tmp.cfg
	echo "isMC: $3" >> tmp.cfg
	echo "verbose: 1" >> tmp.cfg
	
	if [[ $mode == "get"  ]]; then
		root -b -q "systematics.c(\"tmp.cfg\")"
	fi
	root -b -q "draw_sys_err.c(\"tmp.cfg\")"

	rm tmp.cfg
}

mode=$1

run ChPS PbPb 0
run RDpT x 0
run ChPS pp 0


# root -b -q "draw_conf_plots.c(\"tmp.cfg\")"
