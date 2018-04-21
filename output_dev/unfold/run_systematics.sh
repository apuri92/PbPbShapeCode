function run {
	echo "mode: $1" > tmp.cfg
	echo "dataset_type: $2" >> tmp.cfg
	echo "isMC: $3" >> tmp.cfg
	echo "verbose: 1" >> tmp.cfg
	root -b -q "systematics.c(\"tmp.cfg\")"
	# root -b -q "draw_conf_plots.c(\"tmp.cfg\")"
}

run RDpT x 0
run ChPS pp 0
run ChPS PbPb 0


root -b -q "draw_conf_plots.c(\"tmp.cfg\")"
