function auto {
	echo "dataset_type: $1" > auto.cfg
	echo "isMC: $2" >> auto.cfg
	if [[ $mode == "unfold" ]]; then
		./analysis auto.cfg
	fi
	root -b -q "draw_ChPS.c(\"auto.cfg\")"
	rm auto.cfg
}

mode=$1 #if running full unfolding as well

# auto PbPb 1
auto PbPb 0
# auto pp 1
auto pp 0
