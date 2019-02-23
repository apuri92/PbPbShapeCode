echo "pbpb/pp	data/mc	perf/ff	config"
read input_data input_did input_type config


if [[ $input_type == "ff" ]]; then
	
	if [[ $input_data == "pbpb" && $input_did == "data" ]]; then
		file_name="FF_data_out_histo_*_5p02_r001.root"
		data_type="data_ff_pptight"
	fi

	if [[ $input_data == "pbpb" && $input_did == "mc" ]]; then
		file_name="FF_MC_out_histo_*_5p02_r001.root"
		data_type="mc_ff_pptight"
	fi

	if [[ $input_data == "pp" && $input_did == "data" ]]; then
		file_name="FF_data_out_histo_*_5p02_r001.root"
		data_type="pp_data_ff_pptight"
	fi

	if [[ $input_data == "pp" && $input_did == "mc" ]]; then
		file_name="FF_MC_out_histo_*_5p02_r001.root"
		data_type="pp_mc_ff_pptight"
	fi


elif [[ $input_type == "perf" ]]; then

	if [[ $input_data == "pbpb" && $input_did == "mc" ]]; then
		file_name="Perf_MC_JZ*_out_histo_*_5p02_r001.root"
		data_type="mc_perf_pptight"
	fi

	if [[ $input_data == "pp" && $input_did == "mc" ]]; then
		file_name="Perf_MC_JZ*_out_histo_*_5p02_r001.root"
		data_type="pp_mc_perf_pptight"
	fi

fi


files=/usatlas/u/apuri/Work/output/PbPbFFShape/processed_$data_type/config$config/$file_name
mkdir -p /Users/Akshat/Box\ Sync/Research/Analysis/jetFragmentation/PbPbShapeCode/output_dev/raw_results/sys$config/
local_area=/Users/Akshat/Box\ Sync/Research/Analysis/jetFragmentation/PbPbShapeCode/output_dev/raw_results/sys$config/
echo $files

scp $scp_bnl:$files "$local_area"



