# echo "pbpb/pp	data/mc	perf/ff	config"
# read input_data input_did input_type config

input_data=$1
input_did=$2
input_type=$3
config=$4

# echo "data or mc?"
# read input_did

# echo "perf or ff?"
# read input_type

# echo "config?"
# read config

# input_data=$1 #pbpb or pp
# input_did=$2 #data or MC
# input_type=$3 #perf of ff
# config=$4


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

# echo $file_name
# echo $data_type
# echo $config

mkdir -p /Users/Akshat/Box\ Sync/Research/Analysis/jetFragmentation/PbPbShapeCode/output_dev/raw_results/c$config

files=/atlasgpfs01/usatlas/data/apuri/output/PbPbFFShape//processed_$data_type/config$config/$file_name
local_area=/Users/Akshat/Box\ Sync/Research/Analysis/jetFragmentation/PbPbShapeCode/output_dev/raw_results/c$config
echo -e "$files \n -> \n $local_area"

# scp $scp_bnl:$files "$local_area"
scp apuri@acas0003:$files "$local_area"
# scp $scp_bnl:$processed_pp_mc_perf_pptight/config$config/Perf_MC_JZ*_out_histo_pp_5p02_r001.root output_dev/raw_results/
# scp $scp_bnl:$processed_mc_ff_pptight/config4/FF_MC_out_histo_PbPb_5p02_r001.root output_dev/raw_results/
# scp $scp_bnl:$processed"_pp_mc_ff_pptight/config"$config"/FF_MC_out_histo_pp_5p02_r001.root" ./output_dev/raw_results/x.root





