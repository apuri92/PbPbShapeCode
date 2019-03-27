
systematics=(1005) #10 1011 1014 1015)
for i in ${systematics[@]}; do
	echo "************************************"
	source get_full_data.sh pbpb data ff $i
	source get_full_data.sh pbpb mc ff $i

	source get_full_data.sh pp data ff $i
	source get_full_data.sh pp mc ff $i

	# cp output_dev/raw_results/c1000/FF_data_out_histo_pp_5p02_r001.root output_dev/raw_results/c$i/FF_data_out_histo_pp_5p02_r001.root 
	# cp output_dev/raw_results/c1000/FF_MC_out_histo_pp_5p02_r001.root output_dev/raw_results/c$i/FF_MC_out_histo_pp_5p02_r001.root 

done
