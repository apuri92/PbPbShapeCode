
systematics=(2000 2001 2006 2007 2008 2009 2010 2011 2014 2015 2016 2017 2018 2019 2020 2021 2022 2023)
for i in ${systematics[@]}; do
	echo "************************************"
#	source get_full_data.sh pbpb data ff $i
	source get_full_data.sh pbpb mc ff $i

#	source get_full_data.sh pp data ff $i
	source get_full_data.sh pp mc ff $i

	# cp output_dev/raw_results/c1000/FF_data_out_histo_pp_5p02_r001.root output_dev/raw_results/c$i/FF_data_out_histo_pp_5p02_r001.root 
	# cp output_dev/raw_results/c1000/FF_MC_out_histo_pp_5p02_r001.root output_dev/raw_results/c$i/FF_MC_out_histo_pp_5p02_r001.root 

done
