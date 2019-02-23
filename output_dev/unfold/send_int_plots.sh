local_area="/Users/Akshat/Box Sync/Research/Analysis/jetFragmentation/PbPbShapeCode"
intnote_area="/Users/Akshat/Dropbox/trackjet_corr_intnote/"

corrections="/Users/Akshat/Dropbox/trackjet_corr_intnote/figures_corrections/"
general="/Users/Akshat/Dropbox/trackjet_corr_intnote/figures_general/"
performance="/Users/Akshat/Dropbox/trackjet_corr_intnote/figures_performance/"
results="/Users/Akshat/Dropbox/trackjet_corr_intnote/figures_results/"
systematics="/Users/Akshat/Dropbox/trackjet_corr_intnote/figures_systematics/"
UE="/Users/Akshat/Dropbox/trackjet_corr_intnote/figures_UE/"


#event acceptances
cp "$local_area"/output_dev/output_pdf/EventAccept_pp.pdf $general
cp "$local_area"/output_dev/output_pdf/EventAccept_PbPb.pdf $general


#shape response (pos cor factors)
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/ShapeResponse2D_PbPb.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/ShapeResponse2D_pp.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/pos_corr_factors_pp.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/pos_corr_factors_PbPb.pdf $corrections
#resp matrices
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/resp_matrix_ChPS_PbPb_MC.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/resp_matrix_ChPS_pp_MC.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/resp_matrix_jet_PbPb_MC.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/resp_matrix_jet_pp_MC.pdf $corrections
#jet spectra closure
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/spect_closure_pp_MC.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/spect_closure_PbPb_MC.pdf $corrections
#evolution plots
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/evol_PbPb_MC.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/evol_PbPb_data.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/evol_pp_MC.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/evol_pp_data.pdf $corrections
#final chps as function of track pt and r
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/ChPS_final_PbPb_MC.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/ChPS_final_pp_MC.pdf $corrections

cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/RatioProj_PbPb.pdf $corrections
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/RatioProj_pp.pdf $corrections





#final injet as function of track pt and r
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/ChPS_final_inJet_pp_data.pdf $results
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/ChPS_final_inJet_PbPb_data.pdf $results
#UE plots
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/ChPS_dR_UE_PbPb_data.pdf $results
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/ChPS_dR_UE_PbPb_MC.pdf $results
#RDptR
cp "$local_area"/output_dev/unfold/output_pdf_nominal/ChPS_final_ratio_dR_data.pdf $results
cp "$local_area"/output_dev/unfold/output_pdf_nominal/pp/ChPS_final_inJet_pp_MC.pdf $results
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/ChPS_final_inJet_PbPb_MC.pdf $results


cp "$local_area"/output_dev/unfold/output_pdf_nominal/conf/ChPS_final_dR_CONF_DpT_data_jet9_cent*.pdf $results
cp "$local_area"/output_dev/unfold/output_pdf_nominal/conf/ChPS_final_dR_CONF_DpT_data_jet7_cent*.pdf $results
cp "$local_area"/output_dev/unfold/output_pdf_nominal/conf/RDpT_final_ratio_dR_CONF_data_trk_cent*.pdf $results
cp "$local_area"/output_dev/unfold/output_pdf_nominal/conf/RDpT_final_ratio_dR_CONF_data_jet9_cent*.pdf $results
cp "$local_area"/output_dev/unfold/output_pdf_nominal/conf/RDpT_final_ratio_dR_CONF_data_jet7_cent*.pdf $results


#systematics
cp "$local_area"/output_dev/unfold/output_pdf_nominal/systematics/ChPS_dR_sys_PbPb_error.pdf $systematics
cp "$local_area"/output_dev/unfold/output_pdf_nominal/systematics/ChPS_dR_sys_pp_error.pdf $systematics
cp "$local_area"/output_dev/unfold/output_pdf_nominal/systematics/RDpT_dR_sys_error.pdf $systematics

#UE factors
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/UE_factors.pdf $UE
cp "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/UE_factors_r.pdf $UE
#B2S
cp  "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/ChPS_B2S_PbPb_data.pdf $UE
cp  "$local_area"/output_dev/unfold/output_pdf_nominal/PbPb/ChPS_B2S_PbPb_MC.pdf $UE




