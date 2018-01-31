#root -b -q "combine_samples.c(isPbPb,cut)"
#root -b -q "draw_eff_trketa.c(isPbPb,cut)"


root -b -q "combine_samples.c(1, \"ppTight\")"
# root -b -q "draw_eff_trketa.c(1, \"ppTight\")"
root -b -q "draw_eff_jetpt_jety.c(\"ppTight\")"
