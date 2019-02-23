# package="$scp_bnl:/direct/usatlas+u/apuri/Work/pPb_analysis/jetFragmentation_pPb_xAOD/pPbFragmentation/"
package="apuri@acas0003:/direct/usatlas+u/apuri/Work/pPb_analysis/jetFragmentation_pPb_xAOD/pPbFragmentation/"

#package="$cern:/afs/cern.ch/user/a/apuri/private/tmp/jetFragmentation_pPb_xAOD/pPbFragmentation/"
echo $package
scp $package/../testing_$1/hist*root .
#scp $package/../JZ*/hist*root .
