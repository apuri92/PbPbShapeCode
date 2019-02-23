#!/bin/bash

# package="$scp_bnl:/direct/usatlas+u/apuri/Work/pPb_analysis/jetFragmentation_pPb_xAOD/pPbFragmentation"
package="apuri@acas0003:/direct/usatlas+u/apuri/Work/pPb_analysis/jetFragmentation_pPb_xAOD/pPbFragmentation"

#package="apuri@lxplus.cern.ch:/afs/cern.ch/user/a/apuri/private/jetFragmentation_pPb_xAOD/pPbFragmentation/"
#package="apuri@lxplus.cern.ch:/afs/cern.ch/user/a/apuri/private/tmp/jetFragmentation_pPb_xAOD/pPbFragmentation/"

echo $package
if [[ $1 = *"cxx" ]] ; then
	echo "sending $1 to Root"
	scp ./pPbFragmentation/Root/$1 $package/Root/
elif [[ $1 = *"h" ]] ; then
	echo "sending $1 to pPbFragmentation"
	scp ./pPbFragmentation/pPbFragmentation/$1 $package/pPbFragmentation/
elif [[ $1 = "util" ]] ; then
	echo "sending $1 to Util"
	scp ./pPbFragmentation/util/* $package/util/
elif [[ $1 = "cfg" ]] ; then
	echo "sending $1 to analysis"
	scp *sh $package/../
	scp *cfg $package/../
elif [[ $1 = "Link" ]] ; then
	echo "sending $1 to analysis"
	scp ./pPbFragmentation/Root/LinkDef.h $package/Root/
elif [[ $1 = "all" ]] ; then
	source send.sh "*.cxx"
	# source send.sh "Link"
	source send.sh "*.h"
	# source send.sh "util"
	# source send.sh "cfg"

else
	echo "Specify file to send"
fi
