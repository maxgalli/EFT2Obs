#!/usr/bin/env bash
EFT2OBS_DIR=/afs/cern.ch/work/g/gallim/MGStudies/EFT2Obs_MattSetup3
export LHAPDF_CONFIG_PATH="${EFT2OBS_DIR}/lhapdf/bin/lhapdf-config"
export PYTHONPATH="${EFT2OBS_DIR}/lhapdf/lib64/python2.7/site-packages:${PYTHONPATH}"
export PYTHONPATH="${EFT2OBS_DIR}/scripts:${PYTHONPATH}"
export PYTHONPATH="${EFT2OBS_DIR}/additional_scripts:${PYTHONPATH}"
export LD_LIBRARY_PATH="${EFT2OBS_DIR}/lhapdf/lib:${LD_LIBRARY_PATH}"
export RIVET_ANALYSIS_PATH=${EFT2OBS_DIR}/RivetPlugins
export MG_DIR="MG5_aMC_v2_6_7"
export MG_TARBALL="MG5_aMC_v2.6.7.tar.gz"
export RIVET_VERSION="3.0.1"
export DEBUG_SCRIPTS=0

if [ -f "${EFT2OBS_DIR}/local/rivetenv.sh" ]; then
	source ${EFT2OBS_DIR}/local/rivetenv.sh
fi

if [ "$DEBUG_SCRIPTS" -eq "1" ]; then
	set -x
fi

#[[ ":$PYTHONPATH:" != *":$EFT2OBS_DIR/${MG_DIR}:"* ]] && PYTHONPATH="$EFT2OBS_DIR/${MG_DIR}:${PYTHONPATH}"
