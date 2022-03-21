#!/usr/bin/env bash

set -e

proc=$1
rw_num=$2
nevents=$3
ncores=$4

#make mgrunscript
pushd /afs/cern.ch/work/m/mknight/private/EFT/topU3l_prelim/EFT2Obs
        #where launch 0001 is
	start_line=$(grep -n '0001' cards/${proc}/reweight_card.dat | sed 's/:.*//')

	start=0001
	end=0002
	n_par=$(sed -n "/${start}/, /${end}/{ /${start}/! { /${end}/! p } }" cards/${proc}/reweight_card.dat | wc -l)

	#where the 'set SMEFT 2 0.1' commands for this reweighting point start
	start=$(echo "$start_line + 1 + ($rw_num - 1) * ($n_par + 1)" | bc)
	end=$(echo "$start_line - 1 + $rw_num * ($n_par + 1)" | bc)

	{
	  echo "shower=OFF"
	  echo "reweight=OFF"
	  echo "done"
	  echo "set gridpack True"
	  echo "set nevents $nevents"
	} > ${TMPDIR}/mgrunscript
	sed -n "${start},${end}p" cards/${proc}/reweight_card.dat >> ${TMPDIR}/mgrunscript
	echo "done" >> ${TMPDIR}/mgrunscript
popd

RUNLABEL="rw_${rw_num}_nevents_${nevents}"

mkdir -p working_${RUNLABEL}
cp /afs/cern.ch/work/m/mknight/private/EFT/topU3l_prelim/EFT2Obs/proc_dir_${proc}.tar.gz working_${RUNLABEL}/

pushd working_${RUNLABEL}
  tar -xzf proc_dir_${proc}.tar.gz
  pushd ${proc}
    mkdir -p Events/${RUNLABEL}

    if [ "$ncores" -gt "0" ]; then
  	./bin/generate_events ${RUNLABEL} --nb_core="${ncores}" < ${TMPDIR}/mgrunscript > Events/${RUNLABEL}/generation.log
    else
	./bin/generate_events ${RUNLABEL} < ${TMPDIR}/mgrunscript > Events/${RUNLABEL}/generation.log
    fi

    rm Events/${RUNLABEL}/unweighted_events.lhe.gz
    
    mkdir -p /afs/cern.ch/work/m/mknight/private/EFT/topU3l_prelim/EFT2Obs/directResults/${proc}
    rm -rf /afs/cern.ch/work/m/mknight/private/EFT/topU3l_prelim/EFT2Obs/directResults/${proc}/${RUNLABEL}
    cp -r Events/${RUNLABEL} /afs/cern.ch/work/m/mknight/private/EFT/topU3l_prelim/EFT2Obs/directResults/${proc}/${RUNLABEL}
  popd
popd

rm -r working_${RUNLABEL}
