cd cards
for i in H_*_SMEFT;
do
    cp -r ${i} ${i}sim-3pars
    rm ${i}sim-3pars/reweight_card.dat ${i}sim-3pars/param_card.dat
    sed -i 's/SMEFTsim_A_U35_MwScheme_UFO_v2_1/SMEFTsim_U35_MwScheme_UFO/g' ${i}sim-3pars/proc_card.dat
    sed -i "s/$i/${i}sim-3pars/g" ${i}sim-3pars/proc_card.dat
    sed -i "s/JetMatching:merge = on/JetMatching:merge = off/g" ${i}sim-3pars/pythia8_card.dat
done
