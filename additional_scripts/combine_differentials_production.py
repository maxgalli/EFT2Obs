"""
Run from /afs/cern.ch/work/g/gallim/MGStudies/EFT2Obs_MattSetup3 using python3
"""
import numpy as np
import json
from copy import copy

# Values for 125.9 reported in the xlsx file reported here https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG#Higgs_cross_sections_and_decay_b
cross_sections = {
        "ggH_SMEFTatNLO": 4.796e01,
        "ggH_SMEFTsim_topU3l": 4.796e01,
        "qqH_SMEFTsim_topU3l": 1.347,
        "bbH_SMEFTsim_topU3l": 4.757e-01,
        "ttH_SMEFTsim_topU3l": 4.971e-01,
        "WH_lep_tau_SMEFTsim_topU3l": 5.845e-02 + 9.220e-02,
        "ZH_lep_tau_SMEFTsim_topU3l": 2.924e-02 + 1.742e-01,
        "tHq_SMEFTsim_topU3l": 7.370e-02,
        "tHW_SMEFTsim_topU3l": 2.835e-03
}

obs_dc_dict = {
    "pt_h": ["Hgg", "HZZ", "HWW", "Htt", "HttBoosted", "HbbVBF"],
    #"pt_h": ["HbbVBF"],
    #"deltaphijj": ["Hgg", "HZZ"]
}

histograms_base_dir = "condor/pass2/"

eq_base_dir = "ConvertedEquations/"


for obs in obs_dc_dict:
    print("\n\nObservable: {}".format(obs))
    for dc in obs_dc_dict[obs]:
        print("\nDecay channel: {}".format(dc))
        output_dct = {} # keys are edges

        with open("{}ggH_SMEFTatNLO_{}_combined_{}.json".format(eq_base_dir, dc, obs), "r") as f:
            dct = json.load(f)
        bins = list(dct.keys())
        bins = [float(b) for b in bins]
        bins.sort()
        bins = [str(b) for b in bins]
        print(bins)

        xs_pieces = {}

        # this has to be done because, due to che chg-chgtil fits, we have to merge both ggH@NLO and ggH-SMEFTsim with everything else
        xs_list = [copy(list(cross_sections.keys())), copy(list(cross_sections.keys()))]
        xs_list[0].remove("ggH_SMEFTsim_topU3l")
        xs_list[1].remove("ggH_SMEFTatNLO")

        for cross_sections_names in xs_list:
            #print("Cross section names: {}".format(cross_sections_names))
            for i_bin, edge in enumerate(bins):
                print(i_bin, edge)
                output_dct[edge] = {}
                # Make denominator
                #print("Making denominator")
                den = 0
                factors = {} # dictionary of factors for a certain bin i_bin, keys are the production modes
                obs_json = obs
                xs_pieces[edge] = {}
                for xs_name in cross_sections_names:
                    if obs == "pt_h":
                        if dc == "Hgg" and xs_name != "ggH_SMEFTsim_topU3l":
                            obs_json = "pt_gg"
                    rivet_json_name = "{}_{}/Rivet.json".format(xs_name, dc)
                    if xs_name == "ggH_SMEFTatNLO":
                        rivet_json_name = "ggH_SMEFTatNLO_{}_loop/Rivet.json".format(dc)
                    #print("Will open file inside {}".format(rivet_json_name))
                    to_open = "{}{}".format(histograms_base_dir, rivet_json_name)
                    #print("Will open file {}".format(to_open))
                    with open(to_open) as f:
                        rivet_dct = json.load(f)
                    try:
                        #print(rivet_dct["{}_active_bins".format(obs_json)])
                        active_left_edges = [str(tpl[0]) for tpl in rivet_dct["{}_active_bins".format(obs_json)]]
                        #print(active_left_edges)
                        i_active_edge = active_left_edges.index(edge)
                        #print(i_active_edge)
                        try:
                            sigma_ij_mg = rivet_dct["{}[rw0000]".format(obs_json)][i_active_edge][0]
                        except KeyError:
                            print("{}[rw0000] not found in {}".format(obs_json, to_open))
                            sigma_ij_mg = 0.
                    except ValueError:
                        print("{} does not seem to be in the list of active bins!".format(edge))
                        sigma_ij_mg = 0.
                    #print("sigma_ij_mg = {}".format(sigma_ij_mg))
                    sigma_yr = cross_sections[xs_name]
                    #print("sigma_yr = {}".format(sigma_yr))
                    try:
                        sigma_mg = rivet_dct["h_sigma[rw0000]"][0][0]
                        #print("sigma_mg = {}".format(sigma_mg))
                    except:
                        print("h_sigma was not found in file {}!".format(to_open))
                        sigma_mg = 0.0000001 
                        print("Converting sigma_ij_mg to 0 to avoid un-realistic cases")
                        sigma_ij_mg = 0.
                    to_add = sigma_ij_mg * sigma_yr / sigma_mg
                    #print("Add {} to denominator from xs {}".format(to_add, xs_name))
                    den += to_add
                    xs_pieces[edge][xs_name] = {
                            "sigma_ij_mg": sigma_ij_mg,
                            "sigma_yr": sigma_yr,
                            "sigma_mg": sigma_mg
                            }
                for xs_name in cross_sections_names:
                    sigma_ij_mg = xs_pieces[edge][xs_name]["sigma_ij_mg"]
                    sigma_yr = xs_pieces[edge][xs_name]["sigma_yr"]
                    sigma_mg = xs_pieces[edge][xs_name]["sigma_mg"]
                    factors[xs_name] = (sigma_ij_mg * sigma_yr / sigma_mg) / den
                    #if i_bin in (0, 1) and xs_name.startswith('ggH'): 
                    #    print(sigma_ij_mg)
                    #    print(sigma_yr)
                    #    print(sigma_mg)
                    #    print(den)
                    #    print(factors[xs_name])
                print("Factors for dc {}, bin {}, edge {}: {}".format(dc, i_bin, edge, factors))
                print("Sum of factors = {}".format(np.sum(list(factors.values()))))
                
                for xs_name in cross_sections_names:
                    conv_equations_file_name = "{}{}_{}_{}.json".format(eq_base_dir, xs_name, dc, obs)
                    if xs_name == "ggH_SMEFTatNLO":
                        conv_equations_file_name = "{}{}_{}_combined_{}.json".format(eq_base_dir, xs_name, dc, obs)
                    #print("Will open file {}".format(conv_equations_file_name))
                    with open(conv_equations_file_name, "r") as f:
                        try:
                            as_and_bs = json.load(f)[edge]
                        except KeyError:
                            as_and_bs = []
                    for key in as_and_bs:
                        if key.startswith(("A", "B")):
                            try:
                                output_dct[edge][key] += as_and_bs[key] * factors[xs_name]
                            except KeyError:
                                output_dct[edge][key] = as_and_bs[key] * factors[xs_name]

            # Write 
            out_file = "{}FullProduction{}_{}_{}.json".format(eq_base_dir, "WithSMEFTsim" if "ggH_SMEFTsim_topU3l" in cross_sections_names else "", dc, obs)
            print("Dumping output file {}".format(out_file))
            with open(out_file, "w") as f:
                json.dump(output_dct, f, indent=4)
