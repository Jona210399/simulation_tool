import numpy as np


def loguniform(low = 0, high = 1, size=None):
    return np.exp(np.random.uniform(np.log(low), np.log(high), size))

def sample_random_parameters(session_path, session_path_old, index, bandgap_ev):

    randomized_settings = {}
    
    ### hierarchy: From active layer to l1/l3, and then to the electrodes
    # KF: l1 and l3 seem to be switched compared to the online version (check that is just the naming)
    # l2: active layer
    l2_E_c = np.random.randn()*0.5 + 3.0
    # KF absolute values do not matter, so E_c can also be constant if the bandgap is randomized
    l2_E_v = l2_E_c + bandgap_ev
    # l1: l1_E_c above l2_E_c: don't block holes
    l1_E_c = np.random.randn()*0.1 + l2_E_c + 0.1 
    # KF: I would allow energy level offsets of up to 0.5 eV (for all cases). They cause S-shapes that are found in measurements  
    # l1: l1_E_v above l2_E_v: block electrons
    l1_E_v = np.random.randn()*0.1 + l2_E_v + 0.1
    # l3: l3_E_c below l2_E_c: block holes
    l3_E_c = np.random.randn()*0.1 + l2_E_c - 0.1
    # l3: l3_E_v below l2_E_v: don't block electrons
    l3_E_v = np.random.randn()*0.1 + l2_E_v - 0.1
    # electrodes: W_L equal or above l1_E_c: don't block holes --> hardcoded in simss
    # KF: electrodes can be aligned to the transport layer if energy level offsets between active layers and transport layers are allowed
    # if randomized, I would go for a smaller offset (np.random.randn()*0.5)
    W_L = np.random.randn()*0.5 + l1_E_c
    if W_L < l1_E_c:
        W_L = l1_E_c
    # electrodes: W_R equal or below l3_E_v: don't block electrons --> hardcoded in simss
    W_R = np.random.randn()*0.5 + l3_E_v
    if W_R > l3_E_v:
        W_R = l3_E_v
    
    ### General settings
    # Contacts
    randomized_settings["W_L"] = W_L # default 3
    randomized_settings["W_R"] = W_R # default 5
    randomized_settings["S_n_L"] = -1e-7 # default -1e-7
    randomized_settings["S_p_L"] = -1e-7 # default -1e-7
    randomized_settings["S_n_R"] = -1e-7 # default -1e-7
    randomized_settings["S_p_R"] = -1e-7 # default -1e-7
    randomized_settings["R_shunt"] = -1 # default -1
    # KF: should be randomized, something like loguniform(low = 1e-4, high = 1e-2)
    randomized_settings["R_series"] = 0 # default 0
    # KF: should be randomized, something like loguniform(low = 1e-3, high = 1e1)
    # KF: series and parallel resistance are not part of the drift diffusion calculation. if computation time is an issue,
    # they can also be added using Kirchhoff's law
    
    # Optics
    randomized_settings["G_frac"] = 1 # np.clip(np.random.randn()*0.3+1.0, a_min=0.0, a_max=1.0) # default 1
    # KF: I would not randomize
    randomized_settings["genProfile"] = "calc" # 
    # KF: if you set this to 'none', you can vary G_ehp to direcltly generate a large range of short-circuit currents (which are 
    # of course not longer related to the stack)
    randomized_settings["L_TCO"] = loguniform(low = 1e-8, high = 3e-7) # default: 1.1e-7, if 0 it is not used
    # KF: I suggest loguniform(low = 1e-8, high = 4e-7)
    randomized_settings["L_BE"] = loguniform(low = 1e-8, high = 3e-7) # default 2e-7, has to be >0
    # KF: I suggest loguniform(low = 1e-8, high = 1e-7). Silver becomes completely reflective at ~70 nm

    # Optics
    randomized_settings["Vstep"] = 0.025 # default 0.025
        
    # get the missing default parameters and add them here
    simss_config = get_simss_config(session_path, session_path_old, index)
    for key in simss_config.keys():
        if key not in randomized_settings.keys():
            randomized_settings[key] = simss_config[key]

    ### Layer settings
    # General
    randomized_settings["l1.L"] = loguniform(low = 3e-9, high = 3e-7) # default 3E-8
    # I have never seen a transport layer thicker than 50 nm, but for foundation model purposes this is ok
    layerfactor_l1 = (np.log(randomized_settings["l1.L"]) - np.log(3e-9))/(np.log(3e-7) - np.log(3e-9))
    # KF: I know I talked about the interdependece bewtween thickness and mobility, but I would skip the layer factor because it is a bit arbitrary
    # recombination rates should be considered as well, and the model might learn a certain bias
    randomized_settings["l1.eps_r"] = np.random.randn()*1.0 + 4.0 # default 5
    randomized_settings["l1.E_c"] = l1_E_c # default 3
    randomized_settings["l1.E_v"] = l1_E_v # default 5.5
    randomized_settings["l1.N_c"] = loguniform(low = 1e16, high = 1e27) # default 1E26
    # KF: not more than +- one order of magnitude. e.g. loguniform(low = 1e25, high = 1e27)
    randomized_settings["l1.N_D"] = loguniform(low = 1e16, high = 1e21) # default 0
    randomized_settings["l1.N_A"] = loguniform(low = 1e16, high = 1e21) # default 0
    # KF: I have never seen an influence of these parameters, but that's most likely we only use thin transport layers

    randomized_settings["l2.L"] = loguniform(low = 5e-9, high = 7.5e-7) # default 2.4E-7
    # KF: OPV 50nm - 200nm, perovskite 300nm - ~1µm
    layerfactor_l2 = (np.log(randomized_settings["l2.L"]) - np.log(5e-9))/(np.log(7.5e-7) - np.log(5e-9))
    randomized_settings["l2.eps_r"] = np.random.randn()*1.0 + 4.0 # default 4
    randomized_settings["l2.E_c"] = l2_E_c # default 3
    randomized_settings["l2.E_v"] = l2_E_v #  # default 5
    randomized_settings["l2.N_c"] = loguniform(low = 1e24, high = 1e27) # default 1E26
    randomized_settings["l2.N_D"] = loguniform(low = 1e16, high = 1e21) # default 0
    randomized_settings["l2.N_A"] = loguniform(low = 1e16, high = 1e21) # default 0

    randomized_settings["l3.L"] = loguniform(low = 3e-9, high = 3e-7) # default 3E-8
    layerfactor_l3 = (np.log(randomized_settings["l3.L"]) - np.log(3e-9))/(np.log(3e-7) - np.log(3e-9))
    randomized_settings["l3.eps_r"] = np.random.randn()*1.0 + 4.0 # default 3.5
    randomized_settings["l3.E_c"] = l3_E_c # default 2.5
    randomized_settings["l3.E_v"] = l3_E_v # default 5
    randomized_settings["l3.N_c"] = loguniform(low = 1e24, high = 1e27) # default 1E26
    randomized_settings["l3.N_D"] = loguniform(low = 1e16, high = 1e21) # default 0
    randomized_settings["l3.N_A"] = loguniform(low = 1e16, high = 1e21) # default 0
    
    # Mobilities
    randomized_settings["l1.mu_n"] = 10.0**(-np.random.randn()*2.0 - 7.0 + layerfactor_l1*2.0) # loguniform(low = 1e-9, high = 1e-2) # default 1e-6
    # KF: I would use loguniform, (low = 1e-9, high = 1e-0)
    randomized_settings["l1.mu_p"] = 10.0**(-np.random.randn()*2.0 - 7.0 + layerfactor_l1*2.0) # loguniform(low = 1e-9, high = 1e-2) # default 1e-6
    # KF: treat electron and hole mobility in the same way
    randomized_settings["l1.mobnDep"] = 0 # not randomized
    randomized_settings["l1.mobpDep"] = 0 # not randomized
    randomized_settings["l1.gamma_n"] = 0 # not randomized
    randomized_settings["l1.gamma_p"] = 0 # not randomized
    
    randomized_settings["l2.mu_n"] = 10.0**(-np.random.randn()*2 - 4.0 + layerfactor_l2*2.0) # loguniform(low = 1e-6, high = 1e0) # default 1e-2
    # KF: I would use loguniform, OPV: loguniform(low = 5e-9, high = 1e-6), pero and in general loguniform(low = 5e-9, high = 1e0)
    randomized_settings["l2.mu_p"] = 10.0**(-np.random.randn()*2 - 4.0 + layerfactor_l2*2.0) # loguniform(low = 1e-6, high = 1e0) # default 1e-2
    # KF: treat electron and hole mobility in the same way
    randomized_settings["l2.mobnDep"] = 0 # not randomized
    randomized_settings["l2.mobpDep"] = 0 # not randomized
    randomized_settings["l2.gamma_n"] = 0 # default 0.0001, not randomized
    randomized_settings["l2.gamma_p"] = 0 # not randomized

    randomized_settings["l3.mu_n"] = 10.0**(-np.random.randn()*2.0 - 7.0 + layerfactor_l3*2.0) # loguniform(low = 1e-9, high = 1e-2) # default 1e-8
    # KF: see above
    randomized_settings["l3.mu_p"] = 10.0**(-np.random.randn()*2.0 - 7.0 + layerfactor_l3*2.0) # loguniform(low = 1e-9, high = 1e-2) # default 1e-8
    # KF: see above
    randomized_settings["l3.mobnDep"] = 0 # not randomized
    randomized_settings["l3.mobpDep"] = 0 # not randomized
    randomized_settings["l3.gamma_n"] = 0 # not randomized
    randomized_settings["l3.gamma_p"] = 0 # not randomized
    
    # Interface-layer-to-right
    randomized_settings["l1.nu_int_n"] = 1E3 # default 1e3, m/s, interface transfer velocity of electrons, to layer to the right
    randomized_settings["l1.nu_int_p"] = 1E3 # default 1e3, m/s, interface transfer velocity of holes, to layer to the right
    #l1_trap_probabilitly = 0.1
    #if l1_int_trap_probabilitly > np.random.random():
    #    l1_int_N_t_int = loguniform(low = 1e-12, high = 1e-8)
    #else:
    #    l1_N_t_int = 0
    randomized_settings["l1.N_t_int"] = 0 # default 0, m^-2, trap density at interface with layer to the right
    # KF: needs to be larger than 0 (it does not make sense to randomize the capture coefficients if trap density is zero)
    # KF: interface traps should be randomized for l1 and l2 (at the interfaces between transport layers and absorber, not at the interfaces between transport
    # layers and electrodes)
    #max_E_c_12 = max(randomized_settings["l1.E_c"], randomized_settings["l2.E_c"])
    #min_E_v_12 = min(randomized_settings["l1.E_v"], randomized_settings["l2.E_v"])
    #l1_E_t_int = np.random.uniform(low = max_E_c_12 + 0.01, high = min_E_v_12 - 0.01)
    randomized_settings["l1.E_t_int"] = 4 # default 4, eV, energy level of traps at interface
    # KF: all trap energies should be randomized to cover the whole band gap (as it is done for the active layer)
    randomized_settings["l1.intTrapFile"] = "none" # not randomized, name of file with interface trap energy profile (or 'none'). If specified, overrides E_t_int
    randomized_settings["l1.intTrapType"] = -1 # not randomized, Trap type for the right interface: -1: acceptor, 0: neutral, 1: donor
    # KF: I have never changed this parameter, but I would randomize it to see what happens
    randomized_settings["l1.C_n_int"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13, m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)
    randomized_settings["l1.C_p_int"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13, m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band)
    # all trap-assisted recombination in layers and at interfaces is determined by the product of trap density and capture coefficient, so only either N_t_int or C_p/n_int needs to be randomized

    randomized_settings["l2.nu_int_n"] = 1E3 # default 1e3
    randomized_settings["l2.nu_int_p"] = 1E3 # default 1e3
    #l2_int_trap_probabilitly = 0.1
    #if l2_int_trap_probabilitly > np.random.random():
    #    l2_N_t_int = loguniform(low = 1e-12, high = 1e-8)
    #else:
    #    l2_N_t_int = 0
    randomized_settings["l2.N_t_int"] = 0 # default 0
    #max_E_c_23 = max(randomized_settings["l2.E_c"], randomized_settings["l3.E_c"])
    #min_E_v_23 = min(randomized_settings["l2.E_v"], randomized_settings["l3.E_v"])
    #l2_E_t_int = np.random.uniform(low = max_E_c_23 + 0.01, high = min_E_v_23 - 0.01)
    randomized_settings["l2.E_t_int"] = 4 # default 4
    randomized_settings["l2.intTrapFile"] = "none" # not randomized
    randomized_settings["l2.intTrapType"] = 1 # not randomized
    randomized_settings["l2.C_n_int"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13
    randomized_settings["l2.C_p_int"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13
    # KF: here same comments as for the interface traps at the right side of l1
    
    randomized_settings["l3.nu_int_n"] = 1E3 # default 1e3
    randomized_settings["l3.nu_int_p"] = 1E3 # default 1e3
    #l3_int_trap_probabilitly = 0.1
    #if l3_int_trap_probabilitly > np.random.random():
    #    l3_N_t_int = loguniform(low = 1e-12, high = 1e-8)
    #else:
    #    l3_N_t_int = 0
    randomized_settings["l3.N_t_int"] = 0 # default 0
    randomized_settings["l3.E_t_int"] = 4 # default 4
    randomized_settings["l3.intTrapFile"] = "none" # not randomized
    randomized_settings["l3.intTrapType"] = 1 # not randomized
    randomized_settings["l3.C_n_int"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13
    randomized_settings["l3.C_p_int"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13
    # if this is the interface to the electrode (as in the online version), no interface trapping needs to be assumed
    
    # Ions --> not randomized here
    randomized_settings["l1.N_anion"] = 0 # not randomized
    randomized_settings["l1.N_cation"] = 0 # not randomized
    randomized_settings["l1.mu_anion"] = 1E-14 # not randomized
    randomized_settings["l1.mu_cation"] = 1E-14 # not randomized
    randomized_settings["l1.ionsMayEnter"] = 0 # not randomized
    # KF: yes, no ions in the transport layer
    
    randomized_settings["l2.N_anion"] = 0 # not randomized
    randomized_settings["l2.N_cation"] = 0 # not randomized
    # KF: these values should be randomized for perovskites, but increasing the values by less than an order of magnitude leads to convergence issues. I am not sure 
    # how to handle this in a general way (we can ask the development team)
    randomized_settings["l2.mu_anion"] = 1E-14 # not randomized
    randomized_settings["l2.mu_cation"] = 1E-14 # not randomized
    randomized_settings["l2.ionsMayEnter"] = 1 # not randomized
    
    randomized_settings["l3.N_anion"] = 0 # not randomized
    randomized_settings["l3.N_cation"] = 0 # not randomized
    randomized_settings["l3.mu_anion"] = 1E-14 # not randomized
    randomized_settings["l3.mu_cation"] = 1E-14 # not randomized
    randomized_settings["l3.ionsMayEnter"] = 0 # not randomized
    # KF: yes, no ions in the transport layer
    
    # Generation and recombination
    randomized_settings["l1.G_ehp"] = 0 # not randomized
    randomized_settings["l1.layerGen"] = 0 # not randomized
    randomized_settings["l1.nkLayer"] = f"{session_path_old}/../Data/nk_PEDOT.txt" # not randomized here
    randomized_settings["l1.fieldDepG"] = 0 # not randomized
    randomized_settings["l1.P0"] = 0 # not randomized
    randomized_settings["l1.a"] = 1E-9 # not randomized
    randomized_settings["l1.thermLengDist"] = 2 # not randomized
    randomized_settings["l1.k_f"] = 1E6 # not randomized
    randomized_settings["l1.k_direct"] = 1E-18 # not randomized
    randomized_settings["l1.preLangevin"] = 1 # not randomized
    randomized_settings["l1.useLangevin"] = 1 # not randomized
    # KF: set useLangevin = 0 for ALL the layers. For useLangvin = 1, band-to-band recombination is calculated from the mobilities (and not from k_direct)
    # by a formula that does not correspond to observation (that's why the empirical pre-Langevin factor is introduced)

    randomized_settings["l2.G_ehp"] = 7.25E27 # default 7.25E27, m^-3 s^-1, generation rate of electron-hole pairs in this layer
    # KF: could be randomized if wen want to disregard the layer stack
    randomized_settings["l2.layerGen"] = 1 # not randomized
    randomized_settings["l2.nkLayer"] = f"{session_path}/nk_{index}.txt" # randomized before #  this also determines l2.bandgap
    randomized_settings["l2.fieldDepG"] = 0 # indirectly randomized - KF: only if the parameter is set to 1
    # the following parameters until k_f all refer to the Onsager-Braun model, which is only relevant for OPV and only if fieldDepG = 1
    # I have never used this model, and according to Jan Anton it is mainly here for historic reasons (I have never used it)
    randomized_settings["l2.P0"] = np.clip(np.random.randn()*0.2, a_min=0.0, a_max=1.0) # default 0, 0<=P0<1, fraction of quenched excitons that direcltly yield free carriers
    if randomized_settings["l2.P0"] > 0.0:
        randomized_settings["l2.fieldDepG"] = 1
    randomized_settings["l2.a"] = loguniform(low = 0.5e-9, high = 1e-8) # m, charge separation distance, Braun model used
    randomized_settings["l2.thermLengDist"] = 2 # distribution of a, 1 for delta function, 2 for Gaussian, 3 for exponential and 4 for r^2 exponential 5 for r^4 Gaussian
    randomized_settings["l2.k_f"] = loguniform(low = 1e4, high = 1e8) # default 1e6, 1/s, decay rate
    randomized_settings["l2.k_direct"] = loguniform(low = 1e-19, high = 1e-17) # default 1e-18, m3/s, direct (band-to-band, bimolecular) recombination rate
    # KF: I suggest logunifrom(low = 1e-18 high = 1e-10)
    randomized_settings["l2.preLangevin"] = 1 # not randomized, Langevin recombination prefactor
    randomized_settings["l2.useLangevin"] = 1 # not randomized, (1) use Langevin to calc. recombination or not (<>1, kdirect is used)

    randomized_settings["l3.G_ehp"] = 0 # not randomized
    randomized_settings["l3.layerGen"] = 0 # not randomized
    randomized_settings["l3.nkLayer"] = f"{session_path_old}/../Data/nk_Ca.txt" # not randomized here
    randomized_settings["l3.fieldDepG"] = 0 # not randomized
    randomized_settings["l3.P0"] = 0 # not randomized
    randomized_settings["l3.a"] = 1E-9 # not randomized
    randomized_settings["l3.thermLengDist"] = 2 # not randomized
    randomized_settings["l3.k_f"] = 1E6 # not randomized
    randomized_settings["l3.k_direct"] = 1E-18 # not randomized
    randomized_settings["l3.preLangevin"] = 1 # not randomized
    randomized_settings["l3.useLangevin"] = 1 # not randomized
    
    # Bulk trapping
    randomized_settings["l1.N_t_bulk"] = loguniform(low = 1e17, high = 1e23) # default 0, m^-3, trap density (in bulk)
    # KF: only the product of trap density and capture coefficient is relenvant, so in principle only one parameter needs to be randomized 
    randomized_settings["l1.C_n_bulk"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13, m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)
    randomized_settings["l1.C_p_bulk"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13, m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band)
    randomized_settings["l1.E_t_bulk"] = np.random.uniform(low = randomized_settings["l1.E_c"]+0.01, high = randomized_settings["l1.E_v"]-0.01) # default 4, eV, energy level of all traps
    randomized_settings["l1.bulkTrapFile"] = "none" # not randomized here, name of file with bulk trap energy profile (or 'none'). If specified, overrides E_t_bulk
    randomized_settings["l1.bulkTrapType"] = -1 # not randomized here, Trap type of bulk traps: -1: acceptor, 0: neutral, 1: donor
    # KF: does not need to be randomized if 
    
    randomized_settings["l2.N_t_bulk"] = loguniform(low = 1e17, high = 1e23) # default 4.42E20
    randomized_settings["l2.C_n_bulk"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13
    randomized_settings["l2.C_p_bulk"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13
    randomized_settings["l2.E_t_bulk"] = np.random.uniform(low = randomized_settings["l2.E_c"]+0.01, high = randomized_settings["l2.E_v"]-0.01) # default 4 
    randomized_settings["l2.bulkTrapFile"] = "none" # not randomized here
    randomized_settings["l2.bulkTrapType"] = -1 # not randomized here
    
    randomized_settings["l3.N_t_bulk"] = loguniform(low = 1e17, high = 1e23) # default 0
    randomized_settings["l3.C_n_bulk"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13
    randomized_settings["l3.C_p_bulk"] = loguniform(low = 1e-15, high = 1e-11) # default 1e-13
    randomized_settings["l3.E_t_bulk"] = np.random.uniform(low = randomized_settings["l3.E_c"]+0.01, high = randomized_settings["l3.E_v"]-0.01) # default 4 
    randomized_settings["l3.bulkTrapFile"] = "none" # not randomized here
    randomized_settings["l3.bulkTrapType"] = -1 # not randomized here
    
    return(randomized_settings)
    
    

def get_simss_config(session_path, session_path_old, index):
    simss_config = {
        # General
        "T": 295,
        
        # Layers
        "l1": f"{session_path}/L1_parameters_{index}.txt",
        "l2": f"{session_path}/L2_parameters_{index}.txt",
        "l3": f"{session_path}/L3_parameters_{index}.txt",
        
        # Contacts
        "W_L": 3.0,
        "W_R": 5.0,
        "S_n_L": -1e-7,
        "S_p_L": -1e-7,
        "S_n_R": -1e-7,
        "S_p_R": -1e-7,
        "R_shunt": -1,
        "R_series": 0, #np.random.random()*0.1,
        
        # Optics
        "G_frac": 1.0,
        "genProfile": "calc",
        "L_TCO": 1.1e-7,
        "L_BE": 2e-7,
        "nkSubstrate": f"{session_path_old}/../Data/nk_SiO2.txt",
        "nkTCO": f"{session_path_old}/../Data/nk_ITO.txt",
        "nkBE": f"{session_path_old}/../Data/nk_Al.txt",
        "spectrum": f"{session_path_old}/../Data/AM15G.txt",
        "lambda_min": 3.5e-7,
        "lambda_max": 8e-7,
        
        # Numerical Parameters
        "NP": 400,
        "tolPois": 0.00001,
        "maxDelV": 10,
        "maxItPois": 1000,
        "maxItSS": 1000,
        "currDiffInt": 2,
        "tolDens": 1e-8,
        "couplePC": 4,
        "minAcc": 0.5,
        "maxAcc": 1.0,
        "ignoreNegDens": 1,
        "failureMode": 0,
        "grad": 1,
        
        # Voltage range
        "Vdist": 1,
        "preCond": 0,
        "Vpre": 0,
        "fixIons": 0,
        "Vscan": 1,
        "Vmin": -0.5,
        "Vmax": 1.7,
        "Vstep": 0.025, # default 0.025
        "Vacc": 0,
        "NJV": 100,
        "untilVoc": 0,
        
        # User interface
        "timeout": -1,
        "pauseAtEnd": 0,
        "autoTidy": 1,
        "useExpData": 0,
        "expJV": "expJV.csv",
        "fitMode": "lin",
        "fitThreshold": 0.8,
        "JVFile": "JV.dat",
        "varFile": "Var.dat",
        "limitDigits": 1,
        "outputRatio": 1,
        "scParsFile": "scPars.dat",
        "logFile": "log.txt"
    }
    return(simss_config)

def get_layer1_config(session_path, session_path_old, index):
    layer1_config = {
        # General
        "l1.L": 3E-8,
        "l1.eps_r": 5,
        "l1.E_c": 3.0,
        "l1.E_v": 5.5,
        "l1.N_c": 1E26,
        "l1.N_D": 0,
        "l1.N_A": 0,
        
        # Mobilities
        "l1.mu_n": 1e-6,
        "l1.mu_p": 1e-6,
        "l1.mobnDep": 0,
        "l1.mobpDep": 0,
        "l1.gamma_n": 0,
        "l1.gamma_p": 0,
        
        # Interface-layer-to-right
        "l1.nu_int_n": 1E3,
        "l1.nu_int_p": 1E3,
        "l1.N_t_int": 0,
        "l1.E_t_int": 4,
        "l1.intTrapFile": "none",
        "l1.intTrapType": -1,
        "l1.C_n_int": 1E-13,
        "l1.C_p_int": 1E-13,
        
        # Ions
        "l1.N_anion": 0,
        "l1.N_cation": 0,
        "l1.mu_anion": 1E-14,
        "l1.mu_cation": 1E-14,
        "l1.ionsMayEnter": 0,
        
        # Generation and recombination
        "l1.G_ehp": 0,
        "l1.layerGen": 0,
        "l1.nkLayer": f"{session_path_old}/../Data/nk_PEDOT.txt",
        "l1.fieldDepG": 0,
        "l1.P0": 0,
        "l1.a": 1E-9,
        "l1.thermLengDist": 2,
        "l1.k_f": 1E6,
        "l1.k_direct": 1E-18,
        "l1.preLangevin": 1,
        "l1.useLangevin": 1,
        
        # Bulk trapping
        "l1.N_t_bulk": 0,
        "l1.C_n_bulk": 1E-13,
        "l1.C_p_bulk": 1E-13,
        "l1.E_t_bulk": 4,
        "l1.bulkTrapFile": "none",
        "l1.bulkTrapType": -1
    }
    return(layer1_config)
    
def get_layer2_config(session_path, session_path_old, index, bandgap_ev):
    print(bandgap_ev)
    layer2_config = {
        # General
        "l2.L": 5.4E-8,
        "l2.eps_r": 4,
        "l2.E_c": 3.0,
        "l2.E_v": 3 + bandgap_ev,
        "l2.N_c": 1E25,
        "l2.N_D": 0,
        "l2.N_A": 0,
        
        # Mobilities
        "l2.mu_n": 0.01,
        "l2.mu_p": 0.01,
        "l2.mobnDep": 0,
        "l2.mobpDep": 0,
        "l2.gamma_n": 0.0001,
        "l2.gamma_p": 0,
        
        # Interface-layer-to-right
        "l2.nu_int_n": 1E3,
        "l2.nu_int_p": 1E3,
        "l2.N_t_int": 0,
        "l2.E_t_int": 4,
        "l2.intTrapFile": "none",
        "l2.intTrapType": 1,
        "l2.C_n_int": 1E-13,
        "l2.C_p_int": 1E-13,
        
        # Ions
        "l2.N_anion": 0,
        "l2.N_cation": 0,
        "l2.mu_anion": 1E-14,
        "l2.mu_cation": 1E-14,
        "l2.ionsMayEnter": 1,
        
        # Generation and recombination
        "l2.G_ehp": 7.25E27,
        "l2.layerGen": 1,
        "l2.nkLayer": f"{session_path}/nk_{index}.txt",
        "l2.fieldDepG": 0,
        "l2.P0": 0,
        "l2.a": 1E-9,
        "l2.thermLengDist": 2,
        "l2.k_f": 1E6,
        "l2.k_direct": 1E-18,
        "l2.preLangevin": 1,
        "l2.useLangevin": 1,
        
        # Bulk trapping
        "l2.N_t_bulk": 4.42E21,
        "l2.C_n_bulk": 1E-13,
        "l2.C_p_bulk": 1E-13,
        "l2.E_t_bulk": 4,
        "l2.bulkTrapFile": "none",
        "l2.bulkTrapType": -1
    }
    return(layer2_config)

def get_layer3_config(session_path, session_path_old, index):
    layer3_config = {
        # General
        "l3.L": 3E-8,
        "l3.eps_r": 3.5,
        "l3.E_c": 2.5,
        "l3.E_v": 5,
        "l3.N_c": 1E26,
        "l3.N_D": 0,
        "l3.N_A": 0,
        
        # Mobilities
        "l3.mu_n": 1e-8,
        "l3.mu_p": 1e-8,
        "l3.mobnDep": 0,
        "l3.mobpDep": 0,
        "l3.gamma_n": 0,
        "l3.gamma_p": 0,
        
        # Interface-layer-to-right
        "l3.nu_int_n": 1E3,
        "l3.nu_int_p": 1E3,
        "l3.N_t_int": 0,
        "l3.E_t_int": 4,
        "l3.intTrapFile": "none",
        "l3.intTrapType": 1,
        "l3.C_n_int": 1E-13,
        "l3.C_p_int": 1E-13,
        
        # Ions
        "l3.N_anion": 0,
        "l3.N_cation": 0,
        "l3.mu_anion": 1E-14,
        "l3.mu_cation": 1E-14,
        "l3.ionsMayEnter": 0,
        
        # Generation and recombination
        "l3.G_ehp": 0,
        "l3.layerGen": 0,
        "l3.nkLayer": f"{session_path_old}/../Data/nk_Ca.txt",
        "l3.fieldDepG": 0,
        "l3.P0": 0,
        "l3.a": 1E-9,
        "l3.thermLengDist": 2,
        "l3.k_f": 1E6,
        "l3.k_direct": 1E-18,
        "l3.preLangevin": 1,
        "l3.useLangevin": 1,
        
        # Bulk trapping
        "l3.N_t_bulk": 0,
        "l3.C_n_bulk": 1E-13,
        "l3.C_p_bulk": 1E-13,
        "l3.E_t_bulk": 4,
        "l3.bulkTrapFile": "none",
        "l3.bulkTrapType": -1
    }
    return(layer3_config)
    
    
    
def get_simss_config_text(simss_config):
    simss_config_text = f"""** SimSS Simulation Setup:
** Don't change the order of the parameters, comments can be added anywhere,
** but only after an '*'. Use '**' if you want your comment to be left-justified.
** version: 5.18

**General***************************************************************************
T = {simss_config["T"]}                           * K, absolute temperature

**Layers****************************************************************************
l1 = {simss_config["l1"]}            * parameter file for layer 1, mandatory
l2 = {simss_config["l2"]}            * parameter file for layer 2
l3 = {simss_config["l3"]}            * parameter file for layer 3

**Contacts**************************************************************************
W_L = {simss_config["W_L"]}                           * eV, work function left electrode (= cathode), or 'sfb'
W_R = {simss_config["W_R"]}                           * eV, work function right electrode (= cathode), or 'sfb'
S_n_L = {simss_config["S_n_L"]}                     * m/s, surface recombination of electrons at the left electrode
S_p_L = {simss_config["S_p_L"]}                     * m/s, surface recombination of holes at the left electrode
S_n_R = {simss_config["S_n_R"]}                     * m/s, surface recombination of electrons at the right electrode
S_p_R = {simss_config["S_p_R"]}                     * m/s, surface recombination of holes at the right electrode
                                  * nb: use negative values if Sn/pR/L should be infinite
R_shunt = {simss_config["R_shunt"]}                      * Ohms m2, shunt resistance. Use negative value for infinite R_shunt
R_series = {simss_config["R_series"]}                      * Ohms m2, series resistance.

**Optics****************************************************************************
G_frac = {simss_config["G_frac"]}                        * fraction of Gmax used in solar cell
genProfile = {simss_config["genProfile"]}                 * name of file generation profile (or 'none' or 'calc')
L_TCO = {simss_config["L_TCO"]}                    * m, thickness of the TCO. Set to 0 if layer is not used
L_BE = {simss_config["L_BE"]}                       * m, thickness of back electrode, must be >0
nkSubstrate = {simss_config["nkSubstrate"]} * name of file with n,k values of substrate
nkTCO = {simss_config["nkTCO"]}        * name of file with n,k values of TCO
nkBE = {simss_config["nkBE"]}          * name of file with n,k values of back electrode
spectrum = {simss_config["spectrum"]}      * name of file that contains the spectrum
lambda_min = {simss_config["lambda_min"]}               * m, lower bound wavelength
lambda_max = {simss_config["lambda_max"]}                 * m, upper bound wavelength

**Numerical Parameters**************************************************************
NP = {simss_config["NP"]}                          * integer, number of grid points, must be at least 5 per layer.
tolPois = {simss_config["tolPois"]}                 * V, abs. tolerance of iterative Poisson solver
maxDelV = {simss_config["maxDelV"]}                      * maximum change (in Vt) of the potential per loop
maxItPois = {simss_config["maxItPois"]}                  * max. number it. Poisson loop
maxItSS = {simss_config["maxItSS"]}                    * max. number it. main loop
currDiffInt = {simss_config["currDiffInt"]}                   * Calc. current from differential (1) or integral (2) expression
tolDens = {simss_config["tolDens"]}                    * relative tolerance of density solver
couplePC = {simss_config["couplePC"]}                      * >= 0, coupling between Poisson equation and continuity equations
minAcc = {simss_config["minAcc"]}                      * >0, min. acceleration parameter
maxAcc = {simss_config["maxAcc"]}                      * <2, max. acceleration parameter
ignoreNegDens = {simss_config["ignoreNegDens"]}                 * whethe 1) or not(<>1) to ignore negative densities
failureMode = {simss_config["failureMode"]}                   * how treat failed (t,V,G) points: 0: stop, 1: ignore, 2: skip
grad = {simss_config["grad"]}                          * determines shape of exp. grid, increase grad for smaller h[1]

**Voltage range of simulation*******************************************************
Vdist = {simss_config["Vdist"]}                         * 1 for uniform (specified by Vstep), 2 for logarithmic (specified by Vacc and NJV)
preCond = {simss_config["preCond"]}                       * pre-conditioning, yes(1)/no(0)
Vpre = {simss_config["Vpre"]}                          * V, pre-conditioned voltage
fixIons = {simss_config["fixIons"]}                       * fix ions at first applied voltage? yes(1) or no (0).
Vscan = {simss_config["Vscan"]}                         * integer, 1 for forward sweep direction, -1 for reverse sweep
Vmin = {simss_config["Vmin"]}                          * V
Vmax = {simss_config["Vmax"]}                        * V
Vstep = {simss_config["Vstep"]}                     * V
Vacc = {simss_config["Vacc"]}                         * V, point of accumulation of row of V's, note: Vacc should be
                                  * slightly larger than Vmax or slightly lower than Vmin
NJV = {simss_config["NJV"]}                         * number of JV points in logarithmic distribution
untilVoc = {simss_config["untilVoc"]}                      * if 1 then SimSS will stop at Voc

**User interface*******************************************************************
timeout = {simss_config["timeout"]}                      * s, max run time, use negative value for unlimited run time.
pauseAtEnd = {simss_config["pauseAtEnd"]}                    * pause at the end of the simulation yes(1) or no (0)
autoTidy = {simss_config["autoTidy"]}                      * if 1, then the program will always tidy up this file
useExpData = {simss_config["useExpData"]}                    * if 1, SimSS will try to read JV_Exp and use it
expJV = {simss_config["expJV"]}                 * name of file with experimental JV characteristics
fitMode = {simss_config["fitMode"]}                     * lin or log: use J or log(J) in calc. of fit error
fitThreshold = {simss_config["fitThreshold"]}                * threshold of fraction converged points in calc. fit error
JVFile = {simss_config["JVFile"]}                   * name of the file with simulated JV characteristics
varFile = {simss_config["varFile"]}                 * name of the file with (x,V,n,p,Jn,etc) or none for no file.
limitDigits = {simss_config["limitDigits"]}                   * if 1, then number of digits in output is limited
outputRatio = {simss_config["outputRatio"]}                   * Output to varFile every outputRatio voltages
scParsFile = {simss_config["scParsFile"]}           * name of file with solar cell parameters
logFile = {simss_config["logFile"]}                 * name of log file
"""
    return(simss_config_text)


def get_layer_config_text(layer_config):
    layer_config_text = f"""** SIMsalabim Layer parameters:
** Don't change the order of the parameters, comments can be added anywhere,
** but only after an '*'. Use '**' if you want your comment to be left-justified.
** version: 5.18

**General**************************************************************************
L = {layer_config["L"]}                       * m, device length/thickness
eps_r = {layer_config["eps_r"]}                      * relative dielectric constant
E_c = {layer_config["E_c"]}                        * eV, conduction band edge
E_v = {layer_config["E_v"]}                      * eV, valence band edge
N_c = {layer_config["N_c"]}                     * m^-3, DOS of conduction and valence bands
N_D = {layer_config["N_D"]}                        * m^-3, ionised n-doping
N_A = {layer_config["N_A"]}                        * m^-3, ionised p-doping

**Mobilities************************************************************************
mu_n = {layer_config["mu_n"]}                    * m^2/Vs, zero field mobility
mu_p = {layer_config["mu_p"]}                    * m^2/Vs, zero field mobility
mobnDep = {layer_config["mobnDep"]}                    * 0 : const. mob, 1 : field-dependent
mobpDep = {layer_config["mobpDep"]}                    * 0 : const. mob, 1 : field-dependent
gamma_n = {layer_config["gamma_n"]}                    * (m/V)^0.5, field dependence of mob, Poole-Frenkel form
gamma_p = {layer_config["gamma_p"]}                    * (m/V)^0.5, field dependence of mob, Poole-Frenkel form

**Interface-layer-to-right**********************************************************
nu_int_n = {layer_config["nu_int_n"]}                 * m/s, interface transfer velocity of electrons, to layer to the right
nu_int_p = {layer_config["nu_int_p"]}                 * m/s, interface transfer velocity of holes, to layer to the right
N_t_int = {layer_config["N_t_int"]}                    * m^-2, trap density at interface with layer to the right
E_t_int = {layer_config["E_t_int"]}                    * eV, energy level of traps at interface
intTrapFile = {layer_config["intTrapFile"]}             * name of file with interface trap energy profile (or 'none'). If specified, overrides E_t_int
intTrapType = {layer_config["intTrapType"]}               * Trap type for the right interface: -1: acceptor, 0: neutral, 1: donor
C_n_int = {layer_config["C_n_int"]}                * m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)
C_p_int = {layer_config["C_p_int"]}                * m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band)

**Ions******************************************************************************
N_anion = {layer_config["N_anion"]}                    * m^-3, concentration of negative ions
N_cation = {layer_config["N_cation"]}                   * m^-3, concentration of positive ions
mu_anion = {layer_config["mu_anion"]}               * m^2/Vs, mobility of negative ions (take 0 if they don't move)
mu_cation = {layer_config["mu_cation"]}              * m^2/Vs, mobility of positive ions (take 0 if they don't move)
ionsMayEnter = {layer_config["ionsMayEnter"]}               * may ions enter from other layers? yes(1) or no(<>1)

**Generation and recombination******************************************************
G_ehp = {layer_config["G_ehp"]}                      * m^-3 s^-1, generation rate of electron-hole pairs in this layer
layerGen = {layer_config["layerGen"]}                   * does this layer generate electron/hole pairs? yes(1) or no (0)
nkLayer = {layer_config["nkLayer"]} * name of file with n,k values of this layer
fieldDepG = {layer_config["fieldDepG"]}                  * field dependent generation yes (1) or no (0)
P0 = {layer_config["P0"]}                         * 0<=P0<1, fraction of quenched excitons that direcltly yield free carriers
a = {layer_config["a"]}                       * m, charge separation distance, Braun model used
thermLengDist = {layer_config["thermLengDist"]}              * distribution of a, 1 for delta function, 2 for Gaussian
                               * 3 for exponential and 4 for r^2 exponential 5 for r^4 Gaussian
k_f = {layer_config["k_f"]}                      * 1/s, decay rate
k_direct = {layer_config["k_direct"]}               * m3/s, direct (band-to-band, bimolecular) recombination rate
preLangevin = {layer_config["preLangevin"]}                * Langevin recombination prefactor
useLangevin = {layer_config["useLangevin"]}                * (1) use Langevin to calc. recombination or not (<>1, kdirect is used)

**Bulk trapping**************************************************************************
N_t_bulk = {layer_config["N_t_bulk"]}                   * m^-3, trap density (in bulk)
C_n_bulk = {layer_config["C_n_bulk"]}               * m^3/s, capture coefficient for electrons (put to 0 to exclude capture from and emission to the conduction band)
C_p_bulk = {layer_config["C_p_bulk"]}               * m^3/s, capture coefficient for holes (put to 0 to exclude capture from and emission to the valence band)
E_t_bulk = {layer_config["E_t_bulk"]}                   * eV, energy level of all traps
bulkTrapFile = {layer_config["bulkTrapFile"]}            * name o file with bulk trap energy profile (or 'none'). If specified, overrides E_t_bulk
bulkTrapType = {layer_config["bulkTrapType"]}              * Trap type of bulk traps: -1: acceptor, 0: neutral, 1: donor
"""
    return(layer_config_text)
