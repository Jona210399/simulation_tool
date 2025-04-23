import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------
# Read results.txt and record successful simulations
# -------------------------------
def read_successful_simulations(results_file):
    """
    Reads the results.txt file and returns a list of simulation IDs
    that completed successfully (i.e. lines that contain 'DONE').
    """
    successful_ids = []
    with open(results_file, 'r') as f:
        for line in f:
            if "DONE" in line:
                # Assume the first token is the simulation ID
                sim_id = line.split()[0]
                successful_ids.append(sim_id)
    return successful_ids

def read_config_file(fn, prefix = ""):
    param_dict = {}
    for line in open(fn, "r"):
        line = line.replace(" ", "").replace("\n", "")
        if len(line)>0:
            if line[0] != "*":
                line = line.split("*")[0].split("=")
                param_dict[f"{prefix}{line[0]}"] = line[1]
    return(param_dict)

def read_settings(sim_dir, sim_id):
    param_dict_main = read_config_file(f"{sim_dir}/simulation_setup_custom_{sim_id}.txt")
    param_dict_l1 = read_config_file(f"{sim_dir}/L1_parameters_{sim_id}.txt", prefix="l1.")
    param_dict_l2 = read_config_file(f"{sim_dir}/L2_parameters_{sim_id}.txt", prefix="l2.")
    param_dict_l3 = read_config_file(f"{sim_dir}/L3_parameters_{sim_id}.txt", prefix="l3.")
    param_dict_final = param_dict_main | param_dict_l1 | param_dict_l2 | param_dict_l3
    nk_data = []
    for lineidx, line in enumerate(open(f"{sim_dir}/nk_{sim_id}.txt", "r")):
        if lineidx > 0:
            nk_data.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    nk_data = np.array(nk_data)
    return(param_dict_final, nk_data)
            

# -------------------------------
# Read data from a given simulation run folder
# -------------------------------
def read_simulation_data(sim_dir, sim_id):
    """
    Reads jV, device characteristics, EQE, and UV-vis data from the simulation directory.
    
    Assumes that within sim_dir the following files exist:
      - jV_0.txt
      - device_characteristics_0.txt
      - EQE_0.dat
      - AbsorptionSpectrum_0.txt
      - reflection_transmission_spectrum_0.txt
    
    Returns a dictionary with the following keys:
      'jV': { 'V': np.array, 'J': np.array }
      'dc': { 'jsc': float, 'voc': float, 'ff': float }
      'EQE': { 'lambda': np.array, 'EQE': np.array }
      'uvvis': { 'wavelengths': np.array, 'A_prime': np.array, 'R': np.array }
    """
    
    param_dict, nk_data = read_settings(sim_dir, sim_id)
    
    data_dict = {}

    # jV data
    try:
        jv_file = os.path.join(sim_dir, f"jV_{sim_id}.txt")
        df_jv = pd.read_csv(jv_file, sep='\s+')
        data_dict['jV'] = {'V': df_jv['Vext'].to_numpy(), 'J': df_jv['Jext'].to_numpy()}
    except Exception as e:
        print(f"Error reading {jv_file}: {e}")
        return None

    # Device characteristics (extract second line; expected order: jsc, ..., voc, ..., ff, ... as in your code)
    try:
        dc_file = os.path.join(sim_dir, f"device_characteristics_{sim_id}.txt")
        with open(dc_file, "r") as f:
            lines = f.readlines()
        # Assuming the second line (index 1) contains the relevant numbers
        values = [float(x) for x in lines[1].split()]
        # According to your code: jsc = values[0], voc = values[2], ff = values[8]
        data_dict['dc'] = {'jsc': values[0], 'voc': values[2], 'ff': values[8]}
    except Exception as e:
        print(f"Error reading {dc_file}: {e}")
        return None

    # EQE data
    try:
        eqe_file = os.path.join(sim_dir, f"EQE_{sim_id}.dat")
        df_eqe = pd.read_csv(eqe_file, sep='\s+')
        # Expecting columns named 'lambda' and 'EQE'
        data_dict['EQE'] = {'lambda': df_eqe['lambda'].to_numpy(), 'EQE': df_eqe['EQE'].to_numpy()}
    except Exception as e:
        print(f"Error reading {eqe_file}: {e}")
        return None

    # UV-vis absorption and reflection
    try:
        uvvis_file = os.path.join(sim_dir, f"AbsorptionSpectrum_{sim_id}.txt")
        uvvis_list = []
        with open(uvvis_file, "r") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2:
                    uvvis_list.append([float(parts[0]), float(parts[1])])
        uvvis_data = np.array(uvvis_list)
        # Convert wavelength from meters to nm
        wavelengths = uvvis_data[:, 0] * 1e9
        A = uvvis_data[:, 1]
        # Compute effective absorption: A' = 1 - exp(-A)
        A_prime = 1 - np.exp(-A)
    except Exception as e:
        print(f"Error reading {uvvis_file}: {e}")
        return None

    try:
        rt_file = os.path.join(sim_dir, f"reflection_transmission_spectrum_{sim_id}.txt")
        rt_list = []
        with open(rt_file, "r") as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 3:
                    rt_list.append([float(parts[0]), float(parts[1]), float(parts[2])])
        rt_data = np.array(rt_list)
        # Reflection is taken as the second column.
        R = rt_data[:, 1]
    except Exception as e:
        print(f"Error reading {rt_file}: {e}")
        return None

    data_dict['uvvis'] = {'wavelengths': wavelengths, 'A_prime': A_prime, 'R': R}

    return data_dict, param_dict, nk_data

def isnumber(x):
    try:
        x = float(x)
    except:
        return(False)
    return(True)

# -------------------------------
# (c) Process all data: calculate PCE and record device parameters and curves
# -------------------------------
def process_all_simulations(results_file, simulation_runs_dir):
    successful_ids = read_successful_simulations(results_file)
    print(f"Found {len(successful_ids)} successful simulations.")

    # Lists to store device parameters for histogram/distribution plots
    pce_list = []
    jsc_list = []
    voc_list = []
    ff_list = []

    # Lists to store curves (each element is a dict with curve data)
    jV_curves = []
    EQE_curves = []
    uvvis_curves = []  # contains wavelengths, absorption (A_prime), and reflection
    param_dicts_all = []
    param_lists_all = []
    nk_data_all = []
    
    # Loop over each successful simulation folder
    for idx, sim_id in enumerate(successful_ids):
        sim_dir = os.path.join(simulation_runs_dir, f"simulation_run_{sim_id}")
        if not os.path.isdir(sim_dir):
            print(f"Directory {sim_dir} not found. Skipping simulation {sim_id}.")
            continue

        sim_data, param_dict, nk_data = read_simulation_data(sim_dir, sim_id)
        if sim_data is None:
            print(f"Data for simulation {sim_id} could not be read.")
            continue

        # Device characteristics
        dc = sim_data['dc']
        jsc = dc['jsc']      # in A/m²
        voc = dc['voc']      # in V
        ff = dc['ff']        # as fraction (e.g., 0.7)
        jsc_list.append(jsc)
        voc_list.append(voc)
        ff_list.append(ff)
        param_dicts_all.append(param_dict)
        param_list = []
        for key in param_dict.keys():
            if isnumber(param_dict[key]):
                param_list.append(float(param_dict[key]))
        param_lists_all.append(param_list)
        nk_data_all.append(nk_data)

        # Calculate PCE in percentage: (Jsc [A/m²] * Voc [V] * FF) gives power density in W/m².
        # Divide by 1000 (1000 W/m²) and multiply by 100 to get percentage.
        pce = -1.0 * (jsc * voc * ff) * 100.0 / 1000.0  # equivalently, pce = (jsc * voc * ff) / 10.
        pce_list.append(pce)

        # jV curve data
        jV_curves.append(sim_data['jV'])

        # EQE curve data
        EQE_curves.append(sim_data['EQE'])

        # UV-vis curve data
        uvvis_curves.append(sim_data['uvvis'])
        if (idx+1)%1000 == 0 or (idx+1) == len(successful_ids):
            print(f"{idx+1} out of {len(successful_ids)} results read")
            
    # Convert lists to numpy arrays where applicable
    device_characteristics = {
        'PCE': np.array(pce_list),
        'Jsc': np.array(jsc_list),
        'Voc': np.array(voc_list),
        'FF': np.array(ff_list)
    }

    return device_characteristics, jV_curves, EQE_curves, uvvis_curves, param_dicts_all, param_lists_all, nk_data_all, successful_ids

# -------------------------------
# (d) Plot distributions and all curves
# -------------------------------
def plot_distributions(device_characteristics, filename = "device_characteristics.png"):
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)
    plt.hist(device_characteristics['PCE'], bins=20, color='C0', edgecolor='black')
    plt.xlabel("PCE [%]")
    plt.ylabel("Counts")
    plt.title("PCE distribution")

    plt.subplot(2, 2, 2)
    plt.hist(device_characteristics['Jsc'], bins=20, color='C1', edgecolor='black')
    plt.xlabel("Jsc [A/m²]")
    plt.ylabel("Counts")
    plt.title("Jsc distribution")

    plt.subplot(2, 2, 3)
    plt.hist(device_characteristics['Voc'], bins=20, color='C2', edgecolor='black')
    plt.xlabel("Voc [V]")
    plt.ylabel("Counts")
    plt.title("Voc distribution")

    plt.subplot(2, 2, 4)
    plt.hist(device_characteristics['FF']*100, bins=20, color='C3', edgecolor='black')
    plt.xlabel("Fill Factor [%]")
    plt.ylabel("Counts")
    plt.title("FF distribution")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print("plotted device characteristics distributions")

def plot_jV_curves(jV_curves, filename = "jV.png", selection = [], labels = []):
    plt.figure(figsize=(8, 6))
    if len(selection) == 0:
        jV_curves_todo = jV_curves
        linewidth = 0.02
        for curve in jV_curves_todo:
            plt.plot(curve['V'], curve['J'], linewidth=linewidth, alpha=0.5)
    else:
        jV_curves_todo = [jV_curves[idx] for idx in selection]
        linewidth = 1
        for idx, curve in enumerate(jV_curves_todo):
            plt.plot(curve['V'], curve['J'], linewidth=linewidth, alpha=0.5, label = labels[idx])
        plt.legend(loc="upper left")
    plt.xlabel("Voltage [V]")
    plt.ylabel("Current density [A/m²]")
    plt.axhline(0, color='k', linestyle='--')
    plt.axvline(0, color='k', linestyle='--')
    plt.ylim(-200, 200)
    plt.savefig(filename, dpi=300)
    plt.close()
    print("plotted jV curves")

def plot_EQE_curves(EQE_curves, filename = "EQE.png", selection = [], labels = []):
    plt.figure(figsize=(8, 6))
    if len(selection) == 0:
        EQE_curves_todo = EQE_curves
        linewidth = 0.02
        for curve in EQE_curves_todo:
            plt.plot(curve['lambda']*1e9, curve['EQE'], linewidth=linewidth, alpha=0.5)
    else:
        EQE_curves_todo = [EQE_curves[idx] for idx in selection]
        linewidth = 1
        for idx, curve in enumerate(EQE_curves_todo):
            plt.plot(curve['lambda']*1e9, curve['EQE'], linewidth=linewidth, alpha=0.5, label = labels[idx])
        plt.legend(loc="lower left")
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("EQE [%]")
    plt.ylim(0, 1)
    plt.savefig(filename, dpi=300)
    plt.close()
    print("plotted EQE spectra")


def plot_absorption_curves(uvvis_curves, filename = "absorption.png", selection = [], labels = []):
    plt.figure(figsize=(8, 6))
    if len(selection) == 0:
        uvvis_curves_todo = uvvis_curves
        linewidth = 0.02
        for curve in uvvis_curves_todo:
            plt.plot(curve['wavelengths'], curve['A_prime'], linewidth=linewidth, alpha=0.5)
    else:
        uvvis_curves_todo = [uvvis_curves[idx] for idx in selection]
        linewidth = 1
        for idx, curve in enumerate(uvvis_curves_todo):
            plt.plot(curve['wavelengths'], curve['A_prime'], linewidth=linewidth, alpha=0.5, label = labels[idx])
        plt.legend(loc="lower left")
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Effective absorption [%]")
    plt.savefig(filename, dpi=300)
    plt.close()
    print("plotted absorption spectra")

def plot_reflection_curves(uvvis_curves, filename = "reflection.png", selection = [], labels = []):
    plt.figure(figsize=(8, 6))
    if len(selection) == 0:
        uvvis_curves_todo = uvvis_curves
        linewidth = 0.02
        for curve in uvvis_curves_todo:
            plt.plot(curve['wavelengths'], curve['R'], linewidth=linewidth, alpha=0.5)
    else:
        uvvis_curves_todo = [uvvis_curves[idx] for idx in selection]
        linewidth = 1
        for idx, curve in enumerate(uvvis_curves_todo):
            plt.plot(curve['wavelengths'], curve['R'], linewidth=linewidth, alpha=0.5, label = labels[idx])
        plt.legend(loc="upper left")
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Reflection")
    plt.savefig(filename, dpi=300)
    plt.close()
    print("plotted reflection spectra")

# -------------------------------
# Main script
# -------------------------------
def main():
    # Path to the master results file and simulation runs directory.
    results_file = "results.txt"
    simulation_runs_dir = "simulation_runs"

    # Process all simulations.
    device_characteristics, jV_curves, EQE_curves, uvvis_curves, param_dicts_all, param_lists_all, nk_data_all, successful_ids = process_all_simulations(results_file, simulation_runs_dir)

    # save all the data
    jv_savedata = [curve['J'] for curve in jV_curves]
    np.savetxt("analysis/combined_data/jv_j.txt", jv_savedata)
    np.savetxt("analysis/combined_data/jv_v.txt", jV_curves[0]['V'])
    EQE_savedata = [curve['EQE'] for curve in EQE_curves]
    np.savetxt("analysis/combined_data/EQE.txt", EQE_savedata)
    np.savetxt("analysis/combined_data/EQE_wavelengths.txt", EQE_curves[0]['lambda'])
    A_savedata = [curve['A_prime'] for curve in uvvis_curves]
    np.savetxt("analysis/combined_data/absorption_A.txt", A_savedata)
    np.savetxt("analysis/combined_data/absorption_wavelengths.txt", uvvis_curves[0]['wavelengths'])
    R_savedata = [curve['R'] for curve in uvvis_curves]
    np.savetxt("analysis/combined_data/reflection_R.txt", R_savedata)
    np.savetxt("analysis/combined_data/reflection_wavelengths.txt", uvvis_curves[0]['wavelengths'])
    np.savetxt("analysis/combined_data/PCEs.txt", device_characteristics['PCE'])
    np.savetxt("analysis/combined_data/Jscs.txt", device_characteristics['Jsc'])
    np.savetxt("analysis/combined_data/Vocs.txt", device_characteristics['Voc'])
    np.savetxt("analysis/combined_data/FFs.txt", device_characteristics['FF'])
    np.savetxt("analysis/combined_data/param_lists.txt", param_lists_all)
    n_savedata = [curve[:,1] for curve in nk_data_all]
    k_savedata = [curve[:,2] for curve in nk_data_all]
    np.savetxt("analysis/combined_data/nk_n.txt", n_savedata)
    np.savetxt("analysis/combined_data/nk_k.txt", k_savedata)
    np.savetxt("analysis/combined_data/nk_lambda.txt", nk_data_all[0][:,0])

    # Plot distributions of PCE, Jsc, Voc, and FF.
    plot_distributions(device_characteristics, filename = "analysis/device_characteristics_all.png")
    # Plot all jV curves, all EQE curves, and all UV-vis absorption/reflection curves.
    plot_jV_curves(jV_curves, filename = "analysis/jV_all.png")
    plot_EQE_curves(EQE_curves, filename = "analysis/EQE_all.png")
    plot_absorption_curves(uvvis_curves, filename = "analysis/absorption_all.png")
    plot_reflection_curves(uvvis_curves, filename = "analysis/reflection_all.png")

    # Plot the best 5 jV, EQE curves, and UV-vis absorption/reflection curves.
    best_idxs = np.argsort(device_characteristics["PCE"])[-5:][::-1]
    print(best_idxs)
    labels = []
    for idx in best_idxs:
        #print([(key, device_characteristics[key][idx]) for key in device_characteristics.keys()])
        jsc = device_characteristics["Jsc"][idx]
        voc = device_characteristics["Voc"][idx]
        pce = device_characteristics["PCE"][idx]
        ff = device_characteristics["FF"][idx]*100.0
        jv_label = f"PCE={pce:.1f}%, jsc={jsc:.1f}A/m², Voc={voc:.2f}V, FF={ff:.1f}% ({successful_ids[idx]})"
        labels.append(jv_label)
        print(jv_label)
        os.system(f"cp -r simulation_runs/simulation_run_{successful_ids[idx]}/*png analysis/best/.")
    plot_jV_curves(jV_curves, filename = "analysis/jV_best.png", selection = best_idxs, labels = labels)
    plot_EQE_curves(EQE_curves, filename = "analysis/EQE_best.png", selection = best_idxs, labels = labels)
    plot_absorption_curves(uvvis_curves, filename = "analysis/absorption_best.png", selection = best_idxs, labels = labels)
    plot_reflection_curves(uvvis_curves, filename = "analysis/reflection_best.png", selection = best_idxs, labels = labels)

if __name__ == "__main__":
    main()

