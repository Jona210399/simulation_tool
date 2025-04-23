import os
import sys
import uuid

import joblib
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pySIMsalabim as sim
from pySIMsalabim.experiments.EQE import run_EQE
from pySIMsalabim.experiments.JV_steady_state import run_SS_JV

import simulation_tool.templates as templates

np.random.seed(0)
if len(sys.argv) >= 2:
    try:
        np.random.seed(int(sys.argv[1]))
        print(f"replaced random seed with {sys.argv[1]}")
    except:
        pass

SESSION_PATH = "/home/jona/Promotion/Repositories/simsalabim/pySIMsalabim/SIMsalabim/SimSS"

SPECTRUM = "/home/jona/Promotion/Repositories/simsalabim/pySIMsalabim/SIMsalabim/Data/AM15G.txt"


        
def prep_optical_data(plot = False):
    optical_data = []
    for lineidx, line in enumerate(open("params/nk_P3HTPCBM_BHJ_original.txt", "r")):
        if lineidx>0:
            optical_data.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
    optical_data = np.array(optical_data)

    if plot:
        plt.figure()
        plt.plot(optical_data[:,0]*1e9, optical_data[:,1], "k-", label="n")
        plt.plot(optical_data[:,0]*1e9, optical_data[:,2], "b-", label="k")
        plt.legend()
        plt.xlabel("Wavelength [nm]")
        plt.ylabel("n/k [a.u.]")
        plt.savefig("optical_data_nk_P3HTPCBM_BHJ_original.png", dpi=120)
        plt.close()
    wavelengths = optical_data[:, 0]
    return(wavelengths)

def extract_lx(settings_dict, lx):
    new_settings_dict = {}
    for key in settings_dict.keys():
        if f"{lx}." in key:
            new_key = key.replace(f"{lx}.", "")
            new_settings_dict[new_key] = settings_dict[key]
    return(new_settings_dict)
    
def prep_random_simulation(index, manual_randomization = False):
    dpi = 60
    session_path_old = SESSION_PATH
    session_path = f"{os.getcwd()}/simulation_runs/simulation_run_{index}"
    if not os.path.exists("simulation_runs"):
        os.makedirs("simulation_runs")
    if not os.path.exists(session_path):
        os.makedirs(session_path)
        
    os.system(f"cp {session_path_old}/simss {session_path}/.")

    # FAKE N/K DATA
    ###############
    
    wavelengths = np.linspace(3e-7, 8e-7, 501)
    #wavelengths2 = prep_optical_data(plot = False)
    #print(np.max(np.abs(wavelengths - wavelengths2)))
    
    n_fake = generate_fake_n(wavelengths*1e9)
    k_fake, bandgap_nm = generate_fake_k(wavelengths*1e9)
    
    h = 4.135667696e-15  # Planck's constant in eV*s
    c = 2.99792458e8     # Speed of light in m/s
    wavelength_m = bandgap_nm * 1e-9  # Convert nm to meters
    bandgap_ev = (h * c) / wavelength_m  # Energy in eV
    
    plt.figure()
    plt.plot(wavelengths*1e9, n_fake, "k-", label="n")
    plt.plot(wavelengths*1e9, k_fake, "b-", label="k")
    plt.legend()
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("n/k [a.u.]")
    plt.savefig(f"{session_path}/optical_data_fake_{index}.png", dpi=dpi)
    plt.close()

    outfile=open(f"{session_path}/nk_{index}.txt", "w")
    outfile.write("lambda n k\n")
    for idx in range(len(wavelengths)):
        outfile.write(f"{wavelengths[idx]} {n_fake[idx]} {k_fake[idx]}\n")
    outfile.close()
    
    
    # CONFIG FILES
    ##############
    if manual_randomization:
        simss_config = templates.get_simss_config(session_path, session_path_old, index)
        layer1_config = templates.get_layer1_config(session_path, session_path_old, index)
        layer2_config = templates.get_layer2_config(session_path, session_path_old, index, bandgap_ev)
        layer3_config = templates.get_layer3_config(session_path, session_path_old, index)
        simss_config_text = templates.get_simss_config_text(simss_config)
        layer1_config_text = templates.get_layer_config_text(extract_lx(layer1_config, "l1"))
        layer2_config_text = templates.get_layer_config_text(extract_lx(layer2_config, "l2"))
        layer3_config_text = templates.get_layer_config_text(extract_lx(layer3_config, "l3"))
        randomized_settings = {}
        for key in simss_config.keys():
            randomized_settings[key] = simss_config[key]
        for key in layer1_config.keys():
            randomized_settings[key] = layer1_config[key]
        for key in layer2_config.keys():
            randomized_settings[key] = layer2_config[key]
        for key in layer3_config.keys():
            randomized_settings[key] = layer3_config[key]
    else:
        randomized_settings = templates.sample_random_parameters(session_path, session_path_old, index, bandgap_ev)
        simss_config_text = templates.get_simss_config_text(randomized_settings)
        layer1_config_text = templates.get_layer_config_text(extract_lx(randomized_settings, "l1"))
        layer2_config_text = templates.get_layer_config_text(extract_lx(randomized_settings, "l2"))
        layer3_config_text = templates.get_layer_config_text(extract_lx(randomized_settings, "l3"))

    
    simss_device_parameters_filename = f'{session_path}/simulation_setup_custom_{index}.txt'
    outfile=open(simss_device_parameters_filename, "w")
    outfile.write(simss_config_text)
    outfile.close()
        
    simss_device_parameters_filename = f'{session_path}/L1_parameters_{index}.txt'
    outfile=open(simss_device_parameters_filename, "w")
    outfile.write(layer1_config_text)
    outfile.close()
    simss_device_parameters_filename = f'{session_path}/L2_parameters_{index}.txt'
    outfile=open(simss_device_parameters_filename, "w")
    outfile.write(layer2_config_text)
    outfile.close()
    simss_device_parameters_filename = f'{session_path}/L3_parameters_{index}.txt'
    outfile=open(simss_device_parameters_filename, "w")
    outfile.write(layer3_config_text)
    outfile.close()

    # plot the band diagram
    x_position = 0  # starting position for the first layer
    W_L = randomized_settings["W_L"]
    W_R = randomized_settings["W_R"]
    L_TCO = randomized_settings["L_TCO"]*1e9
    L_BE = randomized_settings["L_BE"]*1e9
    Ls = [randomized_settings[f"l{index+1}.L"]*1e9 for index in range(3)]
    E_cs = [randomized_settings[f"l{index+1}.E_c"] for index in range(3)]
    E_vs = [randomized_settings[f"l{index+1}.E_v"] for index in range(3)]
    
    # Plot each layer as a box.
    fig, ax = plt.subplots(figsize=(10, 6))
    # layer 0
    width = L_TCO
    height = W_L
    color = f"C0"
    rect = patches.Rectangle((x_position, 0), width, height, linewidth=1, edgecolor="black", facecolor=color, alpha=0.3)
    ax.add_patch(rect)
    ax.plot([x_position, x_position + L_TCO], [height, height], color=color, linestyle="--")
    x_position += L_TCO
    # layers 1, 2, and 3
    for i in range(3):
        L = Ls[i]
        E_c = E_cs[i]
        E_v = E_vs[i]
        name = f"Layer {i+1}"
        width = L
        height1 = E_c
        height2 = 10.0 - E_v
        color = f"C{i+1}"
        # Create a rectangle representing the band gap of the layer.
        rect1 = patches.Rectangle((x_position, 0), width, height1, linewidth=1, edgecolor="black", facecolor=color, alpha=0.3)
        rect2 = patches.Rectangle((x_position, E_v), width, height2, linewidth=1, edgecolor="black", facecolor=color, alpha=0.3)
        ax.add_patch(rect1)
        ax.add_patch(rect2)
        # Plot dashed lines for the conduction and valence band edges.
        ax.plot([x_position, x_position + L], [E_c, E_c], color=color, linestyle="--")
        ax.plot([x_position, x_position + L], [E_v, E_v], color=color, linestyle="--")
        # Annotate the layer name in the middle of the box.
        ax.text(x_position + width/2, E_c + height/2, name, ha="center", va="center", fontsize=12, fontweight="bold")
        # Update the x-position for the next layer.
        x_position += L
    # layer 4
    width = L_BE
    height = W_R
    color = f"C0"
    rect = patches.Rectangle((x_position, 0), width, height, linewidth=1, edgecolor="black", facecolor=color, alpha=0.3)
    ax.add_patch(rect)
    ax.plot([x_position, x_position + L_BE], [height, height], color=color, linestyle="--")
    x_position += L_BE
    total_thickness = x_position
    
    # Set plot labels and title.
    ax.set_xlabel("Distance (m)")
    ax.set_ylabel("Energy (eV)")
    ax.set_title("Solar Cell Band Diagram")
    ax.set_xlim(-0.1*total_thickness, total_thickness+0.1*total_thickness)
    # Determine the overall energy range from all layers and contacts.
    #energy_min = min(min(E_cs), W_L, W_R) - 0.5
    #energy_max = max(max(E_vs), W_L, W_R) + 0.5
    ax.set_ylim(1.0, 7.0)
    
    plt.tight_layout()
    plt.savefig(f"{session_path}/band_diagram_{index}.png", dpi=dpi)
    plt.close()
    return(session_path)
    
    
def run_simulation(index, session_path, do_EQE = True):
    dpi = 60
    print(f"Run simulation {index}")
    # Set the JV parameters
    Gfrac = 1.0 # Fractions of the generation rate to simulate
    UUID = str(uuid.uuid4())

    # Run the JV simulation
    simss_device_parameters_filename = f'{session_path}/simulation_setup_custom_{index}.txt'
    try:
        ret, mess = run_SS_JV(simss_device_parameters_filename, session_path, JV_file_name = 'JV.dat', varFile= 'Var.dat', G_fracs = [Gfrac], parallel = False, max_jobs = 3, UUID=UUID) #, cmd_pars=[{'par': 'l2.L', 'val': '400e-9'}])
    except:
        print(f"ERROR: jV {index} failed")
        return(f"ERROR: jV {index} failed")
    if ret != 0:
        print(f"ERROR: jV {index} failed: {ret}, {mess}")
        return(f"ERROR: jV {index} failed: {ret}, {mess}")
    print(f"jV simulation {index} done")

    log_filename_old = f'{session_path}/log_Gfrac_1.0_{UUID}.txt'
    log_filename_new = f'{session_path}/log_jV_{index}.txt'
    os.system(f"mv {log_filename_old} {log_filename_new}")
    
    detailed_output_filename_old = f'{session_path}/Var_Gfrac_1.0_{UUID}.dat'
    os.system(f"rm {detailed_output_filename_old}")
    
    jV_filename_old = f'{session_path}/JV_Gfrac_1.0_{UUID}.dat'
    jV_filename_new = f'{session_path}/jV_{index}.txt'
    os.system(f"mv {jV_filename_old} {jV_filename_new}")
    data = pd.read_csv(jV_filename_new, sep=r'\s+')
    device_characteristics_filename_old = f'{session_path}/scPars_Gfrac_1.0_{UUID}.txt'
    device_characteristics_filename_new = f'{session_path}/device_characteristics_{index}.txt'
    if os.path.exists(device_characteristics_filename_old):
        os.system(f"mv {device_characteristics_filename_old} {device_characteristics_filename_new}")
        data2 = [float(x) for x in open(device_characteristics_filename_new, "r").readlines()[1].split()]
        jsc, voc, ff, vmpp, mpp = data2[0], data2[2], data2[8], data2[4], data2[6]
    else:
        print(f"ERROR: jV output {index} not found")
        return(f"ERROR: jV output {index} not found")
        
    plt.figure()
    plt.plot(data['Vext'], data['Jext'], label=f"Voc={voc:.2f}, Jsc={jsc:.1f}, FF={ff*100.0:.1f}%")
    #plt.plot(data['Vext'], data['Jext'])
    plt.xlabel('V$_{ext}$ [V]')
    plt.ylabel('J$_{ext}$ [A m$^{-2}$]')
    plt.xlim([-0.5, 1.7])
    plt.ylim([-200, 100])
    plt.plot([-10, 10], [0, 0], "k--")
    plt.plot([0, 0], [-1000.0, 1000.0], "k--")
    plt.legend(loc="upper left")
    plt.xlabel("Voltage [V]")
    plt.ylabel("Current density [A/m^2]")
    plt.savefig(f"{session_path}/jV_{index}.png", dpi=dpi)
    plt.close()
    
    uvvis_filename_old = f'{session_path}/AbsorptionSpectrum.txt'
    uvvis_filename_new = f'{session_path}/AbsorptionSpectrum_{index}.txt'
    os.system(f"mv {uvvis_filename_old} {uvvis_filename_new}")
    RT_filename_old = f'{session_path}/reflection_transmission_spectrum.txt'
    RT_filename_new = f'{session_path}/reflection_transmission_spectrum_{index}.txt'
    os.system(f"mv {RT_filename_old} {RT_filename_new}")
    uvvis_data = []
    for lineidx, line in enumerate(open(uvvis_filename_new, "r")):
        uvvis_data.append([float(line.split()[0]), float(line.split()[1])])
    uvvis_data = np.array(uvvis_data)
    RT_data = []
    for lineidx, line in enumerate(open(RT_filename_new, "r")):
        #RT_data.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
        RT_data.append([float(line.split()[0]), float(line.split()[1])])
    RT_data = np.array(RT_data)
    wavelengths = uvvis_data[:,0]*1e9
    A = uvvis_data[:,1]
    R = RT_data[:,1]
    A_prime = 1.0 - np.exp(-A)
    T = 1.0 - R - A_prime
    plt.figure()
    plt.plot(wavelengths, A_prime, "-", color="C0", label="absorption")
    plt.plot(wavelengths, R, "-", color="C1", label="reflection")
    plt.plot(wavelengths, T, "-", color="C2", label="transmission")
    #plt.plot(wavelengths, checksum, "-", color="C3", label="T+R+A")
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Absorption [a.u.]")
    plt.legend(loc="best")
    plt.savefig(f"{session_path}/uvvis_absorption_spectrum_{index}.png", dpi=dpi)
    plt.close()

    e_of_x_filename_old = f'{session_path}/E_of_x.txt'
    e_of_x_filename_new = f'{session_path}/E_of_x_{index}.txt'
    os.system(f"mv {e_of_x_filename_old} {e_of_x_filename_new}")
    e_of_x_data = []
    for lineidx, line in enumerate(open(e_of_x_filename_new, "r")):
        e_of_x_data.append([float(line.split()[0]), float(line.split()[1])])
    e_of_x_data = np.array(e_of_x_data)
    plt.figure()
    plt.plot(e_of_x_data[:,0]*1e9, e_of_x_data[:,1], "k-")
    plt.xlabel("Position [nm]")
    plt.ylabel("Squared electrical field [a.u.]")
    plt.savefig(f"{session_path}/E_of_x_{index}.png", dpi=60)
    plt.close()

    alpha_of_x_filename_old = f'{session_path}/alpha_of_x.txt'
    alpha_of_x_filename_new = f'{session_path}/alpha_of_x{index}.txt'
    os.system(f"mv {alpha_of_x_filename_old} {alpha_of_x_filename_new}")
    alpha_of_x_data = []
    for lineidx, line in enumerate(open(alpha_of_x_filename_new, "r")):
        alpha_of_x_data.append([float(line.split()[0]), float(line.split()[1])])
    alpha_of_x_data = np.array(alpha_of_x_data)
    plt.figure()
    plt.plot(alpha_of_x_data[:,0]*1e9, alpha_of_x_data[:,1]*1e-9, "k-")
    plt.xlabel("Position [nm]")
    plt.ylabel("Absorption coefficient [1/nm]")
    plt.savefig(f"{session_path}/alpha_of_x_{index}.png", dpi=dpi)
    plt.close()

    if not do_EQE:
        os.system(f'rm {session_path}/simss')
        return(f"{index} DONE")
    
    lambda_min = 350
    lambda_max = 900
    lambda_step = 10
    Vext = 0
    EQE_outfile_name = 'EQE.dat'
    spectrum = SPECTRUM
    try:
        ret_eqe, mess_eqe = run_EQE(simss_device_parameters_filename, session_path, spectrum, lambda_min, lambda_max, lambda_step, Vext, outfile_name = EQE_outfile_name, JV_file_name = 'JV.dat', run_mode = False, parallel = False, force_multithreading = False, threadsafe=True)
    except:
        print(f"ERROR: EQE {index} failed")
        return(f"ERROR: EQE {index} failed")
    if ret != 0:
        print(f"ERROR: EQE {index} failed: {ret}, {mess}")
        return(f"ERROR: EQE {index} failed: {ret}, {mess}")
        
    print(f"EQE simulation {index} done")
    
    EQE_filename_old = f'{session_path}/EQE.dat'
    EQE_filename_new = f'{session_path}/EQE_{index}.dat'
    os.system(f"mv {EQE_filename_old} {EQE_filename_new}")
    if os.path.exists(EQE_filename_new):
        df = pd.read_csv(EQE_filename_new, sep = r'\s+')
    else:
        print(f"ERROR: EQE output {index} not found")
        return(f"ERROR: EQE output {index} not found")
    plt.plot(df['lambda']*1e9,df['EQE'])
    plt.xlabel('Wavelength [nm]')
    plt.ylabel('EQE')
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("EQE [%%]")
    plt.savefig(f"{session_path}/EQE_{index}.png", dpi=dpi)
    plt.close()

    # Clean up the output files (comment out if you want to keep the output files)
    #sim.clean_all_output(session_path)
    sim.delete_folders('tmp', session_path)
    os.system(f'rm {session_path}/simss')
    return(f"{index} DONE")

def run_random_simulation(index, manual_randomization = False, do_EQE = True):
    session_path = prep_random_simulation(index, manual_randomization)
    success = run_simulation(index, session_path, do_EQE)
    return(success)
    
def gaussian(x, mu, s, alpha):
    """
    Computes the Gaussian function with a given mean (mu), standard deviation (s), and prefactor (alpha).
    
    Parameters:
        x (float or array-like): The input value(s) where the Gaussian is evaluated.
        mu (float): The mean of the Gaussian.
        s (float): The standard deviation of the Gaussian.
        alpha (float): The prefactor (amplitude) of the Gaussian.
    
    Returns:
        float or np.ndarray: The value of the Gaussian function at x.
    """
    return alpha * np.exp(-((x - mu) ** 2) / (2 * s ** 2))
    
def generate_fake_n(wavelengths):
    num_gaussians = 10
    means = np.random.random(num_gaussians)*600+200
    sigmas = np.random.random(num_gaussians)*60+10
    alphas = np.random.random(num_gaussians)*0.5
    results = np.zeros((len(wavelengths)))+1.5
    for idx in range(num_gaussians):
        results += gaussian(wavelengths, means[idx], sigmas[idx], alphas[idx])
    return(results)

def generate_fake_k(wavelengths):
    num_gaussians = 10
    bandgap = np.random.random()*200+550
    means = np.random.random(num_gaussians)*(bandgap-200)+200
    sigmas = np.random.random(num_gaussians)*60+10
    alphas = np.random.random(num_gaussians)*0.2
    results = np.zeros((len(wavelengths)))
    for idx in range(num_gaussians):
        results += gaussian(wavelengths, means[idx], sigmas[idx], alphas[idx])
    return(results, np.max(means))


run_random_simulation(0, manual_randomization = False, do_EQE = True)

exit()



parallel = joblib.Parallel(8, return_as="generator")

startindex = 0
num_todo = 8
jobs = []
for index in range(startindex, startindex+num_todo):
    j = joblib.delayed(run_random_simulation)(index)
    #run_random_simulation(index)
    jobs.append(j)

results = list(parallel(jobs))
outfile=open("results.txt", "w")
for idx, x in enumerate(results):
    outfile.write(f"{idx} {x}\n")
outfile.close()




