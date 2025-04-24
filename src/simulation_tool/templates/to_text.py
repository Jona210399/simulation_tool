def extract_lx(settings_dict: dict[str, str | int | float], lx: str):
    new_settings_dict = {}
    for key in settings_dict.keys():
        if f"{lx}." in key:
            new_key = key.replace(f"{lx}.", "")
            new_settings_dict[new_key] = settings_dict[key]
    return new_settings_dict


def get_simss_config(simss_config: dict[str, str | int | float]) -> str:
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
    return simss_config_text


def get_layer_config(
    layer_config: dict[str, str | int | float],
    layer_name: str,
) -> str:
    layer_config = extract_lx(layer_config, layer_name)
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
    return layer_config_text
