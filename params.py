#!/usr/bin/env python

__author__ = "Alex Goldsack"

"""
Various parameters controlling SKReact, kept here for cleanliness
"""

import pandas as pd
import numpy as np
import ROOT



# Logo name
LOGO_FILE = "skreact_logo.png"

# Reactor .xls info filepath
REACT_DIR = "./react_p/"
REACT_PICKLE = "reactors_main.pkl"

# Force SKReact to import info from .xls files
FORCE_XLS_IMPORT = True

# Prints reactor names if true. Prints filenames regardless.
VERBOSE_IMPORT = False

# Prints errors when importing from .xls file if true.
VERBOSE_IMPORT_ERR = False

# Roughly select reactors within this range using
# approx. long/lat distance, 1 deg = 111 km
R_THRESH = 1000000.0  # km
R_THRESH_DEG = R_THRESH / 111.0

# Geoneutrinos luminosity filename
GEO_FILE = "geoneutrino-luminosity.knt"

# Smearing information
WIT_SMEAR_FILE = "smear_main.csv"

# SKReact Parameters
# WIN_X,Y are not used, get screen resolution instead.
WIN_SCALE = 1
WIN_X = 900
WIN_Y = 1100
# Matplotlib def figure sizes
FIG_SCALE = WIN_SCALE
FIG_X = FIG_SCALE*4
FIG_Y = FIG_SCALE*2.5

# Super-Kamiokande Info
SK_LAT = 36.4267  # deg
SK_LONG = 137.3104  # deg
SK_ALT = 0.370  # km
SK_R = 1690  # radius /cm
SK_HH = 1810  # half height /cm

# E Spectrum Hyper-parameters
# +1 on bins because it's inclusive to the last E
E_MIN = 0.0005  # MeV
E_MAX = 15.9995  # MeV
E_BINS = 16000  # Use tidy numbers
E_INTERVAL = (E_MAX - E_MIN + 0.001) / (E_BINS)
# List of energies to calculate spectrum at
_energies = np.linspace(E_MIN, E_MAX, E_BINS)
# Linspace has rounding errors
ENERGIES = np.array([float("%.4f" % energy) for energy in _energies])

# Smearing is quite intensive, so give option to reduce
# number of bins
SMEAR_BINS = E_BINS  # Number of bins to show smeared spec
SMEAR_INTERVAL = (E_MAX - E_MIN + 0.001) / (SMEAR_BINS)

# List of energies to calculate smeared spectrum at
_smear_energies = np.linspace(E_MIN, E_MAX, SMEAR_BINS)
SMEAR_ENERGIES = np.array([float("%.4f" % energy) for energy in _smear_energies])

# Scaling factor to make SKReact flux match others
# Shouldn't be needed, but temporary fix
FLUX_SCALE = 1
# FLUX_SCALE = 1339/731 # Geoneutrinos.org
# FLUX_SCALE = 2106/1785 # VERY ROUGH KamLAND


# INTERACTION =================================================================
IBD_MIN = 1.806  # MeV
N_P_O = 8  # n protons in O
SK_FM = 22.5  # kt
SK_ID_M = 32.0  # kt

# PDG 2022 constant
M_E = 0.51099895000  # MeV (mass of e)
DEL_NP = 1.293  # MeV (Diff between n and p)
M_P = 938.27208816  # MeV (mass of p)
M_O_U = 15.999  # u
U = 1.66053906660e-27  # kg


# CHOOSE ID_M or FM BASED ON WHAT THE SMEARING INFO WAS CALCULATED USING
M_WATER = (2 + M_O_U) * U  # kg
N_WATER_SK = SK_ID_M * 1e6 / M_WATER
SK_N_P = N_WATER_SK * 2  # Free protons TODO: Look into O interactions


# List of offset energies to use as index for offset spectra
UP_ENERGIES = np.array([float("%.4f" % (energy + DEL_NP)) for energy in _energies])
DOWN_ENERGIES = np.array([float("%.4f" % (energy - DEL_NP)) for energy in _energies])
OFFSET_UP_DICT = dict(zip(ENERGIES, UP_ENERGIES))
OFFSET_DOWN_DICT = dict(zip(ENERGIES, DOWN_ENERGIES))
UNDO_OFFSET_UP_DICT = dict(zip(UP_ENERGIES, ENERGIES))
UNDO_OFFSET_DOWN_DICT = dict(zip(DOWN_ENERGIES, ENERGIES))
# Get which elements overlap between ENERGIES and DOWN ENERGIES
DOWN_OVERLAP_ENERGIES, ENERGIES_OVERLAP_I_DOWN, DOWN_OVERLAP_I = np.intersect1d(
    ENERGIES, DOWN_ENERGIES, return_indices=True
)
# Smearing matrix corresponds to ENERGIES
# Mask out energies in DOWN_ENERGIES not in ENERGIES
DOWN_MASK = DOWN_ENERGIES > ENERGIES[0]

# REACTOR NU PRODUCTION =======================================================

# Fuel factors for different reactors (fission/power)
# TODO: Which values are best? Also, MOX is different for different reactors
# TODO: Get better values for FBRs rather than copying PWR
# NOTE : HTGR is a newly added reactor type in SKReact-v2.0
core_types = ["PWR", "BWR", "LWGR", "GCR", "HTGR", "PHWR", "MOX", "FBR"]

#Based on DB fission fraction
__fuel_makeup_data = {
    "U_235":  [0.564, 0.564, 0.564, 0.564, 0.564, 0.543, 0.000, 0.564],
    "Pu_239": [0.304, 0.304, 0.304, 0.304, 0.304, 0.411, 0.081, 0.304],
    "U_238":  [0.076, 0.076, 0.076, 0.076, 0.076, 0.022, 0.708, 0.076],
    "Pu_241": [0.056, 0.056, 0.056, 0.056, 0.056, 0.024, 0.212, 0.056],
}

#Old SKReact value = 1st value in Baldoncini table
# __fuel_makeup_data = {
#     "U_235":  [0.538, 0.538, 0.538, 0.538, 0.538, 0.543, 0.000, 0.538],
#     "Pu_239": [0.328, 0.328, 0.328, 0.328, 0.328, 0.411, 0.081, 0.328],
#     "U_238":  [0.078, 0.078, 0.078, 0.078, 0.078, 0.022, 0.708, 0.078],
#     "Pu_241": [0.056, 0.056, 0.056, 0.056, 0.056, 0.024, 0.212, 0.056],
# }

#Average value based on Baldoncini = Arithmetic mean (no weight considered here)
# __fuel_makeup_data = {
#     "U_235":  [0.573, 0.573, 0.573, 0.573, 0.573, 0.543, 0.000, 0.573],
#     "Pu_239": [0.299, 0.299, 0.299, 0.299, 0.299, 0.411, 0.081, 0.299],
#     "U_238":  [0.075, 0.075, 0.075, 0.075, 0.075, 0.022, 0.708, 0.075],
#     "Pu_241": [0.053, 0.053, 0.053, 0.053, 0.053, 0.024, 0.212, 0.053],
# }

#Lowest in 235U
# __fuel_makeup_data = {
#     "U_235":  [0.480, 0.480, 0.480, 0.480, 0.480, 0.543, 0.000, 0.480],
#     "Pu_239": [0.370, 0.370, 0.370, 0.370, 0.370, 0.411, 0.081, 0.370],
#     "U_238":  [0.070, 0.070, 0.070, 0.070, 0.070, 0.022, 0.708, 0.070],
#     "Pu_241": [0.080, 0.080, 0.080, 0.080, 0.080, 0.024, 0.212, 0.080],
# }

#Highest in 235U
# __fuel_makeup_data = {
#     "U_235":  [0.650, 0.650, 0.650, 0.650, 0.650, 0.543, 0.000, 0.650],
#     "Pu_239": [0.240, 0.240, 0.240, 0.240, 0.240, 0.411, 0.081, 0.240],
#     "U_238":  [0.070, 0.070, 0.070, 0.070, 0.070, 0.022, 0.708, 0.070],
#     "Pu_241": [0.040, 0.040, 0.040, 0.040, 0.040, 0.024, 0.212, 0.040],
# }

#Test
# __fuel_makeup_data = {
#     "U_235":  [0.564, 0.564, 0.564, 0.564, 0.564, 0.564, 0.564, 0.564],
#     "Pu_239": [0.304, 0.304, 0.304, 0.304, 0.304, 0.304, 0.304, 0.304],
#     "U_238":  [0.076, 0.076, 0.076, 0.076, 0.076, 0.076, 0.076, 0.076],
#     "Pu_241": [0.056, 0.056, 0.056, 0.056, 0.056, 0.056, 0.056, 0.056],
# }



FUEL_MAKEUP = pd.DataFrame(__fuel_makeup_data, index=core_types)

# E per fission (denom of fuel factor's num)
U_235_Q   = 202.36
U_235_DQ  = 0.26
U_238_Q   = 205.99
U_238_DQ  = 0.52
PU_239_Q  = 211.12
PU_239_DQ = 0.34
PU_241_Q  = 214.26
PU_241_DQ = 0.33

# 5th order polynomial coeffs and error for E_spectra 1 thru 6
# TODO: Put into dataframe
E_SPEC_N_ORDER = 5
U_235_A = [3.217, -3.111, 1.395, -3.690e-1, 4.445e-2, -2.053e-3]
U_235_DA = [4.09e-2, 2.34e-2, 4.88e-3, 6.08e-4, 7.77e-5, 6.79e-6]

PU_239_A = [6.413, -7.432, 3.535, -8.820e-1, 1.025e-1, -4.550e-3]
PU_239_DA = [4.57e-2, 2.85e-2, 6.44e-3, 9.11e-4, 1.38e-4, 1.29e-5]

U_238_A = [4.833e-1, 1.927e-1, -1.283e-1, -6.762e-3, 2.233e-3, -1.536e-4]
U_238_DA = [1.24e-1, 5.86e-2, 1.11e-2, 1.92e-3, 2.84e-4, 2.86e-5]

PU_241_A = [3.251, -3.204, 1.428, -3.675e-1, 4.254e-2, -1.896e-3]
PU_241_DA = [4.37e-2, 2.60e-2, 5.66e-3, 7.49e-4, 1.02e-4, 9.03e-6]


# OSCILLATION =================================================================
# Osc. Params
# S_XX = sin^2(theta_XX), DM_blah = delta(m_blah^2)
# S_2_XX = sin^2(2*theta_XX)
# NH = Normal Hierarchy, IH = Inverted

# From NuFit Jul 2019 (with SK atm constraints)  IBD rate = 8.89896
# DM_21 = 7.39e-5  # eV^2
# DM_31 = 2.528e-3  # eV^2
# DM_23 = 2.510e-3  # eV^2
# S_12 = 0.310
# THET_12 = 33.82
# S_23_NH = 0.563
# S_23_IH = 0.565
# THET_23 = 48.6
# S_13_NH = 0.02237
# S_13_IH = 0.02259
# THET_13 = 8.6
# S_13_IH = 0.02223
# S_23_IH = 0.569
# DM_23 = DM_31 - DM_21  # eV^2




# From NuFit November 2022 (with SK atm constraints)  IBD rate = 8.9693
# DM_21 = 7.41e-5  # eV^2
# DM_31 = 2.507e-3  # eV^2
# DM_23 = DM_31 - DM_21  # eV^2
# S_12  = 0.303
# THET_12 = 33.41
# S_23_NH = 0.451
# S_23_IH = 0.539
# S_13_NH = 0.02225
# S_13_IH = 0.0220
# THET_13 = 8.58
# THET_23 = 42.2


# From PDG 2022 Prog.Theor.Exp.Phys. 2022, 083C01 (2022)
DM_21 = 7.53e-5  # eV^2
DM_23 = 2.453e-3  # eV^2 (NH)
DM_31 = DM_21 + DM_23
S_12 = 0.307
THET_12 = 33.64708 # calculated from S_12
S_23_NH = 0.546
S_23_IH = 0.539
S_13_NH = 0.0220
S_13_IH = 0.0220
THET_13 = 8.52981 # calculated from S_13_NH/IH


# +1sigma -> decrease rate  IBD rate = 8.78933
# DM_21 = 7.62e-5  # eV^2       #-> decrease rate  IBD rate = 8.92358
# DM_31 = 2.533e-3  # eV^2      #-> no impact      IBD rate = 8.9693
# S_12 = 0.315                  #-> decrease rate  IBD rate = 8.84577
# THET_12 = 34.16
# S_23_NH = 0.470               #-> no impact  IBD rate = 8.9693
# THET_23 = 43.3
# S_13_NH = 0.02281             #-> decrease rate  IBD rate = 8.95928
# THET_13 = 8.69

# -1sigma -> increase rate IBD rate = 9.16304
# DM_21 = 7.21e-5  # eV^2         #-> increase rate  IBD rate = 9.01216
# DM_31 = 2.480e-3  # eV^2        #-> no impact      IBD rate = 8.9693
# S_12 = 0.291                    #-> increase rate  IBD rate = 9.11055
# THET_12 = 32.69
# S_23_NH = 0.435                 #-> no impact  IBD rate = 8.9693
# THET_23 = 41.3
# S_13_NH = 0.02166               #-> increase rate  IBD rate = 8.97987
# THET_13 = 8.47


C_12 = 1 - S_12
C_23_NH = 1 - S_23_NH
C_23_IH = 1 - S_23_IH
C_13_NH = 1 - S_13_NH
C_13_IH = 1 - S_13_IH
S_2_12 = 4 * S_12 * C_12
S_2_13 = 4 * S_13_NH * C_13_NH


# FITTING ====================================================================

# Directory where the files to fit are
FIT_DIR = "./"

# Number of steps for each parameter per cycle
N_STEPS = 10
# Number of cycles with finer granularity
N_CYCLES = 3
# Factor of n cycle parameter space to cover for n+1 cycle
CYCLE_FACTOR = 0.1
# Search within +- FIT_RANGE_SCALE*DM_21_FIT
FIT_RANGE_SCALE = 0.2

# Defining (half) range to fit in, comments after are the 1 sig errors from nufit
DM_21_FIT = DM_21
DM_21_RANGE = FIT_RANGE_SCALE * DM_21_FIT
THET_12_FIT = THET_12
THET_12_RANGE = FIT_RANGE_SCALE * THET_12_FIT
DM_31_FIT = DM_31
DM_31_RANGE = FIT_RANGE_SCALE * DM_31_FIT
THET_13_FIT = THET_13
THET_13_RANGE = FIT_RANGE_SCALE * THET_13_FIT


# Misc.============================================================================

EARTH_R = 6371  # km
EARTH_R_POLAR = 6356  # km
EARTH_R_EQUATOR = 6378  # km

EV_J = 1.60218e-19  # = 1 eV in J

POSITRON_PDG = -11  # pdg code for positron


#============================================================================
# Get precomputed PWR spectrum

inFileName = "/mnt/c/Users/lolos/Documents/Soft/skreact/SKReact/data/ReactorModel_SK.root"
inFile = ROOT.TFile.Open( inFileName , "UPDATE" )

spectrum = inFile.Get( "nspec_Tot" )
spectrum.SetDirectory(0)

covariance = inFile.Get( "ncov_Tot" )
covariance.SetDirectory(0)

hbin = spectrum.GetNbinsX()
emin = spectrum.GetBinLowEdge(1)
emax = spectrum.GetBinLowEdge(hbin+1)

hTot_SK = ROOT.TH1D( 'nspec_SK', 'nspec_SK', hbin, emin, emax )
oscillation_prefactor = ROOT.TH1D( 'oscillation_prefactor_SK', 'oscillation_prefactor_SK', hbin, emin, emax )
hTot_SK.SetDirectory(0)
oscillation_prefactor.SetDirectory(0)

IsTFileUpdated = False


