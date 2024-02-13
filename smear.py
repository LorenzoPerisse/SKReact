#!/usr/bin/env python

__author__ = "Alex Goldsack"

""" 
    Produces a smearing matrix based on WIT energy reconstruction resolution
    for the SKReact app
"""

import matplotlib.pyplot as plt
import matplotlib
from params import *
import pandas as pd
import numpy as np
import math


def gaussian(x, mu, sig, c=1):
    return c * np.exp(-(x - mu) ** 2 / (2 * (sig ** 2)))


"""
    Contains the smearing information needed for folding MC/unfolding data
    calculated from a given set of gaussian parameters describing detectors's
    energy reconstruction for given true energies
"""


class Smear:
    def __init__(self, filename):
        wit_dat = pd.read_csv(filename, index_col="e")
        # Create series of energies
        energies_df = pd.DataFrame(
            np.nan, index=ENERGIES, columns=["c", "mu", "sig", "eff"]
        )
        # Get rid of points which overlap with WIT points
        energies_df = energies_df[~energies_df.index.isin(wit_dat.index)]
        # Concat with wit dat
        full_dat = pd.concat([wit_dat, energies_df])
        full_dat.sort_index(inplace=True)
        # Interpolate to fill the new values between and AFTER WIT points
        full_dat.interpolate(limit_direction="forward", inplace=True)

        # Get rid of WIT points we don't want
        full_dat = full_dat[full_dat.index.isin(ENERGIES)]
        # Check the shape of the matrix works
        if full_dat.shape[0] != E_BINS:
            print("SHAPE ISSUE IN SMEARING MATRIX")
            print("Check for floating point errors in smearing data.")
            print("Matrix shape:")
            print(full_dat.shape)
            exit()

        # Calc gaussian for each row with SMEAR_BINS bins
        # Modify area according to efficiency
        gauss_list = []
        gauss_ints = []
        for row in full_dat.itertuples():
            # Check if it is below the WIT smear data
            # Assume it won't be detected at all if so
            # print(row)
            if math.isnan(row.mu):
                smear_gauss = [0] * SMEAR_BINS
            else:
                # Need to multiply by new bin interval
                # to get right frequency density
                smear_gauss = [
                    gaussian(energy, row.mu, row.sig, SMEAR_INTERVAL * row.eff * row.c)
                    for energy in SMEAR_ENERGIES
                ]
            gauss_ints.append(np.trapz(smear_gauss, x=SMEAR_ENERGIES))

            gauss_list.append(smear_gauss)

        # import matplotlib.pyplot as plt
        # plt.plot(wit_dat.index.tolist(), wit_ints)
        # plt.plot(SMEAR_ENERGIES, gauss_ints)
        # wit_dat["eff"].plot()
        # plt.show()
        # For looking at example smearing gauss
        # Could use cycler to make colours prettier
        # for i in range(len(gauss_list)):
        #     if (i%(100) == 0 and (sum(gauss_list[i]) != 0)):
        #     # if (ENERGIES[i]%1 < 0.01 and (sum(gauss_list[i]) != 0)):
        #         # if(DOWN_ENERGIES[i] < 1-DEL_NP): continue
        #         # if(DOWN_ENERGIES[i] > 9-DEL_NP): break 
        #         plt.plot(DOWN_ENERGIES,gauss_list[i],
        #             label = "%i MeV" % ENERGIES[i],
        #             color = "C%i" % (i/100))
        #         plt.vlines(x=ENERGIES[i],
        #             ymin=0,
        #             ymax=0.006,
        #             color = "C%i" % (i/100),
        #             linestyles="--")
        #         plt.legend(loc="upper left")
        #         # print(i)
        #         # print(ENERGIES[i])
        #         # print(np.trapz(gauss_list[i],x=SMEAR_ENERGIES)/SMEAR_INTERVAL)
        #         # print()
        # plt.xlim(1-DEL_NP,E_MAX-DEL_NP)
        # plt.xlabel("Positron Energy")
        # plt.ylabel("Arbitrary Units")
        # plt.show()
        # exit()

        # Now just fill all NaNs with 0s
        full_dat.fillna(0, inplace=True)

        self.smear_mat = np.vstack(gauss_list)
        # self.inverse_smear = np.linalg.inv(self.smear_mat)
        self.effs = full_dat["eff"]
        return

    """
        Takes a pandas Series input int_spec_type spectrum, offsets it to the e+
        spectrum if needed and multiplies by the smearing matrix to produce a
        detected e+ spectrum of same length. Assumes the spectrum is vanishing 
        at the higher end.
    """

    def smear(self, int_spec):
        # Has to offset neutrino pos_spectrum to positron pos_spectrum
        int_spec = int_spec[DOWN_MASK]
        # And append zeroes to the spectrum
        int_spec = np.append(int_spec, np.zeros(E_BINS - int_spec.size))

        # print(np.trapz(np.multiply(int_spec,self.effs.to_numpy()),
        #     dx = E_INTERVAL))

        # fig, ax1 = plt.subplots()
        # ax1.plot(ENERGIES, int_spec)
        # ax2 = ax1.twinx()
        # ax2.plot(ENERGIES, self.effs)
        # plt.plot(ENERGIES, np.multiply(int_spec, self.effs.to_numpy()))
        # plt.show()
        # exit()

        # plt.plot(ENERGIES, int_spec)
        # plt.show()

        # The proof for this is left as an exercise to the reader
        smeared_np = np.matmul(int_spec, self.smear_mat)
        return smeared_np

    """
        Calculates inverse smearing matrix for unfolding from the already
        calculated matrix
    """

    def inverse_smear(self, spec):
        return np.matmul(spec.to_numpy(), self.inverse_mat)
