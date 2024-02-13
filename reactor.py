#!/usr/bin/env python

__author__ = "Alex Goldsack"

""" 
The portion of SKReact dealing with neutrino production in
nuclear reactors and subsequent oscillation.
"""

import params
from params import *
import pandas
from math import sin, cos, tan, sqrt, radians
from calendar import monthrange
import numpy as np
import math

# Calculating xsec for each energy
e_e = lambda e: e - DEL_NP
p_e = lambda e: math.sqrt(e_e(e) ** 2 - M_E * M_E)
e_exp = lambda e: e ** (
    -0.07056 + 0.02018 * math.log(e) - 0.001953 * (math.log(e)) ** 3
)
xsec = lambda e: 1e-43 * p_e(e) * e_e(e) * e_exp(e)  # cm^2

# Set up list of xsecs for set of ENERGIES
xsecs = np.array([xsec(e) if e > IBD_MIN else 0 for e in ENERGIES])



class Reactor:

    # Initialiser
    def __init__(
        self,
        country,
        name,
        latitude,
        longitude,
        elevation,
        core_type,
        mox,
        p_th,
        lf_monthly,
        default=True,
        calc_spec=False,
    ):

        self.country = country
        self.name = name
        self.latitude  = latitude
        self.longitude = longitude
        self.elevation = elevation
        # core_type is checked later, need to remove whitespace
        self.core_type = core_type.rstrip()
        self.mox = mox
        self.p_th = p_th  # MW, series of yearly reference power
        self.dist_to_sk = self._dist_to_sk()
        self.lf_monthly = lf_monthly  # Pandas series
        self.p_monthly = self._p_monthly()
        self.p_r_monthly = self._p_r_monthly()
        self.n_ints_monthly = []  # Will be set to series later
        self.default = default  # If the reactor came from the xls
        self.current_flux = -1.0
        if calc_spec:
            self.prod_spec          = self._prod_spec()  # Produced
            self.def_osc_pre_factor = self.osc_pre_factor()  # Total oscillated prefactor at SK
            self.def_osc_spec       = self.osc_spec()  # Oscillated
            self.def_int_spec       = self.int_spec(self.def_osc_spec)  # Interacted
            self.dir_flux_at_sk     = self._dir_flux_at_sk()  # Direction of nu flux at SK

    # Monthly power output calculate from load factor and p_th
    def _p_monthly(self):
        # Same format as lf
        index = self.lf_monthly.index
        # p_list = [self.p_th * lf/100 for lf in self.lf_monthly.tolist()]
        p_list = []
        for date, lf in self.lf_monthly.items():
            p_list.append(self.p_th[date[:4]] * lf / 100)
        return pd.Series(p_list, index=index)

    # Monthly power/r^2 output calculate from p_monthly and dist_to_sk
    def _p_r_monthly(self):
        # Same format as lf
        index = self.p_monthly.index
        p_r_list = [p / (self.dist_to_sk ** 2) for p in self.p_monthly.tolist()]
        return pd.Series(p_r_list, index=index)

    # Need an explicit function rather than just append,
    def add_to_lf(self, date, lf):
        # print(self.name)
        # print(lf)
        self.lf_monthly.loc[date] = lf
        self.p_monthly.loc[date] = lf *self.p_th[date[:4]]
        self.p_r_monthly[date] = lf * self.p_th[date[:4]] / (self.dist_to_sk ** 2)
        # self.lf_monthly.set_value(date, lf)
        # self.p_monthly.set_value(date, lf * self.p_th[date[:4]])
        # self.p_r_monthly.set_value(
        #     date, lf * self.p_th[date[:4]] / (self.dist_to_sk ** 2)
        # )
        return

    # Sets on sets on sets
    def set_country(self, country):
        self.country = country
        return

    def set_name(self, name):
        self.name = name
        return

    # Need to recalculate dist to sk when changing pos
    def set_latitude(self, latitude):
        self.latitude = latitude
        self.dist_to_sk = self._dist_to_sk()
        return

    def set_longitude(self, longitude):
        self.longitude = longitude
        self.dist_to_sk = self._dist_to_sk()
        return
    
    def set_elevation(self, elevation):
        self.elevation = elevation
        self.dist_to_sk = self._dist_to_sk()
        return

    # Need to reproduce E spec if core type changes
    def set_core_type(self, core_type):
        self.core_type = core_type
        self.prod_spec = self._prod_spec()
        return

    def set_mox(self, mox):
        self.mox = mox
        self.prod_spec = self._prod_spec()
        return

    def set_p_th(self, p_th):
        self.p_th = p_th
        return

    # Should be pd Series, maybe I should assert?
    def set_lf_monthly(self, lf_monthly):
        self.lf_monthly = lf_monthly
        self.p_monthly = self._p_monthly()
        self.p_r_monthly = self._p_r_monthly()
        return

    def set_all_spec(self):
        self.set_prod_spec()
        self.set_osc_pre_factor()
        self.set_osc_spec()
        self.set_int_spec()
        return

    # This depends on number of energy bins so need to recalculate
    # when number of bins change
    def set_prod_spec(self):
        self.prod_spec = self._prod_spec()

    def set_osc_spec(self):
        self.def_osc_spec = self.osc_spec(period="/s")

    def set_osc_pre_factor(self):
        self.def_osc_pre_factor = self.osc_pre_factor(period="/s")

    # Interacted spec takes osc spec as argument anyway
    def set_int_spec(self):
        self.def_int_spec = self.int_spec(self.def_osc_spec)

    # Directionality of nu flux at SK
    def set_dir_flux_at_sk(self):
        self.dir_flux_at_sk = self._dir_flux_at_sk()

    # Calculate the number of neutrinos produced in given period
    # CURRENT STATE IS DEPRICATED, CANNOT GUARANTEE IT PRODUCES GOOD NUMBERS
    # TODO: Move the common calcs outside the if statement
    def n_nu(self, period="Max"):
        # Pre-calculating the nu per second at reference power for self
        nu_per_s = self.p_th * NU_PER_MW
        if period == "Max" or period == "max":  # Yearly at reference P
            return 365 * 24 * 60 * 60 * nu_per_s
        elif len(period) == 15:  # Inclusive period YYYY/MM-YYYY/MM
            year_start = int(period[:4])
            month_start = int(period[5:7])
            year_end = int(period[8:12])
            month_end = int(period[13:])

            # Cycle through all months calculating nu per month
            month_range_start = month_start
            month_range_end = 13
            n_nu_tot = 0
            for year in range(year_start, year_end + 1):
                # Start from Jan after first year
                if year != year_start:
                    month_range_start = 1
                # Only go up to end of period in final year
                if year == year_end:
                    month_range_end = month_end + 1  # For inclusivity
                for month in range(month_range_start, month_range_end):
                    n_days_in_month = monthrange(year, month)[1]
                    # Query the specific month from the LF series
                    lf_month = self.lf_monthly["%i/%02i" % (year, month)]
                    lf_month /= 100  # To be a factor, not %age
                    n_nu_month = n_days_in_month * 24 * 60 * 60
                    n_nu_month *= lf_month * nu_per_s

                    n_nu_tot += n_nu_month

            return n_nu_tot
        elif len(period) == 7:  # Specific month YYYY/MM
            year = int(period[:4])
            month = int(period[5:])
            n_days_in_month = monthrange(year, month)[1]
            lf_month = self.lf_monthly["%i/%02i" % (year, month)]
            lf_month /= 100
            n_nu_month = n_days_in_month * 24 * 60 * 60
            n_nu_month *= lf_month * nu_per_s
            return n_nu_month
        elif len(period) == 4:  # Specific year YYYY
            year = int(period)

            # Cycle through all months calculating nu per month
            n_nu_tot = 0
            for month in range(1, 13):
                n_days_in_month = monthrange(year, month)[1]
                # Query the specific month from the LF series
                lf_month = self.lf_monthly["%i/%02i" % (year, month)]
                lf_month /= 100

                n_nu_month = n_days_in_month * 24 * 60 * 60
                n_nu_month *= lf_month * nu_per_s

                n_nu_tot += n_nu_month

            return n_nu_tot
        else:
            print(
                'reactor.n_nu() requires either YYYY/MM, YYYY or "Max" '
                "(per year) for period of nu production."
            )
            exit()

    """ 
    Earth bulges a the equator, this gives distance in km to
    centre of the Earth as a function of latitude
    """

    def _dist_to_earth_centre(self, latitude):
        a = EARTH_R_EQUATOR ** 2 * cos(latitude)
        b = EARTH_R_POLAR ** 2 * sin(latitude)
        c = EARTH_R_EQUATOR * cos(latitude)
        d = EARTH_R_POLAR * sin(latitude)

        r = sqrt((a * a + b * b) / (c * c + d * d))

        return r

    """
    Returns sin of geocentric latitude from geodetic latitude
    """

    def _sin_geocentric(self, latitude):
        tan_a = EARTH_R_POLAR * tan(latitude) / EARTH_R_EQUATOR
        sin_a = tan_a / sqrt(1 + tan_a * tan_a)
        return sin_a

    """
    Returns cos of geocentric latitude from geodetic latitude
    """

    def _cos_geocentric(self, latitude):
        tan_a = EARTH_R_POLAR * tan(latitude) / EARTH_R_EQUATOR
        cos_a = 1 / sqrt(1 + tan_a * tan_a)
        return cos_a

    """
    Use Lat and Long info to calc distance to SK in km
    Assume reactors are at sea level, very reasonable
    assumption, given most are on coastline
    """

    def _dist_to_sk(self):
        lat_react_rad = radians(self.latitude)
        long_react_rad = radians(self.longitude)
        elev_react = self.elevation * 0.001     # elevation initially given in meter, needed in km

        lat_sk_rad = radians(SK_LAT)
        long_sk_rad = radians(SK_LONG)

        r_react = self._dist_to_earth_centre(lat_react_rad) + elev_react    #in km
        r_sk = self._dist_to_earth_centre(lat_sk_rad) + SK_ALT              #in km

        x_react = r_react * self._cos_geocentric(lat_react_rad) * cos(long_react_rad)
        y_react = r_react * self._cos_geocentric(lat_react_rad) * sin(long_react_rad)
        z_react = r_react * self._sin_geocentric(lat_react_rad)

        x_sk = r_sk * self._cos_geocentric(lat_sk_rad) * cos(long_sk_rad)
        y_sk = r_sk * self._cos_geocentric(lat_sk_rad) * sin(long_sk_rad)
        z_sk = r_sk * self._sin_geocentric(lat_sk_rad)

        dist = sqrt(
            (x_react - x_sk) ** 2 + (y_react - y_sk) ** 2 + (z_react - z_sk) ** 2
        )

        return dist


    """
    Getting E spectrum from text file
    """
    def _f_from_TFile(self, energy):
        bin = spectrum.FindBin(energy)
        flux = spectrum.GetBinContent(bin)
        return flux
    

    """
    Getting E spectrum from text file
    """
    def _ferr_from_TFile(self, energy):
        bin = spectrum.FindBin(energy)
        err = np.sqrt( covariance.GetBinContent(bin, bin) )
        return err
    

    """
    Use precomputed E spectrum produced by PWR reactors,
    NOTE: The spectrum produced is PER SECOND at reference power p_th
    """
    def _prod_spec(self):
        core_type = self.core_type
        if self.mox:
            core_type = "MOX"

        # Fuel fractions for this type of core (p_i in literature)
        u_235_frac  = FUEL_MAKEUP.loc[core_type]["U_235"]
        pu_239_frac = FUEL_MAKEUP.loc[core_type]["Pu_239"]
        u_238_frac  = FUEL_MAKEUP.loc[core_type]["U_238"]
        pu_241_frac = FUEL_MAKEUP.loc[core_type]["Pu_241"]

        # E release per fission in MJ
        # P is in MW Q is in MeV, so change Q to MJ
        E_fission = (u_235_frac*U_235_Q + u_238_frac*U_238_Q + pu_239_frac*PU_239_Q + pu_241_frac*PU_241_Q) * EV_J

        tot_spectrum = [
            (1./E_fission) * self._f_from_TFile(energy)
            for energy in ENERGIES
        ]

        spectrum_data = {
            "Total": tot_spectrum,
        }

        prod_spec_dat = list(
            zip(
                tot_spectrum,
            )
        )

        prod_spec = np.array(
            prod_spec_dat,
            dtype=[
                ("Total", "f4"),
            ],
        )

        return prod_spec

    """
    Produces tuple of maximum and minimum e spectra, calculated
    by finding the max and min coeffs
    """
    def _prod_spec_err(self):
        core_type = self.core_type
        if self.mox:
            core_type = "MOX"

        # Fuel fractions for this type of core (p_i in literature)
        u_235_frac  = FUEL_MAKEUP.loc[core_type]["U_235"]
        pu_239_frac = FUEL_MAKEUP.loc[core_type]["Pu_239"]
        u_238_frac  = FUEL_MAKEUP.loc[core_type]["U_238"]
        pu_241_frac = FUEL_MAKEUP.loc[core_type]["Pu_241"]

        # E release per fission in MJ
        # P is in MW Q is in MeV, so change Q to MJ
        E_fission = (u_235_frac*U_235_Q + u_238_frac*U_238_Q + pu_239_frac*PU_239_Q + pu_241_frac*PU_241_Q) * EV_J

        tot_spec_err = [
            (1./E_fission) * self._ferr_from_TFile(energy)
            for energy in ENERGIES
        ]

        spec_err_data = {
            "Total": tot_spec_err,
        }

        e_spec_err_tot = pd.DataFrame(spec_err_data, index=ENERGIES)

        return e_spec_err_tot, e_spec_err_tot


    """
    Return oscillation probability for given E at dist_to_sk
    """
    def p_ee(
        self,
        e,
        dm_21=DM_21,
        dm_23=DM_23,
        dm_31=DM_31,
        s_12=S_12,
        s_13=S_13_NH,
        c_12=C_12,
        c_13=C_13_NH,
        s_2_12=S_2_12,
        s_2_13=S_2_13
    ):
        l = self.dist_to_sk
        if e > IBD_MIN:
            # The terms from the propagator which will go in trigs
            # prop_31 = 1.267 * dm_31 * l * 1e3 / e
            # prop_32 = 1.267 * dm_32 * l * 1e3 / e
            # prop_21 = 1.267 * dm_21 * l * 1e3 / e
            p_31 = math.sin(1.267 * dm_31 * l * 1e3 / e) ** 2
            p_23 = math.sin(1.267 * dm_23 * l * 1e3 / e) ** 2
            p_21 = math.sin(1.267 * dm_21 * l * 1e3 / e) ** 2

            # p = 1 - 4 * s_12 * c_13 * c_13 * c_12 * (math.sin(prop_21)) ** 2
            # p -= 4 * s_13 * (math.sin(prop_31)) ** 2
            p = s_2_13 * (c_12 * p_31 + s_12 * p_23)
            p = 1 - s_2_12 * c_13 * c_13 * p_21 - p
            return max(p, 0)
        else:
            return 0


    """
    Calculating the oscillation pre factor that can be used for uncertainty propagation of ALL oscillated nu E at SK (flux [/cm^-2])
    """
    # TODO: Add in hierarchy support (I think it barely changes it)
    def osc_pre_factor(
        self,
        dm_21=DM_21,
        dm_31=DM_31,
        s_12=S_12,
        s_13=S_13_NH,
        c_12=C_12,
        c_13=C_13_NH,
        period="/s",
    ):
        core_type = self.core_type
        if self.mox:
            core_type = "MOX"

        # Fuel fractions for this type of core (p_i in literature)
        u_235_frac  = FUEL_MAKEUP.loc[core_type]["U_235"]
        pu_239_frac = FUEL_MAKEUP.loc[core_type]["Pu_239"]
        u_238_frac  = FUEL_MAKEUP.loc[core_type]["U_238"]
        pu_241_frac = FUEL_MAKEUP.loc[core_type]["Pu_241"]

        #l in km
        l = self.dist_to_sk

       # #E release per fission in MJ
        E_fission = (u_235_frac*U_235_Q + u_238_frac*U_238_Q + pu_239_frac*PU_239_Q + pu_241_frac*PU_241_Q) * EV_J

        # Osc spec per second
        if period == "/s":
            ps = [self.p_ee(e, dm_21, dm_31, s_12, s_13, c_12, c_13) for e in ENERGIES]
            ps = [p / (4 * math.pi * (l * 1e5) ** 2) / E_fission for p in ps]
        # Calculate based off of load factor
        else:
            # Finding total load factor
            year_start = int(period[:4])
            month_start = int(period[5:7])
            year_end = int(period[8:12])
            month_end = int(period[13:])

            # Cycle through all months summing load factor*t
            lf_sum = 0
            month_range_start = month_start
            month_range_end = 13
            n_nu_tot = 0
            # Make list of p_th for each month, find avg
            p_ths = []
            for year in range(year_start, year_end + 1):
                # Start from Jan after first year
                if year != year_start:
                    month_range_start = 1
                # Only go up to end of period in final year
                if year == year_end:
                    month_range_end = month_end + 1  # For inclusivity
                for month in range(month_range_start, month_range_end):
                    n_days_in_month = monthrange(year, month)[1]
                    # Query the specific month from the LF series
                    # print(self.lf_monthly)
                    lf_month = float(self.lf_monthly["%i/%02i" % (year, month)])
                    lf_month /= 100  # To be a factor, not %age
                    lf_sum += lf_month * n_days_in_month
                    p_ths.append(self.p_th[str(year)])

            avg_p_th = sum(p_ths) / len(p_ths)

            # lf_sum is sum of monthly load factors, so
            # p_th*lf_sum*(seconds in month) is integrated power
            # months had to do in sum cause months are stupid
            spec_pre_factor = avg_p_th * lf_sum * 24 * 60 * 60 

            if (
                # If the osc params are unchanged, don't recalculate
                math.isclose(dm_21, DM_21, rel_tol=1e-4)
                and math.isclose(c_13, C_13_NH, rel_tol=1e-4)
                and math.isclose(s_12, S_12, rel_tol=1e-4)
                and math.isclose(s_13, S_13_NH, rel_tol=1e-4)
            ):
                # Don't need to recalculate, just scale
                return self.def_osc_pre_factor * spec_pre_factor
            else:
                # From PHYSICAL REVIEW D 91, 065002 (2015)
                # E in MeV, l in km
                # Calculate the factor for the incoming spectrum to convert to flux in cm-2
                ps = [
                    self.p_ee(e, dm_21, dm_31, s_12, s_13, c_12, c_13) for e in ENERGIES
                ]
                ps = [p * spec_pre_factor / (4 * math.pi * (l * 1e5) ** 2) / E_fission for p in ps]

        osc_pre_factor = np.ones(E_BINS)
        osc_pre_factor = np.multiply(osc_pre_factor, ps)
        return osc_pre_factor
    

    """
    Calculating the spectrum of ALL oscillated nu E at SK (flux [/cm^-2])
    """
    # TODO: Add in hierarchy support (I think it barely changes it)
    def osc_spec(
        self,
        dm_21=DM_21,
        dm_31=DM_31,
        s_12=S_12,
        s_13=S_13_NH,
        c_12=C_12,
        c_13=C_13_NH,
        period="/s",
    ):

        l = self.dist_to_sk
        # Osc spec per second
        if period == "/s":
            ps = [self.p_ee(e, dm_21, dm_31, s_12, s_13, c_12, c_13) for e in ENERGIES]
            ps = [p / (4 * math.pi * (l * 1e5) ** 2) for p in ps]
        # Calculate based off of load factor
        else:
            # Finding total load factor
            year_start = int(period[:4])
            month_start= int(period[5:7])
            year_end   = int(period[8:12])
            month_end  = int(period[13:])

            # Cycle through all months summing load factor*t
            lf_sum = 0
            month_range_start = month_start
            month_range_end = 13
            n_nu_tot = 0
            # Make list of p_th for each month, find avg
            p_ths = []
            for year in range(year_start, year_end + 1):
                # Start from Jan after first year
                if year != year_start:
                    month_range_start = 1
                # Only go up to end of period in final year
                if year == year_end:
                    month_range_end = month_end + 1  # For inclusivity
                for month in range(month_range_start, month_range_end):
                    n_days_in_month = monthrange(year, month)[1]
                    # Query the specific month from the LF series
                    # print(self.lf_monthly)
                    lf_month = float(self.lf_monthly["%i/%02i" % (year, month)])
                    lf_month /= 100  # To be a factor, not %age
                    lf_sum += lf_month * n_days_in_month
                    p_ths.append(self.p_th[str(year)])

            avg_p_th = sum(p_ths) / len(p_ths)

            # lf_sum is sum of monthly load factors, so
            # p_th*lf_sum*(seconds in month) is integrated power
            # months had to do in sum cause months are stupid
            spec_pre_factor = avg_p_th * lf_sum * 24 * 60 * 60

            if (
                # If the osc params are unchanged, don't recalculate
                math.isclose(dm_21, DM_21, rel_tol=1e-4)
                and math.isclose(c_13, C_13_NH, rel_tol=1e-4)
                and math.isclose(s_12, S_12, rel_tol=1e-4)
                and math.isclose(s_13, S_13_NH, rel_tol=1e-4)
            ):
                # Don't need to recalculate, just scale
                return self.def_osc_spec * spec_pre_factor
            else:
                # From PHYSICAL REVIEW D 91, 065002 (2015)
                # E in MeV, l in km
                l = self.dist_to_sk

                # Calculate the factor for the incoming spectrum to convert to flux
                ps = [
                    self.p_ee(e, dm_21, dm_31, s_12, s_13, c_12, c_13) for e in ENERGIES
                ]
                ps = [p * spec_pre_factor / (4 * math.pi * (l * 1e5) ** 2) for p in ps]

        osc_spec = np.multiply(self.prod_spec["Total"], ps)

        return osc_spec

    """
    Spectrum of INTERACTED oscillated nu E at SK
    Takes oscillated spec as list and multiplies by xsec
    """
    def int_spec(self, osc_spec):
        # From PHYSICAL REVIEW D 91, 065002 (2015)
        int_spec = np.multiply(osc_spec, (SK_N_P * xsecs))
        return int_spec
    

    """
    Use Lat and Long info to calc flux directionality at SK
    Result is a vector with direction (x,y,z) and with norm osc_spec()
    Assume reactors are at sea level, very reasonable
    assumption, given most are on coastline
    """
    def _dir_flux_at_sk(self):
        lat_react_rad = radians(self.latitude)
        long_react_rad = radians(self.longitude)
        elev_react = self.elevation * 0.001     # elevation initially given in meter, needed in km

        lat_sk_rad = radians(SK_LAT)
        long_sk_rad = radians(SK_LONG)

        r_react = self._dist_to_earth_centre(lat_react_rad) + elev_react    # in km
        r_sk = self._dist_to_earth_centre(lat_sk_rad) + SK_ALT              # in km

        x_react = r_react * self._cos_geocentric(lat_react_rad) * cos(long_react_rad)
        y_react = r_react * self._cos_geocentric(lat_react_rad) * sin(long_react_rad)
        z_react = r_react * self._sin_geocentric(lat_react_rad)

        x_sk = r_sk * self._cos_geocentric(lat_sk_rad) * cos(long_sk_rad)
        y_sk = r_sk * self._cos_geocentric(lat_sk_rad) * sin(long_sk_rad)
        z_sk = r_sk * self._sin_geocentric(lat_sk_rad)

        dist = sqrt((x_react - x_sk) ** 2 + (y_react - y_sk) ** 2 + (z_react - z_sk) ** 2)

        x_dir = (x_react - x_sk) / dist
        y_dir = (y_react - y_sk) / dist
        z_dir = (z_react - z_sk) / dist
        

        # if self.mox:
        #     print(self.country, ' / ', self.name, ' / ', self.core_type, '- MOX / ', self.p_th[0], ' / ', dist, " / ", x_dir, ' / ', y_dir, ' / ', z_dir, ' / ', E_INTERVAL*sum(self.def_osc_spec) )
        # else:
        #     print(self.country, ' / ', self.name, ' / ', self.core_type,      ' / ', self.p_th[0], ' / ', dist, " / ", x_dir, ' / ', y_dir, ' / ', z_dir, ' / ', E_INTERVAL*sum(self.def_osc_spec) )
        
        # print(self.country, ' / ', self.name, ' / ', self.core_type,      ' / ', self.p_th[0], ' / ', r_sk, " / ", x_sk, ' / ', y_sk, ' / ', z_sk, ' / ', E_INTERVAL*sum(self.def_osc_spec) )

        return [x_dir, y_dir, z_dir]