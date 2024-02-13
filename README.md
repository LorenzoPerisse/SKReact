# SKReact v2.0
A reactor neutrino simulator for Super-Kamiokande, simulating production in
reactors, oscillation to SK, cross section for IBD on free protons, and detector
smearing. Built referencing 
[An *et al.* (2022)](https://arxiv.org/abs/2203.06686) (for reactor data-driven model),
[Perisse *et al.* (2023)](https://arxiv.org/abs/2304.14992v2) (for summation model),
[M. Baldoncini *et al.*, *"Reference worldwide model for antineutrinos from reactors"*](https://arxiv.org/abs/1411.6475) (for reactor survey).
Note that the IBD cross section implemented in SKReact is using a simplified formula corresponding to the Naive+ one from [Strumia and Vissani (2003)](https://arxiv.org/abs/astro-ph/0302055).



## List of main updates

### v1.3
* Single reactor spectrum modeling: conversion method based on [Mueller *et al.* (2011)](https://arxiv.org/abs/1101.2663).
* Oscillation parameters: NuFit, 3-flavor fit, July 2019.

### v2.0
* Single reactor spectrum modeling: reactor data-driven method based on [An *et al.* (2022)](https://arxiv.org/abs/2203.06686) + summation method based on [Perisse *et al.* (2023)](https://arxiv.org/abs/2304.14992v2).
* ROOT file used as input and output files.
* Oscillation parameters: NuFit, 3-flavor fit, November 2022.
* Oscillation prefactor computation, used for covariance matrix computation.
* Covariance matrix computation for reactor specturm at SK.
* Reactor elevation added to reactor data (taken from IAEA data).
* Use of Daya Bay fission fraction from 2017 (previously using value from M. Baldoncini *et al.*).
* 3 default data sheets `DB2020.xlsx`,`DB2021.xlsx` and `DB2022.xlsx` provided in `react_p/`. They are associated to the full SK-VI data taking period with the same average load factor for each month (and null load factors outside of the SK-VI period).
* Spectrum binning is 1keV and range up to 16 MeV.
* Constant updated with PDG 2022 data.



## Installation and Running
SKReact uses python 3(.7.0), it does not support python 2.
Install the required modules using pip: 

`$ pip install -r requirements.txt`.

By default SKReact uses `.xls/.xlsx` data sheets (FORCE_XLS_IMPORT = True in params.py) summarizing ractor activities.
Data are taken from the [IAEA-PRIS database](https://pris.iaea.org/PRIS/home.aspx).
If are not a IAEA-PRIS member, you can stil get PRIS data kindly compiled by [INFN](https://www.fe.infn.it/radioactivity/antineutrino/index.html#download) or through [https://reactors.geoneutrinos.org/](https://reactors.geoneutrinos.org/).

Place the `.xls/xlsx` files into `react_p/` in the SKReact directory, with their names unchanged (i.e. they must be labeled DBXXXX where XXXX is the associated year).
Simply create the directory `react_p/` if it does not already exist, or change the
`REACT_DIR` variable in `params.py` if you want to pull the data from somewhere else.
You will compute and saved into the ROOT file by default the spectrum associated to the summed period of all `.xls/.xlsx` files located in `react_p/`.

You can also use a pickle file such as `reactors_main.pkl`.
If you don't have a `.pkl` file, you can generate your own based on the `.xls/.xlsx` files and by symply running SKReact with FORCE_XLS_IMPORT = False.
The first time running without a `reactors_main.pkl` file will be very slow to
start as the relevant information needs to be pulled and compiled from the
`.xls/.xlsx` files. SKReact tries to handle all inconcistencies across the files and
will print any changes in reactor information as it compiles, if you want
verbose import/import errors, change the `VERBOSE_IMPORT` and
`VERBOSE_IMPORT_ERR` values in `params.py`.
Running from then on with a `reactors_main.pkl` file will be much faster on
startup. There is a chance you will need to recompile from `.xls/.xlsx` files between
SKReact releases.

To run the program, go to `SKReact/` and simply run `python3 skreact.py`.

You will need a modeling of a reactor antineutrino spectrum for a single reactor (typically a PWR spectrum).
If you are a member of Super-Kamiokande, you can find this file on the collaboration cluster: /home/lperisse/Work/Reactor/model/ReactorModel_SK.root
Otherwise you can send me a request at [lorenzo.perisse@cnrs.fr](lorenzo.perisse@cnrs.fr).
If you want to use your own reactor prediction, make sure the values in params.py are consistent between your reactor prediction and E_MIN, E_MAX, E_BINS and E_INTERVAL.


## Features

### Reactor List
![reactor_list](../assets/reactor_list.png?raw=true)

This is a list of the reactors contributing to the total spectrum at SK i.e. the
reactors imported from the PRIS. It is ordered by largest *flux* contribution at
Super-K. Select a reactor to highlight it on the spectra plots. The last
selected reactor's information will appear in the reactor info pane.

### Reactor Info
![reactor_info](../assets/reactor_info.png?raw=true)

Here you can see the relevant information about the last selected reactor. You
can edit the information as you like then click the *Update* button to
re-generate the spectra.

To edit a load factor or power information, select a value in the listbox, type
the desired value in the entry box below and hit return.

### Period Selection/Load Factor Plot
![lf](../assets/lf.png?raw=true)

This plot lets you visualise the expected reactor signal/reactor power
information, by default showing the number of interactions each month (for
the total (blue) and highlighted reactors). All values dynamically update
with reactor info changes **except** the number of interactions, which is
calculated on startup. You can select the period to focus on using the drop
boxes. NOTE: the selection is **inclusive**, i.e. the end month selected is
included in the simulation, so to inspect one month, you set both start and
end to the same month (will likely print errors with plotting to terminal, will
be fixed in future update).

### Spectrogram
![spectro](../assets/spectro.png?raw=true)

A spectrogram (calculated on startup) to visualise how the spectrum changes over
time. Currently "for show" but will make more use of this in the future.

### Produced Spectrum
![prod_spec](../assets/prod_spec.png?raw=true)

Shows the produced spectrum (/s) and n produces for the last selected reactor
(varies with reactor type). Click the checkboxes to show individual fuel
contributions (assumes fixed fuel fractions).

### Oscillated Spectrum
![osc_spec](../assets/osc_spec.png?raw=true)

The total oscillated spectrum at SK with some numbers. Can save the plot
(with its current size, will update to produce full sized plot in future), as
well as produce a `.csv` of the spectrum by selecting `.csv` as the extension
when saving. Can vary the oscillation parameters using the sliders:

![osc_slider](../assets/osc_slider.png?raw=true)

This is slow (the spectra for each reactor needs to be calculated and summed on
the fly), but is nice for visualising the main parameters' effects on the
spectra.

The oscillated spectrum is also saved into a ROOT file defined in `params.py` at l.337.
Note that this action is performed only once and for the full period considered when SKReact is launched.
At the same time, an energy dependant factor called *oscillation prefactor* is saved in the ROOT file and is used to compute the covariance matrix in a latter step.

### Interacted Spectrum
![int_spec](../assets/int_spec.png?raw=true)

The final spectrum after multiplying the oscillated spectrum by the IBD cross
section (Naive+ formula from [Strumia and Vissani (2003)](https://arxiv.org/abs/astro-ph/0302055)). 
You can view the neutrino or positron spectrum (i.e. offset by 1.806 MeV) by selecting the radio button for each in the options.

![int_spec_options](../assets/int_spec_options.png?raw=true)

From this you can save the plot (also as `.csv` as with the oscillated
spectrum) or generate a `.nuance` file from the spectrum for input into detector
simulators (will always generate positron spectrum). 

### Smearing and Fitting

On startup, if `smear_main.csv` is found (change `WIT_SMEAR_FILE` in params.py
for other filenames) SKReact will generate a smearing matrix to multiply by the
interacted spectrum, giving a "detected" spectrum.

`smear_main.csv` contains the gaussian parameters for the reconstructed energy
from fixed energy in the detector and must be of the form (with headers):

`e,c,mu,sig,eff`

Where `e` is the energy reconstructed, `c,mu,sig` are the gaussian parameters
and `eff` is the efficiency of detection at that energy. You can generate this
yourself using your detector simulator by generating samples at fixed energies,
reconstructing them and fitting a gaussian. If you are an SK collaborator,
contact [lorenzo.perisse@cnrs.fr](lorenzo.perisse@cnrs.fr) or
[alexander.goldsack@physics.ox.ac.uk](alexander.goldsack@physics.ox.ac.uk) for a
WIT smearing file.

SKReact also has a (crude) fitter, which will vary selected oscillation
parameters to fit a given detected spectrum with the calculated smeared
spectrum. Click "Import data to fit" on the interacted spectrum options frame
and select a `.csv` file of format `e,bin content` (without headers). Then
select which parameters you'd like to fit, and run. The fitter cyclicly finds a
minima in a range defined under `# FITTING` in `params.py`.

![fit_win](../assets/fit_win.png?raw=true)

NOTE: The fitter in SKReact is very crude and requires some tuning and
optimisation/parallelisation.



## Covariance matrix 

A covariance matrix for the oscillated spectrum at SK can be computed by using the C++ macro `make_final_spec.C` located in the `macros/` directory.
This macro needs 3 SKReact output files to compute the covariance matrix associated to the oscillation parameters: one output based on the parameter central values, one output for the parameters at +1 standard deviation and one for the values at -1 standard  deviation.
These 3 files are obtained by running 3 times SKReact, each time with different oscillation parameters that can be commented/uncommented in params.py at lines 230 to 258.
Rename each SKReact output file as ReactorModel_SKReact_v2p0.root, ReactorModel_SKReact_v2p0_oscillationUp.root and ReactorModel_SKReact_v2p0_oscillationDown.root (or modify the names at lines 338, 339 and 340 of the macro).
You can then run the macro by launching ROOT in a terminal and by compiling the macro directly in ROOT (you might want to modify the variable path_SKReact at l.332).

`
.L make_final_spec.C
macro()
`



## Contact

For any issues, suggestions, enhancements or bugs, please post an issue on
the repository's issues page. To contact me directly, please email
[lorenzo.perisse@cnrs.fr](lorenzo.perisse@cnrs.fr) or
[alexander.goldsack@physics.ox.ac.uk](alexander.goldsack@physics.ox.ac.uk),
however please be aware my support for SKReact is primarily limited to my own
requirements and so I will have to weigh up the time/personal usefulness for any
modifications.

Please enjoy SKReact!
