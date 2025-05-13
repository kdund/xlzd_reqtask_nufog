import numpy as np
import matplotlib.pyplot as plt
import inference_interface as ii
from tqdm import tqdm

import pickle as pkl
from multihist import Histdd

import flamedisx as fd

import os



# Script based on Robert James' generation notebook 

def generate_wimp_template(fname = "deleteme.h5",
                           detector_performance = "good",
                           N_hist = int(1e8)):
    if detector_performance == "good":
        cs1_bins = np.linspace(5,100,101)
        cs2_bins = np.linspace(3.1, 4.1, 101)
        detector_configuration = dict(
        configuration_name = "60t",
        drift_field_v_cm = 80,
        gas_field_kV_cm = 7.5,
        elife_ns = 10e6,
        g1 = 0.27,
        )
    elif detector_performance =="bad":
        cs1_bins = np.linspace(5,100,101)
        cs2_bins = np.linspace(2.8, 3.7, 101)
        detector_configuration = dict(
        configuration_name = "60t",
        drift_field_v_cm = 25,
        gas_field_kV_cm = 6,
        elife_ns = 10e6,
        g1 = 0.27,
        )
    else:
        raise NotImplementedError("detector_performance must be good or bad")

    sources = dict()
    source_WIMP = fd.xlzd.XLZDWIMPSource
    wimp_masses = [11,40., 5000]
    sources = {k:i(**detector_configuration) for k,i in sources.items()}
    sources_wimp = {"WIMP{:.0f}".format(wimp_mass):source_WIMP(wimp_mass=wimp_mass,**detector_configuration) for wimp_mass in wimp_masses }
    sources.update(**sources_wimp)
    binning = dict(
        bins = [cs1_bins, cs2_bins],
        axis_names = ["cS1", "log10_cS2"]
            )

    N_hist_batch = 100
    N_mu = int(1e6)

    pdfs = dict()

    for n,source in sources.items():
        print("filling ",n)
        hist = Histdd(**binning)
        for i in tqdm(range(N_hist_batch)):
            data = source.simulate(int(N_hist/N_hist_batch))
            hist.add(data["cs1"], np.log10(data["cs2"]))
        mu = source.estimate_mu(n_trials=N_mu)
        hist.histogram *= mu / hist.n
        #plt.clf()
        #hist.plot()
        #plt.title(n)
        #plt.show()
        pdfs[n] = hist

    histogram_names = sorted(pdfs.keys())
    ii.multihist_to_template(
        [pdfs[k] for k in histogram_names], 
        fname, 
        histogram_names = histogram_names, 
        )


if __name__ == "__main__":

    #r0.2 WIMP templates: 
    #generate_wimp_template(detector_performance = "good", N_hist = int(1e8), fname = "../data/pdfs_wimp_60t_r0.2_good.h5")
    #generate_wimp_template(detector_performance = "bad", N_hist = int(1e8), fname = "../data/pdfs_wimp_60t_r0.2_bad.h5")
    generate_wimp_template(detector_performance = "good", N_hist = int(1e8), fname = "../data/pdfs_wimp_60t_r0.3_good.h5")
    generate_wimp_template(detector_performance = "bad", N_hist = int(1e8), fname = "../data/pdfs_wimp_60t_r0.3_bad.h5")






