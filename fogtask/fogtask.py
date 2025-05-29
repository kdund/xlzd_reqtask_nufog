import numpy as np
import scipy.stats as sps
import yaml
from importlib_resources import files
import inference_interface as ii
import pickle as pkl
from flamedisx.xlzd import XLZDvERSource, XLZDPb214Source, XLZDKr85Source, XLZDXe124Source, XLZDXe136Source, XLZDvNRSolarSource, XLZDvNROtherLNGSSource, XLZDWIMPSource, XLZDNeutronSource
from tqdm import tqdm
from multihist import Histdd
from itertools import product as iterproduct
from copy import deepcopy

default_version = "v0.5"

def test_util():
    print("I'm helping!")

def product_dict(**kwargs):
    """
        returns product of the list of each element of the dictionary as a dictionary
        from https://stackoverflow.com/users/118160/seth-johnson on 
        https://stackoverflow.com/questions/5228158/cartesian-product-of-a-dictionary-of-lists
    """
    keys = kwargs.keys()
    for instance in iterproduct(*kwargs.values()):
        yield dict(zip(keys, instance))

def get_parameters(version=default_version):
    fname = files("detector_parameters").joinpath("parameters_{version:s}.yaml".format(version=version))
    with open(fname,"r") as file:
        ret = yaml.safe_load(file)
    return ret

def get_detector_parameters(version=default_version):
    return get_parameters(version)["parameters"]

def get_template_parameters(version=default_version):
    parameters = get_parameters(version)["parameters"]
    parameters.update(get_parameters(version)["background_rates"])
    ret_fix = dict()
    ret_iter = dict()
    nominal_parameters = dict()
    template_format_string = ""
    for k,i in sorted(parameters.items()):
        if "range" in i.keys():
            ret_iter[k] = i["range"]
            nominal_parameters[k] = i.get("nominal_value", i["value"])
            template_format_string += k+"_{"+k+":"+i.get("format",".2f")+"}"
        else:
            ret_fix[k] = i["value"]
            nominal_parameters[k] = i["value"]
    return ret_fix, ret_iter, nominal_parameters, template_format_string 

def save_dict_to_pickle(ret, file_name, signal_pattern = "WIMP"):
    pdfs = [dict(), dict()]
    for k,h in ret.items():
        if signal_pattern in k:
            pdfs[1][k] = h
        else:
            pdfs[0][k] = h
    pkl.dump(pdfs, open(file_name, "wb"))

def save_dict_to_ii(ret, file_name):
    histogram_names = sorted(ret.keys())
    ii.multihist_to_template(
        [ret[k] for k in histogram_names], 
        file_name, 
        histogram_names = histogram_names, 
        )


def generate_template_set(parameters, analysis_parameters, n_samples = int(1e7), file_name = None, use_radius = False):
    """
    For a fixed set of parameters, generate all templates used in the XLZD flamefit template generation. If file_name is
    not None, store the templates both as pickle and inference_interface files
    """
    common_pass_parameters = dict(
    cS1_min = analysis_parameters['cs1_range']['value'][0],
    cS1_max = analysis_parameters['cs1_range']['value'][-1],
    log10_cS2_min = analysis_parameters['cs2_range']['value'][0],
    log10_cS2_max = analysis_parameters['cs2_range']['value'][-1],
    s2_thr = analysis_parameters['s2_threshold']['value'],
    coin_level = analysis_parameters['coincidence_threshold']['value'],
    drift_field_V_cm = parameters['drift_field'],
    gas_field_kV_cm = parameters['gas_field'],
    elife_ns = parameters['electron_livetime'] * 1e6,
    g1 = parameters['PMT_quantum_efficiency'],
    temperature_K = parameters['temperature'],
    pressure_bar = parameters['pressure'],
    num_pmts = parameters['n_pmts'],
    double_pe_fraction = parameters['p_dpe'],
    g1_gas = parameters['g1_gas'],
    s2Fano = parameters['s2_fano'],
    spe_res = parameters['spe_resolution'],
    spe_thr = parameters['spe_threshold'],
    spe_eff = parameters['spe_efficiency'],
    configuration = parameters['lce_configuration'],
    )

    fd_sources = dict()
    fd_sources["SolarER"] = XLZDvERSource(                                         **common_pass_parameters)
    fd_sources["Pb214"]   = XLZDPb214Source(activity_mBq_kg = parameters["Pb214"], **common_pass_parameters)
    fd_sources["Kr85"]    = XLZDKr85Source(activity_ppt = parameters["Kr85"],      **common_pass_parameters)
    fd_sources["Xe136"]   = XLZDXe136Source(                                       **common_pass_parameters)
    fd_sources["Xe124"]   = XLZDXe124Source(                                       **common_pass_parameters)
    fd_sources["CEvNS_solar"] = XLZDvNRSolarSource(                                **common_pass_parameters)
    fd_sources["CEvNS_other_LNGS"] =XLZDvNROtherLNGSSource(                        **common_pass_parameters)
    fd_sources["neutrons"]= XLZDNeutronSource(                                     **common_pass_parameters)

    wimp_masses =  analysis_parameters["wimp_mass"]["value"]
    if type(wimp_masses) != list: 
        wimp_masses = [wimp_masses]
    for wimp_mass in wimp_masses:
        fd_sources["WIMP_{:.1f}".format(wimp_mass)]= XLZDWIMPSource(wimp_mass = wimp_mass, **common_pass_parameters)

    cs1_bins = np.linspace(analysis_parameters['cs1_range']['value'][0],
                           analysis_parameters['cs1_range']['value'][-1],
                           analysis_parameters["cs1_bins"]["value"])
    log10cs2_bins = np.linspace(analysis_parameters['cs2_range']['value'][0],
                                analysis_parameters['cs2_range']['value'][-1],
                                analysis_parameters["cs2_bins"]["value"])
    rsq_bins = np.linspace(analysis_parameters['rsq_range']['value'][0],
                           analysis_parameters['rsq_range']['value'][-1],
                           analysis_parameters["rsq_bins"]["value"])

    if use_radius:
        hist_args = dict(
                bins = (cs1_bins, log10cs2_bins, rsq_bins),
                axis_names = ['cS1', 'log10_cS2', 'rsq'])
    else:
        hist_args = dict(
                bins = (cs1_bins, log10cs2_bins),
                axis_names = ['cS1', 'log10_cS2'])

    mus = dict()
    templates = dict()
    for k, source in tqdm(fd_sources.items()):
        hist = Histdd(**hist_args)

        data = source.simulate(n_samples)
        cS1 = data['cs1'].values
        log10_cS2 = np.log10(data['cs2'].values)
        rsq = (data['r'].values)**2

        if use_radius:
            hist.add(cS1, log10_cS2, rsq)
        else:
            hist.add(cS1, log10_cS2)

        mus[k] = source.estimate_mu(n_trials = n_samples)
        hist.histogram = mus[k] * (hist.histogram / hist.n)

        templates[k] = hist

    # Different normalisation procedure
    templates['neutrons'] = templates['neutrons'] / templates['neutrons'].n * mus['CEvNS_other_LNGS'] * parameters['neutron']

    for sname, template in templates.items():
        print(f'{sname}: {template.n}')

    if file_name is not None:
        fname = file_name.format(**parameters)

        save_dict_to_pickle(templates, fname+".pkl", "WIMP")
        save_dict_to_ii(templates, fname+".ii.h5")

    return templates

def generate_all_wimp_templates(
        version = default_version,
        n_samples = int(1e7),
        file_name_pattern = "{parameter_string:s}",
        nominal_only = True,
        use_radius = False,
        ):

    all_parameters = get_parameters(version=version)
    analysis_parameters = all_parameters["wimp_analysis_parameters"]

    ret_fix, ret_iter, nominal_parameters, template_format_string = get_template_parameters(version=version)

    parameters = deepcopy(nominal_parameters)

    if nominal_only:
        parameter_string = template_format_string.format(**parameters)
        file_name = file_name_pattern.format(parameter_string=parameter_string)

        generate_template_set(parameters=parameters,
                              analysis_parameters = analysis_parameters,
                              n_samples = n_samples,
                              file_name = file_name,
                              use_radius = use_radius)
    else:
        for pars in product_dict(**ret_iter):
            parameters.update(pars)
            parameter_string = template_format_string.format(**parameters)
            file_name = file_name_pattern.format(parameter_string=parameter_string)

            generate_template_set(parameters=parameters,
                                  analysis_parameters = analysis_parameters,
                                  n_samples = n_samples,
                                  file_name = file_name,
                                  use_radius = use_radius)
