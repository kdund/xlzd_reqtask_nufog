import numpy as np
import scipy.stats as sps
import yaml
from importlib.resources import files
import inference_interface as ii
import pickle as pkl
from os.path import isfile

from flamedisx.xlzd import XLZDvERSource, XLZDPb214Source, XLZDKr85Source, XLZDXe124Source, XLZDXe136Source
from flamedisx.xlzd import XLZDvNRSolarSource, XLZDvNROtherLNGSSource, XLZDvNROtherSURFSource, XLZDNeutronSource
from flamedisx.xlzd  import XLZDWIMPSource, XLZDEFTScalarO6Source
from flamedisx.xlzd  import XLZDALPGalacticDMSource, XLZDHiddenPhotonSource

from tqdm import tqdm
from multihist import Histdd
from itertools import product as iterproduct
from copy import deepcopy

default_version = "v1.0"

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

def cross_dict(nominal_parameters=dict(),**kwargs):
    """
    return a list of dicts where only one parameter at a time is varied
    """
    ret = []
    keys = kwargs.keys()
    
    prototype = {k:nominal_parameters[k] for k in keys}
    ret.append(prototype)

    for k in keys:
        for par in kwargs[k]:
            if prototype[k] != par:
                rret = deepcopy(prototype)
                rret[k] = par
                ret.append(rret)
    return ret


def get_parameters(mode, version=default_version):
    fname = files("detector_parameters").joinpath(f'parameters_{mode}_{version}.yaml')
    with open(fname,"r") as file:
        ret = yaml.safe_load(file)
    return ret

def get_template_parameters(mode, version=default_version):
    parameters = get_parameters(mode, version)["parameters"]
    parameters.update(get_parameters(mode, version)["background_rates"])
    ret_fix = dict()
    ret_iter = dict()
    nominal_parameters = dict()
    template_format_string = ""
    for k,i in sorted(parameters.items()):
        if "range" in i.keys():
            ret_iter[k] = i["range"]
            nominal_parameters[k] = i.get("nominal_value", i["value"])
            template_format_string += '_'+k+"_{"+k+":"+i.get("format",".2f")+"}"
        else:
            ret_fix[k] = i["value"]
            nominal_parameters[k] = i["value"]
    return ret_fix, ret_iter, nominal_parameters, template_format_string 

def save_dict_to_pickle(ret, file_name, signal_type="WIMP"):
    pdfs = [dict(), dict()]
    for k,h in ret.items():
        if signal_type in k:
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


def generate_template_set(mode, signal_type, parameters, analysis_parameters, n_samples = int(1e7), file_name = None, use_radius = False):
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
    if mode in ['LENR', 'HENR']:
        fd_sources["SolarER"] = XLZDvERSource(                                         **common_pass_parameters)
        fd_sources["Pb214"]   = XLZDPb214Source(activity_uBq_kg = parameters["Pb214"], **common_pass_parameters)
        fd_sources["Kr85"]    = XLZDKr85Source(activity_ppt = parameters["Kr85"],      **common_pass_parameters)
        fd_sources["Xe136"]   = XLZDXe136Source(                                       **common_pass_parameters)
        fd_sources["Xe124"]   = XLZDXe124Source(                                       **common_pass_parameters)
        fd_sources["CEvNS_solar"] = XLZDvNRSolarSource(                                **common_pass_parameters)
        fd_sources["CEvNS_other_LNGS"] =XLZDvNROtherLNGSSource(                        **common_pass_parameters)
        fd_sources["CEvNS_other_SURF"] =XLZDvNROtherSURFSource(                        **common_pass_parameters)
        fd_sources["neutrons"]= XLZDNeutronSource(                                     **common_pass_parameters)
    elif mode == 'LEER':
        fd_sources["SolarER"] = XLZDvERSource(                                         **common_pass_parameters)
        fd_sources["Pb214"]   = XLZDPb214Source(activity_uBq_kg = parameters["Pb214"], **common_pass_parameters)
        fd_sources["Kr85"]    = XLZDKr85Source(activity_ppt = parameters["Kr85"],      **common_pass_parameters)
        fd_sources["Xe136"]   = XLZDXe136Source(                                       **common_pass_parameters)
        fd_sources["Xe124"]   = XLZDXe124Source(                                       **common_pass_parameters)
    else:
        raise ValueError(f'Invalid mode {mode}.')

    if signal_type == 'WIMP':
        if mode == 'LENR':
            signal_source = XLZDWIMPSource
            signal_parameter = 'wimp_mass'
        elif mode == 'HENR':
            signal_source = XLZDEFTScalarO6Source
            signal_parameter = 'mass_GeV'
        else:
             raise ValueError(f'Invalid mode {mode}.')
    elif signal_type == 'ALP':
        signal_source = XLZDALPGalacticDMSource
        signal_parameter = 'mass_keV'
    elif signal_type == 'HP':
        signal_source = XLZDHiddenPhotonSource
        signal_parameter = 'mass_keV'
    else:
        raise ValueError(f'Invalid signal type {signal_type}.')
    
    masses =  analysis_parameters["mass"]["value"]
    if type(masses) != list:
        masses = [masses]
    for mass in masses:
        signal_dict = {signal_parameter: mass}
        fd_sources[f'{signal_type}{mass:.0f}']= signal_source(**signal_dict, **common_pass_parameters)
    if mode == 'LENR':
        fd_sources["WIMP"]= XLZDWIMPSource(wimp_mass = analysis_parameters["wimp_mass_benchmark"]["value"], **common_pass_parameters)

    if mode in ['LENR', 'HENR']:
        cs1_bins = np.linspace(analysis_parameters['cs1_range']['value'][0],
                              analysis_parameters['cs1_range']['value'][-1],
                              analysis_parameters['cs1_bins']["value"])
        log10cs2_bins = np.linspace(analysis_parameters['cs2_range']['value'][0],
                                    analysis_parameters['cs2_range']['value'][-1],
                                    analysis_parameters['cs2_bins']["value"])
        rsq_bins = np.linspace(analysis_parameters['rsq_range']['value'][0],
                               analysis_parameters['rsq_range']['value'][-1],
                               analysis_parameters['rsq_bins']["value"])
        if use_radius:
            hist_args = dict(
                    bins = (cs1_bins, log10cs2_bins, rsq_bins),
                    axis_names = ['cS1', 'log10_cS2', 'rsq'])
        else:
            hist_args = dict(
                    bins = (cs1_bins, log10cs2_bins),
                    axis_names = ['cS1', 'log10_cS2'])
    elif mode == 'LEER':
        recoE_bins = np.linspace(analysis_parameters['recoE_range']['value'][0],
                                 analysis_parameters['recoE_range']['value'][-1],
                                 analysis_parameters["recoE_bins"]['value'])
        hist_args = dict(
                bins = (recoE_bins,),
                axis_names = ['recoE'])
    else:
        raise ValueError(f'Invalid mode {mode}.')


    mus = dict()
    templates = dict()

    for k, source in tqdm(fd_sources.items()):
        data = source.simulate(n_samples)
        hist = Histdd(**hist_args)

        if mode in ['LENR', 'HENR']:
            cS1 = data['cs1'].values
            log10_cS2 = np.log10(data['cs2'].values)
            rsq = (data['r'].values)**2

            hist_2d_args = dict(
                    bins = (cs1_bins, log10cs2_bins),
                    axis_names = ['cS1', 'log10_cS2'])
            hist_2d = Histdd(**hist_2d_args)
            hist_2d.add(cS1, log10_cS2)

            if use_radius:
                hist_1d_args = dict(
                    bins = (rsq_bins,),
                    axis_names = ['rsq'])
                hist_1d = Histdd(**hist_1d_args)
                hist_1d.add(rsq)
                hist_1d = hist_1d / hist_1d.n

                hist.histogram = np.array([hist_2d.histogram * rsq_scaling for rsq_scaling in hist_1d.histogram])
                hist.histogram = np.transpose(hist.histogram, [1, 2, 0])
            else:
                hist = hist_2d

        elif mode == 'LEER':
            import matplotlib.pyplot as plt
            hist.add(data['ces_er_equivalent'].values)
            plt.figure()
            hist.plot()

        else:
            raise ValueError(f'Invalid mode {mode}.')

        mus[k] = source.estimate_mu(n_trials = n_samples)
        hist.histogram = mus[k] * (hist.histogram / hist.n)

        templates[k] = hist

    # Different normalisation procedures
    if 'neutrons' in templates:
        templates['neutrons'] = templates['neutrons'] / templates['neutrons'].n * mus['CEvNS_other_LNGS'] * parameters['neutron']
    if 'WIMP' in templates:
        templates['WIMP'] = templates['WIMP'] * (analysis_parameters["wimp_cross-section_benchmark"]["value"] / 1e-45)

    for sname, template in templates.items():
        print(f'{sname}: {template.n}')

    if file_name is not None:
        fname = file_name.format(**parameters)

        save_dict_to_pickle(templates, fname + mode + "_" + ".pkl", signal_type)
        save_dict_to_ii(templates, fname+".ii.h5")

    return templates

def generate_templates(
        mode = 'LENR',
        signal_type = 'WIMP',
        version = default_version,
        n_samples = int(1e7),
        file_name_pattern = "{version}{parameter_string}",
        scan_parameters = False, #use only nominal if False, if "product" get all variations, if "cross" vary only one parameter at a time
        use_radius = False,
        skip_generated = True,
        ):

    all_parameters = get_parameters(mode=mode, version=version)
    analysis_parameters = all_parameters["wimp_analysis_parameters"]

    ret_fix, ret_iter, nominal_parameters, template_format_string = get_template_parameters(mode=mode, version=version)

    parameters = deepcopy(nominal_parameters)

    if not scan_parameters:
        parameters_array = [parameters]
    elif scan_parameters == "product":
        parameters_array = [p for p in product_dict(**ret_iter)]
    elif scan_parameters == "cross":
        parameters_array = [p for p in cross_dict(nominal_parameters=nominal_parameters, **ret_iter)]
    else:
        raise NotImplementedError("scan_parameters must be False, 'cross' or 'product'!")

    n_pars = len(parameters_array)

    for pars in tqdm(parameters_array, desc="Generating several templates", total=n_pars):
        parameters.update(pars)
        parameter_string = template_format_string.format(**parameters)
        file_name = file_name_pattern.format(version=version, parameter_string=parameter_string)


        if isfile(file_name+".ii.h5") and skip_generated:
            pass
        else:
            generate_template_set(mode=mode, signal_type=signal_type,
                          parameters=parameters,
                          analysis_parameters = analysis_parameters,
                          n_samples = n_samples,
                          file_name = file_name,
                          use_radius = use_radius)
