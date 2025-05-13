import numpy as np
import scipy.stats as sps
import yaml
from importlib_resources import files

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
    ret_fix = dict()
    ret_iter = dict()
    nominal_parameters = dict()
    template_format_string = ""
    for k,i in sorted(parameters.items()):
        if type(i["value"]) == list:
            ret_iter[k] = i["value"]
            nominal_parameters[k] = i.get("nominal_value", i["value"][0])
            template_format_string += k+"_{"+k+":"+i.get("format",".2f")+"}"
        else:
            ret_fix[k] = i["value"]
            nominal_parameters[k] = i["value"]
    return ret_fix, ret_iter, nominal_parameters, template_format_string
