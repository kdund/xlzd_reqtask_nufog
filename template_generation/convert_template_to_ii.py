import multihist as mh
import numpy as np
import pickle as pk
import inference_interface as ii


def convert_file(fname):
    file = open(fname, "rb")
    objects = pk.load(file)

    hs = []
    hns = []

    for i in range(len(objects)):
        for k,h in objects[i].items():
            hs.append(h)
            hns.append(k)
    ii.multihist_to_template(hs, fname.replace(".pkl",".h5"), hns)

if __name__ == "__main__":
    convert_file("../data/pdfs_CEvNS_disco_60t_r0.2_good.pkl")
    convert_file("../data/pdfs_CEvNS_disco_60t_r0.2_bad.pkl")
