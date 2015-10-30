__author__ = 'peeyush'
import subprocess, os
import tarfile


def create_odir():
    odir = 'results'
    indir = ['alignedLane', 'cache', 'peaks']
    cdpath = '/ps/imt/e/HL60_Christene/further_analysis'##os.getcwd() when run script from the folder of interest
    if not os.path.exists(os.path.join(cdpath, odir)): os.mkdir(os.path.join(cdpath, odir))
    for i in indir:
        print os.path.join(cdpath, odir, i)
        if not os.path.exists(os.path.join(cdpath, odir, i)):
            os.mkdir(os.path.join(cdpath, odir, i))
    print os.path.join(cdpath, odir)
    return os.path.join(cdpath, odir)


def ensure_path(path):
    if not os.path.exists(path):
        os.mkdir(path)
