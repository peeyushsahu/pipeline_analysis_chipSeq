import subprocess as sp
import os
__author__ = 'peeyush'


class MACS:
    '''
    Peak calling for model-based analysis of chip seq data
    '''
    def __init__(self):
        self.macs_dir = "/home/sahu/Documents/MACS-1.3.7.1/bin"
        self.macs_cmd = os.path.join(self.macs_dir, 'macs')

    def get_version(self):
        if not hasattr(self, 'version'):
            cmd = ['python', self.macs_cmd, '--version']
            p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE)
            stdout, stderr = p.communicate()
            self.version = stdout
        print('MACS version:', self.version)