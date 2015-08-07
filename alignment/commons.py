__author__ = 'peeyush'
import subprocess, os
import tarfile


def create_odir():
    odir = 'results'
    indir = ['alignedLane', 'cache', 'peaks']
    cdpath = os.getcwd()
    if not os.path.exists(os.path.join(cdpath, odir)): os.mkdir(os.path.join(cdpath, odir))
    for i in indir:
        print os.path.join(cdpath, odir, i)
        if not os.path.exists(os.path.join(cdpath, odir, i)):
            os.mkdir(os.path.join(cdpath, odir, i))
    print os.path.join(cdpath, odir)
    return os.path.join(cdpath, odir)


def create_sample_dir(path):
    if not os.path.exists(path):
        os.mkdir(path)




def read_file_blocked(filename, block_size = 1 * 1024 * 1024):
    """Read a (possibly compressed) file block by uncompressed block, yielding the blocks"""
    p = None
    zf = None
    if filename.endswith('.tar.bz2'):
        cmd = ['lbzip2 -d | tar -xO']
        stdin = open(filename,'rb')
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr = subprocess.PIPE, stdin=stdin, shell=True)
    elif filename.endswith('.tgz') or filename.endswith('.tar.gz'):
        if not os.path.exists(filename):
            raise IOError("[Errno 2] No such file or directory: '%s" % filename)
        cmd = ["tar",'-xOf', filename] #as cool as python is, 5x speedup by calling tar directly instead of using tarfile is *not* funny.
        p = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                             stderr = subprocess.PIPE,
                             stdin=subprocess.PIPE,
                            bufsize=0)
    elif filename.endswith('.zip'):
        import zipfile
        zf = zipfile.ZipFile(filename,'r')
        names = zf.namelist()
        if len(names) > 1:
            raise ValueError("This zipfile contains more than one compressed sequence file. That's unexpected, augment read_file_blocked")
        op = zf.open(names[0])
    elif filename.endswith('.tar'):
        #Note delegation instead of using the common control flow... since it generates multiple pipes
        for block in _read_file_blocked_tarfile(filename, block_size):
            yield block
        return #no need for all the rest...
    elif filename.endswith('.fastq.gz'):
        cmd = "tar -xvOf %s | gzip -cd" % (filename)
        print cmd
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=0, shell=True)
        try:
            while True:
                block = p.stdout.read(block_size)
                if not block:
                    break
                else:
                    yield block
                    del block
        finally:
            p.wait()
            p.stdout.close()
            p.stderr.close()
    #else:
        #op = exptools.common.open_file(filename, 'rb')
    if not p is None:
        input_pipe = p.stdout
        if not p.stdin is None:
            p.stdin.close()
    else:
        raise ValueError('File type does not match')
    try:
        while True:
            block = input_pipe.read(block_size)
            if not block:
                break
            else:
                yield block
                del block
    finally:
        input_pipe.close()
        if not p is None:
            p.wait()
            p.stderr.close()
        if not zf is None:
            zf.close()

def _read_file_blocked_tarfile(filename, block_size):
 #assume for now that it's a tar of multiple fastq.gz, and possibly a csv that's being ignored, and will explode on other files... Such data is for example produced by Braunschweig after their 1.8 Illumina upgrade
    try:
        tf = tarfile.open(filename)
        for tarinfo in tf.getmembers():
            print len(tarfile)
            if tarinfo.size == 0: #ignore empty entries (directories)
                continue
            if tarinfo.name.endswith('.csv'):
                continue
            if tarinfo.name.endswith('.fastq.gz'):
                cmd = "tar -xvOf %s %s | gzip -cd" % (filename, tarinfo.name)
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=0, shell=True)
                try:
                    while True:
                        block = p.stdout.read(block_size)
                        if not block:
                            break
                        else:
                            yield block
                            del block
                finally:
                    p.wait()
                    p.stdout.close()
                    p.stderr.close()
            else:
                raise ValueError("chipseq.common.read_file_blocked does not know how to handle this file (%s) in tar archive %s" % (tarinfo.name, filename))
    finally:
        tf.close()