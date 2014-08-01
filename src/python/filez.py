import gzip
import os.path
import subprocess

###################################################################################################

def open(filename, mode="rb", compresslevel=9):
    """
    Function that allows transparent usage of dictzip, gzip and ordinary files
    """
    if mode[0] == 'r':

        if os.path.splitext(filename)[1].lower() in (".gz",".dz"):
            return gzip.GzipFile( filename, mode, compresslevel )

        for ext in (".dz",".DZ",".gz",".GZ"):
            if os.path.exists( filename + ext ):
                return gzip.GzipFile( filename + ext, mode, compresslevel )

    return file(filename, mode)

###################################################################################################

def openurl(url, command = "scp -q %s /dev/stdout"):
    """
    Function that opens a pipe to a URL with a user-specified command.  It will also uncompress .gz/.bz2 files on the fly
    """

    # shortcuts -- to avoid having to quote commands...
    if command == "scp":
        command = "scp -q %s /dev/stdout"

    if command == "wget":
        command = "wget -q -O - %s"

    # use underscore as synonym to a space -- solves quoting headaches
    command = command.replace('_',' ')
    p1 = subprocess.Popen( (command % url).split(), stdout=subprocess.PIPE )

    if url.lower().endswith('.gz'):
        p2 = subprocess.Popen( ["zcat"], stdin=p1.stdout, stdout=subprocess.PIPE )
    elif url.lower().endswith('.bz2'):
        p2 = subprocess.Popen( ["bzcat"], stdin=p1.stdout, stdout=subprocess.PIPE )
    else:
        p2 = p1

    return p2.stdout

###################################################################################################
