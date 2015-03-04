#!/usr/bin/python

"""
Top-level interface to the Platypus tool-set
"""
from __future__ import division

import sys
import runner

###################################################################################################

def notImplemented(x):
    """
    Print an error message and raise an exception for commands that are not yet
    implemented.
    """
    print "\n\n%s command is not yet implemented\n\n" %(x)
    sys.exit(0)

###################################################################################################

if __name__ == "__main__":

    possCommands = {}
    possCommands["callVariants"] = runner.callVariants
    possCommands["continueCalling"] = runner.continueCalling

    if len(sys.argv) == 1 or sys.argv[1] not in possCommands.keys():
        print "\n\n"
        print "Invalid usage: use Platypus as follows:"
        print "\n"
        for k in possCommands.keys():
            print "python Platypus.py", k, "[Options]"
        print "\n"
        print "For a list of possible options for a specific command, type 'python Platypus.py Command -h'"
        print "\n"
        sys.exit(0)
    if sys.version_info[0] != 2 or sys.version_info[1] < 6:
        print "\n\n"
        print "Platypus works only with Python versions 2.6 and greater. Python 3.X is not yet supported."
        print "\n\n"
        sys.exit(0)
    else:
        command = sys.argv[1]
        possCommands[command](sys.argv[2:])
