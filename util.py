# Richard Darst, September 2008

import sys

from saiga12.io import io_open

if __name__ == "__main__":
    if sys.argv[1] in ('--help', '-h', 'help'):
        print """usage: util.py <command> <arguments>...

        commands are:
        hash -- print command line hashes of all filename arguments
        """
    
    elif sys.argv[1] == 'hash':
        for name in sys.argv[2:]:
            S = io_open(name)
            print name, S.hash()

