#/usr/bin/python

# MLTreeMap Imagemaker TSN v. 0.0

# Include the necessary libraries
try:

except:


def os_type():
    """Return the operating system of the user."""
    x = sys.platform
    if x:
        hits = re.search(r'darwin', x, re.I)
        if hits:
            return 'mac'
        hits = re.search(r'win', x, re.I)
        if hits:
            return 'win'
        hits = re.search(r'linux', x, re.I)
        if hits:
            return 'linux'


def pathDelim():
    """Return the path deliminator based on the user's operating system."""
    ostype = os_type()
    if ostype == 'win':
        return '\\'
    if ostype in ['linux', 'mac']:
        return '/'


PATHDELIM = str(pathDelim())


def getParser():
    """TK"""
    parser = argparse.ArgumentParser(description='example input:\n\
                      python tsn_imagemaker.py -i example_input/\n\n\
                      mandatory input option:\n\
                      -i: either a) a path to a directory with input files \
                      or b) to a specific input file.\n\
                      note that in case a) all input files of the same type \
                      (i.e. p, g, n, h, r, m, d) will be concatenated.\n\n\
                      optional input options:\n\
                      -o: output directory (default: output/)\n\
                      -d: use different colors for different datasets \
                      (0: don\'t use this mode (default),\
                       1: use this mode,\
                       11: use this mode but concatenate percent values \
                      over datasets). \
                      See the user guide for a more detailed description.\n\
                      -t: display trees inclusive text and labels \
                      (1, default) \
                      or without text (0)\n\
                      -b: parameter for the size of the placement bubbles \
                      (default: 1)\n\
                      -r: display 16s and 18s rRNA hits in different trees \
                      (2, default) \
                      or one tree of life (1)\n\
                      -h: high quality placement bubbles (1, default) or \
                      low quality bubbles (0). \
                      High quality bubbles work perfectly with inkscape \
                      but can cause trouble with Adobe Illustrator.')
    parser.add_argument('-i', '--input', required=True, help='either \
                        a) a path to a directory with input files or \
                        b) to a specific input file')
    parser.add_argument('-o', '--output', default='output', help='output \
                        directory (default: output/)')
    parser.add_argument('-d', '--colors', default=0, choices=[0,1,11], \
                        help='use different colors for different datasets \
                        (0: don\'t use this mode (default),\
                         1: use this mode,\
                         11: use this mode but concatenate percent values \
                             over datasets)')
    parser.add_argument('-t', '--text', default=1, choices=[0,1], \
                        help='display trees inclusive text and labels \
                        (1, default) or without text(0)')
    parser.add_argument('-b', '--bubble_size', default=1, type=int, \
                        help='parameter for the size of the placement \
                        bubbles (default: 1)')
    parser.add_argument('-r', '--rRNA', default=2, choices=[1,2], \
                        help='display 16s and 18s rRNA hits in different \
                        trees (2, default) or one tree of life (1)')
    parser.add_argument('-h', '--hq_bubbles', default=1, choices=[0,1], \
                        help='high quality bubbles work perfectly with \
                        inkscape but can cause trouble with Adobe \
                        Illustrator')
    parser.add_argument('-i', '--input', required=True, default=0, type=int, choices=[0, 's'], help='')
   return parser

def read_user_input(parser):
    """TK"""
    args = parser.parse_args()

    # Remove any pre-existing path delim on the output directory
    while re.search(r'/\Z', args.output) or re.search(r'\\\Z', args.output):
        args.output = args.output[:-1]

    # Add the path delim to the output directory
    args.output += PATHDELIM

    # Check to see if output directory already exists; 
    # if not, create output directory
    while os.path.isdir(args.output):
        print 'WARNING: Your output directory ' + str(args.output) + \
              'already exists! Overwrite [1], quit [2], or \
              change directory [3]?'
        answer = int(raw_input())

        while not answer in [1,2,3]:
            answer = int(raw_input('Invalid input. Please choose \
                     1, 2, or 3.\n'))

        if answer == 1:
            print 'Do you really want to overwrite the old output directory?'
            print 'All data in it will be lost!'
            answer2 = raw_input('Yes [y] or no [n]?\n')

            while not answer2 in ['y', 'n']:
                answer2 = raw_input('Invalid input. Please choose y or n.\n')

            if answer2 == 'y':
                shutil.rmtree(args.output)
            else:
                sys.exit('Exit MLTreeMap Imagemaker\n')
        elif answer == 2:
            sys.exit('Exit MLTreeMap Imagemaker\n')
        else:
            args.output = raw_input('Please enter the path to the new \
                          directory.\n')

            while re.search(r'/\Z', args.output) \
                  or re.search(r'\\\Z', args.output):
                args.output = args.output[:-1]

            args.output += PATHDELIM

    os.makedirs(args.output)
    return args

def main(argv):
    print 'Start MLTreeMap Imagemaker TSN v. 0.0'

if __name__ == "__main__":
    main(sys.argv[1:])
