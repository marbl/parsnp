#!/usr/bin/env python
# See the LICENSE file included with this software for license information.

import os, sys, string, getopt, random,subprocess, time, glob,operator, math, datetime,numpy #pysam
from collections import defaultdict
import csv
import tempfile
import shutil
import re
import logging
import argparse
import signal
import inspect
from multiprocessing import *

__version__ = "1.2"
reroot_tree = True #use --midpoint-reroot

try:
    import dendropy
except ImportError:
    reroot_tree = False

#check for sane file names
ALIGNER_TO_IDX = {
        "mafft": "1",
        "muscle": "2",
        "fsa": "3",
        "prank": "4"
}

VERBOSE = 0
VERSION = "v1.5"
PHI_WINDOWSIZE = 1000
TOTSEQS=0
PARSNP_DIR = sys.path[0]


############################################# Logging ##############################################
BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE = range(8)
#These are the sequences need to get colored ouput
CSI=""#"\x1B["
BOLDME = ""#CSI+'\033[1m'
STATUS_BLUE = ""#CSI+'\033[94m'
OK_GREEN = ""#CSI+'\033[92m'#'32m'
SKIP_GRAY = ""#CSI+'\033[37m'
WARNING_YELLOW = ""#CSI+'\033[93m'
ERROR_RED = ""#CSI+'\033[91m'
ENDC = ""#CSI+'0m'
RESET_SEQ = "\033[0m"
COLOR_SEQ = "\033[1;%dm"
BOLD_SEQ = "\033[1m"

def formatter_message(message, use_color = True):
    if use_color:
        message = message.replace("$RESET", RESET_SEQ).replace("$BOLD", BOLD_SEQ)
    else:
        message = message.replace("$RESET", "").replace("$BOLD", "")
    return message

COLORS = {
    'DEBUG': BLUE,
    'INFO': WHITE,
    'WARNING': YELLOW,
    'ERROR': RED,
    'CRITICAL': RED
}

class ColoredFormatter(logging.Formatter):
    def __init__(self, msg, datefmt = None, use_color = True):
        logging.Formatter.__init__(self, fmt=msg, datefmt=datefmt)
        self.use_color = use_color

    def format(self, record):
        levelname = record.levelname
        if self.use_color and levelname in COLORS:
            levelname_color = COLOR_SEQ % (30 + COLORS[levelname]) + levelname + RESET_SEQ
            record.levelname = levelname_color
        return logging.Formatter.format(self, record)
####################################################################################################


########################################### Environment ############################################
try:
    os.environ["PARSNPDIR"]
    PARSNP_DIR = os.environ["PARSNPDIR"]
except KeyError:
    PARSNP_DIR = sys.path[0]
SIGINT = False

try:
    os.environ["PYTHONPATH"] = PARSNP_DIR + os.pathsep + os.environ["PYTHONPATH"]
except KeyError:
    os.environ["PYTHONPATH"] = PARSNP_DIR + os.pathsep

frozenbinary = True
application_path = ""
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)
    frozenbinary = False

if frozenbinary:
   utilPath = PARSNP_DIR
   libPath = os.path.abspath(os.path.join(utilPath, "..", "lib"))
   if os.path.exists(libPath):
      oldLDPath = ""
      needToAdd = True
      if "LD_LIBRARY_PATH" in os.environ:
          oldLDPath = os.environ["LD_LIBRARY_PATH"]
          if libPath in oldLDPath:
              needToAdd = False
      elif "DYLD_FALLBACK_LIBRARY_PATH" in os.environ:
         oldLDPath = os.environ["DYLD_FALLBACK_LIBRARY_PATH"]
         if libPath in oldLDPath:
            needToAdd = False
      if needToAdd:
         os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = libPath + os.pathsep + oldLDPath
         os.environ["LD_LIBRARY_PATH"] = libPath + os.pathsep + oldLDPath

OSTYPE="linux"
p = subprocess.Popen("echo `uname`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(checkStdout, checkStderr) = p.communicate()
if checkStderr != b"":
    sys.stderr.write(WARNING_YELLOW+"Warning: Cannot determine OS, defaulting to %s\n"%(OSTYPE)+ENDC)
else:
    OSTYPE = checkStdout.decode('utf-8').strip()

binary_type = "linux"
if OSTYPE == "Darwin":
    binary_type = "osx"
else:
    binary_type = "linux"


# Should save parsnp alias for the main entry point 
# if not os.path.lexists("%s/bin/parsnp"%(PARSNP_DIR)):
    # os.system("ln -s %s/bin/parsnp %s/bin/parsnp"%(PARSNP_DIR, PARSNP_DIR))
if not os.path.lexists("%s/bin/harvest"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/harvest_%s %s/bin/harvest"%(PARSNP_DIR,binary_type,PARSNP_DIR))
if not os.path.lexists("%s/bin/ft"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/fasttree_%s %s/bin/ft"%(PARSNP_DIR,binary_type,PARSNP_DIR))
if not os.path.lexists("%s/bin/phiprofile"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/Profile_%s %s/bin/phiprofile"%(PARSNP_DIR,binary_type,PARSNP_DIR))
####################################################################################################


######################################## Utility Functions #########################################
def get_os():
    p = subprocess.Popen(
            "echo `uname`",
            shell=True,
            stdin=None,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    (checkStdout, checkStderr) = p.communicate()
    if checkStderr != b'':
        OSTYPE = "Linux"
        logger.warning("Cannot determine OS, defaulting to %s"%(OSTYPE))
    else:
        OSTYPE = checkStdout.decode('utf-8').strip()
    if OSTYPE == "Darwin":
        binary_type = "osx"
    else:
        binary_type = "linux"
    return OSTYPE, binary_type


def handler(signum, frame):
    global SIGINT
    SIGINT = True
    logger.critical('Caught request to terminate by user (CTRL+C), exiting now, bye')
    sys.exit(128)

signal.signal(signal.SIGINT, handler)


#TODO Merge run fns
def run_phipack(query,seqlen,workingdir):
    currdir = os.getcwd()
    os.chdir(workingdir)
    command = "%s/bin/phiprofile -o -v -n %d -w 100 -m 100 -f %s > %s.out"%(PARSNP_DIR,seqlen,query,query)
    run_command(command,1)
    os.chdir(currdir)

def run_fasttree(query,workingdir,recombination_sites):
    currdir = os.getcwd()
    os.chdir(workingdir)
    command = "%s/bin/ft -nt -quote -gamma -slow -boot 100 seq.fna > out.tree"%(PARSNP_DIR)
    run_command(command,1)
    os.chdir(currdir)


def run_bng(query,workingdir):
    currdir = os.getcwd()
    os.chdir(workingdir)
    command = "%s/run_fasta"%()
    run_command(command, 1)
    command = "%s/run_fasta"%()
    run_command(command, 1)
    command = "%s/run_fasta"%()
    run_command(command, 1)
    os.chdir(currdir)


#TODO Merge wrappers
def parallelWrapper(params):
   try:
        jobID = params["jobID"]
        result = {}
        result["jobID"] = jobID
        result["status"] = 0
        run_mummer(params["ref"], params["query"], params["prefix"])
        result["status"] = 1
        return result
   except KeyboardInterrupt:
        result["status"] = 0
        logger.info("Keyboard interrupt in thread %d, quitting\n"%(jobID))
        return result
   except Exception:
        result["status"] = 0
        logger.info("Other error in thread %d, quitting\n"%(jobID))
        return result


def parallelFtWrapper(params):
   try:
        jobID = params["jobID"]
        result = {}
        result["jobID"] = jobID
        result["status"] = 0
        run_fasttree(params["query"], params["dir"], params["recombination"])
        result["status"] = 1
        return result
   except KeyboardInterrupt:
        result["status"] = 0
        logger.info("Keyboard interrupt in thread %d, quitting\n"%(jobID))
        return result
   except Exception:
        result["status"] = 0
        logger.info("Other error in thread %d, quitting\n"%(jobID))
        return result


def parallelPhiWrapper(params):
   try:
        jobID = params["jobID"]
        result = {}
        result["jobID"] = jobID
        result["status"] = 0
        if params["seqlen"] >= 1000:
            run_phipack(params["query"],params["seqlen"],params["dir"])
            result["status"] = 1
        else:
            result["status"] = 2
        return result
   except KeyboardInterrupt:
        result["status"] = 0
        logger.info("Keyboard interrupt in thread %d, quitting\n"%(jobID))
        return result
   except Exception:
        result["status"] = 0
        logger.info("Other error in thread %d, quitting\n"%(jobID))
        return result


def run_command(command,ignorerc=0):
   global SIGINT
   logger.debug(command)
   p = subprocess.Popen(command, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,close_fds=True,executable="/bin/bash")
   fstdout,fstderr = p.communicate()
   rc = p.returncode
   if VERBOSE:
      logging.debug(fstderr)

   if rc != 0 and not SIGINT and not ignorerc and "rm " not in command and "ls " not in command and "unlink " not in command and "ln " not in command and "mkdir " not in command and "mv " not in command:
      logger.error("""The following command failed:
      >>$ {}
      Please veryify input data and restart Parsnp.
      If the problem persists please contact the Parsnp development team.""".format(command))
      sys.exit(rc)


def is_valid_file_path(parser, arg):
    if not os.path.exists(arg) and arg != "!" and arg != None and arg != "":
        logger.critical("The file %s does not exist!" % arg)
    else:
        return arg


def is_valid_dir(parser, arg):
    if not os.path.exists(arg):
        parser.error( "The directory %s does not exist\n" % (arg))
    if len(glob.glob("%s/*"%(arg))) == 0:
        parser.error("The director %s is empty"%(arg))


def parse_args():
    parser = argparse.ArgumentParser(description="""
    Parsnp quick start for three example scenarios:
    1) With reference & genbank file:
    python Parsnp.py -g <reference_genbank_file1 reference_genbank_file2 ...> -d <seq_file1 seq_file2 ...>  -p <threads>

    2) With reference but without genbank file:
    python Parsnp.py -r <reference_genome> -d <seq_file1 seq_file2 ...> -p <threads>

    3) Autorecruit reference to a draft assembly:
    python Parsnp.py -q <draft_assembly> -d <seq_file1 seq_file2 ...> -p <threads>
    """, formatter_class=argparse.RawTextHelpFormatter)
    #TODO Use lambda to check files and directories
    input_output_args = parser.add_argument_group(title="Input/Output")
    input_output_args.add_argument(
        "-c",
        "--curated",
        action = "store_true",
        help = "(c)urated genome directory, use all genomes in dir and ignore MUMi?")
    input_output_args.add_argument(
        "-d",
        "--sequences",
        type = str,
        nargs = '+',
        required = True,
        help = "A list of files containing genomes/contigs/scaffolds")
    input_output_args.add_argument(
        "-r",
        "--reference",
        type = lambda fname: is_valid_file_path(parser, fname),
        default = "",
        help = "(r)eference genome (set to ! to pick random one from sequence dir)")
    #TODO Accept as space-separated input and parse automatically w/ argparse
    input_output_args.add_argument(
        "-g",
        "--genbank",
        nargs = '*',
        help = "A list of Genbank file(s) (gbk)")
    input_output_args.add_argument(
        "-o",
        "--output-dir",
        type = str,
        default = "[P_CURRDATE_CURRTIME]")
    input_output_args.add_argument(
        "-q",
        "--query",
        type = str,
        help = "Specify (assembled) query genome to use, in addition to genomes found in genome dir")

    MUMi_args = parser.add_argument_group(title="MUMi")
    MUMi_mutex_args = MUMi_args.add_mutually_exclusive_group()
    #TODO whats the default?
    MUMi_mutex_args.add_argument(
        "-U",
        "--max-mumi-distr-dist",
        "--MUMi",
        type = float,
        default = 0.5,
        help = "Max MUMi distance value for MUMi distribution")
    #TODO Not parsed in current parsnp version and had a duplicate -i flag. Is this no longer used?
    MUMi_mutex_args.add_argument(
        "-mmd",
        "--max-mumi-distance",
        type = float,
        help = "Max MUMi distance (default: autocutoff based on distribution of MUMi values)")
    MUMi_args.add_argument(
        "-F",
        "--fastmum",
        action = "store_true",
        help = "Fast MUMi calculation")
    MUMi_args.add_argument(
        "-M",
        "--mumi_only",
        "--onlymumi",
        action = "store_true",
        help = "Calculate MUMi and exit? overrides all other choices!")

    MUMi_rec_prog = MUMi_args.add_mutually_exclusive_group()
    MUMi_rec_prog.add_argument(
        "--use-mummer-mumi",
        action = "store_true",
        help = "Use mummer for MUMi distance genome recruitment")
    MUMi_rec_prog.add_argument(
        "--use-ani",
        action = "store_true",
        help = "Use ani for genome recruitment")
    MUMi_args.add_argument(
        "--min-ani",
        type = float,
        default = 90,
        help = "Min ANI value to allow for genome recruitment.")
    MUMi_rec_prog.add_argument(
        "--use-mash",
        action = "store_true",
        help = "Use mash for genome recruitment")
    MUMi_args.add_argument(
        "--max-mash-dist",
        type = float,
        default = .1,
        help = "Max mash distance.")

    MUM_search_args = parser.add_argument_group(title="MUM search")
    #new, default to lower, 12-17
    MUM_search_args.add_argument(
        "-a",
        "--min-anchor-length",
        "--anchorlength",
        type = str,
        default = "1.1*(Log(S))",
        help = "Min (a)NCHOR length (default = 1.1*(Log(S)))")
    MUM_search_args.add_argument(
        "-C",
        "--max-cluster-d",
        "--clusterD",
        type = int,
        default = 300,
        help = "Maximal cluster D value")
    MUM_search_args.add_argument(
        "-z",
        "--min-cluster-size",
        "--minclustersize",
        type = int,
        default = 21,
        help = "Minimum cluster size")
    #TODO -z was a duplicate flag but no longer parsed as min-lcb-size in the current parsnp version
    # MUM_search_args.add_argument(
            # "-z",
            # "--min-lcb-size",
            # type = int,
            # default = 25,
            # help = "Min LCB si(z)e")

    LCB_alignment_args = parser.add_argument_group(title="LCB alignment")
    LCB_alignment_args.add_argument(
        "-D",
        "--max-diagonal-difference",
        "--DiagonalDiff",
        metavar = "MAX_DIAG_DIFF",
        type = str,
        default="0.12",
        help = "Maximal diagonal difference. Either percentage (e.g. 0.2) or bp (e.g. 100bp)")
    LCB_alignment_args.add_argument(
        "-n",
        "--alignment-program",
        "--alignmentprog",
        type = str,
        choices = list(ALIGNER_TO_IDX.keys()),
        default = "muscle",
        help = "Alignment program to use")
    LCB_alignment_args.add_argument(
        "-u",
        "--unaligned",
        action = "store_true",
        help = "Ouput unaligned regions")

    recombination_args = parser.add_argument_group("Recombination filtration")
    #TODO -x was a duplicate flag but no longer parsed as filter-phipack-snps in the current parsnp version
    # recombination_args.add_argument(
            # "-x",
            # "--filter-phipack-snps",
            # action = "store_true",
            # help = "Enable filtering of SNPs located in PhiPack identified regions of recombination")

    probe_design_args = parser.add_argument_group("Probe design")
    probe_design_args.add_argument(
        "-b",
        "--probe",
        action = "store_true",
        help = "Remove genome length constraints to search for MUMs in concatenated sequences much larger than reference")

    misc_args = parser.add_argument_group("Misc")
    misc_args.add_argument(
        "-p",
        "--threads",
        type = int,
        default = 1,
        help = "Number of threads to use")
    misc_args.add_argument(
        "-P",
        "--max-partition-size",
        type = int,
        default = 15000000,
        help = "Max partition size (limits memory usage)")
    misc_args.add_argument(
        "-v",
        "--verbose",
        action = "store_true",
        help = "Verbose output")
    misc_args.add_argument(
        "-V",
        "--version",
        action = "version",
        version = "%(prog)s " + __version__)

    todo_args = parser.add_argument_group("Need to be placed in a group")
    todo_args.add_argument(
        "-e",
        "--extend",
        action = "store_true")
    todo_args.add_argument(
        "-l",
        "--layout",
        action = "store_true")
    todo_args.add_argument(
        "-x",
        "--xtrafast",
        action = "store_true")
    todo_args.add_argument(
        "-s",
        "--split",
        action = "store_true",
        help = "Split genomes by n's")
    todo_args.add_argument(
        "-i",
        "--inifile",
        "--ini-file",
        type = str)
    todo_args.add_argument(
        "-m",
        "--mum-length",
        "--mumlength",
        type = str,
        default = "1.1*(Log(S))",
        help = "TODO!!!")
    return parser.parse_args()
####################################################################################################
#print("-g = <bool>: auto-launch (g)ingr? (default = NO)"


if __name__ == "__main__":
    t1 = time.time()
    sys.stderr.write( BOLDME+"|--Parsnp %s--|\n"%(VERSION))
    sys.stderr.write( BOLDME+"For detailed documentation please see --> http://harvest.readthedocs.org/en/latest\n")


    parsnp_dir= sys.path[0]
    #print parsnp_dir
    #PARSNP_DIR = parsnp_dir
    opts = []
    args = []

    OSTYPE, BINARY_TYPE = get_os()
    args = parse_args()
    currdir = os.getcwd()
    logging_level = logging.DEBUG if args.verbose else logging.INFO
    ref = args.reference
    if ref == '!':
        randomly_selected_ref = True
    input_files = args.sequences
    query = args.query
    anchor = args.min_anchor_length
    #TODO I'm guessing mummer_mumi was intended to be an option?
    use_mummer_mumi = args.use_mummer_mumi
    use_ani = args.use_ani
    use_mash = args.use_mash
    use_parsnp_mumi = not (use_mash or use_mummer_mumi or use_ani)
    mum = args.mum_length
    maxpartition = args.max_partition_size
    fastmum = args.fastmum
    cluster = args.max_cluster_d
    curated = args.curated
    aligner = ALIGNER_TO_IDX[args.alignment_program.lower()]
    threads = args.threads
    unaligned = "0" if not args.unaligned else "1"
    mincluster = args.min_cluster_size
    diagdiff = args.max_diagonal_difference
    splitseq = args.split
    extend = args.extend
    layout = args.layout
    xtrafast = args.xtrafast
    inifile = args.inifile
    inifile_exists = args.inifile is not None
    mumi_only = args.mumi_only
    mumidistance = args.max_mumi_distr_dist
    max_mash_dist = args.max_mash_dist
    min_ani_cutoff = args.min_ani
    outputDir = args.output_dir
    probe = args.probe
    genbank_file = ""
    genbank_files = []
    genbank_files_cat = ""
    genbank_ref = ""
    reflen = 0
    use_gingr = ""
    filtreps = False

    repfile = ""
    multifasta = False
    ref_seqs = {}

    # Instantiate logger
    logger = logging.getLogger("Parsnp")
    logger.setLevel(logging_level)
    ch = logging.StreamHandler()
    ch.setLevel(logging_level)
    # create formatter and add it to the handlers
    formatter = ColoredFormatter('%(asctime)s - %(levelname)s - %(message)s',
                                 datefmt="%H:%M:%S")
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)

    # Create output dir
    if outputDir == "." or outputDir == "./" or outputDir == "/":
        logger.critical("Specified output dir is current working dir or root dir! will clobber any parsnp.* results")
        sys.exit(1)
    elif outputDir == "[P_CURRDATE_CURRTIME]":
        today = datetime.datetime.now()
        timestamp = "P_"+today.isoformat().replace("-","_").replace(".","").replace(":","").replace("T","_")
        outputDir = os.getcwd()+os.sep+timestamp
    os.makedirs(outputDir, exist_ok=True)
    shutil.rmtree(os.path.join(outputDir, "tmp"), ignore_errors=True)
    os.makedirs(os.path.join(outputDir, "tmp"), exist_ok=True)

    input_files_processed = []
    for input_f in input_files:
        if os.path.isdir(input_f):
            for f in os.listdir(input_f):
                f = os.path.join(input_f, f)
                if os.path.isfile(f):
                    input_files_processed.append(f)
        elif os.path.isfile(input_f):
            input_files_processed.append(input_f)
    input_files = input_files_processed

    # Parse reference if necessary
    if ref and ref != "!":
        try:
            rf = open(ref, 'r')
            rfd = rf.read()
            refseqs = rfd.split(">")[1:]
            currpos = 0
            if len(refseqs) > 1:
                multifasta = True
                for seqnum, seq in enumerate(refseqs):
                    seq = seq.split('\n', 1)[1]
                    fastalen = len(seq) - seq.count('\n')
                    ref_seqs[currpos + fastalen] = seqnum
                    currpos += fastalen
            rf.close()
        except IOError as e:
            logger.critical(" Reference genome file %s not found\n"%(ref))
            sys.exit(1)

    # Validate genbank files
    #TODO Make this a function
    # return genbank_ref
    if args.genbank:
        genbank_files = args.genbank
        genbank_files_processed = []
        for genbank_f in genbank_files:
            if os.path.isdir(genbank_f):
                for f in os.listdir(genbank_f):
                    f = os.path.join(genbank_f, f)
                    if os.path.isfile(f):
                        genbank_files_processed.append(f)
            elif os.path.isfile(genbank_f):
                genbank_files_processed.append(genbank_f)
        genbank_files = genbank_files_processed
        ctcmd = "cat "

        first = True
        #genbank_ref = ""
        for genbank_file in genbank_files:
            if len(genbank_file) <= 1:
                continue
            ctcmd += genbank_file + " "
            genbank_ref = os.path.join(outputDir, "tmp", os.path.basename(genbank_file)+".fna")
            try:
                #parse out reference, starts at ORIGIN ends at //, remove numbers,
                rf = open(genbank_file,'r')
                genbank_ref_d = open(genbank_ref, "a+")
                while True:
                    giline = rf.readline()
                    if "VERSION" and "GI" in giline:
                        break
                    elif giline == None or giline == "":
                        logger.critical("Genbank file %s malformatted \n"%(genbank_file))
                        sys.exit(1)
                if len(giline) <= 2:
                    logger.critical("Genbank file %s malformatted \n"%(genbank_file))
                    sys.exit(1)
                genbank_ref_d.write(">gi|"+giline.split("GI:")[-1])
                first = False
                ntdata = False
                data = ""
                for line in rf:
                    if ntdata:
                        if "//" in line:
                            ntdata = False
                            break
                        data += line[9:].replace(" ","")
                    if "ORIGIN" in line:
                         ntdata = True

                rf.close()
                if len(data) < 10:
                      logger.critical("Genbank file %s contains no sequence data\n"%(genbank_file))
                      sys.exit(1)
                genbank_ref_d.write(data.upper())
                genbank_ref_d.close()
            except IOError as e:
                logger.critical("Genbank file %s not found\n"%(genbank_file))
                sys.exit(1)



    sortem = True
    ref_string = ref
    genome_string = ""
    if len(input_files) > 1:
        genome_string = "\n\t"
        if len(input_files) > 4:
            genome_string += "\n\t".join(input_files[:2])
            genome_string += "\n\t...{} more file(s)...\n\t".format(len(input_files) - 4)
            genome_string += "\n\t".join(input_files[-2:])
        else:
            genome_string += "\n\t".join(input_files)
    else:
        genome_string = input_files[0]
    if len(ref) == 0 and len(genbank_ref) != 0:
        #we are parsing from genbank, set ref to genbank_ref && turn off sorting
        ref = genbank_ref
        if len(genbank_files) > 1:
            ref_string = "\n\t"
            if len(genbank_files) > 4:
                ref_string += "\n\t".join(genbank_files[:2])
                ref_string += "\n\t...{} more file(s)...\n\t".format(len(genbank_files) - 4)
                ref_string += "\n\t".join(genbank_files[-2:])
            else:
                ref_string = "\n\t".join(genbank_files)
        else:
            ref_string += genbank_files[0]

        sortem = False

    autopick_ref = False
    if (not ref and not query) or not input_files:
        logger.critical("No seqs provided, yet required. exit!")
        sys.exit(0)  # TODO Should this exit value be 0?
    elif not ref and query:
        logger.warning("No reference genome specified, going to autopick from %s as closest to %s\n"%(seqdir, query))
        autopick_ref = True
        ref = query

    logger.info("""
{}
SETTINGS:
|-refgenome:\t{}
|-genomes:\t{}
|-aligner:\t{}
|-outdir:\t{}
|-OS:\t{}
|-threads:\t{}
{}
    """.format(
        (len(outputDir)+17)*"*",
        "autopick" if ref == '!' else ref_string,
        genome_string,
        args.alignment_program,
        outputDir,
        OSTYPE,
        threads,
        (len(outputDir)+17)*"*"))



    logger.info("<<Parsnp started>>")

    #1)read fasta files (contigs/scaffolds/finished/DBs/dirs)
    # logger.info("Reading Genbank file(s) for reference (.gbk) %s"%("\t".join(genbank_files)))
    if len(genbank_file) == 0:
        logger.info("No genbank file provided for reference annotations, skipping..")

    allfiles = []
    fnaf_sizes = {}
    allfile_dict = {} 
    reflen = 0
    fnafiles = []
    if ref == "!":
        ref = random.choice(input_files)

    # Check if reference genome is aligned
    with open(ref, 'r') as ff:
        hdr = ff.readline()
        seq = ff.read()
        if hdr[0] != ">":
            logger.critical("Reference {} has improperly formatted header.".format(ref))
            sys.exit(1)
        for line in seq.split('\n'):
            if '-' in line and line[0] != ">":
                logger.warning("Reference genome sequence %s has '-' in the sequence!"%((ref)))
        reflen = len(seq) - seq.count('\n')

    for input_file in input_files:
        ff = open(input_file, 'r')
        hdr = ff.readline()
        seq = ff.read()
        name_flag = True
        seqlen = len(seq) - seq.count('\n')
        if hdr[0] != ">":
            logger.error("{} has improperly formatted header. Skip!".format(input_file))
            continue
        elif '-' in seq:
            logger.warning("Genome sequence %s seems to aligned! Skip!"%((input_file)))
            # continue
        elif seqlen <= 20:
            logger.error("File %s is less than or equal to 20bp in length. Skip!"%(input_file))
            continue
        sizediff = float(reflen)/float(seqlen)

        # EDITED THIS TO CHANGE GENOME THRESHOLD
        # WILL NOW CONSIDER CONCATENATED GENOMES THAT ARE MUCH BIGGER THAN THE REFERENCE
        if not args.probe and sizediff <= 0.6:
                logger.warning("File %s is too long compared to reference!"%(input_file))
                continue
        else:
            if sizediff >= 1.4:
                logger.warning("File %s is too short compared to reference genome!"%(input_file))
                continue
        fnafiles.append(input_file)
        fnaf_sizes[input_file] = seqlen
        ff.close()

    # if ref in fnafiles:
        # fnafiles.remove(ref)

    #sort reference by largest replicon to smallest
    if sortem and os.path.exists(ref) and not autopick_ref:
        logger.debug("Sorting reference replicons")
        ff = open(ref, 'r')
        seqs = ff.read().split(">")[1:]
        seq_dict = {}
        seq_len = {}
        for seq in seqs:
            try:
                hdr, seq = seq.split("\n",1)
            except ValueError:
                # TODO Why do we ignore when theres a header but no sequence?
                continue
            seq_dict[hdr] = seq
            seq_len[hdr] = len(seq) - seq.count('\n')
        seq_len_sort = sorted(iter(seq_len.items()), key=operator.itemgetter(1), reverse=True)
        ref = os.path.join(outputDir, os.path.basename(ref)+".ref")
        ffo = open(ref, 'w')
        for hdr, seq in seq_len_sort:
            ffo.write(">%s\n"%(hdr))
            ffo.write("%s"%(seq_dict[hdr]))
        ff.close()
        ffo.close()
    else:
        ref = genbank_ref

    # TODO stray comment: remove any query sequences 30% diff in length
    allfiles = [os.path.basename(ref)]
    #write INI file
    if not inifile_exists:
        logger.debug("Writing .ini file")
        if xtrafast or 1:
            extend = False

        inifiled = open("%s/template.ini"%(PARSNP_DIR), 'r').read()
        inifiled = inifiled.replace("$REF", ref)
        inifiled = inifiled.replace("$EXTEND", "%d"%(extend))
        inifiled = inifiled.replace("$ANCHORS", str(anchor))
        inifiled = inifiled.replace("$MUMS", str(mum))
        inifiled = inifiled.replace("$MINCLUSTER", str(mincluster))
        inifiled = inifiled.replace("$CLUSTERD", str(cluster))
        inifiled = inifiled.replace("$THREADS", str(threads))
        inifiled = inifiled.replace("$ALIGNER", str(aligner))
        inifiled = inifiled.replace("$DIAGDIFF", str(diagdiff))
        inifiled = inifiled.replace("$RECOMBFILT", "%d"%(xtrafast))
        inifiled = inifiled.replace("$OUTDIR", outputDir)
        if fastmum:
            inifiled = inifiled.replace("$PARTPOS","%d"%(0.2*reflen))
        else:
            inifiled = inifiled.replace("$PARTPOS","%s"%(maxpartition))

        file_string = ""
        for cnt, fna_file in enumerate(fnafiles, 1):
            file_string += "file%d=%s\n"%(cnt, fna_file)
            file_string += "reverse%d=0\n"%(cnt)
        inifiled_mumi = inifiled.replace("$FILES\n", file_string)
        inifiled_mumi = inifiled_mumi.replace("calcmumi=0","calcmumi=1")
        inifile_mumi = open(os.path.join(outputDir, "all_mumi.ini"), 'w')
        inifile_mumi.write(inifiled_mumi)
        inifile_mumi.close()

    #2)get near neighbors (mumi distance)
    if os.path.exists(os.path.join(outputDir, "alltogether.fasta")):
        os.remove(os.path.join(outputDir, "alltogether.fasta"))
    if os.path.exists(os.path.join(outputDir, "blocks/b1")):
        ftrm = glob.glob(os.path.join(outputDir, "blocks/b*"))
        for f in ftrm:
            shutil.rmtree(f)

    fileidx = -1

    hit_dict = {}
    qry_hit_dict = {}
    hdr_dict = {}
    length_dict = {}

    TOTSEQS= len(fnafiles) + 1
    seqids_list = []

    if len(fnafiles) < 1 or ref == "":
        logger.critical("Parsnp requires 2 or more genomes to run, exiting")
        logger.debug("Only files found are: {}\n{} ".format(fnafiles, ref))
        sys.exit(1)

    mumi_dict = {}
    if not curated:
        logger.info("Recruiting genomes...")
        finalfiles = []
        auto_ref = ""
        if use_parsnp_mumi:
            if not inifile_exists:
                command = "%s/bin/parsnp %sall_mumi.ini"%(PARSNP_DIR,outputDir+os.sep)
            else:
                # TODO why are we editing the suffix of a provided file?
                command = "%s/bin/parsnp %s"%(PARSNP_DIR,inifile.replace(".ini","_mumi.ini"))
            run_command(command)
            try:
                mumif = open(os.path.join(outputDir, "all.mumi"),'r')
                for line in mumif:
                    line = line.rstrip('\n')
                    idx, mi = line.split(":")
                    mumi_dict[int(idx)-1] = float(mi)
            except IOError:
                logger.error("MUMi file generation failed... use all?")
                for i, _ in enumerate(fnafiles):
                    mumi_dict[i] = 1
            lowest_mumi = 100

            if autopick_ref:
                for idx in list(mumi_dict.keys()):
                    #TODO is there a way to organize these via dict rather than list? Seems error prone
                    if mumi_dict[idx] < lowest_mumi:
                        auto_ref = fnafiles[idx]
                        ref = auto_ref
                        lowest_mumi = mumi_dict[idx]
            mumi_f = ""
            if mumi_only and not curated:
                mumi_f = open(os.path.join(outputDir, "recruited_genomes.lst"),'w')


            sorted_x = sorted(iter(mumi_dict.items()), key=operator.itemgetter(1))
            mumivals = []
            for scnt, item in enumerate(sorted_x):
                if scnt > 100 or scnt >= len(sorted_x):
                    break
                if float(item[1]) < float(mumidistance):
                    mumivals.append(float(item[1]))
            minv = minv = numpy.percentile(mumivals, 0) if len(mumivals) > 0 else 1.0
            dvals = mumivals

            stdv = 0
            hpv = 0
            if len(dvals) > 0:
                stdv = numpy.std(dvals)
                hpv = minv + (3*stdv)

            for idx in mumi_dict.keys():
                if mumi_dict[idx] < (float(mumidistance)) or curated:
                    if fastmum and mumi_dict[idx] > hpv:
                        continue
                    #TODO if 1, why is this?
                    if 1 or auto_ref != fnafiles[idx]:
                        if mumi_only:
                            mumi_f.write(os.path.abspath(fnafiles[idx])+",%f"%(mumi_dict[idx])+"\n")
                        finalfiles.append(fnafiles[idx])
                        allfiles.append(fnafiles[idx])

        else:
            try:
                # tmp_dir = tempfile.mkdtemp()
                tmp_dir = outputDir
                all_genomes_fname = os.path.join(tmp_dir, "genomes.lst")
                # with open(all_genomes_fname, 'w') as all_genomes_file:
                    # all_genomes_file.write("\n".join(fnafiles))
                for f in fnafiles:
                    print(f)
                if use_mash:
                    if randomly_selected_ref:
                        logger.warning("You are using a randomly selected genome to recruit genomes from your input...")
                    mash_out = subprocess.check_output([
                            "mash", "dist", "-t", 
                            "-d", str(max_mash_dist), 
                            "-p", str(threads), 
                            ref, 
                            "-l", all_genomes_fname],
                        stderr=open(os.path.join(outputDir, "mash.err"), 'w')).decode('utf-8')
                    finalfiles = [line.split('\t')[0] for line in mash_out.split('\n')[1:] if line != '']
                elif use_ani:
                    if randomly_selected_ref:
                        print(" ".join([
                                "fastANI", 
                                "--ql", all_genomes_fname, 
                                "--rl", all_genomes_fname, 
                                "-t", str(threads),
                                "-o", os.path.join(outputDir, "fastANI.tsv")]))
                        subprocess.check_call([
                                "fastANI", 
                                "--ql", all_genomes_fname, 
                                "--rl", all_genomes_fname, 
                                "-t", str(threads),
                                "-o", os.path.join(outputDir, "fastANI.tsv")],
                            stderr=open(os.path.join(outputDir, "fastANI.err"), 'w'))
                    else:
                        subprocess.check_call([
                                "fastANI", 
                                "--q", ref, 
                                "--rl", all_genomes_fname, 
                                "-t", str(threads),
                                "-o", os.path.join(outputDir, "fastANI.tsv")],
                            stderr=open(os.path.join(outputDir, "fastANI.err"), 'w'))
                    genome_to_genomes = defaultdict(set)
                    with open(os.path.join(outputDir, "fastANI.tsv")) as results:
                        print(len(list(results)))
                    with open(os.path.join(outputDir, "fastANI.tsv")) as results:
                        for line in results:
                            # FastANI results file -> Query, Ref, ANI val, extra stuff,,,
                            line = line.split('\t')
                            # if float(line[2]) >= min_ani_cutoff:
                            genome_to_genomes[line[0]].add(line[1])
                        
                        # for g in genome_to_genomes:
                            # print(len(g))
                        ani_ref = max(genome_to_genomes, key=(lambda key: len(genome_to_genomes[key])))
                        if autopick_ref:
                            auto_ref = ani_ref
                        finalfiles = list(genome_to_genomes[ani_ref])
                        print(set(fnafiles).difference(set(finalfiles)))
                        
                # shutil.rmtree(tmp_dir)
            except subprocess.CalledProcessError as e:
                logger.critical(
                    "Recruitment failed with exception {}. More details may be found in the *.err output log".format(str(e))) 
                # shutil.rmtree(tmp_dir)
            allfiles.extend(finalfiles)

    if curated:
        for f in fnafiles:
            if f not in finalfiles:
                finalfiles.append(f)
            if f not in allfiles:
                allfiles.append(f)

    if mumi_only:
        mumi_f.close()
        sys.exit(1)

    orig_auto_ref = auto_ref
    if os.path.exists(auto_ref) and autopick_ref:
        #TODO This code block is duplicated
        ff = open(auto_ref, 'r')
        seqs = ff.read().split(">")[1:]
        seq_dict = {}
        seq_len = {}
        for seq in seqs:
            try:
                hdr, seq = seq.split("\n",1)
            except ValueError:
                continue
            seq_dict[hdr] = nt
            seq_len[hdr] = len(seq) - seq.count('\n')
        seq_len_sort = sorted(seq_len.iteritems(), key=operator.itemgetter(1))
        seq_len_sort.reverse()
        auto_ref = os.path.join(outputDir, os.path.basename(auto_ref)+".ref")
        ffo = open(ref, 'w')
        for item in seq_len_sort:
            ffo.write(">%s\n"%(item[0]))
            ffo.write(seq_dict[item[0]])
        ff.close()
        ffo.close()
        ref = auto_ref

    inifiled_closest = inifiled
    #TODO This code is duplicated
    if not inifile_exists:
        if len(finalfiles) < 1 or ref == "":
            logger.critical("Parsnp requires 2 or more genomes to run, exiting\n")
            sys.exit(1)

        file_string = ""
        cnt = 1
        file_string_closest = ""
        #TODO whats the point of iterating over one file?
        for cnt, f in enumerate(finalfiles[0:1], 1):
            file_string_closest+="file%d=%s\n"%(cnt, f)
            file_string_closest+="reverse%d=0\n"%(cnt)
        for cnt, f in enumerate(finalfiles, 1):
            file_string+="file%d=%s\n"%(cnt, f)
            file_string+="reverse%d=0\n"%(cnt)
        inifiled = inifiled.replace("$FILES\n",file_string)
        #new, output unaligned regions
        inifiled = inifiled.replace("$UNALIGNED",unaligned)
        inifiled_closest = inifiled.replace("$FILES\n",file_string_closest)

        if fastmum:
            inifiled = inifiled.replace("p=%d"%(0.2*reflen),"p=%s"%(maxpartition))
            inifiled_closest = inifiled.replace("p=%d"%(0.2*reflen),"p=%s"%(maxpartition))
        if autopick_ref:
            inifiled = inifiled.replace(orig_auto_ref,auto_ref)
            inifiled = inifiled.replace(auto_ref,"tmp_"+auto_ref)
            inifiled = inifiled.replace(query,auto_ref)
            inifiled = inifiled.replace("tmp_"+auto_ref,query)
            inifiled_closest = inifiled_closest.replace(auto_ref,"tmp_"+auto_ref)
            inifiled_closest = inifiled_closest.replace(query,auto_ref)
            inifiled_closest = inifiled_closest.replace("tmp_"+auto_ref,query)

        inifile = open(outputDir+os.sep+"parsnpAligner.ini",'w')
        inifile.write(inifiled)
        inifile.close()
        inifile_closest = open(outputDir+os.sep+"psnn.ini",'w')
        inifile_closest.write(inifiled_closest)
        inifile_closest.close()


    #3)run parsnp (cores, grid?)
    logger.info("Running Parsnp multi-MUM search and libMUSCLE aligner...")
    blocks_dir = os.path.join(outputDir, "blocks")
    if not os.path.exists(blocks_dir):
        os.mkdir(blocks_dir)
    command = ""
    run_parsnp = 1
    if run_parsnp:
        successful_run = False
        maxruns = 2
        runcnt = 0
        while not successful_run:
            if not inifile_exists:
                if command == "" and xtrafast and 0:
                    command = "%s/parsnpA_fast %sparsnpAligner.ini"%(PARSNP_DIR,outputDir+os.sep)
                elif command == "":
                    command = "%s/bin/parsnp %sparsnpAligner.ini"%(PARSNP_DIR,outputDir+os.sep)
                else:
                    command = "%s/bin/parsnp %spsnn.ini"%(PARSNP_DIR,outputDir+os.sep)
            else:
                if not os.path.exists(inifile):
                    logger.error("ini file %s does not exist!\n"%(inifile))
                    sys.exit(1)
                command = "%s/bin/parsnp %s"%(PARSNP_DIR,inifile)
            run_command(command)

            if not os.path.exists(os.path.join(outputDir, "parsnpAligner.xmfa")):
                successful_run = False
                runcnt += 1
                if runcnt >= 2:
                    logger.critical("Set of recruited genomes are too divergent for parsnp, please reduce MUMi (%f) and relaunch\n"%(float(mumidistance)))
                    sys.exit(1)
            else:
                successful_run = True
                runcnt += 1
                break
        shutil.move(
                os.path.join(outputDir, "parsnpAligner.xmfa"),
                os.path.join(outputDir, "parsnp.xmfa"))
    xmfafile = open(os.path.join(outputDir, "parsnp.xmfa"),'r')

    file2hdr_dict = {}
    fileid = ""
    blockfiles = []

    #get coverage
    coverage = 0
    totlength = 0
    totseqs = 0
    try:
        cf = open(os.path.join(outputDir, "parsnpAligner.log"))
        for line in cf:
            if "Total coverage among all sequences:" in line:
                coverage = line.split(":",1)[-1].replace("\n","")
                coverage = float(coverage.replace("%",""))/100.0
            elif "Length:" in line:
                totlength += int(line.split(":",1)[-1].replace("\n","").split("bps")[0])
                totseqs += 1
    except IOError:
        logger.critical("ParsnpAligner.log missing, parsnpAligner failed.")
        sys.exit(1)

    #update thresholds
    if coverage < 0.1 and not args.probe:
        if coverage <= 0.01:
            logger.critical("""Aligned regions cover less than 1% of reference genome, something is not right
Adjust params and rerun. If issue persists please submit a GitHub issue""")
            sys.exit(1)
        else:
            logger.warning("""Aligned regions cover less than 10% of reference genome!
Please verify recruited genomes are all strain of interest""")
    else:
        pass
    t2 = time.time()
    elapsed = float(t2)-float(t1)
    #print("-->Getting list of LCBs.."
    allbfiles = glob.glob(os.path.join(blocks_dir, "b*/*"))
    blockfiles = []
    block_startpos = []
    block_dict = {}
    for f in allbfiles:
        if os.path.isfile(f):
            if "seq.fna" in f:
                blockfiles.append(f)
                lf = open(f, 'r')
                header = lf.readline()
                if header[0] != ">":
                    logger.error("Error with LCB: %s\n"%(f))
                    continue

                inf = header.split("+",1)[0]

                rseq = ""
                while 1:
                    lff = lf.readline()
                    if lff[0] == ">":
                        break
                    rseq += lff.replace("\n","")

                spos,epos = inf.split(":",1)[-1].split("-",1)
                block_startpos.append(int(spos))
                block_dict[f] = [int(spos),int(epos), rseq]
                lf.close()
    run_repeat_filter = filtreps

    #initiate parallelPhiPack tasks
    run_recomb_filter = 0

    if xtrafast:
        run_recomb_filter = 1
    else:
        run_recomb_filter = 0

    recombination_sites = {}
    bedfile = ""
    bedfile_dict = {}
    if run_recomb_filter and len(blockfiles) > 0:
        logger.info("Running PhiPack on LCBs to detect recombination...")
        bedfile = open(os.path.join(outputDir, "parsnp.rec"), 'w')
        tasks = []
        processed = []
        for icnt, f in enumerate(blockfiles):
            seq1 = ""
            try:
                bf = open(f, 'r')
                seq1 = bf.read().split(">")[1].split("\n",1)[-1]
                seq1 = seq1.replace("\n","")
                bf.close()
            except IOError:
                pass

            processed.append(f)
            params = {}
            path, f = f.rsplit(os.path.sep,1)
            params["jobID"] = len(tasks)
            params["query"] = f
            params["seqlen"] = len(seq1)
            params["spos"] = block_startpos[icnt]
            params["dir"] = path
            params["output"] = os.path.join(path, "Profile.csv")
            tasks.append(params)

        #run parallelPhiPack
        pool = Pool(processes=int(threads))
        result = pool.map_async(parallelPhiWrapper,tasks).get(sys.maxsize)

        for i in result:
            if (i["status"] == 1):
                #process output
                recregions = ""
                block_spos = tasks[i["jobID"]]["spos"]
                try:
                    recregions = open(tasks[i["jobID"]]["output"],'r').read()
                except IOError:
                    logger.error("File %s doesn't exist, no rec regions or error in PhiPack\n"%(tasks[i["jobID"]]["output"]))
                    continue
                reclines = recregions.split("\n")
                prevpos = 0

                for line in reclines:
                    try:
                        pos,eval = line.split(",")
                    except ValueError:
                        continue
                    pos = int(pos)
                    eval = float("%.5f"%(float(eval)))
                    if eval < 0.01 and eval >= 0:
                        idx = 0
                        srpos = 0
                        if pos-50 > 0:
                            srpos = (pos-50)+block_spos
                        else:
                            srpos = block_spos
                        eval = abs(eval)
                        if not multifasta:
                            bedfile_dict[srpos] = "1\t%s\t%s\tREC\t%.3f\t+\n"%(srpos,pos+50+block_spos,eval)
                        else:
                            chrnum = 1
                            for cs in ref_seqs:
                                if block_spos < chr_spos:
                                    chrnum = ref_seqs[cs]
                            bedfile_dict[srpos] = "%d\t%s\t%s\tREC\t%.3f\t+\n"%(chrnum,srpos,pos+50+block_spos,eval)

                qfile = tasks[i["jobID"]]["query"]

            elif i["status"] != 2:
                logger.critical("Parallel phipack job %d failed\n"%(i["jobID"]))
                raise IOError

        pool.close()
        pool.join()
        brkeys = list(bedfile_dict.keys())
        brkeys.sort()
        for key in brkeys:
            bedfile.write(bedfile_dict[key])
        bedfile.close()

    run_lcb_trees = 0

    annotation_dict = {}
    #TODO always using xtrafast?
    if xtrafast or 1:
        #add genbank here, if present
        if len(genbank_ref) != 0:
            rnc = "%s/bin/harvest -q -o %s/parsnp.ggr -x "%(PARSNP_DIR,outputDir)+outputDir+os.sep+"parsnp.xmfa"
            for file in genbank_files:
                rnc += " -g %s " %(file)
            run_command(rnc)
        else:
            run_command("%s/bin/harvest -q -o %s/parsnp.ggr -f %s -x "%(PARSNP_DIR,outputDir,ref)+outputDir+os.sep+"parsnp.xmfa")

        if run_recomb_filter:
            run_command("%s/bin/harvest -q -b %s/parsnp.rec,REC,\"PhiPack\" -o %s/parsnp.ggr -i %s/parsnp.ggr"%(PARSNP_DIR,outputDir,outputDir,outputDir))
        if run_repeat_filter:
            run_command("%s/bin/harvest -q -b %s,REP,\"Intragenomic repeats > 100bp\" -o %s/parsnp.ggr -i %s/parsnp.ggr"%(PARSNP_DIR,repfile,outputDir,outputDir))

        run_command("%s/bin/harvest -q -i %s/parsnp.ggr -S "%(PARSNP_DIR,outputDir)+outputDir+os.sep+"parsnp.snps.mblocks")

    command = "%s/bin/ft -nt -quote -gamma -slow -boot 100 "%(PARSNP_DIR)+outputDir+os.sep+"parsnp.snps.mblocks > "+outputDir+os.sep+"parsnp.tree"
    logger.info("Reconstructing core genome phylogeny...")
    run_command(command)
    #7)reroot to midpoint
    if os.path.exists("outtree"):
        os.remove("outtree")

    if reroot_tree and len(finalfiles) > 1:
        try:
            mtree = open("%sparsnp.tree"%(outputDir+os.sep), 'r')
            mtreedata = mtree.read()
            mtreedata = mtreedata.replace("\n","")
            tree = dendropy.Tree.get_from_string(mtreedata,"newick")
            tree.reroot_at_midpoint(update_bipartitions=False)
            mftreef = tree.as_string('newick').split(" ",1)[1]
            #print mftreef
            mtreef = open(outputDir+os.sep+"parsnp.final.tree",'w')
            mtreef.write(mftreef)
            mtreef.close()
            os.system("mv %s %s"%(outputDir+os.sep+"parsnp.final.tree",outputDir+os.sep+"parsnp.tree"))
        except IOError:
            logger.error("Cannot process fasttree output, skipping midpoint reroot..\n")


    if len(use_gingr) > 0:
        logger.info("Creating Gingr input file..")
        if xtrafast or 1:
            #if newick available, add
            #new flag to update branch lengths
            run_command("%s/bin/harvest --midpoint-reroot -u -q -i "%(PARSNP_DIR)+outputDir+os.sep+"parsnp.ggr -o "+outputDir+os.sep+"parsnp.ggr -n %s"%(outputDir+os.sep+"parsnp.tree "))


    if float(elapsed)/float(60.0) > 60:
        logger.info("Aligned %d genomes in %.2f hours"%(totseqs,float(elapsed)/float(3600.0)))
    elif float(elapsed) > 60:
        #TODO just format the time to get rid of the above formatting
        logger.info("Aligned %d genomes in %.2f minutes"%(totseqs,float(elapsed)/float(60.0)))
    else:
        logger.info("Aligned %d genomes in %.2f seconds"%(totseqs,float(elapsed)))
    #cleanup
    rmfiles = glob.glob(os.path.join(outputDir, "*.aln"))
    #rmfiles2 = glob.glob(outputDir+os.sep+"blocks/b*/*")
    rmfiles3 = glob.glob(os.path.join(outputDir, "blocks/b*"))
    for f in rmfiles:
        os.remove(f)
    for f in rmfiles3:
        shutil.rmtree(f)

    filepres = 0
    logger.info("Parsnp finished! All output available in %s"%(outputDir))
    logger.debug("Validating output directory contents")
    print("\t1)parsnp.tree:\t\tnewick format tree", end =' ')
    if os.path.exists("%sparsnp.tree"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.tree"%(outputDir+os.sep)) > 0:
        print("\t\t\t["+OK_GREEN+"OK"+"]")
        filepres+=1
    else:
        print("\t|->"+ERROR_RED+"MISSING")
    print(BOLDME+"\t2)parsnp.ggr:\t\tharvest input file for gingr (GUI)", end =' ')
    if os.path.exists("%sparsnp.ggr"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.ggr"%(outputDir+os.sep)) > 0:
        print("\t["+OK_GREEN+"OK"+"]")
        filepres+=1
    else:
        print("\t|->"+ERROR_RED+"MISSING")
    print(BOLDME+"\t3)parsnp.xmfa:\t\tXMFA formatted multi-alignment", end = ' ')
    if os.path.exists("%sparsnp.xmfa"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.xmfa"%(outputDir+os.sep)) > 0:
        print("\t\t["+OK_GREEN+"OK"+"]")
        filepres+=1
    else:
        print("\t|->"+ERROR_RED+"MISSING")
    if filepres == 3:
        pass

    else:
        logger.error("Output files missing, something went wrong. Check logs and relaunch or contact developers for assistance"+"]")
    if os.path.exists("%sblocks"%(outputDir+os.sep)):
        os.rmdir("%sblocks"%(outputDir+os.sep))
    if os.path.exists("allmums.out"):
        os.remove("allmums.out")

    if not VERBOSE and os.path.exists("parsnpAligner.ini"):
        os.remove("parsnpAligner.ini")

    prefix = os.path.join(outputDir, os.path.splitext(os.path.basename(ref))[0])
    if not VERBOSE and os.path.exists("%s.coords"%(prefix)):
        os.remove("%s.coords"%(prefix))

    if not VERBOSE and os.path.exists("%s.delta"%(prefix)):
        os.remove("%s.delta"%(prefix))

    for f in glob.glob(os.path.join(outputDir,"*.reps")):
        if not VERBOSE and os.path.exists(f):
            os.remove(f)

    for f in glob.glob(os.path.join(outputDir, "*.ref")):
        if not VERBOSE and os.path.exists(f):
            os.remove(f)

    if not VERBOSE and os.path.exists("%s/psnn.ini"%(outputDir)):
        os.remove("%s/psnn.ini"%(outputDir))

    if not VERBOSE and os.path.exists("%s/all_mumi.ini"%(outputDir)):
        os.remove("%s/all_mumi.ini"%(outputDir))


    if os.path.exists("%s/parsnp.snps.mblocks"%(outputDir)):
        os.remove("%s/parsnp.snps.mblocks"%(outputDir))

    if not VERBOSE and os.path.exists("%s/all.mumi"%(outputDir)):
        os.remove("%s/all.mumi"%(outputDir))

    if os.path.exists(use_gingr):
        #check if available first
        rc = 0
        if binary_type == "osx":
            logger.info("Launching gingr..")
            os.system("open -n %s --args %s/parsnp.ggr"%(use_gingr,outputDir))

