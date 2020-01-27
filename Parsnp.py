# See the LICENSE file included with this software for license information.


import os, sys, string, getopt, random,subprocess, time, glob,operator, math, datetime,numpy #pysam
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
special_chars = [",","[","]","{","}","(",")","!","\'","\"","*","\%","\<" ,"\>", "|", " ", "`"]
ALIGNER_TO_IDX = {
        "mafft": "1",
        "muscle": "2",
        "fsa": "3",
        "prank": "4"
}
CSI=""#"\x1B["
reset=""#CSI+"m"
BOLDME = ""#CSI+'\033[1m'
STATUS_BLUE = ""#CSI+'\033[94m'
OK_GREEN = ""#CSI+'\033[92m'#'32m'
SKIP_GRAY = ""#CSI+'\033[37m'
WARNING_YELLOW = ""#CSI+'\033[93m'
ERROR_RED = ""#CSI+'\033[91m'
ENDC = ""#CSI+'0m'
PARSNP_DIR = sys.path[0]
try:
    os.environ["PARSNPDIR"]
    PARSNP_DIR = os.environ["PARSNPDIR"]
except KeyError:
    PARSNP_DIR = sys.path[0]
SIGINT = False

try:
    os.environ["PYTHONPATH"] = PARSNP_DIR+os.pathsep+os.environ["PYTHONPATH"]
except KeyError:
    os.environ["PYTHONPATH"] = PARSNP_DIR+os.pathsep

frozenbinary = True
application_path = ""
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)
    frozenbinary = False

if frozenbinary:
   utilPath = PARSNP_DIR
   libPath = os.path.abspath(utilPath + os.sep + ".." + os.sep + "lib")
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

VERBOSE = 0
VERSION = "v1.2"
PHI_WINDOWSIZE = 1000
TOTSEQS=0


#template
#read input
#0) init
## Get start time                                                                                                                                                                                                                                                                  
t1 = time.time()

p = subprocess.Popen(
        "echo `uname`", 
        shell=True, 
        stdin=None, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE)
(checkStdout, checkStderr) = p.communicate()
if checkStderr != b'':
    OSTYPE = "Linux"
    sys.stderr.write(WARNING_YELLOW+"Warning: Cannot determine OS, defaulting to %s\n"%(OSTYPE)+ENDC)
else:
    OSTYPE = checkStdout.decode('utf-8').strip()


if OSTYPE == "Darwin":
    binary_type = "osx"
else:
    binary_type = "linux"

def handler(signum, frame):
    global SIGINT
    SIGINT = True
    print('Caught request to terminate by user (CTRL+C), exiting now, bye')
    sys.exit(128)

signal.signal(signal.SIGINT, handler)


#TODO Merge run fns
def run_phipack(query,seqlen,workingdir):
    currdir = os.getcwd()
    os.chdir(workingdir)
    command = "%s/phiprofile -o -v -n %d -w 100 -m 100 -f %s > %s.out"%(PARSNP_DIR,seqlen,query,query)
    run_command(command,1)
    os.chdir(currdir)

def run_fasttree(query,workingdir,recombination_sites):
    currdir = os.getcwd()
    os.chdir(workingdir)
    command = "%s/ft -nt -quote -gamma -slow -boot 100 seq.fna > out.tree"%(PARSNP_DIR)
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
        sys.stderr.write("Keyboard interrupt in thread %d, quitting\n"%(jobID))
        return result
   except Exception:
        result["status"] = 0
        sys.stderr.write( "Other error in thread %d, quitting\n"%(jobID))
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
        sys.stderr.write( "Keyboard error in thread %d, quitting\n"%(jobID))
        return result
   except Exception:
        result["status"] = 0
        sys.stderr.write( "Other error in thread %d, quitting\n"%(jobID))
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
        sys.stderr.write( "Keyboard error in thread %d, quitting\n"%(jobID))
        return result
   except Exception:
        result["status"] = 0
        sys.stderr.write( "Other error in thread %d, quitting\n"%(jobID))
        return result

                                                                                 


def run_command(command,ignorerc=0):
   global SIGINT
   p = subprocess.Popen(command, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,close_fds=True,executable="/bin/bash")
   fstdout,fstderr = p.communicate()
   rc = p.returncode
   if VERBOSE:
      sys.stderr.write( fstderr)

   if rc != 0 and not SIGINT and not ignorerc and "rm " not in command and "ls " not in command and "unlink " not in command and "ln " not in command and "mkdir " not in command and "mv " not in command:
      sys.stderr.write(ERROR_RED+"**ERROR**"+ENDC+"\n")
      sys.stderr.write( "The following command failed:\n")
      sys.stderr.write( ">>%s\n"%(command))
      sys.stderr.write( "Please veryify input data and restart Parsnp. If the problem persists please contact the Parsnp development team.\n")
      sys.stderr.write(ERROR_RED+"**ERROR**"+ENDC+"\n")
      sys.stderr.write( "\n")
      sys.stderr.write( "\n")
      sys.exit(rc)

sys.stderr.write( BOLDME+"|--Parsnp %s--|\n"%(VERSION)+ENDC)
sys.stderr.write( BOLDME+"For detailed documentation please see --> http://harvest.readthedocs.org/en/latest\n"+ENDC)


if not os.path.lexists("%s/parsnp"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/parsnp %s/parsnp"%(PARSNP_DIR, PARSNP_DIR))
if not os.path.lexists("%s/harvest"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/harvest_%s %s/harvest"%(PARSNP_DIR,binary_type,PARSNP_DIR))
if not os.path.lexists("%s/ft"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/fasttree_%s %s/ft"%(PARSNP_DIR,binary_type,PARSNP_DIR))
if not os.path.lexists("%s/phiprofile"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/Profile_%s %s/phiprofile"%(PARSNP_DIR,binary_type,PARSNP_DIR))

if not os.path.lexists("%s/nucmer"%(PARSNP_DIR)):
    os.system("ln -s %s/MUMmer/nucmer %s/nucmer"%(PARSNP_DIR,PARSNP_DIR))
if not os.path.lexists("%s/show-coords"%(PARSNP_DIR)):
    os.system("ln -s %s/MUMmer/show-coords %s/show-coords"%(PARSNP_DIR,PARSNP_DIR))

#set MUMmer paths
if os.path.exists("%s/MUMmer/nucmer_run"%(PARSNP_DIR)):
    ff = open("%s/MUMmer/nucmer_run"%(PARSNP_DIR))
    ffd = ff.read()
    ff.close()
    ffd = ffd.replace("$MUMMERPATH1",PARSNP_DIR)
    ff = open("%s/MUMmer/nucmer"%(PARSNP_DIR),'w')
    ff.write(ffd)
    ff.close()

def is_valid_file(parser, arg):
    if not os.path.exists(arg) and arg != "!" and arg != None and arg != "":
        parser.error("The file %s does not exist!" % arg)
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
    python Parsnp.py -g <reference_genbank_file1,reference_genbank_file2,..> -d <genome_dir> -p <threads> 
     
    2) With reference but without genbank file:
    python Parsnp.py -r <reference_genome> -d <genome_dir> -p <threads> 

    3) Autorecruit reference to a draft assembly:
    python Parsnp.py -q <draft_assembly> -d <genome_db> -p <threads> 
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
            "--sequence-dir",
            "--sequenceDir",
            type = str, 
            help = "(d)irectory containing genomes/contigs/scaffolds") 
    input_output_args.add_argument(
            "-r", 
            "--reference",
            type = lambda fname: is_valid_file(parser, fname), 
            default = "",
            help = "(r)eference genome (set to ! to pick random one from sequence dir)")
    #TODO Accept as space-separated input and parse automatically w/ argparse
    input_output_args.add_argument(
            "-g", 
            "--genbank",
            type = str, 
            default = "",
            help = "Genbank file(s) (gbk), comma separated list")
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
            "--ini-file",
            "--inifile",
            type = str)
    todo_args.add_argument(
            "-m", 
            "--mum-length",
            "--mumlength",
            type = str, 
            default = "1.1*(Log(S))",
            help = "TODO!!!")
    return parser.parse_args()

#print("-g = <bool>: auto-launch (g)ingr? (default = NO)"


if __name__ == "__main__":
    parsnp_dir= sys.path[0]
    #print parsnp_dir
    #PARSNP_DIR = parsnp_dir
    opts = []
    args = []
    args = parse_args()

    print(args)
    currdir = os.getcwd()
    VERBOSE = args.verbose
    ref = args.reference
    query = args.query
    seqdir = args.sequence_dir
    anchor = args.min_anchor_length
    mum = args.mum_length
    maxpartition = args.max_partition_size
    fastmum = args.fastmum
    cluster = args.max_cluster_d
    curated = args.curated
    try:
        aligner = ALIGNER_TO_IDX[args.alignment_program.lower()]
    except KeyError:
        print("{} not supported as an alignment proram".format(args.alignment_program))
        sys.exit(1)
    threads = args.threads
    unaligned = "0" if not args.unaligned else "1"
    mincluster = args.min_cluster_size
    diagdiff = args.max_diagonal_difference
    splitseq = args.split
    extend = args.extend
    layout = args.layout
    xtrafast = args.xtrafast
    inifile= args.ini_file
    mumi_only = args.mumi_only
    mumidistance = args.max_mumi_distr_dist
    outputDir = args.output_dir
    probe = args.probe
    genbank_file = ""
    genbank_files = []
    genbank_files_str = ""
    genbank_files_cat = ""
    genbank_ref = ""
    reflen = 0
    use_gingr = ""
    inifile_exists = False
    filtreps = False

    repfile = ""
    multifasta = False
    ref_seqs = {}

    # Parse reference if necessary
    if ref and ref != "!":
        try:
            rf = open(ref,'r')
            rfd = rf.read()
            refseqs = rfd.split(">")[1:]
            currpos = 0
            seqnum = 1
            if len(refseqs) > 1:
                multifasta = True
                for seq in refseqs:
                    fastalen = len(seq.split("\n",1)[-1].replace("\n",""))
                    ref_seqs[currpos+fastalen] = seqnum
                    currpos = currpos+fastalen
                    seqnum+=1
            rf.close()
        except IOError:
            sys.stderr.write( "ERROR: Reference genome file %s not found\n"%(ref))
            sys.exit(1)       
    
    # Validate genbank files
    if args.genbank:
        genbank_files_str = args.genbank
        genbank_files = genbank_files_str.split(",")
        ctcmd = "cat "
        
        first = True
        #genbank_ref = ""
        for genbank_file in genbank_files_str.split(","):
            if len(genbank_file) <= 1:
                continue
            ctcmd += genbank_file + " "
            try:
                #parse out reference, starts at ORIGIN ends at //, remove numbers, 
                rf = open(genbank_file,'r')
                if first:
                    genbank_ref = genbank_file+".fna"
                    genbank_ref1 = open(genbank_ref,'w')
                    giline = ""
                    while 1:
                        giline = rf.readline()
                        if "VERSION" and "GI" in giline:
                            break
                        elif giline == None or giline == "":
                            sys.stderr.write( "ERROR: Genbank file %s malformatted \n"%(genbank_file))
                            sys.exit(1)
                    if len(giline) <= 2:
                        sys.stderr.write( "ERROR: Genbank file %s malformatted \n"%(genbank_file))
                        sys.exit(1)
                    genbank_ref1.write(">gi|"+giline.split("GI:")[-1])
                    first = False
                else:
                    genbank_ref1 = open(genbank_ref,'a')
                    giline = ""
                    while 1:
                        giline = rf.readline()
                        if "VERSION" and "GI" in giline:
                            break
                        elif giline == None or giline == "":
                            sys.stderr.write( "ERROR: Genbank file %s malformatted \n"%(genbank_file))
                            sys.exit(1)
                    if len(giline) <= 2:
                        sys.stderr.write( "ERROR: Genbank file %s malformatted \n"%(genbank_file))
                        sys.exit(1)
                    genbank_ref1.write(">gi|"+giline.split("GI:")[-1])
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
                      sys.stderr.write( "ERROR: Genbank file %s contains no sequence data\n"%(genbank_file))
                      sys.exit(1)
                genbank_ref1.write(data.upper())
                genbank_ref1.close()
            except IOError:
                sys.stderr.write( "ERROR: Genbank file %s not found\n"%(genbank_file))
                sys.exit(1)
                
            genbank_files_cat = "%s.cat"%(genbank_files[0])
            os.system(ctcmd+"> "+genbank_files_cat)

    # Create output dir
    if outputDir != "[P_CURRDATE_CURRTIME]":
        if outputDir == "." or outputDir == "./" or outputDir == "/":
            sys.stderr.write( WARNING_YELLOW+"Warning: specified output dir is current working dir or root dir! will clobber any parsnp.* results"+ENDC)
            outputDir = ""
        else:
            os.mkdirs("%s"%(outputDir), exist_ok=True)
    else:
        today = datetime.datetime.now()
        timestamp = "P_"+today.isoformat().replace("-","_").replace(".","").replace(":","").replace("T","_")
        outputDir = os.getcwd()+os.sep+timestamp
        os.mkdir("%s"%(outputDir))

    sortem = True
    if len(ref) == 0 and len(genbank_ref) != 0:
        #we are parsing from genbank, set ref to genbank_ref && turn off sorting
        ref = genbank_ref
        sortem = False

    autopick_ref = False
    if (not ref and not query) or not seqdir:
        sys.stderr.write( ERROR_RED+"ERROR: no seqs provided, yet required. exit!\n"+ENDC)
        sys.exit(0)
    elif not ref and query:
        sys.stderr.write( WARNING_YELLOW+"WARNING: no reference genome specified, going to autopick from %s as closest to %s\n"%(seqdir, query)+ENDC)
        autopick_ref = True
        ref = query

    if 1:
        print((len(outputDir)+17)*"*")
        print(BOLDME+"SETTINGS:"+ENDC)
        if ref != "!":
            print("|-"+BOLDME+"refgenome:\t%s"%(ref)+ENDC)
        else:
            print("|-"+BOLDME+"refgenome:\t%s"%("autopick")+ENDC)
        print("|-"+BOLDME+"aligner:\tlibMUSCLE"+ENDC)
        print("|-"+BOLDME+"seqdir:\t%s"%(seqdir)+ENDC)
        print("|-"+BOLDME+"outdir:\t%s"%(outputDir)+ENDC)
        print("|-"+BOLDME+"OS:\t\t%s"%(OSTYPE)+ENDC)
        print("|-"+BOLDME+"threads:\t%s"%(threads)+ENDC)
        print((len(outputDir)+17)*"*")

    print("\n<<Parsnp started>>\n")

    #1)read fasta files (contigs/scaffolds/finished/DBs/dirs)
    sys.stderr.write( "-->Reading Genome (asm, fasta) files from %s..\n"%(seqdir))
    files = []
    try:
        files1 = os.listdir(seqdir)
        files = []
        for f1 in files1:
            if not os.path.isdir(seqdir+os.sep+f1):
                files.append(f1)
            
        sys.stderr.write( "  |->["+OK_GREEN+"OK"+ENDC+"]\n")
    except IOError:
        sys.stderr.write( ERROR_RED+"ERROR: problem reading files from %s\n"%(seqdir)+ENDC)
        sys.exit(1)
    sys.stderr.write( "-->Reading Genbank file(s) for reference (.gbk) %s..\n"%(genbank_files_str))
    if len(genbank_file) == 0:
        sys.stderr.write( "  |->["+WARNING_YELLOW+"WARNING"+ENDC+"]"+": no genbank file provided for reference annotations, skipping..\n"+ENDC)
    elif not os.path.exists(genbank_file):
        sys.stderr.write( "  |->["+ERROR_RED+"ERROR"+ENDC+"]"+": provided genbank file does not exist, skipping..\n"+ENDC)
    else:
        sys.stderr.write( "  |->["+OK_GREEN+"OK"+ENDC+"]\n")

    fnafiles = []
    allfiles = []
    fnaf_sizes = {}

    allfile_dict = {}
    reflen = 0
    fnafiles1 = []
    for file in files:

       #any file in genome dir will be added..
       if file[0] != "." and file[-1] != "~" and len(file) > 1 and not os.path.isdir(seqdir+os.sep+file):
            ff = open(seqdir+os.sep+file,'r')
            hdr = ff.readline()
            seq = ff.readline()
            if len(seq) > 1 and ("A" in seq.upper() or "G" in seq.upper() or "C" in seq.upper() or "T" in seq.upper() or "N" in seq.upper()) and hdr[0] == ">":
                nameok = True
                for char in special_chars:
                    if char in file:
              
                        print("WARNING: File %s contains a non-supported special character (\'%s\') in file name. Please remove if you'd like to include. For best practices see: http://support.apple.com/en-us/HT202808"%(file,char))
                        nameok = False
                        break
            if nameok:
                fnafiles1.append(file)


    if ref == "!":
        ref = random.choice(fnafiles1)
        ref = seqdir+os.sep+ref
    if 1:
        #ADVAIT EDIT TO FIX HEADER BUG
        ff = open(ref,'r')
        hdr = ff.readline()
        if hdr[0] == ">":
            data = ff.read()
            if data.count("-") > 0:
                sys.stderr.write( "ERROR: ref genome sequence %s seems to aligned! remove and restart \n"%(ref))
                sys.exit(1)
            reflen = len(data) - data.count('\n')
        ff.close()

    for file in files:
       nameok = True

       for char in special_chars:
          if char in file:
              #print("WARNING: File %s contains a non-supported special character (%s) in file name. Please remove if you'd like to include. For best practices see: http://support.apple.com/en-us/HT202808"%(file,char)
              nameok = False
              break
       if not nameok:
           continue
       if file[0] != "." and file[-1] != "~":
            ff = open(seqdir+os.sep+file,'r')
            hdr = ff.readline()
            name_flag = True
            if hdr[0] == ">":
                data = []
                totlen = 0
                for line in ff:
                    if line[0] != ">":
                        data.append(line.replace("\n",""))
                        if "-" in line:
                            sys.stderr.write( "ERROR: genome sequence %s seems to aligned! remove and restart \n"%(file))
                            sys.exit(1)

                        totlen += len(line.replace("\n",""))

                if ref in file or file in ref:
                    reflen = totlen#len(data)
                    continue
                #sorry too small
                if totlen <= 20:
                    sys.stderr.write( "WARNING: File %s is less than or equal to 20bp in length. Skip!\n"%(file))
                    continue
                sizediff = float(reflen)/float(totlen) 
                
                #EDITED THIS TO CHANGE GENOME THRESHOLD, WILL NOW CONSIDER CONCATENATED GENOMES THAT ARE MUCH BIGGER THAN THE REFERENCE
                if not probe:
                    if sizediff <= 0.6 or sizediff >= 1.4:
                        sys.stderr.write( "WARNING: File %s is too short or too long compared to reference. Skip!\n"%(file))
                        continue
                else:
                    if sizediff >= 1.4:
                        sys.stderr.write( "WARNING: File %s is too short compared to reference genome. Skip!\n"%(file))
                        continue

                fnafiles.append(file)
                fnaf_sizes[file] = totlen#len(data)
            ff.close()
    
    #sys.exit(1) TEST SYS EXIT FOR CHECKING THRESHOLD
    if ref in fnafiles:
        sys.stderr.write( "ERROR: reference genome %s also in genome directory, restart and select different reference genome\n"%(ref))
        sys.exit(1)
        
    if ref == "!":
        fnafiles.remove(ref)
    
    #sort reference by largest replicon to smallest
    if sortem and os.path.exists(ref) and not autopick_ref:
        ff = open(ref,'r')
        seqs = ff.read().split(">")[1:]
        seq_dict = {}
        seq_len = {}
        for seq in seqs:
            hdr = ""
            nt = ""
            try:
                hdr,nt = seq.split("\n",1)
            except ValueError:
                continue
            seq_dict[hdr] = nt
            seq_len[hdr] = len(nt.replace("\n",""))
        seq_len_sort = sorted(iter(seq_len.items()), key=operator.itemgetter(1))
        seq_len_sort.reverse()
        ffo = open("%s"%(outputDir+os.sep+ref.split(os.sep)[-1]+".ref"),'w')
        for item in seq_len_sort:
            ffo.write(">%s\n"%(item[0]))
            ffo.write("%s"%(seq_dict[item[0]]))
        ff.close()
        ffo.close()
        ref = outputDir+os.sep+ref.split(os.sep)[-1]+".ref"
    else:
        ref = genbank_ref
    
    #remove any query sequences 30% diff in length
    allfiles = [ref.rsplit(os.sep,1)[-1]]
    #write INI file
    if xtrafast or 1:
        extend = False
    
    inifile1 = open("%s/template.ini"%(PARSNP_DIR),'r')
    inifiled = inifile1.read()
    inifiled = inifiled.replace("$REF",ref)
    inifiled = inifiled.replace("$EXTEND","%d"%(extend))
    inifiled = inifiled.replace("$ANCHORS",str(anchor))
    inifiled = inifiled.replace("$MUMS",str(mum))
    inifiled = inifiled.replace("$MINCLUSTER",str(mincluster))
    inifiled = inifiled.replace("$CLUSTERD",str(cluster))
    inifiled = inifiled.replace("$THREADS",str(threads))
    inifiled = inifiled.replace("$ALIGNER",str(aligner))
    inifiled = inifiled.replace("$DIAGDIFF",str(diagdiff))
    inifiled = inifiled.replace("$RECOMBFILT","%d"%(xtrafast))
    inifiled = inifiled.replace("$OUTDIR",outputDir)
    if fastmum:
        inifiled = inifiled.replace("$PARTPOS","%d"%(0.2*reflen))
    else:
        inifiled = inifiled.replace("$PARTPOS","%s"%(maxpartition))

    finalfiles = []
    #2)get near neighbors (mumi distance)
    if os.path.exists(outputDir+os.sep+"alltogether.fasta"):
        os.system("rm " + outputDir+os.sep+"alltogether.fasta")
    if os.path.exists(outputDir+os.sep+"blocks/b1"):
        ftrm = glob.glob(outputDir+os.sep+"blocks/b*")
        for file in ftrm:
            os.system("rm -rf "+file)
    #processed = []
    #fnafiles = processed
    fileidx = -1

    hit_dict = {}
    qry_hit_dict = {}
    hdr_dict = {}
    length_dict = {}
    
    TOTSEQS=len(fnafiles)+1
    seqids_list = []
    use_mummer_mumi = False
    use_parsnp_mumi = True
    if not inifile_exists:
        if len(fnafiles) < 1 or ref == "":
            sys.stderr.write( "Parsnp requires 2 or more genomes to run, exiting\n")
            print(fnafiles, end =' ') 
            print(ref)
            sys.exit(0)
    
        file_string = ""
        cnt = 1
        for file in fnafiles:
                
            file_string+="file%d=%s\n"%(cnt,seqdir+os.sep+file)
            file_string+="reverse%d=0\n"%(cnt)
            cnt +=1
        inifiled2 = inifiled.replace("$FILES\n",file_string)
        inifiled_mumi = inifiled2
        inifiled_mumi = inifiled_mumi.replace("calcmumi=0","calcmumi=1")
        inifile_mumi = open(outputDir+os.sep+"all_mumi.ini",'w')
        inifile_mumi.write(inifiled_mumi)
        inifile_mumi.close()
    mumi_dict = {}
    if use_parsnp_mumi and not curated:
        sys.stderr.write( "-->Calculating MUMi..\n")
        if not inifile_exists:
            command = "%s/parsnp %sall_mumi.ini"%(PARSNP_DIR,outputDir+os.sep)
        else:
            if not os.path.exists(inifile):
                sys.stderr.write( "Error: ini file %s does not exist!\n"%(inifile))
                sys.exit(1)
            command = "%s/parsnp %s"%(PARSNP_DIR,inifile.replace(".ini","_mumi.ini"))
        run_command(command)
        try:
            mumif = open(outputDir+os.sep+"all.mumi",'r')
            for line in mumif:
                line = line.replace("\n","")
                try:
                    idx,mi = line.split(":")
                    mumi_dict[int(idx)-1] = float(mi)
                except ValueError:
                    pass    
        except IOError:
            #mumi file generation failed, skip.. use all?
            i = 0
            for file in fnafiles:
                mumi_dict[i] = 1
        print("  |->["+OK_GREEN+"OK"+ENDC+"]")
    finalfiles = []
    lowest_mumi = 100
    auto_ref = ""

    if autopick_ref:
        for idx in list(mumi_dict.keys()):
            if mumi_dict[idx] < lowest_mumi:
                auto_ref = seqdir+os.sep+fnafiles[idx]
                ref = auto_ref        
                lowest_mumi = mumi_dict[idx]
    
    mumi_f = ""        
    if mumi_only and not curated:
        if not os.path.exists(outputDir):
            sys.stderr.write("Output directory does not exist! writing to cwd\n")
            outputDir = ""
            mumi_f = open("recruited_genomes.lst",'w')
        else:
            mumi_f = open(outputDir+os.sep+"recruited_genomes.lst",'w')
        
    if VERBOSE:
        print("RECRUITED GENOMES:\n")

    sorted_x = sorted(iter(mumi_dict.items()), key=operator.itemgetter(1))
    scnt = 0
    mumivals = []
    for item in sorted_x:
        if scnt > 100 or scnt >= len(sorted_x):
            break
        if float(item[1]) < float(mumidistance):
            mumivals.append(float(item[1]))
        scnt +=1
    minv=1.0
    if len(mumivals) > 0:
        minv = numpy.percentile(mumivals,0)
    
    dvals = mumivals


    stdv = 0
    hpv = 0
    if len(dvals) > 0:
        stdv = numpy.std(dvals)
        hpv = minv+(3*stdv)


    for idx in mumi_dict.keys():
        if mumi_dict[idx] < (float(mumidistance)) or curated:
            if fastmum and mumi_dict[idx] > hpv:
                continue
            if 1 or auto_ref != fnafiles[idx]:
                if mumi_only:
                    mumi_f.write(os.path.abspath(seqdir+os.sep+fnafiles[idx])+",%f"%(mumi_dict[idx])+"\n")
                if VERBOSE:
                    print("\t"+fnafiles[idx])
                finalfiles.append(fnafiles[idx])
                allfiles.append(fnafiles[idx])
    if VERBOSE:
        print("")

    if curated:
        for file in fnafiles:
            if file not in finalfiles:
                finalfiles.append(file)
            if file not in allfiles:
                allfiles.append(file)

    if mumi_only:
        mumi_f.close()
        sys.exit(1)

    orig_auto_ref = auto_ref
    if os.path.exists(auto_ref) and autopick_ref:
        ff = open(auto_ref,'r')
        seqs = ff.read().split(">")[1:]
        seq_dict = {}
        seq_len = {}
        for seq in seqs:
            hdr = ""
            nt = ""
            try:
                hdr,nt = seq.split("\n",1)
            except ValueError:
                continue
            seq_dict[hdr] = nt
            seq_len[hdr] = len(nt.replace("\n",""))
        seq_len_sort = sorted(seq_len.iteritems(), key=operator.itemgetter(1))
        seq_len_sort.reverse()
        ffo = open("%s"%(outputDir+os.sep+auto_ref.split(os.sep)[-1]+".ref"),'w')
        for item in seq_len_sort:
            ffo.write(">%s\n"%(item[0]))
            ffo.write("%s"%(seq_dict[item[0]]))
        ff.close()
        ffo.close()
        auto_ref = outputDir+os.sep+auto_ref.split(os.sep)[-1]+".ref"
        ref = auto_ref
        #print ref
    #print ref
    inifiled_closest = inifiled
    if not inifile_exists:
        if len(finalfiles) < 1 or ref == "":
            sys.stderr.write( "ERROR: Parsnp requires 2 or more genomes to run, exiting\n")
            sys.exit(0)
    
        file_string = ""
        cnt = 1
        file_string_closest = ""
        for file in finalfiles[0:1]:
            file_string_closest+="file%d=%s\n"%(cnt,seqdir+os.sep+file)
            file_string_closest+="reverse%d=0\n"%(cnt)
            cnt +=1
        cnt = 1
        for file in finalfiles:
            file_string+="file%d=%s\n"%(cnt,seqdir+os.sep+file)
            file_string+="reverse%d=0\n"%(cnt)
            cnt +=1
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
    print("-->Running Parsnp multi-MUM search and libMUSCLE aligner..")
    if not os.path.exists(outputDir+os.sep+"blocks"):
        os.mkdir(outputDir+os.sep+"blocks")
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
                    command = "%s/parsnp %sparsnpAligner.ini"%(PARSNP_DIR,outputDir+os.sep)
                else:
                    command = "%s/parsnp %spsnn.ini"%(PARSNP_DIR,outputDir+os.sep)
            else:
                if not os.path.exists(inifile):
                    sys.stderr.write("Error: ini file %s does not exist!\n"%(inifile))
                    sys.exit(1)
                command = "%s/parsnp %s"%(PARSNP_DIR,inifile)
            run_command(command)
        

            if not os.path.exists(outputDir+os.sep+"parsnpAligner.xmfa"):

                successful_run = False
                runcnt +=1
                if runcnt >= 2:
                    sys.stderr.write("Error: set of recruited genomes are too divergent for parsnp, please reduce MUMi (%f) and relaunch\n"%(float(mumidistance)))
                    sys.exit(1)                
            else:
                successful_run = True
                runcnt +=1
                break
        os.system("mv "+outputDir+os.sep+"parsnpAligner.xmfa "+outputDir+os.sep+"parsnp.xmfa")
    xmfafile = open(outputDir+os.sep+"parsnp.xmfa",'r')
        
    file2hdr_dict = {}
    fileid = ""
    blockfiles = []

    #get coverage
    coverage = 0
    totlength = 0
    totseqs = 0
    try:
        cf = open("%sparsnpAligner.log"%(outputDir+os.sep)) 
        for line in cf:
            if "Total coverage among all sequences:" in line:
                coverage = line.split(":",1)[-1].replace("\n","")
                coverage = float(coverage.replace("%",""))/100.0
            elif "Length:" in line:
                totlength += int(line.split(":",1)[-1].replace("\n","").split("bps")[0])
                totseqs +=1
    except IOError:
        print(ERROR_RED+"parsnpAligner.log missing, parsnpAligner failed, exiting.."+ENDC)
        sys.exit(1)

    #update thresholds
    if coverage <= 0.01:
        sys.stderr.write( "  |->["+ERROR_RED+"ERROR"+ENDC+"]"+": aligned regions cover less than 1% of reference genome, something is not right.. please adjust params and rerun. If problem persists please contact developers (treangen@gmail.com)"+ENDC)
        sys.exit(1)
    elif coverage < 0.1:
        sys.stderr.write( "  |->["+WARNING_YELLOW+"WARNING"+ENDC+"]"+": aligned regions cover less than 10% of reference genome! please verify recruited genomes are all strain of interest"+ENDC)
    else:
        pass
    print("  |->["+OK_GREEN+"OK"+ENDC+"]")
    t2 = time.time()
    elapsed = float(t2)-float(t1)
    #print("-->Getting list of LCBs.."
    allbfiles = glob.glob(outputDir+os.sep+"blocks/b*/*")
    blockfiles = []
    icnt = 0
    block_startpos = []
    block_dict = {}
    for file in allbfiles:
        if os.path.isfile(file):
            if "seq.fna" in file:
                blockfiles.append(file)
                lf = open(file,'r')
                header = lf.readline()
                if header[0] != ">":
                    sys.stderr.write( "Error with LCB: %s\n"%(file))
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
                block_dict[file] = [int(spos),int(epos), rseq]
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
    print("-->Running PhiPack on LCBs to detect recombination..")
    if run_recomb_filter and len(blockfiles) > 0:

        bedfile = open("%s/parsnp.rec"%(outputDir),'w')
        tasks = []
        processed = []
        icnt = 0
        for file in blockfiles:
            seq1 = ""
            try:
                bf = open(file,'r')
                seq1 = bf.read().split(">")[1].split("\n",1)[-1]
                seq1 = seq1.replace("\n","")
                bf.close()
            except IOError:
                pass    

            processed.append(file)
            params = {}
            path,file = file.rsplit(os.sep,1)
            params["jobID"] = len(tasks)
            params["query"] = "%s"%(file)
            params["seqlen"] = len(seq1)
            params["spos"] = block_startpos[icnt]
            params["dir"] = "%s"%(path)
            params["output"] = "%sProfile.csv"%(path+os.sep)#(path+os.sep+file+".out")
            tasks.append(params)
            icnt +=1
    
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
                    if VERBOSE:
                        sys.stderr.write( "File %s doesn't exist, no rec regions or error in PhiPack\n"%(tasks[i["jobID"]]["output"]))
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
                            chr_spos = list(ref_seqs.keys())
                            for cs in chr_spos:
                                if block_spos < chr_spos:
                                    chrnum = ref_seqs[cs]
                            bedfile_dict[srpos] = "%d\t%s\t%s\tREC\t%.3f\t+\n"%(chrnum,srpos,pos+50+block_spos,eval)
                             

                    
                qfile = tasks[i["jobID"]]["query"]

            elif i["status"] != 2:
                sys.stderr.write( "Error: parallel phipack job %d failed\n"%(i["jobID"]))
                raise IOError

        pool.close()
        pool.join()
        brkeys = list(bedfile_dict.keys())
        brkeys.sort()
        for key in brkeys:
            bedfile.write(bedfile_dict[key])
        bedfile.close()

    if run_recomb_filter:
        sys.stderr.write("  |->["+OK_GREEN+"OK"+ENDC+"]\n")
    else:        
        sys.stderr.write("  |->["+SKIP_GRAY+"SKIP"+ENDC+"]\n")
    run_lcb_trees = 0

    annotation_dict = {}
    if xtrafast or 1:
        #add genbank here, if present
        if len(genbank_ref) != 0:
            rnc = "%s/harvest -q -o %s/parsnp.ggr -x "%(PARSNP_DIR,outputDir)+outputDir+os.sep+"parsnp.xmfa"
            for file in genbank_files:
                rnc += " -g %s " %(file)
            run_command(rnc)
        else:
            run_command("%s/harvest -q -o %s/parsnp.ggr -f %s -x "%(PARSNP_DIR,outputDir,ref)+outputDir+os.sep+"parsnp.xmfa")

        if run_recomb_filter:
            run_command("%s/harvest -q -b %s/parsnp.rec,REC,\"PhiPack\" -o %s/parsnp.ggr -i %s/parsnp.ggr"%(PARSNP_DIR,outputDir,outputDir,outputDir))
        if run_repeat_filter:
            run_command("%s/harvest -q -b %s,REP,\"Intragenomic repeats > 100bp\" -o %s/parsnp.ggr -i %s/parsnp.ggr"%(PARSNP_DIR,repfile,outputDir,outputDir))

        run_command("%s/harvest -q -i %s/parsnp.ggr -S "%(PARSNP_DIR,outputDir)+outputDir+os.sep+"parsnp.snps.mblocks")

    command = "%s/ft -nt -quote -gamma -slow -boot 100 "%(PARSNP_DIR)+outputDir+os.sep+"parsnp.snps.mblocks > "+outputDir+os.sep+"parsnp.tree"
    print("-->Reconstructing core genome phylogeny..")
    run_command(command)
    #7)reroot to midpoint
    if os.path.exists("outtree"):
        os.system("rm outtree")
       
    if reroot_tree and len(finalfiles) > 1:
        #print("-->Midpoint reroot.."
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
            sys.stderr.write( "ERROR: cannot process fasttree output, skipping midpoint reroot..\n")
    print("  |->["+OK_GREEN+"OK"+ENDC+"]")


    if 1 or len(use_gingr) > 0:
        print("-->Creating Gingr input file..")
        if xtrafast or 1:
            #if newick available, add
            #new flag to update branch lengths
            run_command("%s/harvest --midpoint-reroot -u -q -i "%(PARSNP_DIR)+outputDir+os.sep+"parsnp.ggr -o "+outputDir+os.sep+"parsnp.ggr -n %s"%(outputDir+os.sep+"parsnp.tree "))
 
    print("  |->["+OK_GREEN+"OK"+ENDC+"]")

    print("-->Calculating wall clock time.. ")
    if float(elapsed)/float(60.0) > 60:
        print("  |->"+BOLDME+"Aligned %d genomes in %.2f hours"%(totseqs,float(elapsed)/float(3600.0))+ENDC)
    elif float(elapsed) > 60:
        print("  |->"+BOLDME+"Aligned %d genomes in %.2f minutes"%(totseqs,float(elapsed)/float(60.0))+ENDC)
    else: 
        print("  |->"+BOLDME+"Aligned %d genomes in %.2f seconds"%(totseqs,float(elapsed))+ENDC)
    #cleanup
    rmfiles = glob.glob(outputDir+os.sep+"*.aln")
    #rmfiles2 = glob.glob(outputDir+os.sep+"blocks/b*/*")
    rmfiles3 = glob.glob(outputDir+os.sep+"blocks/b*")
    for file in rmfiles:
        os.system("rm %s"%(file))
    for file in rmfiles3:
        os.system("rm -rf %s"%(file))

    filepres = 0
    print(BOLDME+"\n<<Parsnp finished! All output available in %s>>"%(outputDir)+ENDC)
    print("")
    print(BOLDME+"Validating output directory contents..."+ENDC)
    print(BOLDME+"\t1)parsnp.tree:\t\tnewick format tree"+ENDC, end =' ')
    if os.path.exists("%sparsnp.tree"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.tree"%(outputDir+os.sep)) > 0:
        print("\t\t\t["+OK_GREEN+"OK"+ENDC+"]")
        filepres+=1
    else:
        print("\t|->"+ERROR_RED+"MISSING"+ENDC)
    print(BOLDME+"\t2)parsnp.ggr:\t\tharvest input file for gingr (GUI)"+ENDC, end =' ')
    if os.path.exists("%sparsnp.ggr"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.ggr"%(outputDir+os.sep)) > 0:
        print("\t["+OK_GREEN+"OK"+ENDC+"]")
        filepres+=1
    else:
        print("\t|->"+ERROR_RED+"MISSING"+ENDC)
    print(BOLDME+"\t3)parsnp.xmfa:\t\tXMFA formatted multi-alignment"+ENDC, end = ' ')
    if os.path.exists("%sparsnp.xmfa"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.xmfa"%(outputDir+os.sep)) > 0:
        print("\t\t["+OK_GREEN+"OK"+ENDC+"]")
        filepres+=1
    else:
        print("\t|->"+ERROR_RED+"MISSING"+ENDC)
    if filepres == 3:
        pass

    else:
        print("\t\t["+ERROR_RED+"Output files missing, something went wrong. Check logs and relaunch or contact developers for assistance"+ENDC+"]")
    print("")
    if os.path.exists("%sblocks"%(outputDir+os.sep)):
        os.rmdir("%sblocks"%(outputDir+os.sep))
    if os.path.exists("allmums.out"):
        os.remove("allmums.out")

    if not VERBOSE and os.path.exists("parsnpAligner.ini"):
        os.remove("parsnpAligner.ini")

    prefix = outputDir+os.sep+ref.rsplit(".",1)[0].rsplit(os.sep)[-1]
    if not VERBOSE and os.path.exists("%s.coords"%(prefix)):
        os.remove("%s.coords"%(prefix))

    if not VERBOSE and os.path.exists("%s.delta"%(prefix)):
        os.remove("%s.delta"%(prefix))

    files = glob.glob("%s/*.reps"%(outputDir))
    for file in files:
        if not VERBOSE and os.path.exists(file):
            os.remove(file)


    files = glob.glob("%s/*.ref"%(outputDir))
    for file in files:
        if not VERBOSE and os.path.exists(file):
            os.remove(file)

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
            print(">>Launching gingr..")
            os.system("open -n %s --args %s/parsnp.ggr"%(use_gingr,outputDir))

