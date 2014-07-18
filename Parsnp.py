# See the LICENSE file included with this software for license information.

import os, sys, string, getopt, random,subprocess, time, glob,operator, math, datetime,numpy #pysam
import signal
import inspect
from multiprocessing import *

reroot_tree = True
try:
    import dendropy
except ImportError:
    reroot_tree = False


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

application_path = ""
if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(sys.executable)
elif __file__:
    application_path = os.path.dirname(__file__)    

#DONT NEED
#import gcb_new


if 1:
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
PHI_WINDOWSIZE = 1000
TOTSEQS=0


#template
#read input
#0) init
## Get start time                                                                                                                                                                                                                                                                  
t1 = time.time()
OSTYPE="linux"

p = subprocess.Popen("echo `uname`", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(checkStdout, checkStderr) = p.communicate()
if checkStderr != "":
    sys.stderr.write(WARNING_YELLOW+"Warning: Cannot determine OS, defaulting to %s\n"%(OSTYPE)+ENDC)
else:
    OSTYPE = checkStdout.strip()


binary_type = "linux"
if OSTYPE == "Darwin":
    binary_type = "osx"
else:
    binary_type = "linux"

def handler(signum, frame):
    global SIGINT
    SIGINT = True
    print 'Caught request to terminate by user (CTRL+C), exiting now, bye'
    sys.exit(128)

signal.signal(signal.SIGINT, handler)
def map_samples_to_reference(samples,ref,threads):
    command = "%s/bowtie2-build -o 2 %s "+outputDir+os.sep+"refidx"%(PARSNP_DIR,ref)
    run_command(command)
    os.mkdir(outputDir+os.sep+"reads") 
    for key in samples.keys():
        ispe = samples[key][0]
        isint = samples[key][1]
        min = samples[key][2]
        max = samples[key][3]
        isfq = samples[key][4]
        files = samples[key][-1]
        #command = "./bowtie2 -p %s --local --no-unal --quiet"%(threads)
        command = "%s/bowtie2 -p %s --no-unal --quiet"%(PARSNP_DIR,threads)
        if isfq:
            command+= " -q "
        else:
            command+= " -f "
        command += " -x "+outputDir+os.sep+"refidx "
        if not ispe:
            command += "-U,"
            for file in files[:-1]:
                command += file+","
            command += files[-1]
        else:
            #punt for now, all Unpaired..
            if isint:
                command += "-U,"
                for file in files[:-1]:
                    command += file+","
                command += files[-1]
            else:
                command += " -1 "+files[0]
                command += " -2 "+files[1]
        command += " -S "+outputDir+os.sep+"reads/%s.sam "%(key)
        run_command(command)
def run_freebayes(samfile,ref, xmfa_list):
    if not os.path.exists(samfile):
        sys.stderr.write( "SAM file %s does not exist!\n"%(samfile))
        sys.exit(1)
    currdir = os.getcwd()
    hdr = open(ref,'r').readline().split(" ")[0].split("\t")[0].strip()[1:]
    command = "%s/samtools view -S -b -h -o %s.bam %s"%(PARSNP_DIR,samfile,samfile)
    run_command(command)
    command = "%s/samtools sort %s.bam %s.srt"%(PARSNP_DIR,samfile,samfile)
    run_command(command)
    command = "%s/samtools index %s.srt.bam"%(PARSNP_DIR,samfile)
    run_command(command)
    command = "%s/samtools index %s"%(PARSNP_DIR,ref)
    run_command(command)
    command = "%s/samtools mpileup -uf %s %s.sam.srt.bam > %s.mpileup"%(PARSNP_DIR,ref,samfile,samfile)
    run_command(command)

    vcf_final = []
    first = 0
    for xmfa in xmfa_list:
        command = "%s/freebayes -j --report-all-haplotype-alleles --failed-alleles %s.vcf.failed -i -r \"%s\":%d..%d -C 0 -G 0 -F 0 -0 -E 0 -X -u -P 1 --fasta-reference %s -b %s.srt.bam -v %s.vcf"%(PARSNP_DIR,samfile,hdr,xmfa[0],xmfa[1],ref,samfile,samfile)
        run_command(command)
        vf = open("%s.vcf"%(samfile),'r')
        for line in vf.xreadlines():
            if line[0:2] == "##" or len (line) < 2:
                continue
            elif line[0:1] == "#" and not first:
                vcf_final.append(line)
                first = True
            elif "#" not in line:
                vcf_final.append(line)
        vf.close()
    ofv = open("%s.final.vcf"%(samfile),'w')
    for line in vcf_final:
        ofv.write(line)
    ofv.close()

def run_mummer(ref,query,prefix):
    command = "%s/mum64 -mum -F -L -c -b -l 20 %s %s  > "%(PARSNP_DIR,ref,query)+outputDir+os.sep+"%s.out"%(prefix)
    run_command(command)

def run_phipack(query,seqlen,workingdir):
    currdir = os.getcwd()
    os.chdir(workingdir)
    command = "%s/phiprofile -o -v -n %d -w 100 -m 100 -f %s > %s.out"%(PARSNP_DIR,seqlen,query,query)
    run_command(command,1)
    os.chdir(currdir)

def run_zorro(query,seqlen,workingdir):
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
    run_command(command,1)
    command = "%s/run_fasta"%()
    run_command(command,1)
    command = "%s/run_fasta"%()
    run_command(command,1)
    os.chdir(currdir)
    
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
        sys.stderr.write("Keyboard error in thread %d, quitting\n"%(jobID))
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
def revcomp(seq):
    seq = seq.replace("\n","")
    seq = string.upper(seq)
    rc = ""

    for char in seq:
        #print char                                                                                                                                                                                                    
        if char == "A":
            rc += "T"
        elif char == "T":
            rc += "A"
        elif char == "C":
            rc += "G"
        elif char == "G":
            rc += "C"
        elif char == "N":
            rc+= "N"
    return rc[::-1]

def findrepsref(ref,outdir):
    #run NUCMER on ref to find repeats and filter
    prefix = outdir+os.sep+ref.rsplit(".",1)[0].rsplit(os.sep)[-1]
    command = "%s/nucmer -p %s --maxmatch --nosimplify -l 30 %s %s"%(PARSNP_DIR,prefix,ref,ref)
    run_command(command)
    command  = "%s/show-coords -T -L 100 -I 90 %s.delta > %s.coords"%(PARSNP_DIR,prefix,prefix)
    run_command(command)
    of1 = open("%s.coords"%(prefix),'r')
    seq_dict = {}
    seq_len = {}
    ref1 = open(ref,'r')
    seqs = ref1.read().split(">")[1:]
    cnt = 1
    for seq in seqs:
        hdr,nt = seq.split("\n",1)
        hdr = hdr.split(" ",1)[0].strip()
        seq_dict[hdr] = cnt
        seq_len[hdr] = len(nt.replace("\n",""))
        cnt+=1
    ref1.close()
    bed1f = "%s.reps"%(prefix)
    bed1 = open(bed1f,'w')
    qryweight = {}
    for line in of1.xreadlines():
        data = line.replace("\n","").split("\t")
        if len(data) < 5:
            continue

        ref = data[-2].strip()
        idx = 0
        try:
            idx = seq_dict[ref]
        except KeyError:
            continue
        slen = 0
        try:
            slen = seq_len[ref]
        except KeyError:
            continue
        spos = int(data[0])
        epos = int(data[1])
        if float(epos-spos)/float(slen) >= 0.3:
            continue
        bed1.write("%d\t%d\t%d\n"%(idx,spos,epos))
    bed1.close()
    return bed1f
def layoutseq(ref,file,genomesdir,split):
    
    qryorder = []
    qryrev = {}
    prefix = file.rsplit(".",1)[0]
    #prefix = prefix.split("_")[0]
    if split:
        data = open("%s"%(genomesdir+os.sep+file),'r')
        of1 = open("%s.split"%(genomesdir+os.sep+file),'w')
        seqs = data.read().split(">")[1:]
        for seq in seqs:
            hdr,inf = seq.split("\n",1)
            inf = inf.upper()
            parts = inf.split("NNN")
            pcnt = 1
            for part in parts:
                of1.write(">"+hdr+"_%d\n"%(pcnt))
                of1.write(part)
        data.close()
        of1.close()
            
    if 1 or not os.path.exists("%s.coords"%(genomesdir+os.sep+prefix)):
        if split:
            #os.system("./nucmer -p %s -l 25 -c 60 %s %s.split"%(prefix,ref,genomesdir+os.sep+file))
            command = "%s/nucmer -p %s -l 25 -c 60 %s %s.split"%(PARSNP_DIR,prefix,ref,genomesdir+os.sep+file)
            run_command(command)
            
        else:
            #os.system("./nucmer -p %s -l 25 -c 60 %s %s"%(prefix,ref,genomesdir+os.sep+file))
            command = "%s/nucmer -p %s -l 25 -c 60 %s %s"%(PARSNP_DIR,prefix,ref,genomesdir+os.sep+file)
            run_command(command)
            #p = subprocess.Popen(command, shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,close_fds=True,executable="/bin/bash")
        
        #os.system("./show-coords -T -r -L 100 -I 95 %s.delta > %s.coords"%(prefix,prefix))
        command  = "%s/show-coords -T -r -L 100 -I 95 %s.delta > %s.coords"%(PARSNP_DIR,prefix,prefix)
        run_command(command)
        #continue                                                                                                                                                                                                                                                                        
    of1 = open("%s.coords"%(prefix),'r')
    qryweight = {}
    for line in of1.xreadlines():
        data = line.replace("\n","").split("\t")
        if len(data) < 5:
            continue

        data2 = data[-2:]#.split("\t",1)                                                                                                                                                                                                                                                 
        if len(data2) < 2:
            continue
        ref = data2[0].strip()
        qry = data2[1].strip()
        if "TAG" in qry:
            continue
        length = int(data[5])

        if qry not in qryorder:
            qryorder.append(qry)
            qryweight[qry] = length
        else:
            if length > qryweight[qry]:
                qryorder.remove(qry)
                qryorder.append(qry)
            else:
                pass
        data4 = []
        qryspos = int(data[2])
        qryepos = int(data[3])
        if qryspos > qryepos:
            qryrev[qry] = 1
        else:
            qryrev[qry] = 0
    qrydict = {}
    if split:
        data = open(genomesdir+os.sep+file+".split",'r').read()
    else:
        data = open(genomesdir+os.sep+file,'r').read()
    seqs = data.split(">")[1:]
    orig_hdr = {}
    idx = 0
    for seq in seqs:
        hdr,fasta = seq.split("\n",1)

        hdrn = hdr.split(" ")[0]
        orig_hdr[hdrn] = hdr
        try:
            ff = fasta
            if qryrev[hdrn] == 1:
                ff = revcomp(ff)
            qrydict[hdrn] = ff
        except KeyError:
            pass#print hdr                                                                                                                                                                                                                                                               
    #print qryorder
    of2 = open(genomesdir+os.sep+file+".sorted",'w')
    for hdr in qryorder:

        try:
            of2.write(">"+orig_hdr[hdr]+"\n")
            of2.write(qrydict[hdr])
        except KeyError:
            sys.exit(0)
            pass#print hdr                                                                                                                                                                                                                                                               
    of2.close()


#end def layout()

sys.stderr.write( BOLDME+"|--Parsnp v1.0--|\n"+ENDC)
sys.stderr.write( BOLDME+"For detailed documentation please see --> http://harvest.readthedocs.org/en/latest\n"+ENDC)


if not os.path.exists("%s/parsnp"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/parsnp %s/parsnp"%(PARSNP_DIR, PARSNP_DIR))
if not os.path.exists("%s/harvest"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/harvest_%s %s/harvest"%(PARSNP_DIR,binary_type,PARSNP_DIR))
if not os.path.exists("%s/ft"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/fasttree_%s %s/ft"%(PARSNP_DIR,binary_type,PARSNP_DIR))
if not os.path.exists("%s/phiprofile"%(PARSNP_DIR)):
    os.system("ln -s %s/bin/Profile_%s %s/phiprofile"%(PARSNP_DIR,binary_type,PARSNP_DIR))
if not os.path.exists("%s/nucmer"%(PARSNP_DIR)):
    os.system("ln -s %s/MUMmer/nucmer %s/nucmer"%(PARSNP_DIR,PARSNP_DIR))
if not os.path.exists("%s/show-coords"%(PARSNP_DIR)):
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

def usage():
    print "usage: parsnp [options] [-g|-r|-q](see below) -d <genome_dir> -p <threads>"
    print ""
    print "Parsnp quick start for three example scenarios: "
    print "1) With reference & genbank file: "
    print " >parsnp -g <reference_genbank_file1,reference_genbank_file2,..> -d <genome_dir> -p <threads> "
    print ""
    print "2) With reference but without genbank file:"
    print " >parsnp -r <reference_genome> -d <genome_dir> -p <threads> "
    print ""
    print "3) Autorecruit reference to a draft assembly:"
    print " >parsnp -q <draft_assembly> -d <genome_db> -p <threads> "
    print ""
    print "[Input parameters]"
    print "<<input/output>>"
    print " -c = <flag>: (c)urated genome directory, use all genomes in dir and ignore MUMi? (default = NO)"
    print " -d = <path>: (d)irectory containing genomes/contigs/scaffolds"
    print " -r = <path>: (r)eference genome (set to ! to pick random one from genome dir)"
    print " -g = <string>: Gen(b)ank file(s) (gbk), comma separated list (default = None)"
    print " -o = <string>: output directory? default [./P_CURRDATE_CURRTIME]"
    print " -q = <path>: (optional) specify (assembled) query genome to use, in addition to genomes found in genome dir (default = NONE)"
    print ""
    print "<<MUMi>>"
    print " -U = <float>: max (M)UMi distance (default: autocutoff based on distribution of MUMi values)"
    print " -M = <flag>: calculate MUMi and exit? overrides all other choices! (default: NO)"
    print ""   
    print "<<MUM search>>"
    print " -a = <int>: min (a)NCHOR length (default = 1.1*Log(S))"
    print " -C = <int>: maximal cluster D value? (default=100)"
    print " -z = <path>: min LCB si(z)e? (default = 25)"
    print ""
    print "<<LCB alignment>>"
    print " -D = <float>: maximal diagonal difference? Either percentage (e.g. 0.2) or bp (e.g. 100bp) (default = 0.12)"    
    print " -e = <flag> greedily extend LCBs? experimental! (default = NO)"
    print " -n = <string>: alignment program (default: libMUSCLE)"        
    print ""
    print "<<SNP filtration>>"
    print " -R = <flag>: disable (R)epeat filtering?"
    print " -x = <flag>: enable recombination filtering? (default: NO)"
    print ""
    print "<<Misc>>"
    print " -h = <flag>: (h)elp: print this message"
    print " -p = <int>: number of threads to use? (default= 1)"
    print " -P = <int>: max partition size? limits memory usage (default= 15000000)"
    print " -v = <flag>: (v)erbose output? (default = NO)"
    print ""

#hidden, not yet supported options


#print "-q = <path>: (optional) specify (assembled) query genome to use, in addition to genomes found in genome dir (default = NONE)"
#print "-s = <flag>: (s)plit genomes by n's (default = NO)"
#print "-t = <string>: gene anno(t)ations file (ptt)"
#print "-z = <path>: min cluster si(z)e? (default = 10)"
#print "options: annotate, stopafter, startafter, fq, fa"
#print "-F = <flag>: fast MUMi calc? (default=NO)"
#print "-f = <file/args>: (f)ile containing samples/reads to include in tree (not-assembled)"
#print "-g = <bool>: auto-launch (g)ingr? (default = NO)"
#print "-h: (h)elp?"
#print "-i = <string>: (use) existing parsnp (i)ni file"
#print "-l = <bool>: layout query seqs w.r.t reference genome"

if __name__ == "__main__":
    parsnp_dir= sys.path[0]
    #print parsnp_dir
    #PARSNP_DIR = parsnp_dir
    opts = []
    args = []
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hxved:C:F:D:i:t:g:m:MU:o:a:cln:p:P:q:r:R:sz:f:", ["help","xtrafast","verbose","extend","sequencedir","clusterD","DiagonalDiff","iniFile","geneAnnotations","genbank","mumlength","onlymumi","MUMi","outputDir","anchorlength","curated","layout","aligNmentprog","threads","max-partition-size","query","reference","nofiltreps","split","minclustersiZe","fastqlibs"])
    except getopt.GetoptError, err:
        # print help information and exit:                                                                                                                                                                                                        
        print str(err) # will print something like "option -a not recognized"                                                                                                                                                                     
        usage()
        sys.exit(2)
    if len(opts) < 2:
        usage()
        sys.exit(2)
    ref = ""
    currdir = os.getcwd()
    seqdir = "./genomes"
    fastq_samples = {}
    anchor = "1.1*(Log(S))"
    mum = "1.3*(Log(S))"
    maxpartition = 15000000
    fastmum = True
    cluster = "300"
    curated = False
    aligner = "2"
    threads = "32"
    mincluster = "21"
    diagdiff = "0.12"
    splitseq = False
    extend = False
    layout = False
    xtrafast = False
    inifile=""
    mumi_only = False
    mumidistance = 0.5
    genbank_file = ""
    genbank_files = []
    genbank_files_str = ""
    genbank_files_cat = ""
    genbank_ref = ""
    outputDir = ""
    pttfile = ""
    query = ""
    reflen = 0
    use_gingr = ""
    inifile_exists = False
    req_params = {}
    req_params["genbank"] = 0
    req_params["refgenome"] = 0
    req_params["genomedir"] = 0
    filtreps = True
    frozenbinary = True
    if not frozenbinary and not os.path.exists("./MUMmer/nucmer"):
        filtreps = False   
    repfile = ""
    multifasta = False
    ref_seqs = {}

    for o, a in opts:
        if o in ("-v","--verbose"):
            VERBOSE = True
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-R","--nofiltreps"):
            filtreps = False
        elif o in ("-r","--reference"):
            ref = a
            if a != "!":
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
                            ref_seqs[currpos+fastalen] = seqnum#seq.split("\n",1)[0]
                            currpos = currpos+fastalen
                            seqnum+=1
                        
                    rf.close()
                except IOError:
                    sys.stderr.write( "ERROR: Reference genome file %s not found\n"%(ref))
                    sys.exit(1)       
                req_params["refgenome"] = 1
                 
        elif o in ("-d","--sequenceDir"):
            seqdir = a
            if not os.path.exists(seqdir):
                sys.stderr.write( "ERROR: genome dir %s does not exist\n"%(seqdir))
                sys.exit(1)
            req_params["genomedir"] = 1
        elif o in ("-t","--geneAnnotations"):
            pttfile = a
            try:
                rf = open(pttfile,'r')
                rf.close()
            except IOError:
                sys.stderr.write( "ERROR: ptt file %s not found\n"%(pttfile))
                sys.exit(1)
        elif o in ("-q","--query"):
            query= a
            try:
                rf = open(query,'r')
                rf.close()
            except IOError:
                sys.stderr.write( "ERROR: optional Query file %s provided but not found\n"%(query))
                sys.exit(1)            
        elif o in ("-g","--genbank"):
            genbank_files_str = a
            genbank_files = genbank_files_str.split(",")
            ctcmd = "cat "
            
            first = True
            genbank_ref = ""
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
                        if len(giline) <= 2:
                            sys.stderr.write( "ERROR: Genbank file %s malformatted \n"%(genbank_file))
                            sys.exit(1)
                        genbank_ref1.write(">gi|"+giline.split("GI:")[-1])
                    ntdata = False
                    data = ""
                    for line in rf.xreadlines():
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
            req_params["genbank"] = 1

        elif o in ("-f","--fastqlibs"):
            #PE format: Name1,min,max,file1,file2#Name2,min,max,file1,file2;
            #SE format: Name1,file1#Name2,file1#Name3,file1
            isfile = False
            if len(a.split(",")) == 1:
                isfile = True
                
            if not isfile and "#" not in a:
                print "ERROR: fastq lib not properly formatted"
                print "PE format: Name1,min,max,fileY1,fileY2#Name2,min,max,fileZ1,fileZ2#"
                print "PE format (interleaved): Name1,min,max,fileY12#Name2,min,max,fileZ12#"
                print "SE format: Name1,format,fileY#Name2,format,fileZ#"
                sys.exit(1)
            samples = []
            if not isfile:
                samples = a.split("#")
                samples = samples[:-1]
            else:
                try:
                    sf = open(a,'r')
                except IOError:
                    print "Read info file not available! please verify %s"%(a)
                    sys.exit(1)
                samples = sf.readlines()
                sf.close()
            if 1:
                for sample in samples:
                    fsfiles = []
                    fqfiles = []
                    data = sample.split(",")
                    #print data
                    ispe = False
                    isint = True
                    min = 0
                    max = 0
                    isfq = True
                    if len(data) >= 4:
                        ispe = True
                    else:
                        ispe = False
                        if 'q' in data[1]:
                            isfq = True
                        else:
                            isfq = False
                    
                    if ispe and len(data) > 4:
                        isinit = False
                    if ispe:
                        fqfiles = data[3:]
                        min = int(data[1])
                        max = int(data[2])
                    else:
                        fqfiles = [data[-1]]
                    try:
                        fastq_samples[data[0]]
                        print "Error: duplicate library name, please rename then restart"
                        sys.exit(1)
                    except KeyError:
                        for fqfile in fqfiles:
                            
                            if fqfile not in os.listdir(os.getcwd()):
                                print "ERROR: fasta/q file %s does not exist (or cannot be accessed)"%(fqfile)
                                sys.exit(1)     
                        fastq_samples[data[0]] = [ispe,isint,min,max,isfq,fqfiles] 
            fastqdir = a
        elif o in ("-a","--anchorlength"):
            anchor = a
        elif o in ("-m","--mumlength"):
            mum = a
        elif o in ("-D","--DiagonalDiff"):
            diagdiff = a
        elif o in ("-o","--outputDir"):
            outputDir = a
        elif o in ("-e","--extend"):
            extend = True
        elif o in ("-M","--onlymumi"):
            mumi_only = True
        elif o in ("-U","--MUMi"):
            mumidistance = a
        elif o in ("-x","--xtrafast"):
            xtrafast = True
        elif o in ("-i","--iniFile"):
            inifile = a
            inifile_exists = True
        elif o in ("-c","--curated"):
            curated = True
        elif o in ("-n","--alignmentprog"):
            aligner = a
            if aligner == "muscle":
                aligner = "2"
            elif aligner == "mafft":
                aligner = "1"
            elif aligner == "fsa":
                aligner = "3"
            elif aligner == "prank":
                aligner = "4"
            else:
                aligner = "2"
        elif o in ("-p","--threads"):
            threads = a
        elif o in ("-P","--max-partition-size"):
            maxpartition = a
        elif o in ("-F","--fastmum"):
            fastmum = True
        elif o in ("-C","--clusterD"):
            cluster = a
        elif o in ("-s","--split"):
            splitseq = True
        elif o in ("-l","--layout"):
            layout = True
        elif o in ("-z","--minclustersize"):
            mincluster = a

    
    if outputDir != "":
        today = datetime.datetime.now()
        timestamp = "P_"+today.isoformat().replace("-","_").replace(".","").replace(":","").replace("T","_")
        outputDir2 = timestamp#outputDir+os.sep+timestamp
        if outputDir == "." or outputDir == "./" or outputDir == "/":
            sys.stderr.write( WARNING_YELLOW+"Warning: specified output dir is current working dir! will clobber any parsnp.* results"+ENDC)
            outputDir = ""

        elif os.path.exists("%s"%(outputDir)):
            pass
        else:
            os.mkdir("%s"%(outputDir))
    else:
        today = datetime.datetime.now()
        timestamp = "P_"+today.isoformat().replace("-","_").replace(".","").replace(":","").replace("T","_")
        outputDir = os.getcwd()+os.sep+timestamp
        os.mkdir("%s"%(outputDir))


    if len(ref) == 0 and len(genbank_ref) != 0:
        #we are parsing from genbank, set ref to genbank_ref
        ref = genbank_ref

    print (len(outputDir)+17)*"*"
    print BOLDME+"SETTINGS:"+ENDC
    print "|-"+BOLDME+"aligner:\tlibMUSCLE"+ENDC
    print "|-"+BOLDME+"seqdir:\t%s"%(seqdir)+ENDC
    print "|-"+BOLDME+"outdir:\t%s"%(outputDir)+ENDC
    print "|-"+BOLDME+"OS:\t\t%s"%(OSTYPE)+ENDC
    print "|-"+BOLDME+"threads:\t%s"%(threads)+ENDC
    print (len(outputDir)+17)*"*"

    autopick_ref = False
    if (len(ref) == 0 and len(query) == 0) or len(seqdir) == "":
        sys.stderr.write( ERROR_RED+"ERROR: no seqs provided, yet required. exit!\n"+ENDC)
        sys.exit(0)
    elif len(ref) == 0 and len(query) != 0:
        sys.stderr.write( WARNING_YELLOW+"WARNING: no reference genome specified, going to autopick from %s as closest to %s\n"%(seqdir, query)+ENDC)
        autopick_ref = True
        ref = query
    #1)read fasta files (contigs/scaffolds/finished/DBs/dirs)
    sys.stderr.write( "-->Reading Genome (asm, fasta) files from %s..\n"%(seqdir))
    files = []
    try:
        files = os.listdir(seqdir)
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

       if file[-3:] == ".fa" or file[-4:] == ".fna" or file[-6:] == ".fasta" or file[-4:] == ".fas" or file[-6:] == ".scafs" or file[-5:] == ".scfs" or file[-8:] == ".contigs" or file[-4:] == ".scf" or file[-4:] == ".ctg" or file[-5:] == ".ctgs":

            ff = open(seqdir+os.sep+file,'r')
            hdr = ff.readline()
            if hdr[0] == ">":
                fnafiles1.append(file)


    if ref == "!":
        ref = random.choice(fnafiles1)
        ref = seqdir+os.sep+ref
    else:
        ff = open(ref,'r')
        hdr = ff.readline()
        if hdr[0] == ">":
    
            data = ff.read()
            data = data.replace("\n","")
            reflen = len(data)
        ff.close()
    for file in files:

       if file[-3:] == ".fa" or file[-4:] == ".fna" or file[-6:] == ".fasta" or file[-4:] == ".fas" or file[-6:] == ".scafs" or file[-5:] == ".scfs" or file[-8:] == ".contigs" or file[-4:] == ".scf" or file[-4:] == ".ctg" or file[-5:] == ".ctgs":
    
            ff = open(seqdir+os.sep+file,'r')
            hdr = ff.readline()
            if hdr[0] == ">":

                data = []
                totlen = 0
                for line in ff.xreadlines():
                    if line[0] != ">":
                        data.append(line.replace("\n",""))
                        totlen += len(line.replace("\n",""))

                if ref in file or file in ref:
                    reflen = totlen#len(data)
                    continue
                sizediff = float(reflen)/float(totlen) 
                if sizediff <= 0.5 or sizediff >= 1.5:
                    continue

                fnafiles.append(file)
                fnaf_sizes[file] = totlen#len(data)
            ff.close()
    
    
    if ref in fnafiles:
        sys.stderr.write( "ERROR: reference genome %s also in genome directory, restart and select different reference genome\n"%(ref))
        sys.exit(1)
        
    if ref == "!":
        fnafiles.remove(ref)

    #sort reference by largest replicon to smallest
    if os.path.exists(ref) and not autopick_ref:
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
        seq_len_sort = sorted(seq_len.iteritems(), key=operator.itemgetter(1))
        seq_len_sort.reverse()
        ffo = open("%s"%(outputDir+os.sep+ref.split(os.sep)[-1]+".srt"),'w')
        for item in seq_len_sort:
            ffo.write(">%s\n"%(item[0]))
            ffo.write("%s"%(seq_dict[item[0]]))
        ff.close()
        ffo.close()
        ref = outputDir+os.sep+ref.split(os.sep)[-1]+".srt"

    #remove any query sequences 30% diff in length
    allfiles = [ref.rsplit(os.sep,1)[-1]]
    #write INI file
    if xtrafast or 1:
        extend = False
        
    inifile1 = open("%s/template.ini"%(PARSNP_DIR),'r')
    inifiled = inifile1.read()
    inifiled = inifiled.replace("$REF",ref)
    inifiled = inifiled.replace("$EXTEND","%d"%(extend))
    inifiled = inifiled.replace("$ANCHORS",anchor)
    inifiled = inifiled.replace("$MUMS",mum)
    inifiled = inifiled.replace("$MINCLUSTER",mincluster)
    inifiled = inifiled.replace("$CLUSTERD",cluster)
    inifiled = inifiled.replace("$THREADS",threads)
    inifiled = inifiled.replace("$ALIGNER",aligner)
    inifiled = inifiled.replace("$DIAGDIFF",diagdiff)
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
        os.system("rm -rf "+outputDir+os.sep+"blocks/b*")
    processed = []
    #initiate parallelMummer tasks
    tasks = []
    refg = open(ref,'r')
    #seqf_data = seqf.readlines()
    if 0:
        oseqf = open(ref+".single",'w')
        oseqf.write(refg.readline())
        for line in refg.xreadlines():
           if line[0:1] != ">":
              oseqf.write(line)
        oseqf.close()
    refg.close()
    for file in fnafiles:
        if (file[-3:] == ".fa" or file[-4:] == ".fna" or file[-6:] == ".fasta" or file[-4:] == ".fas" or file[-6:] == ".scafs" or file[-5:] == ".scfs" or file[-8:] == ".contigs" or file[-4:] == ".scf" or file[-4:] == ".ctg" or file[-5:] == ".ctgs") and file not in processed:
            processed.append(file)
            params = {}
            params["jobID"] = len(tasks)
            params["ref"] = "%s"%(ref)#+".single")
            params["query"] = "%s"%(seqdir+os.sep+file+".single")
            params["queryfile"] = file
            params["prefix"] = "%s"%(file+".mummerout")
            params["output"] = outputDir+os.sep+"%s"%(file+".mummerout.out")
            tasks.append(params)
            

    
    fnafiles = processed
    fileidx = -1
    #print len(seqs)
    hit_dict = {}
    qry_hit_dict = {}
    hdr_dict = {}
    length_dict = {}
    
    TOTSEQS=len(fnafiles)+1
    seqids_list = []
    use_mummer_mumi = False#True#False#False
    use_parsnp_mumi = True#False#True
    if not inifile_exists:
        if len(fnafiles) < 1 or ref == "":
            sys.stderr.write( "Parsnp requires 2 or more genomes to run, exiting\n")
            print fnafiles, ref
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
            for line in mumif.xreadlines():
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
        print "  |->["+OK_GREEN+"OK"+ENDC+"]"
    finalfiles = []
    lowest_mumi = 100
    auto_ref = ""

    if autopick_ref:
        for idx in mumi_dict.keys():
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
        print "RECRUITED GENOMES:\n"

    sorted_x = sorted(mumi_dict.iteritems(), key=operator.itemgetter(1))
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
                    print "\t"+fnafiles[idx]
                finalfiles.append(fnafiles[idx])
                allfiles.append(fnafiles[idx])
    if VERBOSE:
        print 

    if curated:
        for file in fnafiles:
            finalfiles.append(file)
            allfiles.append(file)

    if mumi_only:
        mumi_f.close()
        sys.exit(1)

    if use_mummer_mumi:
        print "-->Running MUMmer.."
        pool = Pool(processes=int(threads))
        result = pool.map_async(parallelWrapper,tasks).get(sys.maxint)
        for i in result:
           if (i["status"] == 1):
               #process output
              seqf = open(tasks[i["jobID"]]["output"],'r')
              seq = seqf.read()
              seqf.close()
              os.remove(tasks[i["jobID"]]["output"])
              qfile = tasks[i["jobID"]]["queryfile"]
              pref,data = seq.split(">",1)
              hdr,datax = data.split("\n",1)
              seqid, length = hdr.split("Len =")
              length = length.split("\n")[0]
              seqid = seqid.split("Reverse")[0]
              seqid = seqid.strip()
              seqid = seqid.replace(" ","")

              nseqid = qfile+"_"+seqid

        
              try:
                hdr_dict[nseqid] 
              except KeyError:
                fileidx +=1
                hdr_dict[nseqid] = fileidx
                seqids_list.append(nseqid)
              length = int(length)
              length_dict[nseqid] = length
              try:
                qry_hit_dict[nseqid]
              except KeyError:
                qry_hit_dict[nseqid] = {}
      
              hitlens = 0
              lines = datax.split("\n")
              for line in lines:
                if len(line) < 1:
                    continue
                if line[0] == ">":
                    continue
                line = line.replace("\n","")
                data = line.split(" ")
                datan = []
                for item in data:
                    if len(item) > 0:
                        datan.append(item)
                if len(datan) < 4:
                    continue
                hitsp = int(datan[1])
                hitlens = int(datan[3])
                z = 0
                while z < hitlens:
                   qry_hit_dict[nseqid][z+hitsp] = 1
                   z+=1
           elif i["status"] != 2:
              sys.stderr.write( "Error: parallel mummer job %d failed\n"%(i["jobID"]))
              raise IOError
        pool.close()
        pool.join()
        tasks = []
        genome_mumi_values = []
        TOTSEQS=1
        
        for key in qry_hit_dict.keys():
            minsize = reflen#min(length_dict[key],reflen)
            if minsize == 0:
               sys.stderr.write( "Something is wrong, reference genome size 0bp\n")
               sys.exit(1)
            mumi = float(len(qry_hit_dict[key].values()))/float(minsize)#float(hit_dict[key])/float(length_dict[key])
            mumi_f = float("%.1f"%(mumi))
            print mumi_f, mumi, key, minsize, len(qry_hit_dict[key].values())
            if mumi_f >= 0.5 and mumi not in genome_mumi_values:# and mumi_f != 1.0:#mumi not in genome_mumi_values and mumi_f != 1.0:
                genome_mumi_values.append(mumi)
                if layout:
                    layoutseq(ref,fnafiles[hdr_dict[key]],seqdir,splitseq)
                    if fnafiles[hdr_dict[key]]+".sorted" not in finalfiles:
                        finalfiles.append(fnafiles[hdr_dict[key]]+".sorted")
                        allfiles.append(fnafiles[hdr_dict[key]]+".sorted")
                        TOTSEQS+=1
                else:
                    if fnafiles[hdr_dict[key]] not in finalfiles:
                        finalfiles.append(fnafiles[hdr_dict[key]])
                        allfiles.append(fnafiles[hdr_dict[key]])
                        TOTSEQS+=1
        qry_hit_dict = {}

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
        ffo = open("%s"%(outputDir+os.sep+auto_ref.split(os.sep)[-1]+".srt"),'w')
        for item in seq_len_sort:
            ffo.write(">%s\n"%(item[0]))
            ffo.write("%s"%(seq_dict[item[0]]))
        ff.close()
        ffo.close()
        auto_ref = outputDir+os.sep+auto_ref.split(os.sep)[-1]+".srt"
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
    print "-->Running Parsnp multi-MUM search and libMUSCLE aligner.."
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
        for line in cf.xreadlines():
            if "Total coverage among all sequences:" in line:
                coverage = line.split(":",1)[-1].replace("\n","")
                coverage = float(coverage.replace("%",""))/100.0
            elif "Length:" in line:
                totlength += int(line.split(":",1)[-1].replace("\n","").split("bps")[0])
                totseqs +=1
    except IOError:
        print ERROR_RED+"parsnpAligner.log missing, parsnpAligner failed, exiting.."+ENDC
        sys.exit(1)

    if coverage < 0.1:
        sys.stderr.write( "  |->["+ERROR_RED+"ERROR"+ENDC+"]"+": aligned regions cover less than 10% of reference genome, something is not right.. please adjust params and rerun. If problem persists please contact developers (treangen@gmail.com)"+ENDC)
        sys.exit(1)
    elif coverage < 0.3:
        sys.stderr.write( "  |->["+WARNING_YELLOW+"WARNING"+ENDC+"]"+": aligned regions cover less than 30% of reference genome! please verify recruited genomes are all strain of interest"+ENDC)
    else:
        pass
    print "  |->["+OK_GREEN+"OK"+ENDC+"]"
    t2 = time.time()
    elapsed = float(t2)-float(t1)
    #print "-->Getting list of LCBs.."
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
    print "-->Determining repetitive regions.."

    run_repeat_filter = filtreps
    if run_repeat_filter:
        repfile = findrepsref(ref,"%s"%(outputDir))
        if os.path.exists("%s.delta"%(os.getcwd()+os.sep+ref.split(os.sep)[-1])):
            os.system("rm %s.delta"%(os.getcwd()+os.sep+ref.split(os.sep)[-1]))
        if os.path.exists("%s.coords"%(os.getcwd()+os.sep+ref.split(os.sep)[-1])):
            os.system("rm %s.coords"%(os.getcwd()+os.sep+ref.split(os.sep)[-1]))

    if run_repeat_filter and len(repfile) > 1:
        sys.stderr.write("  |->["+OK_GREEN+"OK"+ENDC+"]\n")
    if run_repeat_filter and len(repfile) <= 1:
        sys.stderr.write("  |->["+ERROR_RED+"ERROR"+ENDC+"]\n")
    elif not run_repeat_filter:        
        sys.stderr.write("  |->["+SKIP_GRAY+"SKIP"+ENDC+"]\n")
    #initiate parallelPhiPack tasks
    run_recomb_filter = 0
   
    if xtrafast:
        run_recomb_filter = 1
    else:
        run_recomb_filter = 0

    recombination_sites = {}
    bedfile = ""
    bedfile_dict = {}
    print "-->Running PhiPack on LCBs to detect recombination.."
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
        result = pool.map_async(parallelPhiWrapper,tasks).get(sys.maxint)
        
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
                            chr_spos = ref_seqs.keys()
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
        brkeys = bedfile_dict.keys()
        brkeys.sort()
        for key in brkeys:
            bedfile.write(bedfile_dict[key])
        bedfile.close()

    if run_recomb_filter:
        sys.stderr.write("  |->["+OK_GREEN+"OK"+ENDC+"]\n")
    else:        
        sys.stderr.write("  |->["+SKIP_GRAY+"SKIP"+ENDC+"]\n")
    run_lcb_trees = 0
    if run_lcb_trees:
        print "-->Running FastTree on individual LCBs.."
    lcb_trees = []
    bucky_trees =[]
    if run_lcb_trees and len(blockfiles) > 0:
        tasks = []
        processed = []
        icnt = 0
        for file in blockfiles:
            processed.append(file)
            params = {}
            path,file = file.rsplit(os.sep,1)
            params["jobID"] = len(tasks)
            params["query"] = "%s"%(file)
            params["spos"] = block_startpos[icnt]
            params["dir"] = "%s"%(path)
            params["recombination"] = recombination_sites
            params["output"] = "%sout.tree"%(path+os.sep)#(path+os.sep+file+".out")
            params["xmfa"] = "%sseq.fna"%(path+os.sep)#(path+os.sep+file+".out")
            
            tasks.append(params)
            icnt +=1
    
        #run parallelPhiPack
        pool = Pool(processes=int(threads))
        result = pool.map_async(parallelFtWrapper,tasks).get(sys.maxint)
        
        for i in result:
            if (i["status"] == 1):
                #process output
                if os.path.exists(tasks[i["jobID"]]["output"]):
                    lcb_trees.append(tasks[i["jobID"]]["output"])
                    lcbt = open(tasks[i["jobID"]]["output"],'r')
                    newicktree = lcbt.readline().replace("\n","")#.replace("\'","")
                    hdrf = open(tasks[i["jobID"]]["xmfa"],'r')
                    hdrd_dict = {}
                    for xline in hdrf.xreadlines():
                        if ">" in xline:
                            hdrd_dict[xline[1:].replace("\n","")] = xline[1:].split(":",1)[0]
                    for hdrd in hdrd_dict.keys():
                        newicktree = newicktree.replace("\'"+hdrd+"\'",hdrd_dict[hdrd])
                    """
                    nt2 = newicktree.split(",")
                    newtree = ""
                    
                    for node in nt2:
                        ex1 = ""
                        if ")" in node[:-3]:# and ")" not in node[-3:]:
                            ex1 = node.rsplit(")",1)[0].rsplit(":",1)
                        else:
                            ex1 = node.rsplit(":",1)
                        newtree += ex1[0].split(":",1)[0].replace("\'","")#+":"+ex1[-1]+","
                        #for exn in ex1[1:-1]:
                        #    if ")" in exn:
                        #        newtree+=")"
                        #    elif "(" in exn:
                        #        newtree+="("
                        #newtree+=ex1[-1]
                        newtree+=":"+ex1[-1]
                        if ")" in node[:-3]:
                            newtree+=")"
                        newtree+= ","
                    #newtree+=";"
                    
                    newtree = newtree.replace(",;",";")
                    newtree = newtree.replace(";,",";")
                    if newtree.count("(") <> newtree.count(")"):
                        newtree = newtree.replace(");","));")
                    """
                    buckyf = open(tasks[i["jobID"]]["output"].replace(".tree",".t"),'w')
                    
                    buckyf.write("#NEXUS\n")
                    buckyf.write("[ID: 123456789]\n")
                    buckyf.write("begin trees;\n")
                    buckyf.write("   translate\n")
                    lookup = {}
                    cnt = 1
                    
                    for file in allfiles:
                        lookup["\'"+file+"\'"] = cnt
                        if cnt == len(allfiles):
                            buckyf.write((8-len(`cnt`))*" "+`cnt`+" "+file+";\n")
                        else:
                            buckyf.write((8-len(`cnt`))*" "+`cnt`+" "+file+",\n")
                        cnt +=1    
                    buckyf.write("   tree rep.1 = ")

                    buckyf.write(newicktree+"\n")#+" 1\n")
                    
                    buckyf.close()
                    #call mbsum to properly format
                    os.system("%s/mbsum %s -o %s"%(PARSNP_DIR,tasks[i["jobID"]]["output"].replace(".tree",".t"),tasks[i["jobID"]]["output"].replace(".tree",".in")))
                    if len(open(tasks[i["jobID"]]["output"].replace(".tree",".in"),'r').read()) > 1:
                        bucky_trees.append(tasks[i["jobID"]]["output"].replace(".tree",".in"))
                    else:
                        print "mbsum failed.."
                else:
                    sys.stderr.write( "Error: parallel fasttree job %d failed\n"%(i["jobID"]))
                    raise IOError

            else:
                sys.stderr.write( "Error: parallel fasttree job %d failed\n"%(i["jobID"]))
                raise IOError
        recfile.close()
        pool.close()
        pool.join()
            
    #run bucky
    use_bucky_concordance = False

    if len(pttfile) > 0:                    
        print "-->Parsing reference annotations.. (ptt)"
    annotation_dict = {}
    if len(pttfile) > 0:
        
        pttfilex = open(pttfile,'r')
        pttfilex.readline()
        pttfilex.readline()
        pttfilex.readline()
        for line in pttfilex.xreadlines():
            line = line.replace("\n","").strip()
            spos,epos = line.split("\t")[0].split("..")
            description = line.split("\t")[-1].replace("\n","")
            gene = line.split("\t")[4]
            if gene == "-":
                gene = line.split("\t")[5]
            
            pos = 0
            spos = int(spos)
            epos = int(epos)
            while pos+spos<epos:
                annotation_dict[pos+epos] = [gene,description]
                pos+=1
        pttfilex.close()
 
    if xtrafast or 1:
        run_command("%s/harvest -q -o %s/parsnp.ggr -f %s -x "%(PARSNP_DIR,outputDir,ref)+outputDir+os.sep+"parsnp.xmfa")

        if run_recomb_filter:
            os.system("%s/harvest -q -b %s/parsnp.rec,REC,\"PhiPack\" -o %s/parsnp.rec.ggr -i %s/parsnp.ggr"%(PARSNP_DIR,outputDir,outputDir,outputDir))
            os.system("mv %s/parsnp.rec.ggr %s/parsnp.ggr"%(outputDir,outputDir))
        if run_repeat_filter:
            run_command("%s/harvest -q -b %s,REP,\"Intragenomic repeats > 100bp\" -o %s/parsnp.rep.ggr -i %s/parsnp.ggr"%(PARSNP_DIR,repfile,outputDir,outputDir))

            os.system("mv %s/parsnp.rep.ggr %s/parsnp.ggr"%(outputDir,outputDir))
        os.system("%s/harvest -q -i %s/parsnp.ggr -S "%(PARSNP_DIR,outputDir)+outputDir+os.sep+"parsnp.snps.mblocks")
        os.system("rm %s/parsnp.ggr"%(outputDir))
    print "  |->["+OK_GREEN+"OK"+ENDC+"]"
    command = "%s/ft -nt -quote -gamma -slow -boot 100 "%(PARSNP_DIR)+outputDir+os.sep+"parsnp.snps.mblocks > "+outputDir+os.sep+"parsnp.tree"
    print "-->Reconstructing core genome phylogeny.."
    run_command(command)
    #7)reroot to midpoint
    if os.path.exists("outtree"):
        os.system("rm outtree")
       
    if reroot_tree and len(finalfiles) > 1:
        #print "-->Midpoint reroot.."
        try:
            mtree = open("%sparsnp.tree"%(outputDir+os.sep), 'r')
            mtreedata = mtree.read()
            mtreedata.replace("\n","")
            tree = dendropy.Tree.get_from_string(mtreedata,"newick")
            tree.reroot_at_midpoint(update_splits=False)
            mftreef = tree.as_string('newick').split(" ",1)[1]
            #print mftreef
            mtreef = open(outputDir+os.sep+"parsnp.final.tree",'w')
            mtreef.write(mftreef)
            mtreef.close()
            os.system("mv %s %s"%(outputDir+os.sep+"parsnp.final.tree",outputDir+os.sep+"parsnp.tree"))
        except IOError:
            sys.stderr.write( "ERROR: cannot process fasttree output, skipping midpoint reroot..\n")
    print "  |->["+OK_GREEN+"OK"+ENDC+"]"


    if 1 or len(use_gingr) > 0:
        print "-->Creating Gingr input file.."
        if xtrafast or 1:
            if len(genbank_files) > 0:
                if run_recomb_filter:
                    os.system("%s/harvest -q -g %s -o "%(PARSNP_DIR,genbank_files_cat)+outputDir+os.sep+"parsnp.ggr -f %s -n "%(ref)+outputDir+os.sep+"parsnp.tree -x "+outputDir+os.sep+"parsnp.xmfa -V "+outputDir+os.sep+"parsnp.vcf -b "+outputDir+os.sep+"parsnp.rec,REC,\"PhiPack\"")
                    if run_repeat_filter:
                        os.system("%s/harvest -q -b %s,REP,\"SNP within intragenomic repeat > 100bp\" -o %s/parsnp.rep.ggr -i %s/parsnp.ggr"%(PARSNP_DIR,repfile,outputDir,outputDir))
                        os.system("mv %s/parsnp.rep.ggr %s/parsnp.ggr"%(outputDir,outputDir))
                elif run_repeat_filter:
                    os.system("%s/harvest -q -g %s -o "%(PARSNP_DIR,genbank_files_cat)+outputDir+os.sep+"parsnp.rep.ggr -f %s -n "%(ref)+outputDir+os.sep+"parsnp.tree -x "+outputDir+os.sep+"parsnp.xmfa -V "+outputDir+os.sep+"parsnp.vcf -b %s,REP,\"SNP within intragenomic repeat > 100bp\""%(repfile))

                    os.system("mv %s/parsnp.rep.ggr %s/parsnp.ggr"%(outputDir,outputDir))
                else:
                    os.system("%s/harvest -q -g %s -o "%(PARSNP_DIR,genbank_files_cat)+outputDir+os.sep+"parsnp.ggr -f %s -n "%(ref)+outputDir+os.sep+"parsnp.tree -x "+outputDir+os.sep+"parsnp.xmfa -V "+outputDir+os.sep+"parsnp.vcf")
            else:
                if run_recomb_filter:
                    os.system("%s/harvest -q -o "%(PARSNP_DIR)+outputDir+os.sep+"parsnp.ggr -f %s -n "%(ref)+outputDir+os.sep+"parsnp.tree -x "+outputDir+os.sep+"parsnp.xmfa -V "+outputDir+os.sep+"parsnp.vcf -b "+outputDir+os.sep+"parsnp.rec,REC,\"PhiPack\"")
                    if run_repeat_filter:
                        os.system("%s/harvest -q -b %s,REP,\"SNP within intragenomic repeat > 100bp\" -o %s/parsnp.rep.ggr -i %s/parsnp.ggr"%(PARSNP_DIR,repfile,outputDir,outputDir))
                        os.system("mv %s/parsnp.rep.ggr %s/parsnp.ggr"%(outputDir,outputDir))
                elif run_repeat_filter:
                    run_command("%s/harvest -q -o "%(PARSNP_DIR)+outputDir+os.sep+"parsnp.ggr -f %s -n "%(ref)+outputDir+os.sep+"parsnp.tree -x "+outputDir+os.sep+"parsnp.xmfa -V "+outputDir+os.sep+"parsnp.vcf -b %s,REP,\"SNP within intragenomic repeats (100bp)\""%(repfile))
                else:
                    os.system("%s/harvest -q -o "%(PARSNP_DIR)+outputDir+os.sep+"parsnp.ggr -f %s -n "%(ref)+outputDir+os.sep+"parsnp.tree -x "+outputDir+os.sep+"parsnp.xmfa -V "+outputDir+os.sep+"parsnp.vcf")
        else:
            os.system("%s/harvest -q -o "%(PARSNP_DIR)+os.sep+outputDir+os.sep+"parsnp.ggr -f %s -n "%(ref)+outputDir+os.sep+"parsnp.tree -v "+outputDir+os.sep+"parsnp.snps.vcf -x "+outputDir+os.sep+"parsnp.xmfa "%(PARSNP_DIR,ref))

 
    print "  |->["+OK_GREEN+"OK"+ENDC+"]"

    print "-->Calculating wall clock time.. "
    if float(elapsed)/float(60.0) > 60:
        print "  |->"+BOLDME+"Aligned %d genomes in %.2f hours"%(totseqs,float(elapsed)/float(3600.0))+ENDC
    elif float(elapsed) > 60:
        print "  |->"+BOLDME+"Aligned %d genomes in %.2f minutes"%(totseqs,float(elapsed)/float(60.0))+ENDC
    else: 
        print "  |->"+BOLDME+"Aligned %d genomes in %.2f seconds"%(totseqs,float(elapsed))+ENDC
    #cleanup
    rmfiles = glob.glob(outputDir+os.sep+"*.aln")
    rmfiles2 = glob.glob(outputDir+os.sep+"*.mummerout.out")
    rmfiles3 = glob.glob(outputDir+os.sep+"blocks/b*/*")
    rmfiles4 = glob.glob(outputDir+os.sep+"blocks/b*")
    for file in rmfiles:
        os.system("rm %s"%(file))
    for file in rmfiles2:
        os.system("rm %s"%(file))
    for file in rmfiles3:
        os.system("rm %s"%(file))
    for file in rmfiles4:
        os.system("rmdir %s"%(file))



    filepres = 0
    print BOLDME+"\n<<Parsnp finished! All output available in %s>>"%(outputDir)+ENDC
    print
    print BOLDME+"Validating output directory contents..."+ENDC
    print BOLDME+"\t1)parsnp.tree:\t\tnewick format tree"+ENDC,
    if os.path.exists("%sparsnp.tree"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.tree"%(outputDir+os.sep)) > 0:
        print "\t\t\t["+OK_GREEN+"OK"+ENDC+"]"
        filepres+=1
    else:
        print "\t|->"+ERROR_RED+"MISSING"+ENDC
    print BOLDME+"\t2)parsnp.vcf:\t\tVCF format SNP calls"+ENDC,
    if os.path.exists("%sparsnp.vcf"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.vcf"%(outputDir+os.sep)) > 0:
        print "\t\t\t["+OK_GREEN+"OK"+ENDC+"]"
        filepres+=1
    else:
        print "\t|->"+ERROR_RED+"MISSING"+ENDC
    print BOLDME+"\t3)parsnp.ggr:\t\tharvest input file for gingr (GUI)"+ENDC,
    if os.path.exists("%sparsnp.ggr"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.ggr"%(outputDir+os.sep)) > 0:
        print "\t["+OK_GREEN+"OK"+ENDC+"]"
        filepres+=1
    else:
        print "\t|->"+ERROR_RED+"MISSING"+ENDC
    print BOLDME+"\t4)parsnp.xmfa:\t\tXMFA formatted multi-alignment"+ENDC,
    if os.path.exists("%sparsnp.xmfa"%(outputDir+os.sep)) and os.path.getsize("%sparsnp.xmfa"%(outputDir+os.sep)) > 0:
        print "\t\t["+OK_GREEN+"OK"+ENDC+"]"
        filepres+=1
    else:
        print "\t|->"+ERROR_RED+"MISSING"+ENDC
    if filepres == 4:
        pass

    else:
        print "\t\t["+ERROR_RED+"Output files missing, something went wrong. Check logs and relaunch or contact developers for assistance"+ENDC+"]"
    print
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


    files = glob.glob("%s/*.srt"%(outputDir))
    for file in files:
        if not VERBOSE and os.path.exists(file):
            os.remove(file)

    if not VERBOSE and os.path.exists("%s/psnn.ini"%(outputDir)):
        os.remove("%s/psnn.ini"%(outputDir))

    if not VERBOSE and os.path.exists("%s/all_mumi.ini"%(outputDir)):
        os.remove("%s/all_mumi.ini"%(outputDir))

    if not VERBOSE and os.path.exists("%s/parsnpAligner.ini"%(outputDir)):
        os.remove("%s/parsnpAligner.ini"%(outputDir))

    if not VERBOSE and os.path.exists("%s/parsnp.unalign"%(outputDir)):
        os.remove("%s/parsnp.unalign"%(outputDir))
    if os.path.exists("%s/parsnp.snps.mblocks"%(outputDir)):
        os.remove("%s/parsnp.snps.mblocks"%(outputDir))

    if not VERBOSE and os.path.exists("%s/all.mumi"%(outputDir)):
        os.remove("%s/all.mumi"%(outputDir))

    if os.path.exists(use_gingr):
        #check if available first
        rc = 0
        if 0 and binary_type == "linux":
            from subprocess import Popen,PIPE
            p = Popen(['xset','-q'], stdout=PIPE,stderr=PIPE)
            p.communicate()
            rc = p.returncode
        if binary_type == "osx":
            print ">>Launching gingr.."
            os.system("open -n %s --args %s/parsnp.ggr"%(use_gingr,outputDir))

