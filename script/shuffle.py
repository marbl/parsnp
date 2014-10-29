import string
import sys
import random
import os

#RC
def invertSeq(sseq):
    #print sseq
    rseq = ""
    count = 0
    tseq = sseq
    tseq = tseq.replace('\n','')
    tseq = tseq[::-1]
    for base in tseq:
        
            
         base = string.upper(base)
         if base == 'A':
             base = 'T'
             rseq+=base
         elif base == 'T':
             base = 'A'
             rseq+=base
         elif base == 'C':
             base = 'G'
             rseq+=base
         elif base == 'G':
             base = 'C'
             rseq+=base
         else:
             
             #continue
             rseq+=base
         count+=1
         if count == 60:
             rseq+="\n"
             count = 0
         

    return rseq
     
if __name__ == "__main__":


     
     random.seed()
     
     infileName = ""
     outfileName = ""
     rsize = 0
     rnumber = 0
     shuffleperseq = 1
     if len(sys.argv) < 6:
         print "\nUsage: shuffleGenome <FastA input file> <output file> <Rearrangement size> <Number of rearrangements> <Shuffled sequences per seq>"
         sys.exit(1)
     else:
         infileName = sys.argv[1]
         outfileName = sys.argv[2]
         rsize = int(sys.argv[3])
         rnumber = int(sys.argv[4])
         shuffleperseq = int(sys.argv[5])
         

     
     filelist = []
     if infileName == ".":
         tempfiles = os.listdir(".")
         for filen in tempfiles:
             if filen.find("fna") != -1:
                 filelist.append(filen)

     else:
         filelist.append(infileName)

     for seqfile in filelist:
         infile = open(seqfile,'r')
         header = infile.readline()
         infiled = infile.read()
         filelen = len(infiled)
         print "\nInput sequence:%s"%seqfile
         print "Input sequence length: %d\n"%filelen

         #split genome into filelen/rsize parts
         partlist = []
         pos = 0
         if rsize > filelen:
             rsize = filelen
         for part in xrange(0,filelen/rsize):
             partlist.append(infiled[pos:(part+1)*rsize])
             pos +=rsize

         partlist.append(infiled[pos:])


         for shuffleit in xrange(0,shuffleperseq): 
             count = 0
             seq = ""
             parttemp = ""
             while count < rnumber:
                 #choose the 1st part
                 part1 = random.randint(0,len(partlist)-1)
                 #choose the 2nd
                 part2 = part1
                 while part2 == part1:
                     part2 = random.randint(0,len(partlist)-1)

                 #now choose:
                 #1)transposition,
                 #2)inversion,

                 operation = random.randint(1,2)
                 if operation == 1:
                     #transposition, translocation
                     parttemp = partlist[part1]
                     partlist[part1] = partlist[part2]
                     partlist[part2] = parttemp
                     print "Transposition"
                     print "  Positions %d and %d swapped"%(part1,part2)
                 elif operation == 2:
                     #inversion
                     partlist[part1] = invertSeq(partlist[part1])
                     print "Inversion"
                     print " Position %d"%part1


                     
                 count +=1

             #create output
             output = header
             output = ">Shuffled-%d-"%(shuffleit+1) + header[1:]
             
             partstring = string.upper(string.join(partlist))
             partstring = partstring.replace(' ', '')
             output += string.strip(partstring)
             
             #write output to file
             fname = seqfile[:-4]
             fname+= "_Shuffled_%d.fna"%(shuffleit+1)
             print "output: %s"%(fname)
             
             fout = open(fname,'w')
             fout.write(output)
             fout.close()
     
