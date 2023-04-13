'''
This script randomly reorganizes regions of the chromosome drastically.

Notes:

will need to tier the breakpoint selection instead. random.choices: Done
After chromosome selection
function reorientation will need weighting added instead of randint: Done

microhomology mediated assembly
for breakends, during reassembly, randomize in range from negative to positive
    Use this to determine amount removed from end or added on.

create an input parameter file: Version 1 made

Issues:fragile sites and cytoband records do not seem to align properly
with coordinates given.

finish write genome: Done
finish write changes: Done
finish reassemble: Done

p arms 13,14,15,21,22 no longer selectable due to entire sequence being N

To do:
Read in bed file: Fragile sites
Read in bed file: Gaps
Or Create 1 bedfile with weights of all regions, gaps weighted at 0.
--chr coord1 coord2 score
-either way, restructure of select_break required, no more chromosome arms.
Determine best way to use BEDPE for variant recording. One record per breakpoint
Unsure how to ensure that proper length of fragment after breakpoint is seen
'''

import numpy
from Bio.Seq import Seq
from Bio import SeqIO
import random
from pybedtools import BedTool
import os
import argparse

qs_and_ps=[]#a list of chromosomes that both the q and p arm have been selected
deleted={}#dictionary of deleted regions

'''
Input List:
Large Breakpoint Count
Breakpoint Gap
input genome
output tag
parameter file # This will have the weighting values
'''
#Parameter file format
'''
copy numbersie(0,1,2,3,4,5)
copy weightsie(1,2,13,4,5)
fragment inversion weight(Inv,Nor)
internal weights(Inv,Del,Dup)
Internal inversion weights 'tandup n','invdup r','invdup l','dup l','dup r','n'
Internal duplication weights 'tan','inv l','inv r','nor l','nor r'
'''
#argparse input gathering
parser = argparse.ArgumentParser(description="Simulates a complex rearrangement")
parser.add_argument('--breaks', '-b', type=int,
                    default=100,
                    help='count of breakpoints used to generate initial fragments')#Initial breakpoint count
parser.add_argument('--gap','-g', type=int,
                    default=50,
                    help='the minimum distance between breakpoints generated')#Minimum gap between breakpoints
parser.add_argument('--ref','-r', type=str,
                    required=True,
                    help='path to the reference genome')#filepath to ref genome
parser.add_argument('--out','-o', type=str,
                    required=True,
                    help='path to output location and tag that will be the first part of each output file')#output location and tag
parser.add_argument('--chrome','-c', type=str,
                    required=True,
                    help='path to the bed file designating chromosomes and lengths')
parser.add_argument('--gaps','-l', type=str,
                    required=True,
                    help='path to bed file with gaps. Requires a centromere for each chromosome')
parser.add_argument('--frag','-f', type=str,
                    required=True,
                    help='path to bed file of fragile sites')
parser.add_argument('--param','-p',type=str,
                    required=True,
                    help='path to parameter file')#location of parameter file
parser.add_argument('--maxarms','-a', type=int,
                    default=4,
                    help='maximum number of chromosome arms possibly effected')
args=parser.parse_args()


#Read in parameters from parameter file
par_fil=open(args.param,'r')
pars=par_fil.readlines()
par_fil.close()
parsdict={}
for i in range(len(pars)):
    pars[i]=pars[i].strip()
    if pars[i][0]=="#":
        continue
    else:
        pars[i]=pars[i].split(':')
        pars[i][1]=pars[i][1].strip().split(',')
        for x in range(len(pars[i][1])):
            pars[i][1][x]=int(pars[i][1][x])
        parsdict[pars[i][0]]=pars[i][1]
cn=parsdict['Copy Number']#copy numbers
cw=parsdict['Copy Weights']#weights for above copy numbers
inw=parsdict['Inversion Weights']
ib=parsdict['Internal Breaks']
iw=parsdict['Internal Rearrangement Weights']
iiw=parsdict['Internal Inversion Alt Weights']
idw=parsdict['Internal Duplication Alt Weights']
#Designate Gap and Initial Breakpoint
breaktotal=args.breaks
gap=args.gap
inputfile=args.ref
outputtag=args.out


def loadref(inputfile):
    '''
    Reads in a fasta file specified into a list of seqrecord objects
    Each chromosome should be its own seqrecord object
    '''
    seq={}
    for seq_record in SeqIO.parse(open(inputfile,mode='r'), "fasta"):
        seq[seq_record.id]=seq_record
    return seq

def loadweights(bedtool, weights, centromeres):
    '''
    Loads bedfile into a dictionary of lists of lists.
    Key: Chromosome identifier
    List 1: [[start, end]]
    List 2: [weight]
'''
    for line in bedtool:
        tag=''
        if centromeres[line.chrom][0] > int(line.stop):
            tag=line.chrom + ' p'
        else:
            tag=line.chrom + ' q'
        if tag in weights:
            weights[tag][0].append([int(line.start),int(line.stop)])
            weights[tag][1].append(int(line.score))
        else:
            weights[tag]=[[[int(line.start),int(line.stop)]],
                                 [int(line.score)]]

def mergebed():#remove gaps with subtract
    '''
    1. remove gaps from chromes and frags by subtracting gaps
    2. remove frags from chromes by subtracting frags
    3. combine frags and chromes into one BedTool
    4. write Bed file
    5. read in Bed file for weights
    return said weights list

Extract centromeres into a dictionary. Use the centromeres to compare
fragments being loaded to divide into p and q arms.
Determine way to remove arms that are fully gaps
    
'''
    chromes=BedTool(args.chrome)
    frag=BedTool(args.frag)
    gaps=BedTool(args.gaps)
    centromeres={}
    for entry in gaps:
        if entry.name == 'centromere':
            centromeres[entry.chrom]=[int(entry.start),int(entry.stop)]
    chromes=chromes.subtract(gaps)
    frag=frag.subtract(gaps)
    frag=frag.intersect(chromes)
    chromes=chromes.subtract(frag)
    weights={}
    loadweights(chromes,weights, centromeres)
    loadweights(frag,weights, centromeres)
    return weights, BedTool(args.chrome)

def makevcf(bednames):#Maybe one day. Might go BedPE, likely will do this
    '''
    Take in a list of bedfile names. Generates a vcf from the information.
    VCF Columns
    Chrom, position, id(bnd), ref(leave blank), alt(leave blank), qual,
    filter(PASS),info
    bed format:
    Chrom, start, end, strand
    '''
    header='''
##fileformat=VCFv4.1
##source=muddler
##fileDate=20220627
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##ALT=<ID=TRA,Description="Translocation">
##FILTER=<ID=UNRESOLVED,Description="An insertion that is longer than the read and thus we cannot predict the full size.">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=MATEID,Number=.,Type=String,Description="IDs of mate breakends">
##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
'''
    out=open(outputtag+'.vcf','w')
    out.write(header)
    inn=loadref(inputfile)
    counter=1
    for bedname in bednames:
        bed=open(bedname,'r')
        prior=bed.readline().strip().split()
        bedlines=bed.readlines()
        bed.close()
        for lin in bedlines:
            line=lin.strip().split()
            writ1='\n'+prior[0]+'\t'
            writ2='\n'+line[0]+'\t'
            pos1=0
            pos2=0
            if prior[-1]=='+':#poswrit1
                pos1=prior[2]
            else:
                pos1=prior[1]
            if line[-1]=='+':#poswrit2
                pos2=line[1]
            else:
                pos2=line[2]
            ref1=str(inn[prior[0]].seq[int(pos1)])
            ref2=str(inn[line[0]].seq[int(pos2)])
            bnd1='BND_'+str(counter)+'_1'
            bnd2='BND_'+str(counter)+'_2'
            alt1=''
            alt2=''
            if prior[-1]=='-':
                alt2=']'+prior[0]+':'+pos1+']'
            else:
                alt2='['+prior[0]+':'+pos1+'['
            if line[-1]=='-':
                alt1=']'+line[0]+':'+pos2+']'
                alt2=ref2+alt2
            else:
                alt1='['+line[0]+':'+pos2+'['
                alt2=alt2+ref2
            if prior[-1]=='-':
                alt1=alt1+ref1
            else:
                alt1=ref1+alt1
            info='SVTYPE=TRA'
            ##Ideally, add in MATEIDs into Info category. bnd_1_2 could also
            ##connect with bnd_11_1 but how to track that information?
            writ1+=str(pos1)+'\t'+bnd1+'\t'+ref1+'\t'+alt1+'\t\tPASS\t'+info
            writ2+=str(pos2)+'\t'+bnd2+'\t'+ref2+'\t'+alt2+'\t\tPASS\t'+info
            
            out.write(writ1)
            out.write(writ2)
            prior=line
            counter+=1
    out.close()
    return

def breakpoint(chrome, region, breaks):
    '''
    Decides where a breakpoint will be in a given region
    chrome: the chromosome being broken up
    region: the region on the chromosome selected via weighted random
    breaks: a dictionary containing the breakpoints selected so far
    returns: 0 or 1. 0 if unable to find a breakpoint not nearby other breaks
    already. 1 if a new breakpoint was added to the breaks dictionary.
    '''
    closeby=True
    point=0
    attempts=0
    print('Selected Site')
    print(region)
    while closeby:
        if chrome not in breaks.keys():
            breaks[chrome]=[]
            point=random.randrange(region[0],region[1])
            breaks[chrome].append(point)
            closeby=False
        else:
            point=random.randrange(region[0],region[1])
            closeby=False
            for spot in breaks[chrome]:
                if spot <= point+gap and spot >= point-gap:
                    closeby=True
                    break
            if not closeby:
                breaks[chrome].append(point)
        attempts+=1
        if attempts >= 10:
            return 0
    return 1

def select_break_2():
    '''
    Master function for selecting the breaks in the chromosomes.
    returns a dictionary of chromosomes(keys) and the breakpoints.
    '''
    chrcount=random.randrange(1,args.maxarms)#number of effected chromosomes
    print('Number of effected chromosome arms')
    print(chrcount)
    weights, chromes=mergebed()
    chrs=[]
    allchr=list(weights.keys())
    for i in range(chrcount+1):
        chrs.append(allchr.pop(random.randrange(0,len(allchr))))
    print('Effected chromosomes')
    print(chrs)
    for a in chrs:#any p and q chromosome regions on the same chromosome?
        for b in chrs:
            if a.split()[0]==b.split()[0]and a.split()[1]!=b.split()[1]:
                qs_and_ps.append(a.split()[0])
    fsites={}
    for chrome in chrs:
        fsites[chrome]=weights[chrome]
    #chromosomes selected. Now to create the breakpoints.
    breaks={}
    breaksmade=0
    while breaksmade < breaktotal:
        chrome=random.choice(list(fsites.keys()))
        region=random.choices(fsites[chrome][0],fsites[chrome][1])[0]
        breaksmade+=breakpoint(chrome, region, breaks)
    return breaks, chromes

def select_break():
    '''
    Master function for selecting the breaks in the chromosomes.
    returns a dictionary of chromosomes(keys) and the breakpoints.
    '''
    chrcount=random.randrange(1,args.a)#number of effected chromosomes
    print('Number of effected chromosome arms')
    print(chrcount)
    chrs=[]
    allchr=["chr1 q","chr1 p","chr2 q","chr2 p","chr3 q","chr3 p",
            "chr4 q","chr5 q","chr6 q","chr4 p","chr5 p","chr6 p",
            "chr7 q","chr8 q","chr9 q","chr7 p","chr8 p","chr9 p",
            "chr10 q","chr11 q","chr12 q","chr10 p","chr11 p","chr12 p",
            "chr13 q","chr14 q","chr15 q",
            "chr16 q","chr17 q","chr18 q","chr16 p","chr17 p","chr18 p",
            "chr19 q","chr20 q","chr21 q","chr19 p","chr20 p",
            "chr22 q","chr23 q","chrX q","chr23 p","chrX p"]
    for i in range(chrcount+1):
        chrs.append(allchr.pop(random.randrange(0,len(allchr))))
    print('Effected chromosome arms')
    print(chrs)
    for a in chrs:#any p and q chromosome regions on the same chromosome?
        for b in chrs:
            if a.split()[0]==b.split()[0]and a.split()[1]!=b.split()[1]:
                qs_and_ps.append(a.split()[0])
    fsites=load_fsites(chrs)
    #chromosomes selected. Now to create the breakpoints.
    breaks={}
    breaksmade=0
    while breaksmade < breaktotal:
        chrome=random.choice(list(fsites.keys()))
        region=random.choices(fsites[chrome][0],fsites[chrome][1])[0]
        breaksmade+=breakpoint(chrome, region, breaks)
    return breaks

def select_break_test():
    '''
    Master function for selecting the breaks in the chromosomes.
    returns a dictionary of chromosomes(keys) and the breakpoints.
    Simply a test function for making breaks on chr 1 q and p
    '''
    
    chrs=["chr1 q","chr1 p"]
    for a in chrs:#any p and q chromosome regions on the same chromosome?
        for b in chrs:
            if a.split()[0]==b.split()[0]and a.split()[1]!=b.split()[1]:
                qs_and_ps.append(a.split()[0])
    fsites=load_fsites(chrs)
    #chromosomes selected. Now to create the breakpoints.
    breaks={}
    breaksmade=0
    while breaksmade < breaktotal:
        chrome=random.choice(list(fsites.keys()))
        region=random.choices(fsites[chrome][0],fsites[chrome][1])[0]
        breaksmade+=breakpoint(chrome, region, breaks)
    return breaks

def define_regions(breaks):
    '''
    Defines the regions that will be rearranged based on a list of breakpoints
    breaks: a dictionary of lists. Chromosome number key with list of integers
    returns: regions. a dictionary of list of lists, like breaks but with lists
    of regions instead of just integers
    '''
    regions={}
    for chrome in breaks:
        breaks[chrome].sort()
        regions[chrome]=[]
        for i in range(len(breaks[chrome])-1):
            regions[chrome].append([breaks[chrome][i],breaks[chrome][i+1]])
    return regions
    
def reorientation(regions):
    '''
    reorients the various regions of the chromosomes.
    First chooses copy number(duplications and deletions)
    Chooses from inverted or normal orientation for the retained segments
    
    '''
    oriented=[]
    
    for chrome in regions:
        deleted[chrome]=[]
        for region in regions[chrome]:
            copies=random.choices(cn,cw)[0]
            if copies==0:
                deleted[chrome].append(region)
            else:
                for copy in range(copies):
                    #oriented takes in region endpoints, then indicates
                    #inversion or normal orientation. List at end is for
                    #internal rearrangements in a later step
                    oriented.append(
                        [region[0],region[1],
                        random.choices([True,False],inw)[0],chrome.split()[0],[]]
                        )
                        
    return oriented

def inside_segments(oriented):
    '''
    Does internal rearrangements, if any, for the various retained segments
    0, 1, or 2 internal breakpoints.
    Still maintains gap between breakpoints.
    *
    '''
    for fragment in oriented:
        if fragment[1]-fragment[0]> gap*3:
            bps=random.choices([0,1,2],ib)[0]
            if bps==2:
                bp=random.randint(0,1)
                choices=[
                    random.randint(fragment[0]+gap,fragment[1]-(2*gap)),
                    random.randint(fragment[0]+(gap*2),fragment[1]-gap)]
                bp1=choices[bp]
                if bp:
                    bp2=random.randint(fragment[0]+gap,bp1-gap)
                else:
                    bp2=random.randint(bp1+gap,fragment[1]-gap)
                fragment[-1].append(bp1)
                fragment[-1].append(bp2)
                fragment[-1].sort()
                orient=random.choices(['inv','del','dup'],iw)[0]
                fragment[-1].append(orient)
                if orient=='inv':
                    fragment[-1].append(random.choices(
                        ['tandup n','invdup r','invdup l','dup l','dup r','n'],
                        iiw)[0])
                elif orient=='dup':
                    fragment[-1].append(random.choices(['tan','inv l','inv r',
                                                        'nor l','nor r'],
                                                       idw)[0])
            elif bps==1:
                bp1=random.randrange(fragment[0]+gap,fragment[1]-gap)
                bp2=fragment[random.randint(0,1)]
                fragment[-1].append(bp1)
                fragment[-1].append(bp2)
                fragment[-1].sort()
                orient=random.choices(['inv','del','dup'],iw)[0]
                fragment[-1].append(orient)
                if orient=='inv':
                    fragment[-1].append(random.choices(
                        ['tandup n','invdup r','invdup l','dup l','dup r','n'],
                        iiw)[0])
                elif orient=='dup':
                    fragment[-1].append(random.choices(['tan','inv l','inv r',
                                                        'nor l','nor r'],
                                                       idw)[0])
        elif fragment[1]-fragment[0]>gap*2:
            bps=random.choices([0,1],[ib[0],ib[1]])
            if bps:
                bp1=random.randrange(fragment[0]+gap,fragment[1]-gap)
                bp2=fragment[random.randint(0,1)]
                fragment[-1].append(bp1)
                fragment[-1].append(bp2)
                fragment[-1].sort()
                orient=random.choice(['inv','del','dup'])
                fragment[-1].append(orient)
                if orient=='inv':
                    fragment[-1].append(random.choices(
                        ['tandup n','invdup r','invdup l','dup l','dup r','n'],
                        iiw)[0])
                elif orient=='dup':
                    fragment[-1].append(random.choices(['tan','inv l','inv r',
                                                        'nor l','nor r'],
                                                       idw)[0])
                    

    return

def reassemble(oriented,breaks,chromes):
    '''
    devnotes
    Make assembled an intermediate datastructure
    Do a pass and create a new datastructure to handle centromere and bring
    all chromosomes to be in chr1, chr2, ... format, rather than have the
    ps and qs floating about still
    
'''
    assembled={}
    chrms=list(breaks.keys())
    for chrome in chrms:
        assembled[chrome]=[]
        if chrome.split()[1]=='p':
            assembled[chrome].append([0,breaks[chrome][0],
                                      False,chrome.split()[0],[]])
    while oriented:
        i=random.randrange(len(oriented))
        oriented[i],oriented[-1]=oriented[-1],oriented[i]
        piece=oriented.pop()
        chrome=random.choice(chrms)
        assembled[chrome].append(piece)
    for chrome in chrms:#add q end and telomere to q arms
        if chrome.split()[1]=='q':
            chrmend=None
            for entry in chromes:
                if entry.chrom==chrome.split()[0]:
                    chrmend=entry.stop
            assembled[chrome].append([breaks[chrome][-1],chrmend,
                                     False,chrome.split()[0],[]])
    completed={}
    telomeres={}
    coords=open(args.chrome,'r')
    for line in coords:
        lin=line.split()
        coords[lin[0]]=int(lin[2])
    coords.close()
    for chrome in chrms:
        chrm=chrome.split()[0]
        if chrm in qs_and_ps:
            centromere=[[breaks[chrm+' p'][-1],breaks[chrm+' q'][0],
                         random.choice([True,False]),chrm,[]]]
            completed[chrm]=assembled[chrm+' p']+centromere+assembled[chrm+' q']
        elif chrome.split()[1]=='p':
            centro_to_q=[[breaks[chrome][-1],telomeres[chrm],#(once more testing only)
                          False,chrm,[]]]
            completed[chrm]=assembled[chrome]+centro_to_q
        elif chrome.split()[1]=='q':
            p_to_centro=[[0,breaks[chrome][0],#(once more testing only)
                          False,chrm,[]]]
            completed[chrm]=p_to_centro+assembled[chrome]
    return completed

def write_changes(completed):
    '''
    writes out a record of the changes made to the genome.
    '''
    ofil=open(outputtag+'.csv','w')
    ofil.write('start,end,inv,chrm,instart,instop,inchange,inchange2\n')
    for chrome in completed:
        ofil.write(str(chrome)+'\n')
        for change in completed[chrome]:
            print('newostr')
            #start,end,invert,chrm,internal
            #internal: start,end,del/dup/inv
                #for del, just 3 long
                #for inv, 4th is
                #'tandup n','invdup r','invdup l','dup l','dup r','n'
                #for dup, 4th is 'tan','inv l','inv r','nor l','nor r'
            ostr=str(change[0])+','+str(change[1])+','+str(change[2])+','+change[3]
            print(ostr)
            print(change)
            if change[4]:
                ostr+=','+str(change[4][0])+','+str(change[4][1])+','+change[4][2]
                print(ostr)
                if change[4][-1]!=change[4][2]:
                    ostr+=','+str(change[4][3])+'\n'
                    print(ostr)
                else:
                    ostr+=',NA\n'
                    print(ostr)
            else:
                ostr+=',NA,NA,NA,NA\n'
            ofil.write(ostr)
    ofil.close()
    return



def write_to_bed(completed):
    '''
    Writes out the changes to multiple bed files. One per rearranged
    chromosome.
    tag+'_'+chromosome name.bed
    so like messedup1_chrm4.bed
    Also now writes out a bedfile of the deleted segments. Internally deleted
    sections are not recorded in this deletion file
    
    '''
    bedfiles=[]
    for chrome in completed:
        ofil=open(outputtag+'_'+chrome+'.bed','w')
        bedfiles.append(outputtag+'_'+chrome+'.bed')
        for seg in completed[chrome]:
            #start,end,invert,chrm,internal
            #internal: start,end,del/dup/inv
                #for del, just 3 long
                #for inv, 4th is
                #'tandup n','invdup r','invdup l','dup l','dup r','n'
                #for dup, 4th is 'tan','inv l','inv r','nor l','nor r'
                #
            coords=[]
            hold=[]
            skip=False
            if seg[-1]:
                if seg[0] not in seg[-1]:
                    #first to first segment
                    coords.append([seg[0],seg[-1][0],'+'])
                #Internal handling. del is ignored as it doesn't show up
                    #in the output
                if seg[-1][2]=='inv':
                    more=seg[-1][3].split()[0]
                    coords.append([seg[-1][0],seg[-1][1],'-'])
                    if more == 'tandup':
                        coords.append([seg[-1][0],seg[-1][1],'-'])
                    elif more == 'invdup':
                        hold = [seg[-1][0],seg[-1][1],'+']
                        if seg[-1][3].split()[1]=='l':
                            coords=[hold]+coords
                            hold=[]
                    elif more == 'dup':
                        hold = [seg[-1][0],seg[-1][1],'-']
                        if seg[-1][3].split()[1]=='l':
                            coords=[hold]+coords
                            hold=[]
                if seg[-1][2]=='dup':
                    more=seg[-1][3].split()[0]
                    coords.append([seg[-1][0],seg[-1][1],'+'])
                    if more == 'tandup':
                        coords=[[seg[0],seg[-1][1],'+']]
                        coords.append([seg[-1][0],seg[-1][1],'+'])
                    elif more == 'inv':
                        skip=True
                        coords=[[seg[0],seg[1],'+']]
                        hold = [seg[-1][0],seg[-1][1],'-']
                        if seg[-1][3].split()[1]=='l':
                            coords=[hold]+coords
                            hold=[]
                    elif more == 'nor':
                        skip=True
                        coords=[[seg[0],seg[1],'+']]
                        hold = [seg[-1][0],seg[-1][1],'+']
                        if seg[-1][3].split()[1]=='l':
                            coords=[hold]+coords
                            hold=[]
                if seg[1] not in seg[-1] and not skip:
                    #second to second segment
                    coords.append([seg[-1][1],seg[1],'+'])
                if hold:
                    coords.append(hold)
                
            else:
                coords=[[seg[0],seg[1],'+']]
            if seg[2]:
                for i in range(len(coords)):
                    odex=-(i+1)
                    olin=seg[3]+'\t'
                    olin+=str(coords[odex][0])+'\t'+str(coords[odex][1]-1)+'\t\t\t'
                    if coords[odex][2]=='+':
                        olin+='-\n'
                    else:
                        olin+='+\n'
                    ofil.write(olin)
            else:
                for i in range(len(coords)):
                    olin=seg[3]+'\t'
                    if seg == completed[chrome][-1]:
                        olin+=str(coords[i][0])+'\t'+str(coords[i][1])+'\t\t\t'
                    else:
                        olin+=str(coords[i][0])+'\t'+str(coords[i][1]-1)+'\t\t\t'
                    olin+=coords[i][2]+'\n'
                    ofil.write(olin)
                    
        ofil.close()
    ofil=open(outputtag+'_deleted.bed','w')
    for chrome in deleted:
        for seg in deleted[chrome]:
            olin=chrome+'\t'+str(seg[0])+'\t'+str(seg[1])+'\n'
            ofil.write(olin)
    ofil.close()
    return bedfiles

def bed_get_genome(bedfiles):
    '''
    Takes the bed files to create a full genome haplotype
    writes the unmodified chromosomes and modified chromosomes to a fasta
    '''
    endseq=[]
    chrms=[]
    for name in bedfiles:
        temp=outputtag+'_temp.fa'
        bedt=BedTool(name)
        bedt=bedt.sequence(fi=inputfile,fo=temp,s=True)
        seq=''
        for seg in SeqIO.parse(temp,"fasta"):
            print(seg.id)
            if isinstance(seq,str):
                seq=seg
            else:
                seq=seq+seg
        tag=name.strip('.bed').split('_')[-1]
        seq.name=''
        seq.id=tag+'_modified'
        seq.description=''
        os.remove(temp)
        endseq.append(seq)
        chrms.append(tag)
    refseq=loadref(inputfile)
    for seq in refseq:
        if seq not in chrms:
            endseq.append(refseq[seq])
    SeqIO.write(endseq, outputtag+'.fa','fasta')
    return

def main():
    '''
    asaaa
    '''
    breaks, chromes=select_break_2()

    #print(qs_and_ps)
    print(breaks)

    regions=define_regions(breaks)

    print('\nRegions\n')
    print(regions)

    oriented=reorientation(regions)

    #print('\nReorientOutside\n')
    #print(oriented)

    inside_segments(oriented)

    #print('\nInside\n')
    #print(oriented)
    #print('\nReassembly and Deleted\n')

    assembled=reassemble(oriented,breaks,chromes)

    for section in assembled:
        print(section)
        for frag in assembled[section]:
            print(frag)
    #print()
    for section in deleted:
        print(section)
        for frag in deleted[section]:
            print(frag)
    #print('\nCentromere\n')
    #print(str(breaks['chr1 p'][-1])+'\n'+str(breaks['chr1 q'][0]))
    write_changes(assembled)
    bedfiles=write_to_bed(assembled)
    bed_get_genome(bedfiles)
    makevcf(bedfiles)
    return



main()

