'''

This module includes functions composing a motif merging program by determining the offset of two motifs, cutting off excess spaces at the end of the motif, and combining instances. 
Uses functions from author ichaudr


'''

from Bio import motifs, SeqIO
from Bio.motifs import Instances, Motif
from pathlib import Path
import sys
import math
motifone_aligned = []
motiftwo_aligned = []

#### Functions used to create input motif from jaspar/fasta/file ####             
def motif_create(file_name,i):
    '''
    Creates or reads a motif from a fasta, jaspar, or xml file.

    Parameters
    ----------
    file_name : String
        file name with motif or instances (Assumes the file directory is in the motif folder within working directory)

    i : integer
        If reading a jaspar or xml file, specifies the wanted motif. 

    Returns
    -------
    motif : Motif Object
        motif created or read from the file

    '''
    entries = Path('motif/')
    for entry in entries.iterdir():
        if file_name in str(entry):
            if file_name != None and '.fas' in file_name:        
                seqlist = []
                seqs = SeqIO.parse(str(entry),'fasta')
                for seq in seqs:
                    seqlist.append(seq.seq)
                motif = motifs.create(seqlist,'ACGT')
            if file_name != None and '.jaspar' in file_name:
                with open(str(entry)) as q:
                    motif = motifs.parse(q,'jaspar')[i]
            if file_name != None and '.xml' in file_name:
                with open(str(entry)) as w:
                    motif = motifs.parse(w,'meme')[i]
            else:
                print('Unrecognized file format')
    return(motif)


### Functions below used to determine the alignment of two motifs that maximizes information content of said alignments ###
def get_alignment_offset(motif, other):
    '''
    @author ichaudr
    Determines the optimal alignment of two motifs by maximizing the information content (ic) in the aligned regions.
    
    Parameters
    ----------
    motif, other: Motif objects
        The two motifs of interest.
    
    Returns
    -------
    offsets: int
        The offset that results in the maxium ic in the alignment. 
    '''
    max_ic = float('-inf')
    for offset in range(-len(motif) + 1, len(other)):
        if offset < 0:
            ic = ic_at(motif, other, -offset)
        else:
            ic = ic_at(other, motif, offset)

        
        if ic > max_ic:
            max_ic = ic
            max_offset = offset
            
    return max_offset

def ic_at(motif, other, offset):
    '''
    @author ichadur
    Caculates the information content, ic, for a specific alignment. The approach makes a temporary motif object containing the overlapping sequences in the alignemnt and taking the average of the pssm.
    
    Parameters
    ----------
    motif, other: Motif objects
        The motifs of interest
    offset: int
        The offset value that results in the alignment of interest. 
    '''

    #Pull the sequences containined in the aligned region of the motifs from each of the motif instances. 
    alignment_len = min(len(motif)-offset, len(other))
    motif_seqs = [site[offset:alignment_len+offset] for site in motif.instances]
    other_seqs = [site[:alignment_len] for site in other.instances]

    # Create the motif and compute the IC
    amotif = Motif(instances=Instances(motif_seqs+other_seqs))
    amotif.pseudocounts = dict(A=0.25, C=0.25, G=0.25, T=0.25)

    #print('Motif Seqs: ' , motif_seqs)
    #print('Other Seqs: ' , other_seqs)
    #print('Offset ', offset)
    #print('IC: ' , amotif.pssm.mean(), '\n\n')

    return amotif.pssm.mean()

### Functions below prepare the two inputted motifs for merging ###
# Approach:
# 1. In the case of a motif without instances, fake instances are created based on the counts to allow for offset value to be determined
# 2. Use the outputted offset value for maximized information content to align the two motifs
def fake_instances(motif):
    ''' 
    Creates fake instances for motifs without reported instances based on the count values
    
    Parameters 
    ----------
    motif: Motif Object
        Motif object that does not possess instances
    
    Returns 
    ----------
    instanlist: list
        List of generated instances that would produce the same matrix as inputted
    '''
    instanlist = []
    counter = 0
    numinstances = motif.counts[0][0] + motif.counts[1][0] + motif.counts[2][0] + motif.counts[3][0]
# Iterates through matrix length and base wise (number of instances = summation of all base counts at a particular position)
    for i in range(0, motif.length):
        for l in range(0, int(numinstances)):
            if i == 0:
                #Intializes the list by adding the reported base based on the count in the matrix. List length = instance count
                if motif.counts[0][0] > 0:
                    instanlist.append('A')
                    motif.counts[0][0] -= 1
                if motif.counts[1][0] > 0:
                    instanlist.append('C')
                    motif.counts[1][0] -= 1
                if motif.counts[2][0] > 0:
                    instanlist.append('G')
                    motif.counts[2][0] -= 1 
                if motif.counts[3][0] > 0:
                    instanlist.append('T')
                    motif.counts[3][0] -= 1
            if i > 0:
                #At each position, adds to existing list (instanlist) elements based on the the property .counts[base][i]
                #Uses counter variable to keep track of which position in list to add base
                if motif.counts[0][i] > 0:
                    instanlist[counter] = instanlist[counter] + 'A'
                    motif.counts[0][i] -= 1
                    counter += 1
                elif motif.counts[1][i] > 0:
                    instanlist[counter] = instanlist[counter] + 'C'
                    counter += 1
                    motif.counts[1][i] -= 1
                elif motif.counts[2][i] > 0:
                    instanlist[counter] = instanlist[counter] + 'G'
                    counter += 1                
                    motif.counts[2][i] -= 1
                elif motif.counts[3][i] > 0:
                    instanlist[counter] = instanlist[counter] + 'T'
                    counter += 1
                    motif.counts[3][i] -= 1
        counter = 0
    #Returns instanlist after iteration through bases and lengths
    return instanlist
def offset_align(listone,listtwo,offset):
    '''
    Aligns two instance lists based on an inputted offset value
    
    Parameters
    -----------
    listone, listtwo: list
        Lists with instances of two motifs that are going to be aligned
    offset: int
        Integer indicating how the first inputted motif's instances should be moved to maximize information content
    
    Results
    -----------
    merged_motif: list
        list of aligned and combined instances ready to be used to construct a merged motif
    
    Variables
    -----------
    S1 - Beginning of Instance in list 1
    E1 - End of Instance in list 1
    S2 - Beginning of Instance in list 2
    E2 - End of Instance in list 2
    k1 - Beginning of aligned motif in list 1
    p1 - End of aligned motif in list 1
    k2 - Beginning of aligned motif in list 1
    p2 - End of aligned motif in list 2
    '''
    merged_motif = []
    S1 = 0
    E1 = len(listone[0]) - 1
    S2 = 0
    E2 = len(listtwo[0]) - 1
    if offset <= 0:
        #Initializes variables for alignment bounds for each motif
        k2 = S2
        p2 = min(E1+offset,E2) + 1
        k1 = abs(offset)
        p1 = min(E1,E2-offset) + 1
    else:
        k2 = offset
        p2 = min(E2, E1+offset) + 1
        k1 = S1
        p1 = min(E2-offset,E1) + 1
    #Appends the aligned portion of the motif's instances to motifone_aligned and motiftwo_aligned, lists indicating the aligned portion of each motif. 
    for items in listone:
        motifone_aligned.append(items[k1:p1])
    for items in listtwo:
        motiftwo_aligned.append(items[k2:p2])
    merged_motif = motifone_aligned + motiftwo_aligned
    #Combines instances and returns merged_motif
    return merged_motif

### Functions below create the merged motif and calculates euclidean distance and KL divergence between two motifs to determine similarity/dissimilarity ###
def motif_col(motif,col):
    '''
    Creates a dictionary with the column values at any particular inputted column.

    Parameters
    ----------
    motif : Motif Object
        Inputted Motif for column creation
    col : Integer
        Integer value for particular location of dictionary creation (motif_col(motif,0) = dictionary with values of first column). 

    Returns
    -------
    col_dict : dict
        Dictionary where keys resemble the particular base (ACGT) and values resemble the pwm frequency of a particular column (held as a floaer value)
    '''
    col_dict = {letw:motif.pwm[letw][col] for letw in 'ACGT'}
    return col_dict

def calc_euclidean(cola, colb):
    '''
    Calculates the euclidean distance between two columns of a pwm.
    
    Parameters
    ----------
    cola, colb: [float]
        Two columns A and B, repsectively, of a PWM. 
    Returns
    -------
    distance: float
        The calculated distance between the two columns. 
    '''
    distance = math.sqrt(sum((cola[let]-colb[let])**2 for let in 'ACGT'))
    return distance  
def calc_kld_distance(cola, colb):
    '''
    Calculates the KL distance between two columns of a pwm.
    
    Parameters
    ----------
    cola, colb: [float]
        Two columns A and B, repsectively, of a PWM. 
    Returns
    -------
    distance: float
        The calculated distance between the two columns. 
    '''
    safe_log2 = lambda x: math.log(x, 2) if x != 0 else 0.0
    distance = (sum(cola[l] * safe_log2(cola[l] / colb[l]) for l in "ACTG" if colb[l] != 0) + 
            sum(colb[l] * safe_log2(colb[l] / cola[l]) for l in "ACTG" if cola[l] != 0))
    return distance

def eucl_dist(motifone,motiftwo):
    '''
    Calculates average euclidean distance between two motifs
    Parameters
    ----------
    motifone : Motif Object
        First motif to be compared
    motiftwo : Motif Object
        Second motif to be compared

    Returns
    -------
    None.

    '''
    eucl_distance = 0
    for cols in range(0,len(motifone)):
        eucl_distance += calc_euclidean(motif_col(motifone,cols),motif_col(motiftwo,cols))
    eucl_distance = eucl_distance/len(motifone)        
    print('The euclidean distance between motif one and motif two is: ' + str(eucl_distance))
    
def kl_dist(motifone,motiftwo):
    '''
    Calculates average KL divergence between two motifs

    Parameters
    ----------
    motifone : Motif Object
        First motif to be compared
    motiftwo : Motif Object
        Second motif to be compared
        
    Returns
    -------
    None.

    '''
    kl_distance = 0
    for cols in range(0,len(motifone)):
        kl_distance += calc_kld_distance(motif_col(motifone,cols),motif_col(motiftwo,cols))
    kl_distance = kl_distance/len(motifone)
    print('The KL distance between motif one and motif two is: ' + str(kl_distance))
    
def motif_merge(motifsone,motifstwo,distance_wanted = True):
    '''
    Merges two motifs using aligned instances

    Parameters
    ----------
    motifsone : Motif Object
        Motif with or without instances
    motifstwo : Motif Object
        Second Motif with or without instances
    distance_wanted: True or False Statement
        True or False parameter used to calculate EC distance and Kl divergence between component motifs and the combined motif

    Returns
    -------
    mixed_motif : Motif Object
        Motif composed of instances of both motifsone and motifstwo.

    '''
    motifsonelist = []
    motifstwolist = []
    #Used to reference the fake_instances function if instances are not reported in the motif
    if motifsone.instances == None:
        motifsonefakeinsta = fake_instances(motifsone)
        motifsone = motifs.create(motifsonefakeinsta,alphabet = 'ACGT')
    if motifstwo.instances == None:
        motifstwofakeinsta = fake_instances(motifstwo)
        motifstwo = motifs.create(motifstwofakeinsta,alphabet = 'ACGT')
    #With fake instances determined if necessary, get_alignment_offset function is called to determine alignment of instances that maximizes information content
    offset = get_alignment_offset(motifsone, motifstwo)
    for instances in motifsone.instances:
        motifsonelist.append(instances)
    for instancesw in motifstwo.instances:
        motifstwolist.append(instancesw)
    #Calls offset_align function with instances of both inputted motifs and offset value to save the combined aligned instances
    motiflist = offset_align(motifsonelist,motifstwolist,offset)
    #Creates a motif based on the combined, aligned instances
    mixed_motif = motifs.create(motiflist,'ACGT')
    if distance_wanted == True:
        motifone_edited = motifs.create(motifone_aligned,'ACGT')
        motiftwo_edited = motifs.create(motiftwo_aligned,'ACGT') 
        #Calculates eucl_distance and kl_dist between aligned motifs, comparing similar columns
        eucl_dist(motifone_edited,motiftwo_edited)
        kl_dist(motifone_edited,motiftwo_edited)
    return mixed_motif

#Functions to transform script into program

'''
def main():
    if len(sys.argv) <= 2:
        print('Lack of info')
    else:
        function = sys.argv[1]
        if function == "motif_create":
            if len(sys.argv) < 4:
                print(sys.argv[2])
            else:
                file_name = sys.argv[2]
                i = int(sys.argv[3])
                motif_create(file_name,i)
        if function == "motif_merge":
            if len(sys.argv) < 6:
                print(sys.argv[2])
            else:
                motifsone = sys.argv[2]
                i1 = int(sys.argv[3])
                motifstwo = sys.argv[4]
                i2 = int(sys.argv[5])
                print(motif_merge(motif_create(motifsone,i1),motif_create(motifstwo,i2)))
if (__name__ == '__main__'):
    main()
'''
