## RENAME this file YourLastName_OOP_FinalProject_2023.py

##Assignment: Add to the constructor and methods of a parent class and child classes
##            which inherit the base class properties. NOTE: You are not allowed
##            to import any specialized libraries for this project (e.g., no Biopython)
##            The idea is for you to write these methods from scratch.

## Begin with the parent Seq class and the child DNA class we created in lecture below.
## 

### Seq Class
#
#  Constructor:
#  (1) Use the string functions upper and strip to clean up self.sequence.
#  (2) Add a variable self.kmers to the constructor and make it equal to an empty list.

#  Methods:
#  (1) Add a method called make_kmers that makes overlapping kmers of a given length from self.sequence
#      appends these to self.kmers. Default kmer parameter=3.
#  (2) Add a method called fasta that returns a fasta formatted string like this:
#      >species gene
#      AGATTGATAGATAGATAT


### DNA Class: INHERITS Seq class
#   
#  Constructor:
#  Use re.sub to change any non nucleotide characters in self.sequence into an 'N'.
#      re.sub('[^ATGCU]','N',sequence) will change any character that is not a
#      capital A, T, G, C or U into an N. (Seq already uppercases and strips.)

#  Methods:
#  (1) Add a method called print_info that is like print_record, but adds geneid and an
#      empty space to the beginning of the string.
#  (2) Add a method called reverse_complement that returns the reverse complement of
#      self.sequence
#  (3) Add a method called six_frames that returns all 6 frames of self.sequence
#      This include the 3 forward frames, and the 3 reverse complement frames

### RNA Class:  INHERITS DNA class
#  
#  Construtor:
#  Use the super() function (see DNA Class example).
#  (1) Automatically change all Ts to Us in self.sequence. 
#  (2) Add self.codons equals to an empty list

#  Methods:
#  (1) Add make_codons which breaks the self.sequence into 3 letter codons
#      and appends these codons to self.codons unless they are less than 3 letters long.
#  (2) Add translate which uses the Global Variable standard_code below to
#      translate the codons in self.codons and returns a protein sequence.

### Protein Class: INHERITS Seq class
#
#  Construtor:
#  Use the super() function (see DNA Class example).
#  Use re.sub to change any non LETTER characters in self.sequence into an 'X'.

#  Methods:
#  The next 2 methods use a kyte_doolittle and the aa_mol_weights dictionaries.
#  (2) Add total_hydro, which return the sum of the total hydrophobicity of a self.sequence
#  (3) Add mol_weight, which returns the total molecular weight of the protein
#      sequence assigned to the protein object. 


import re
import doctest

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}

aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}


class Seq:

    def __init__(self,sequence,gene,species, kmers=[]):
        self.sequence=sequence
        self.gene=gene
        self.species=species
        self.kmers=[]
        self.sequence=self.sequence.upper() #upper cases all nucleotides
        self.sequence=self.sequence.strip() #strip ends of sequence


    def __str__(self):
        return self.sequence
    
    def make_kmers(self, k=3): #makes kmers 
        """This function returns a the self.sequence with kmer length 5
        >>> s=Seq("  gATATAGGACctttaGGACCAC  ","my_gene","H.sapiens")
        >>> s.make_kmers(5)
        ['GATAT', 'ATATA', 'TATAG', 'ATAGG', 'TAGGA', 'AGGAC', 'GGACC', 'GACCT', 'ACCTT', 'CCTTT', 'CTTTA', 'TTTAG', 'TTAGG', 'TAGGA', 'AGGAC', 'GGACC', 'GACCA', 'ACCAC']
        """
        for i in range(0, len(self.sequence)):
            kmer=self.sequence[i: i+k]
            if len(kmer)==k:
                self.kmers.append(kmer)

        return self.kmers
                
    def fasta(self): #formats sewuences into fasta file
        """This function returns fasta format
        >>> s=Seq("  gATATAGGACctttaGGACCAC  ","my_gene","H.sapiens")
        >>> print(s.fasta())
        >H.sapiens my_gene
        GATATAGGACCTTTAGGACCAC
        """
        return ">"+self.species+" "+self.gene+"\n"+ self.sequence


    def print_record(self):
        #print(self.species + " " + self.gene + ": " + self.sequence)
        """This function returns a the self.sequence as a string
        >>> x=Seq("ATATAG","my_gene","H.sapiens")   
        >>> print(x.print_record())
        ATATAG
        """
        return self.sequence

    
class DNA(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence=self.sequence.upper()
        self.sequence=self.sequence.strip()
        self.geneid=geneid #provides geneid
        self.sequence=re.sub('[^ATGCU]','N',self.sequence) #change any nonnulceotide to 'N'
 
    def analysis(self):
        """This function returns number of GC
        >>> d=DNA("  -tcaaaGCGGCGGATCTCCCaaatga","my_dna","D.terebrans","AX5667")
        >>> print(d.analysis())
        13
        """
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence)) #count all GC
        return gc

    def print_info(self):
        #print(self.geneid + " " +self.species + " " + self.gene + ": " + self.sequence)
        """This function prints geneID, species, gene and sequence
        >>> d=DNA("  -tcaaaGCGGCGGATCTCCCaaatga","my_dna","D.terebrans","AX5667")
        >>> print(d.print_info())
        AX5667 D.terebrans my_dna: NTCAAAGCGGCGGATCTCCCAAATGA
        """

        return self.geneid + " " +self.species + " " + self.gene + ": " + self.sequence
        

    def reverse_complement(self): #finds reverse complement of the sequence
        """This function returns reverse complement of the sequence
        >>> d=DNA("  -tcaaaGCGGCGGATCTCCCaaatga","my_dna","D.terebrans","AX5667")
        >>> rc=d.reverse_complement()
        >>> print(rc)
        TCATTTGGGAGATCCGCCGCTTTGAN
        """
        new_sequence=''
        for i in self.sequence:
            if i == "T":
                new_sequence+="A"
            elif i=="A":
                new_sequence+="T"
            elif i=="G":
                new_sequence+="C"
            elif i=="C":
                new_sequence+="G"
            else:
                new_sequence+="N"
            
        reverse=new_sequence[::-1]
        return reverse

    def six_frames(self): #finds six frames of the sequence
        """This function returns all six frames of the sequence
        >>> d=DNA("  -tcaaaGCGGCGGATCTCCCaaatga","my_dna","D.terebrans","AX5667")
        >>> all_6_frames=d.six_frames()
        >>> print(all_6_frames) 
        ['NTCAAAGCGGCGGATCTCCCAAATGA', 'TCAAAGCGGCGGATCTCCCAAATGA', 'CAAAGCGGCGGATCTCCCAAATGA', 'TCATTTGGGAGATCCGCCGCTTTGAN', 'CATTTGGGAGATCCGCCGCTTTGAN', 'ATTTGGGAGATCCGCCGCTTTGAN']
        """
        frames = []
        for i in range(3):
            frames.append(self.sequence[i:])
        reverse = self.reverse_complement()
        for i in range(3):
            frames.append(reverse[i:])
        return frames


class RNA(DNA):

    def __init__(self, sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species,geneid)
        self.sequence=self.sequence.upper()
        self.sequence=self.sequence.strip()
        self.sequence=self.sequence.replace("T", "U") #replace Ts with Us
        self.codons=[]
        
    def make_codons(self): #make a set of codons from the sequence
        """This function returns all codons of the sequence
        >>> r=RNA("  g?ATATAGGACctttaGGACCAC  ","my_rna","G.gallus","R5990999")
        >>> print(r.make_codons())
        ['GNA', 'UAU', 'AGG', 'ACC', 'UUU', 'AGG', 'ACC']
        """
        for i in range(0, len(self.sequence), 3):
            codon=self.sequence[i: i+3]
            if len(codon)==3:
                self.codons.append(codon)
        return self.codons
    
    def translate(self): #translate codons to amino acids
        """This function returns translation of sequence to amino acids
        >>> r=RNA("  g?ATATAGGACctttaGGACCAC  ","my_rna","G.gallus","R5990999")
        >>> print(r.translate())
        <BLANKLINE>
        """
        protein=''
        for i in self.codons:
            try:
                value=standard_code[i]
                protein+=value
            except KeyError:
                protein+='X'
        return protein
        
        


class Protein(Seq):

    def __init__(self, sequence, gene, species, geneid):
        super().__init__(sequence,gene,species)
        self.sequence=re.sub('[^ARNDCQEGHILKMFPSTWYV]','X',self.sequence)
        self.geneid=geneid

    def total_hydro(self): #calculates total hydrophobicity
        """This function returns total hydrophobicity
        >>> testp=Protein('VIKING','test','unknown',999)
        >>> x=testp.total_hydro()
        >>> print(x)
        5.399999999999999
        """
        total_hydro=0
        for i in self.sequence:
            value=kyte_doolittle[i]
            total_hydro+=value
        return total_hydro
            

    def mol_weight(self): #calculates total molecular weight
        """this function returns molecular weight
        >>> testp=Protein('VIKING','test','unknown',999)
        >>> m=testp.mol_weight()
        >>> print(m)
        732.8699999999999
        """
        total_weight=0
        for i in self.sequence:
            value=aa_mol_weights[i]
            total_weight+=value
        return total_weight

if __name__=="__main__":
    

    doctest.testmod(verbose=True)







