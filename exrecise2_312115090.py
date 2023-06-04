# Alon Luboshitz 312115090
import sys
import csv
'''The cell class holds a name (type of cell) and genome.
genome is list of tuples represnted by dna seq and reading frames accordingly.
cell has a presentation of <Celltype,Number of seq in genome>
it has find_ssr, transcribe, translate functions detailed in the class.'''
class Cell:
    def __init__(self,name,genome) -> None:
        self.__name = name
        self.__genome =  genome
        if len(genome) > 0 :  
            self.__numseq = len(genome[0])
        else: self.__numseq = None
    def __str__(self) -> str:
        return '<' + self.__name + ', ' + str(self.__numseq) + '>'
    '''
    given a dna seq check what is most repetative seq with max size of 6.
    we'll do a naive solution going through every index from 0 - end.
    for each index create 1-6 size substring and check if repeats x times.
    if it does (3 or more) add to dict :
    check if exists in dict alrdy - if bigger exchange.'''
    def __find_ssr__(self,dna_seq):
        ssr_dict = {}
        dna_seq = dna_seq.upper()
        # size of substr 1 - 6
        for j in range(1,6):
            # iterating over dna_seq through offset
            for i in range(0, len(dna_seq), 1):
                #validating index
                if not len(dna_seq[i:i+j])%j == 0:
                    break
                #iterating over dna_seq via substr size
                ssr = dna_seq[i:i+j]
                count = 1
                lenn = len(dna_seq)
                for k in range(i+j,len(dna_seq),j):
                    #validating next index
                    if not len(dna_seq[k:k+j])%j == 0:
                        break
                    next_ssr = dna_seq[k:k+j]
                    if ssr == next_ssr:
                        count += 1
                        if k == (lenn - j) and count > 2:
                            ssr_dict = self.__update_count_dict__(count,ssr_dict,ssr)
                    else:
                        #update dictionary
                        if count > 2:
                            ssr_dict = self.__update_count_dict__(count,ssr_dict,ssr)
                        ssr = next_ssr  
                        count = 1    
    
        return ssr_dict
    '''
    this function converts dict into a string
    represnting the keys sorted separted by ',' then their value then ';'
    '''
    def __convert_dic_repsentation__(self,dict):
        if not dict:
            return 'No simple repeats in DNA sequence'
        else:
            text =''
            myKeys = list(dict.keys())
            myKeys.sort()
            for i in myKeys:
                text = text + i + ',' + str(dict[i]) + ';'
            return text[:-1]
    
    '''
    this function is for ssr function, updating given count to a ssr'''
    def __update_count_dict__(self,count,ssr_dict,ssr):
        dict_count = ssr_dict.get(ssr)
        if dict_count is None or count > dict_count:
            ssr_dict[ssr] = count 
        return ssr_dict
    
    '''
    this function will take 2 args.
    first translate rna_seq to AA letters and assign "Z" as a stop codon
    then determine which seq is the longest by iterating over the letters.
    return the longest seq if exists.
    set temp length_seq = 0, protien_seq =empty.
    iterate in jumps of 3 from the position given.
    map 3 Nuc into AA letter. -> translate all AA - set special letter for STOP
    boolean value for new_seq or not.
    if new seq and M -> create new_seq, append till STOP.
    set length.
    if new_length>old_length set new seq to longest.
    return seq.'''
    def __translate__(self,rna_seq,position):
        if position < 0 or position > 2:
            position = 0
        aa_dict = self.__createAminoAcidsDict__()
        upper_rna = rna_seq.upper()
        aa_seq = ''
        aa_length = 0
        new_seq = True
        temp_length = 0
        temp_seq = ''
        for i in range(position, len(upper_rna), 3):
            codon = upper_rna[i:i+3]
            if not len(codon)%3 == 0:
                break 
            codon_letter = aa_dict[codon]
            if new_seq:
                if codon_letter == 'M':
                    temp_seq += 'M' + ';'
                    new_seq = False
                    temp_length = 1
                else:
                    pass
            elif codon_letter == 'Z':
                #termination of protein
                if temp_length > aa_length:
                    aa_seq = temp_seq
                    aa_length = temp_length
                #new protien shorter then previuos keep old and reset new
                temp_seq = ''
                temp_length = 0
                new_seq = True
            else:
                temp_length += 1
                temp_seq += codon_letter + ';'
        if temp_length > aa_length:
            aa_seq = temp_seq
        return aa_seq[:-1]
            
    def __createAminoAcidsDict__(self):
        # set lists for aa by letter.
        S = ['UCU','UCC','UCA','UCG','AGU','AGC']
        F = ['UUU','UUC']
        L = ['UUA','UUG','CUU','CUC','CUA','CUG']
        I = ['AUU','AUC','AUA']
        V = ['GUU','GUC','GUA','GUG']
        P = ['CCU','CCC','CCA','CCG']
        T = ['ACU','ACC','ACA','ACG']
        A = ['GCU','GCC','GCA','GCG']
        Y = ['UAU','UAC']
        Z = ['UAA','UAG','UGA']
        H = ['CAU','CAC']
        Q = [ 'CAA','CAG']
        N = ['AAU','AAC']
        K = ['AAA','AAG']
        D = [ 'GAU','GAC']
        E = ['GAA','GAG']
        C = ['UGU','UGC']
        R = ['CGU','CGC','CGA','CGG','AGA','AGG']
        G = ['GGU','GGC','GGA','GGG']
        # set empty dic
        aa_dic = {}
        aa_dic.update(self.__convert__(S,'S'))
        aa_dic.update(self.__convert__(F,'F'))
        aa_dic.update(self.__convert__(L,'L'))
        aa_dic.update(self.__convert__(I,'I'))
        aa_dic.update(self.__convert__(V,'V'))
        aa_dic.update(self.__convert__(P,'P'))
        aa_dic.update(self.__convert__(T,'T'))
        aa_dic.update(self.__convert__(A,'A'))
        aa_dic.update(self.__convert__(Y,'Y'))
        aa_dic.update(self.__convert__(Z,'Z'))
        aa_dic.update(self.__convert__(H,'H'))
        aa_dic.update(self.__convert__(Q,'Q'))
        aa_dic.update(self.__convert__(N,'N'))
        aa_dic.update(self.__convert__(K,'K'))
        aa_dic.update(self.__convert__(D,'D'))
        aa_dic.update(self.__convert__(E,'E'))
        aa_dic.update(self.__convert__(C,'C'))
        aa_dic.update(self.__convert__(R,'R'))
        aa_dic.update(self.__convert__(G,'G'))
        aa_dic.update(UGG = 'W')
        aa_dic.update(AUG = 'M')
        return aa_dic
    def __convert__(self,aa,letter):
        dict = {key: letter for key in aa}
        return dict
    ''' Transcribe function - Dna seq into Rna.
    1. Get dna seq
    2. replace corsponding letters using build in dictorionary.
    3. Read the seq backwards 5-3 --> 3-5.
    4. return Rna seq '''
    # check for empty string?
    def __transcribe__(self,dna_seq):
        casefold_seq = dna_seq.casefold()
        # iniating replacment dictoinary
        replace_dic = {'c' : 'G', 't' : 'A', 'a' : 'U', 'g': 'C' }
        for char in replace_dic.keys():
            casefold_seq = casefold_seq.replace(char, replace_dic[char])
        return(casefold_seq[::-1])
    def _find_srr_(self,dna_position):
        #unpack tuple in position
        dna_seq = self.__genome[0][dna_position%self.__numseq]
        dict = self.__find_ssr__(dna_seq)
        return self.__convert_dic_repsentation__(dict)
    def _translate_(self,dna_position):
        #dna_tuple = self.__genome[(dna_position%self.__numseq)]
        rna_seq = self._transcribe_(dna_position)
        reading_frame = self.__genome[1][dna_position%self.__numseq]
        translation =  self.__translate__(rna_seq,int(reading_frame))
        if len(translation) == 0:
            return 'Non-coding RNA'
        else:
            return 'Translation: ' + translation
    def _transcribe_(self,dna_position):
        dna_seq = self.__genome[0][dna_position%self.__numseq]
        return self.__transcribe__(dna_seq)
    '''this function return a list with all the genome ssr and translation for this dna.
    if the genome is empty return empty list'''
    def _repertoire_(self):
        rep_list = []
        if self.__numseq is None:
            return rep_list
        else:
            for i in range(self.__numseq):
                ssr_translate_tuple = (self._find_srr_(i),self._translate_(i))
                rep_list.append(ssr_translate_tuple)
        return rep_list

    def get_genome(self):
        return self.__genome
    def get_name(self):
        return self.__name
    
'''this is a stem cell class
it is the precursor class for diffrenet type of cells.
it can multiple it self, mitosis (split to 2), and differntiate into a difrrenet type of cell.
it has two nested cell classes atm: Muscle and Nerve each with its own functions.'''
class StemCell(Cell):
    def __init__(self):
        pass
    def __init__(self, name, genome) -> None:
        super().__init__(name, genome)
    def __mul__(self, P):
        assert P > 0 ,"Cannot multiple Stem Cell by zero or negative number"
        assert isinstance(P,int),"cannot multiple stemcell by non int number"
        stem_cells = []
        stem_cells.append(self)
        for i in range(P-1):
            stem_cells.append(StemCell(self.get_name(),self.get_genome()))
        return stem_cells
    def mitosis(self):
        return self * 2
    def differentiate(self, cell_name, args):
        if (cell_name == 'Nerve Cell'):
            return self.NerveCell(self,args)
        elif (cell_name == 'Muscle Cell'):
            return self.MuscleCell(self,args)
        else: assert False, 'trying to init an unkown cell, please try again.'
    class MuscleCell(Cell):
        def __init__(self):
            pass
        def __init__(self,Stem_cell, args) -> None:
            assert isinstance(Stem_cell,StemCell),"Cannot create muscle not from Stem cell"
            super().__init__('Muscle Cell', Stem_cell.get_genome())
            self.__treshold = args[1]
            self.__file_path = args[0]
        def receive(self,signal):
            #write to file if signal> treshold
            if float(signal) > float(self.__treshold):
                with open(self.__file_path,"a") as file:
                    line = str(signal) + ', I like to move it\n'
                    file.write(line)
            else: pass
    class NerveCell(Cell):
        def __init__(self):
            pass
        def __init__(self, Stem_cell, P) -> None:
            assert isinstance(Stem_cell,StemCell),"Cannot create nerve cell not from Stem cell"
            super().__init__('Nerve Cell', Stem_cell.get_genome())
            self.__coefficient = P
            self.__signal = 0
        def receive(self, signal):
            self.__signal = signal
        def send(self):
            return self.__signal * self.__coefficient
'''this is a managing class of nerve cells. the network class gets a cell list and split
them into nerve cells and muscle cell.
it can send a signal through the nerve cells and print the cell representation.
'''
class NerveNetwork:
    def __init__(self, cell_list) -> None:
        self.__cell_list = cell_list
        self.__nerve_cells = []
        self.__muscle_cell = []
        self.split_to_cells()
    def send_signal(self,signal):
        for i in self.__nerve_cells:
            i.receive(signal)
            signal = i.send()
        self.__muscle_cell[0].receive(signal)
    def __str__(self) -> str:
        str = ''
        for i in self.__nerve_cells:
            str = str + (i.__str__() + '\n')
        for i in self.__muscle_cell:
            str = str + (i.__str__() + '\n')
        str = str + self.convert_cell_repertoire(self.__muscle_cell[0])
        return str
    def split_to_cells(self):
        for cell in self.__cell_list:
            if cell.get_name() == 'Nerve Cell':
                self.__nerve_cells.append(cell)
            else : self.__muscle_cell.append(cell)
    def get_muscle_cell(self):
        return self.__muscle_cell[0]
    
    '''getting a cell convert its repertoire to desired string
    repertoire --> list of tuples, unpack and condece with new line between them'''
    def convert_cell_repertoire(self, cell):
        str = ''
        for pair in cell._repertoire_():
            (ssr, protein)  = pair
            str = str + ssr + '\n' + protein + '\n'
        return str[:-1]



'''this function gets input_path for an input file and validates its according to instrcutions:
for each currect line append the corssepoding cell to a list.
1. check there are only 4 headline. and that they are from type,dna,frame,param
2. check all cells are from MC//NC type
3. check that all the dna_seq are from A,T,C,G building blocks.
4. check that all reading frame are between 1-3.
5. check that there are matching amount of reading frames to dna_seq number.
6. check the parameter field MC //NC 
'''
def validate_file(input_path):
    with open(input_path, 'r') as tsbfile:
        csv_reader = csv.DictReader(tsbfile,delimiter='\t')
        keys = csv_reader.fieldnames
        valid_keys = check_headlines(keys)
        assert valid_keys[0],"File illegal"
        cells = []
        for row in csv_reader:
            type_check = check_type(row,valid_keys[1][0])
            assert type_check[0],"File illegal"
            dna_check = check_DNA(row,valid_keys[1][1])
            assert dna_check[0],"File illegal"
            frame_check = check_readingframe(row,valid_keys[1][2])
            assert frame_check[0],"File illegal"
            assert int(frame_check[2]) == int(dna_check[2]),"File illegal"
            param_check = check_parameter(row,valid_keys[1][3],type_check[2])
            assert param_check[0],"File illegal"
            create_cell(type_check[1],dna_check[1],frame_check[1],param_check[1],cells)
        return cells    
        
'''validate headline are tpye,dna,reading_frame,parameter no matter their location in the headline.
return list with keys by this order.'''         
def check_headlines(keys):
    look_for = ['type','dna','reading_frames','parameter']
    valid_heads = []
    for x in look_for:
        for y in keys:
            if y.casefold().find(x) != -1:
                valid_heads.append(y)
                break
    if len(valid_heads) == len(look_for):
        return (True,valid_heads)
    else : return (False,)
    
'''All check functions get the corrsepoding context via the key.
key been validated by check_headlines.
each function check by directions represented in the note of validate file'''

def check_type(row, look_for):
    name = row.get(look_for) 
    if name == 'NC':
        return (True,name,1)
    elif name == 'MC':
        return (True,name,2)
    else: return (False,)

def check_DNA(row,look_for):
    dna_seq = row.get(look_for)
    #split the dna by "comma"
    dna_seq = dna_seq.split(",")
    num_seq = len(dna_seq)
    nucleotide_set = {'T','G','A','C'}
    for dna in dna_seq:
        #make set of dna_seq and check if the sets equals
        dna_set = set(dna.upper())
        if not dna_set.issubset(nucleotide_set):
            #dna_seq has other letters then a,t,g,c
            return (False,)
        # all seqs have valid chars
        return (True,dna_seq,num_seq)
    
def check_readingframe(row,look_for):         
    reading_frame = row.get(look_for)
    try:
        reading_frame = [float(x) for x in reading_frame.split(",")]
    except: return (False,)
    else:
        num_frames = len(reading_frame)
        valid_frame = [0,1,2]
        #check frames are integers and are valid - 1,2,3
        for x in reading_frame:
            if not x.is_integer():
                return (False,) 
            elif x not in valid_frame:
                return (False,)
        return (True,reading_frame,num_frames)
'''checks parameter by cell type.
1 = NC, 2 = MC'''
def check_parameter(row,look_for,cell_type):
    param = row.get(look_for)
    param = param.split(",")
    num_of_params = len(param)
    #neuron cell, assuming 1 param - float
    if cell_type == 1:
        if num_of_params != 1:
            return (False,)
        try:
            param = float(param[0])
        except:
            return (False,)
        else: 
            if  param < 0:
                #param == int (param) or
                #not decimal number\ negative
                return (False,)
            else: return (True,param)
    else: #cell type 2 MC
        if num_of_params != 2:
            return (False,) 
        try:
            treshold = float(param[1])
        except:
            return (False,)
        else:
            if treshold < 0:
                return (False,) 
        with open(param[0],"w") as inputfile:
            return (True,param)

'''validates there are two arguments, file_path can be opened.
all signals are int and positive'''
def check_args(args):
    assert len(args) == 3, "Wrong number of arguments"
    #validate file can be opened
    with open(args[1],"r") as file:
        pass
    try:
        signals = [float(i) for i in args[2].split(",")]
    except:
        assert False,"signal argument isnt a number"
    else:
            
        for signal in signals:
            assert signal>0,"Negative signal!"
            assert signal.is_integer(),"Non integer signal"
        return signals

'''this function gets a line from the file and creates the intended cell.
line has laready been validated.'''
def create_cell(name,dna_seq,frames,param,cells): 
    genome = (dna_seq,frames) 
    cell = StemCell('name',genome)
    if name == 'NC':
        #nc
        for each_cell in cell.mitosis():
            cells.append(each_cell.differentiate('Nerve Cell',param))
    else : #must be MC
        cells.append(cell.differentiate('Muscle Cell',param))

'''main function - gets two args, validate them create the nerve network
with cells from the input file. prints the network and send signals to the network.'''          
if __name__ == '__main__':
    signals = check_args(sys.argv)
    cells = validate_file(sys.argv[1])
    network = NerveNetwork(cells)
    print(network)
    for i in signals:
        network.send_signal(i)
    

    
       
   
