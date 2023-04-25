# Alon Luboshitz 312115090
''' Transcribe function - Dna seq into Rna.
1. Get dna seq
2. save length -> save length parameters 5', 3'.
3. read Dna from the end and append current letter to Rna seq.
4. return Rna seq '''
# check for empty string?
def transcribe(dna_seq):
    casefold_seq = dna_seq.casefold()
    # iniating replacment dictoinary
    replace_dic = {'c' : 'G', 't' : 'A', 'a' : 'U', 'g': 'C' }
    for char in replace_dic.keys():
        casefold_seq = casefold_seq.replace(char, replace_dic[char])
    print(casefold_seq[::-1])

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
def translate(rna_seq,position):
    if position < 0 or position > 2:
        position = 0
    aa_dict = createAminoAcidsDict()
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
                temp_seq += 'M'
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
            temp_seq += ';' + codon_letter
    if aa_seq == 0:
        return 'Non-coding RNA'
    else:
        return aa_seq



    
def createAminoAcidsDict():
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
    aa_dic.update(convert(S,'S'))
    aa_dic.update(convert(F,'F'))
    aa_dic.update(convert(L,'L'))
    aa_dic.update(convert(I,'I'))
    aa_dic.update(convert(V,'V'))
    aa_dic.update(convert(P,'P'))
    aa_dic.update(convert(T,'T'))
    aa_dic.update(convert(A,'A'))
    aa_dic.update(convert(Y,'Y'))
    aa_dic.update(convert(Z,'Z'))
    aa_dic.update(convert(H,'H'))
    aa_dic.update(convert(Q,'Q'))
    aa_dic.update(convert(N,'N'))
    aa_dic.update(convert(K,'K'))
    aa_dic.update(convert(D,'D'))
    aa_dic.update(convert(E,'E'))
    aa_dic.update(convert(C,'C'))
    aa_dic.update(convert(R,'R'))
    aa_dic.update(convert(G,'G'))
    aa_dic.update(UGG = 'W')
    aa_dic.update(AUG = 'M')
    return aa_dic
def convert(aa,letter):
    dict = {key: letter for key in aa}
    return dict
translate('AUGUAACGUG',2)