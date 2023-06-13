# Alon Luboshitz 312115090
import sys
''' Transcribe function - Dna seq into Rna.
1. Get dna seq
2. replace corsponding letters using build in dictorionary.
3. Read the seq backwards 5-3 --> 3-5.
4. return Rna seq '''
# check for empty string?
def transcribe(dna_seq):
    casefold_seq = dna_seq.casefold()
    # iniating replacment dictoinary
    replace_dic = {'c' : 'G', 't' : 'A', 'a' : 'U', 'g': 'C' }
    for char in replace_dic.keys():
        casefold_seq = casefold_seq.replace(char, replace_dic[char])
    return(casefold_seq[::-1])

'''
this function is for ssr function, updating given count to a ssr'''
def update_count_dict(count,ssr_dict,ssr):
    dict_count = ssr_dict.get(ssr)
    if dict_count is None or count > dict_count:
        ssr_dict[ssr] = count 
    return ssr_dict
        
'''
given a dna seq check what is most repetative seq with max size of 6.
we'll do a naive solution going through every index from 0 - end.
for each index create 1-6 size substring and check if repeats x times.
if it does (3 or more) add to dict :
check if exists in dict alrdy - if bigger exchange.'''
def find_ssr(dna_seq):
    ssr_dict = {}
    dna_seq = dna_seq.upper()
    # size of substr 1 - 6
    for j in range(1,7):
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
                        ssr_dict = update_count_dict(count,ssr_dict,ssr)
                else:
                    #update dictionary
                    if count > 2:
                        ssr_dict = update_count_dict(count,ssr_dict,ssr)
                    ssr = next_ssr  
                    count = 1    
  
    return ssr_dict
'''
this function converts dict into a string
represnting the keys sorted separted by ',' then their value then ';'
'''
def convert_dic_repsentation(dict):
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
    return aa_seq[:-1]
    
    
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
if __name__ == '__main__':
    # if not len(sys.argv) == 5:
    #     print('ERROR! wrong number of arguments')
    #     exit(0)
    # else:
    #     ssr_dict = find_ssr(sys.argv[1])
    #     print(convert_dic_repsentation(ssr_dict))
    #     print('RNA sequence: ' + transcribe(sys.argv[2]))
    #     mrna= translate(sys.argv[3],int(sys.argv[4]))
    #     if len(mrna) == 0:
    #         print('Non-coding RNA')
    #     else:
    #         print('Translation: ' + mrna) 
    
    ssr_dict = find_ssr("ATGCCGAATGCCGAATGCCGAATGCCGATGCCGATGCCG")
    print(convert_dic_repsentation(ssr_dict))
