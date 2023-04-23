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


