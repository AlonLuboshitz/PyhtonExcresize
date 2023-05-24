# Alon Luboshitz 312115090
import sys
from exrecise1_312115090 import *
class Cell:
    def __init__(self,name,genome) -> None:
        self.__name = name
        self.__genome =  genome
        self.__numseq = len(genome)
    def __str__(self) -> str:
        return '<' + self.__name + ', ' + str(self.__numseq) + '>'
    def _find_srr_(self,dna_position):
        #unpack tuple in position
        dna_seq = self.__genome[(dna_position%self.__numseq)[0]]
        return find_ssr(dna_seq)
    def _translate_(self,dna_position):
        pass
    def _transcribe_(self,dna_position):
        pass
    def _repertoire_(self):
        pass
    def get_genome(self):
        return self.__genome
    def get_name(self):
        return self.__name

class StemCell(Cell):
    def __init__(self, name, genome) -> None:
        super().__init__(name, genome)
    def __mul__(self, P):
        pass
    def mitosis(self):
        return self * 2
    def differentiate(self, cell_name, args):
        return Cell
    class MuscleCell(Cell):
        # init muscle cell from 
        def __init__(self,Stem_cell, args) -> None:
            super().__init__('Muscle Cell', Stem_cell.get_genome())
            self.__treshold = args[1]
            self.__file_path = args[0]
        def receive(signal):
            #write to file if signal> treshold
            pass
    class NerveCell(Cell):
        def __init__(self, Stem_cell, P) -> None:
            super().__init__('Nerve Cell', Stem_cell.get_genome())
            self.__coefficient = P
        def receive(self, signal):
            return signal * self.__coefficient
        def send(self, signal):
            return signal * self.__coefficient
class NerveNetwork:
    def __init__(self, cell_list) -> None:
        self.__cell_list = cell_list
        self.__nerve_cells = None
        self.__muscle_cell = None
    def send_signal(self,signal):
        for i in self.__nerve_cells:
            signal = i.receive(signal)
            signal = i.send(signal)
        self.__muscle_cell.receive(signal)
    def __str__(self) -> str:
        for i in self.__cell_list:
            str = (i + '\n')
        return str


    

        