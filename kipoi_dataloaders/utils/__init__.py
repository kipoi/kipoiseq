# alphabets:
DNA = ["A", "C", "G", "T"]
RNA = ["A", "C", "G", "U"]
AMINO_ACIDS = ["A", "R", "N", "D", "B", "C", "E", "Q", "Z", "G", "H",
               "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


from collections import OrderedDict
from .conversion import onehot, onehot_trafo, ReshapeSeq, get_string_trafo
from .readers import TsvReader, FastaReader
