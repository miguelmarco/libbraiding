#clang C++
#clib braiding

from libcpp.list cimport list

cdef extern from "braiding.h" namespace "Braiding":
    list[int] ConjugatingBraid (int n, list[int] word, list[int] word2)
    
    
def conjugatingbraid(braid1, braid2):
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[int] rop = ConjugatingBraid(nstrands, l1, l2)
    sig_off()
    return rop
