#clang C++
#clib braiding

from libcpp.list cimport list

cdef extern from "braiding.h" namespace "Braiding":
    list[list[int]] ConjugatingBraid (int n, list[int] word, list[int] word2)
    list[list[int]] LeftNormalForm (int n, list[int] word)
    list[list[int]] GreatestCommonDivisor(int n, list[int] word1, list[int] word2)
    list[list[int]] LeastCommonMultiple(int n, list[int] word1, list[int] word2)
    list[list[list[int]]] CentralizerGenerators(int n, list[int] word)
    list[list[list[int]]] SuperSummitSet(int n, list[int] word)
    
def conjugatingbraid(braid1, braid2):
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[list[int]] rop = ConjugatingBraid(nstrands, l1, l2)
    sig_off()
    return rop

def leftnormalform(braid):
    nstrands = braid.parent().strands()
    l1 = braid.Tietze()
    sig_on()
    cdef list[list[int]] rop = LeftNormalForm(nstrands, l1)
    sig_off()
    return rop
    
def greatestcommondivisor(braid1, braid2):
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[list[int]] rop = GreatestCommonDivisor(nstrands, l1, l2)
    sig_off()
    return rop
    
def leastcommonmultiple(braid1, braid2):
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[list[int]] rop = LeastCommonMultiple(nstrands, l1, l2)
    sig_off()
    return rop
    
def centralizer(braid):
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = CentralizerGenerators(nstrands, l)
    sig_off()
    return rop
    
def supersummitset(braid):
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = SuperSummitSet(nstrands, l)
    sig_off()
    return rop