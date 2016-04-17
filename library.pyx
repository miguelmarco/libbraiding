#clang C++
#clib braiding
r"""
Cython wrapper for the libbraiding library.

The libbraiding library is a modification of the braiding program
by Juan Gonzalez-Meneses (https://github.com/jeanluct/cbraid)
to expose the functions as a C++ library instead of an interactive
program.


Braids are returned in left normal form as a list of lists. The 
first list contains only an integer, representing the power of
Delta. The subsequent lists are the Tietze lists of the elementary
permutation braids.
"""

from libcpp.list cimport list


cdef extern from "braiding.h" namespace "Braiding":
    list[list[int]] ConjugatingBraid (int n, list[int] word, list[int] word2)
    list[list[int]] LeftNormalForm (int n, list[int] word)
    list[list[int]] RightNormalForm (int n, list[int] word)
    list[list[int]] GreatestCommonDivisor(int n, list[int] word1, list[int] word2)
    list[list[int]] LeastCommonMultiple(int n, list[int] word1, list[int] word2)
    list[list[list[int]]] CentralizerGenerators(int n, list[int] word)
    list[list[list[int]]] SuperSummitSet(int n, list[int] word)
    
def conjugatingbraid(braid1, braid2):
    r"""
    Return a braid that conjugates braid1 to braid2, if such a braid exists.
    
    INPUT:
    
    - ``braid1`` -- the braid to be conjugated.
    
    - ``braid2`` -- the braid to conjugate to.
    
    OUTPUT:
    
    The list of lists that represent a conjugating braid. If the input braids
    are not conjugate, an empty list is returned.
    """
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[list[int]] rop = ConjugatingBraid(nstrands, l1, l2)
    sig_off()
    return rop

def leftnormalform(braid):
    r"""
    Return the left normal form of a braid.
    
    INPUT:
    
    - ``braid`` -- a braid
    
    OUTPUT:
    
    A list of lists with the left normal form. The first list contains the
    power of delta. The subsequent lists are the elementary permutation braids.
    """
    nstrands = braid.parent().strands()
    l1 = braid.Tietze()
    sig_on()
    cdef list[list[int]] rop = LeftNormalForm(nstrands, l1)
    sig_off()
    return rop

def rightnormalform(braid):
    r"""
    Return the right normal form of a braid.
    
    INPUT:
    
    - ``braid`` -- a braid
    
    OUTPUT:
    
    A list of lists with the right normal form. The first list contains the
    power of delta. The subsequent lists are the elementary permutation braids.
    """
    nstrands = braid.parent().strands()
    l1 = braid.Tietze()
    sig_on()
    cdef list[list[int]] rop = RightNormalForm(nstrands, l1)
    sig_off()
    return rop
    
def greatestcommondivisor(braid1, braid2):
    r"""
    Return the greatest common divisor of two braids.
    """
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[list[int]] rop = GreatestCommonDivisor(nstrands, l1, l2)
    sig_off()
    return rop
    
def leastcommonmultiple(braid1, braid2):
    r"""
    Return the least common multiple of two braids.
    """
    nstrands = max(braid1.parent().strands(), braid2.parent().strands())
    l1 = braid1.Tietze()
    l2 = braid2.Tietze()
    sig_on()
    cdef list[list[int]] rop = LeastCommonMultiple(nstrands, l1, l2)
    sig_off()
    return rop
    
def centralizer(braid):
    r"""
    Return a list of generators of the centralizer of a braid.
    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = CentralizerGenerators(nstrands, l)
    sig_off()
    return rop
    
def supersummitset(braid):
    r"""
    Return a list with the super-summit-set of a braid.
    """
    nstrands = braid.parent().strands()
    l = braid.Tietze()
    sig_on()
    cdef list[list[list[int]]] rop = SuperSummitSet(nstrands, l)
    sig_off()
    return rop