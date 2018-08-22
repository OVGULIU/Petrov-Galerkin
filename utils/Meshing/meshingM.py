# This module provides functions that help with the creation of nutils topologies/domains.
# It may later be extended to real meshing.

"""
Author:           Sasa Lukic
Date of creation: 13 Aug 2018
Python version:   3.6.6
Last update:      14 Aug 2018
"""


# TODO: This function is for CG-FE only. WIll need similar function for DG, IGA,
# or even one function working for all of them.
def feSupport1DRectilinear(index: int, degree: int, n_elems: int, lower: float, span: float) -> (float, float, int, int):
    """
    Assign element and local node number (of topology) to a globally numbered node.
    Then return arguments necessary for nutils.rectilinear to create domain.

    Args
    ----
    index:   Type: int.    Global node index
    degree:  Type: int.    Polynomial degree
    n_elems: Type: int.    Total number of elements
    lower:   Type: float:  Lower boundary of 1D domain
    span:    Type: float.  Length of a 1D element    (This formulation might help with refinement later on)

    Meta:
    ----
    Author:           Sasa Lukic
    Date of creation: 13 Aug 2018
    Python version:   3.6.6
    Last update:      15 Aug 2018

    """

    elem_nr      = index // degree     # Global element number (as long as not parallelized)
    elem_node_nr = index % degree      # Element local node number

    # Check if corner node
    if elem_node_nr == 0 :
        # Check if leftmost or rightmost node
        if elem_nr == 0:
            elems = 1
            leftBord = lower
            rightBord = lower + span
            localNodeNr = elem_node_nr
        elif elem_nr == n_elems:
            elems = 1
            leftBord = lower + span * (elem_nr - 1)
            rightBord = leftBord + span
            localNodeNr = degree       # Rightmost node
        else:
            # Nodes of interest are elem_nr and elem_nr - 1
            elems = 2
            leftBord = lower + span * (elem_nr - 1)
            rightBord = leftBord + 2 * span
            localNodeNr = degree      # 0-th node of current element, but degree_th node of topo (includes previous element)
    else:
        # Middle node
        elems = 1
        leftBord = lower + span * elem_nr
        rightBord = leftBord + span
        localNodeNr = elem_node_nr

    return leftBord, rightBord, elems, localNodeNr


# TODO: This function is for FE only. WIll need similar function for IGA,
# or even one function working for both.
def feNrNodes1DRectilinear(n_elems: int,  degree: int, type: str) -> int:
    """
    Return total number of nodes.

    Args
    ----
    n_elems: Type: int.    Total number of elements
    degree:  Type: int.    Polynomial degree
    type:    Type: str.    CG or DG?

    Meta:
    ----
    Author:           Sasa Lukic
    Date of creation: 14 Aug 2018
    Python version:   3.6.6
    Last update:      14 Aug 2018

    """

    if type == 'CG':
        nr_nodes = 1 + n_elems * degree
    elif type == 'DG':
        nr_nodes = n_elems * (degree + 1)
    else:
        print('The requested discretization type is not yet supported. Alternatively, check for typos, please.')
        nr_nodes = -1
    return nr_nodes
