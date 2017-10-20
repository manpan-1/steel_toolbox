# code to be evaluated for max field output value
# https://gist.github.com/crmccreary/1074551

"""
Module containing tools for Abaqus cae
"""

import os
import odbAccess
from abaqusConstants import *



def GetBlockPosition(model, blockPrefix):
    """
    Find a string and return the block number.

    Method to find a given string on the keywords file of a model and return an integer for the position of the first
    occurrence.

    Parameters
    ----------
    model : class
        Abaqus model to search for keyword
    blockPrefix : string
        String to look for

    Notes
    -----
    Includes code from stackoverflow
    https://stackoverflow.com/questions/35572838/how-to-request-nodal-stress-output-in-abaqus-python-script
    Author: agentp - https://stackoverflow.com/users/1004168/agentp


    References
    ----------

    """
    import string
    pos = 0
    for block in model.keywordBlock.sieBlocks:
        if string.lower(block[0:len(blockPrefix)]) == string.lower(blockPrefix):
            return pos
        pos = pos + 1
    return -1


# Look for the max value in a history output
# TO BE FIXED. THE REFERENCE POINTS (rp1key, ho1key etc.) ARE NOT GENERIC.
# Fetch maximum load, displacement and LPF for a riks analysis.
# The method assumes that a) the the odb is located in the current directory
def history_max(odb_name, step_name):
    """
    Look for the max value in a history output.

    Scan a history output on a step on an abaqus result database and return the maximum value.
    Currently, it returns LPF, load and disp history outputs. To be generalised.

    Parameters
    ----------
    odb_name : class
        Abaqus model containing the history results
    step_name : string
        Name of the step

    Attributes
    ----------

    Notes
    -----
    Incomplete: Work only for specific history output names, needs to be generalised.

    References
    ----------

    """

    myOdb = odbAccess.openOdb(path=odb_name + '.odb')
    RIKSstep = myOdb.steps[step_name]
    rp1key = RIKSstep.historyRegions.keys()[1]
    ho1key = RIKSstep.historyRegions[rp1key].historyOutputs.keys()[0]
    rp2key = RIKSstep.historyRegions.keys()[2]
    ho2key = RIKSstep.historyRegions[rp2key].historyOutputs.keys()[0]
    asskey = RIKSstep.historyRegions.keys()[0]
    hoasse = RIKSstep.historyRegions[asskey].historyOutputs.keys()[-1]
    load_hist = RIKSstep.historyRegions[rp1key].historyOutputs[ho1key].data
    disp_hist = RIKSstep.historyRegions[rp2key].historyOutputs[ho2key].data
    lpf_hist = RIKSstep.historyRegions[asskey].historyOutputs[hoasse].data
    maxpos = load_hist.index(max(load_hist, key=lambda x: x[1]))
    load = load_hist[maxpos][1]
    disp = -disp_hist[maxpos][1]
    lpf = lpf_hist[maxpos][1]
    odbAccess.closeOdb(myOdb)
    return lpf, load, disp


def fetch_eigenv(odb_name, step_name, n_eigen):
    """
    Get eigenvalues.

    Return the eigenvalues of a perturbation buckling analysis from an abaqus database.

    Parameters
    ----------
    odb_name : class
        Abaqus model containing the eigenvalues
    step_name : string
        Name of the step
    n_eigen : int
        Number of eigenvalues to return

    Attributes
    ----------

    Notes
    -----

    References
    ----------

    """

    bckl_odb = odbAccess.openOdb(path=odb_name + '.odb')
    bckl_step = bckl_odb.steps[step_name]

    # Gather the eigenvalues
    eigenvalues = ()
    eigen_string = ""
    for J_eigenvalues in range(1, n_eigen + 1):
        current_eigen = float(bckl_step.frames[J_eigenvalues].description[-11:])
        eigenvalues = eigenvalues + (current_eigen,)
        eigen_string = eigen_string + "%.3E " % (current_eigen)

    # Close the odb
    odbAccess.closeOdb(bckl_odb)

    # Return variables
    return eigenvalues, eigen_string
