# code to be evaluated for max field output value
# https://gist.github.com/crmccreary/1074551

"""
Module containing tools for Abaqus cae
"""

import os
import odbAccess
from abaqusConstants import *



# definition of a method to search for a keyword position
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

    Attributes
    ----------

    Notes
    -----

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


def open_odb(odbPath):
    """
    A more sophisticated open odb function.

    Parameters
    ----------
    odbPath : string
        Path and filename of the database (without the '.odb' extension)

    Attributes
    ----------

    Notes
    -----

    References
    ----------

    """

    base, ext = os.path.splitext(odbPath)
    odbPath = base + '.odb'
    new_odbPath = None
    if odbAccess.isUpgradeRequiredForOdb(upgradeRequiredOdbPath=odbPath):
        print('odb %s needs upgrading' % (odbPath,))
        path, file_name = os.path.split(odbPath)
        file_name = base + "_upgraded.odb"
        new_odbPath = os.path.join(path, file_name)
        odbAccess.upgradeOdb(existingOdbPath=odbPath, upgradedOdbPath=new_odbPath)
        odbPath = new_odbPath
    odb = odbAccess.openOdb(path=odbPath, readOnly=True)
    return odb



def field_max(odb, result):
    """
    Look for the max value in a field output.

    Scan a field output on an abaqus result database and return the maximum value

    Parameters
    ----------
    odb : class
        Abaqus model containing the field results
    result : class
        Field output result to search for max

    Attributes
    ----------

    Notes
    -----

    References
    ----------

    """

    result_field, result_invariant = result
    _max = -1.0e20
    for step in odb.steps.values():
        print('Processing Step:', step.name)
        for frame in step.frames:
            if frame.frameValue > 0.0:
                allFields = frame.fieldOutputs
                if (allFields.has_key(result_field)):
                    stressSet = allFields[result_field]
                    for stressValue in stressSet.values:
                        if result_invariant:
                            if hasattr(stressValue, result_invariant.lower()):
                                val = getattr(stressValue, result_invariant.lower())
                            else:
                                raise ValueError('Field value does not have invariant %s' % (result_invariant,))
                        else:
                            val = stressValue.data
                        if (val > _max):
                            _max = val
                else:
                    raise ValueError('Field output does not have field %s' % (results_field,))
    return _max


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
