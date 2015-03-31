#!/usr/bin/env pythons

"""
rnatk.stats contains several python functions useful for \
statistical analysis of RNA secondary structure
"""

__all__ = ['randomXijDict', 'numRandomList', 'kernelXij', 'searchXijc',
           'fig_plot','Mij', 'equLen', 'equLen_file', 'similarity']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2014, MyRNA Project"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "sivanc7@gmail.com"

import ast, sys, math, numpy as np, matplotlib.pyplot as plt, pylab
from scipy import stats
from math import log

def equLen(str1, str2):
    """
    Check for equal length between the string1 and string2.
    If string1 is different to string 2, the system exits. Otherwise
    the system continues.
    Modules required:
    - sys
    Usage: <string1> <string2>
    """
    if len(str1) == len(str2):
        pass
    else:
        print 'At least one sequence is different in length to others. Press enter to exit.'
        raw_input()
        sys.exit()

def equLen_file(inFile):
    """
    Evaluate whether one sequence from the file is different to the others. The function requieres
    a list with elements Bio.SeqRecord.SeqRecord.
    Module required:
    - equLen
    Usage: <inFile>
    """
    for item_A in range(len(inFile)):
        for item_B in range(len(inFile)):
            if item_B-item_A >= 1:
                equLen(str(inFile[item_A].seq), str(inFile[item_B].seq))

def similarity(str1, str2):
    """
    Returns the similarity index between two string. Both string must have equal lenght.
    Module required:
    - equLen (from rnatk.stats)
    Usage: <string 1> <string 2>
    """
    equLen(str1, str2)
    matches = 0.0
    for index in range(len(str1)):
        if str1[index] == str2[index]:
            matches += 1
    similarity = matches/len(str1)
    return similarity

def Mij(pos1, pos2, logBase, mode='basepair'):
    """
    Calculate the mutual information between two RNA positions
    in an alignment. Each position must be converted into a string
    of RNA. The mode defines if the function must compute only those
    sequences obsered in the basepairs (mode=basepairs) or if the function
    must compute the MI no matters the basepairs (mode=free).
    Modules required:
    - sys
    - math
    - equLen
    Usage: <sequence1> <sequence2> <logBase> <mode (default=basepair)>
    """
    equLen(pos1,pos2)
    MI = 0
    if mode == 'basepair':
        allow = ['CG','GC','AU','UA','GU','UG']
    if mode == 'free':
        allow = ['AA', 'AC', 'AG', 'AU',
                 'CA', 'CC', 'CG', 'CU',
                 'GA', 'GC', 'GG', 'GU',
                 'UA', 'UC', 'UG', 'UU']
    for a in allow:
        Fxy, Fx, Fy = 0, 0, 0
        for b in range(len(pos1)):
            if pos1[b]+pos2[b] == a:
                Fxy += 1
        if Fxy > 0:
            Fx = pos1.count(a[0])/float(len(pos2))
            Fy = pos2.count(a[1])/float(len(pos2))
            if (Fxy/float(len(pos1)))*math.log(((Fxy/float(len(pos1)))/(Fx*Fy)),int(logBase)) <= 0:
                MI += 0
            else:
                MI += (Fxy/float(len(pos1)))*math.log(((Fxy/float(len(pos1)))/(Fx*Fy)),int(logBase))
        else:
            MI += 0
    return MI

def numRandomList(times, lower, upper):
    """
    This fuction creates a list with a determinate amount
    of random float numbers that varying from a lower bound to
    an upper bound.
    Modules required:
    - numpy (as np)
    Usage: <times> <lower> <upper>
    """
    list = []
    for a in range(int(times)):
        list.append(np.random.uniform(float(lower), float(upper)))
    return list

def randomXijDict(lenSeq, upper, bioSense=True):
    """
    This function builds a dictionary of (lenSeq**2) keys. Each key
    is associated with a random number that varies from 0 to an upper
    bound.
    Modules required:
    - numpy (as np)
    Usage: <lenSeq> <upper> 
    """
    dict = {}
    for a in range(int(lenSeq)):
        for b in range(int(lenSeq)):
            if bioSense:
                if b-a>3:
                    dict[a,b] = np.random.uniform(0, float(upper))
            else:
                if b-a>=1:
                    dict[a,b] = np.random.uniform(0, float(upper))
    return dict

def kernelXij(XijList):
    """
    Creates a Kernel Density Function from a list of values.
    Modules required:
    - stats (from scipy)
    - numpy (as np)
    Usage: <list>
    """
    KDE = stats.gaussian_kde(np.array(XijList, dtype=np.float))
    return KDE

def searchXijc(KDE, critical, lower, upper):
    """
    Given a KDE, this function search the integral from <lower> to
    <upper>, where <lower> is fixed, that results in <critical>.
    Module required:
    - numpy (as np)
    - stats (from scipy)
    Usage: <KDE> <critical value> <lower limit> <upper limit>
    """
    stop = True
    search = lower
    while stop:
        if KDE.integrate_box_1d(search+0.00001,upper) > critical:
            search += 0.0001
        else: stop = False
    return search, KDE.integrate_box_1d(search,upper)

def fig_plot(array, lower, upper, fileName, xlabel=''):
    """
    Save a png file for a KDE.
    Module required:
    - numpy (as np)
    - stats (from scipy)
    - matplotlib.plot (as plt)
    - pylab
    Usage: <array> <lower limit> <upper limit> <outFile> <type of covariance (string)>
    """
    x1 = np.array(array, dtype=np.float)
    kde1 = stats.gaussian_kde(x1)
    kde2 = stats.gaussian_kde(x1, bw_method = 'silverman')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x1, np.zeros(x1.shape), 'b+', ms=20)
    x_eval = np.linspace(lower, upper, num=200)
    ax.plot(x_eval, kde1(x_eval), 'k-')
    plt.ylabel('Density')
    plt.xlabel(xlabel)
    pylab.savefig(fileName+'.png')
    plt.clf()
