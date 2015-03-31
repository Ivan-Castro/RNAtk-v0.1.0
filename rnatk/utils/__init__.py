#!/usr/bin/env python
# coding=utf-8

"""
rnatk.utils contains several python functions \
useful for all others rnatk scripts
"""

__all__ = ['locate', 'retrieve_varName', 'get_info_byGB', 'tableByGB', 'fillTableByGB',
           'sortGB_table', 'makeGB_LaTeX_table', 'filterXijc', 'makeFE_table', 'similarity',
           'bestPartner', 'dictXij', 'signXij', 'hintonDiag', 'makeFE_LaTeX_table', 'graph_cova']
__author__ = "Ivan Castro"
__copyright__ = "Copyright 2014, MyRNA Project"
__license__ = "GPL"
__version__ = "0.1.0"
__email__ = "sivanc7@gmail.com"

import rnatk, re, subprocess, inspect, numpy as np, matplotlib.pyplot as plt, networkx as nx
from pylab import axes
from Bio import Entrez, SeqIO
from rnatk.stats import (Mij, equLen)
from rnatk.structure2D import (Cij, getFree_energy)
from rnatk.fasta import openFasta, GBcode

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

def bestPartner(varDict):
    """
    Filter the pair positions and search for the best partner with higgest value.
    Caution: use this function only if the structure is supposed to be pseudoknot-free.
    Script usage: bestPartner <varDict>
    """
    removeList = []
    for a in varDict:
        for b in varDict:
            if b>a:
                if (a[0] in b) or (a[1] in b):
                    if varDict[a] < varDict[b]:
                        removeList.append(a)
                    if varDict[a] > varDict[b]:
                        removeList.append(b)
    for c in removeList:
        if c in varDict:
            varDict.pop(c)
    return varDict

def dictXij(matrixXij, covType):
    """
    This function returns a dictionary that contains all pairs 
    (a,b) posible (with exception b-a>3) from the matrixXij and its 
    associated covariance values.
    Modules required:
    - Mij
    - Cij
    """
    Xij = {}
    for a in range(len(matrixXij)):
        for b in range(len(matrixXij)):
            if b-a>3:
                if covType == 'Mij':
                    Xij[a+1,b+1] = Mij(matrixXij[a],matrixXij[b],2)
                if covType == 'Cij':
                    Xij[a+1,b+1] = Cij(matrixXij[a],matrixXij[b])
    return Xij

def filterXijc(dictXij, critical, criteria='higher'):
    """
    Filter the values of covariance according to a critical value.
    Script usage: filterXijc <dictXij> <critical>
    """
    dictRes = {}
    for a in dictXij.items():
        if criteria == 'higher':
            if a[1] >= float(critical):
                dictRes[a[0][0],a[0][1]]=a[1]
        if criteria == 'lower':
            if a[1] >= float(critical):
                dictRes[a[0][0],a[0][1]]=a[1]
    return dictRes

def sizeAdjust(x):
    """
    Adjusts the size of the square for the hintonDiag function.
    """
    adjust = (79151.44761*(x**-4.471344938))*((x-1)**2.353788038)
    return adjust

def hintonDiag(dictVar, outFile, title, maxValue, color='black'):
    """
    This function save a Hinton Diagram into the specified outFile argument. The input 
    data must be a dictionary where the keys are 2-tuple associated with one value.
    Modules required:
    - numpy (as np)
    - matplotlib.pyplot (as plt)
    - axes (from pylab)
    Usage: <dictionary> <outFile> <title> <max value> <color (default=black)>
    """
    x, y, z = [], [] ,[]
    for elem in dictVar:
        x.append(elem[0])
        y.append(elem[1])
        z.append(dictVar[elem[0],elem[1]])
    N_data = max([max(x), max(y)])
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    plt.title(title)
    plt.xlabel('Positions')
    plt.ylabel('Positions')
    plt.axis([-1, x.max()+1, -1, y.max()+1])
    axes().set_aspect('equal')
    plt.grid(True)
    plt.scatter(x, y, s=z*(sizeAdjust(N_data)/maxValue), c='black', marker='s')
    plt.savefig(outFile)
    plt.clf()

def locate(term):
    """
    Locate a determinate file (path+filename) and return the first result.
    Module required:
    - subprocess
    Script usage: locate <term>
    """
    result = (subprocess.Popen(['locate', str(term)], stdout=subprocess.PIPE).stdout.readlines())[0].strip('\n')
    return result

def retrieve_varName(var):
    """
    Returns the name of a variable in str. Credits for Josh (stackoverflow
    userID = 2602718).
    Module required:
    - inspect
    """
    callers_local_vars = inspect.currentframe().f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

def get_info_byGB(code):
    """
    WARNING: before to use this function identify yourself to NCBI using Entrez.email.
    Return the info associated with the code in the GenBank database. This fucntion
    only return the order taxonomy if it ends in 'formes'. Also, the function only fetch
    the family info if it ends with 'idae' or 'aceae'.
    Module required:
    - Entrez (from Bio)
    - SeqIO (from Bio)
    """
    handle = Entrez.efetch(db="nucleotide", id=code, rettype="gb", retmode="text")
    description = SeqIO.read(handle, "genbank")
    info = {}
    info["GB_accession"] = description.id
    for elem in description.annotations['taxonomy']:
        if elem.endswith("formes"): info["order"] = elem
        if elem.endswith("idae") or elem.endswith('aceae'): info["family"] = elem
    info["id"] = description.annotations["gi"]
    info["Scientific_Name"] = description.annotations["organism"]
    return info
    
def tableByGB(inFile):
    """
    WARNING: before to use this function identify yourself to NCBI using Entrez.email.
    This function build a table based on a list. Each element is the info for each sequence.
    Module required:
    - Entrez (from Bio)
    - SeqIO (from Bio)
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    inFile = GBcode(inFile)
    GB_table = []
    for elem in inFile:
        GB_table.append(get_info_byGB(elem.description))
    return GB_table

def fillTableByGB(GB_table):
    """
    Fill some empty slots of info in the table returned by the function tableByGB.
    Module required:
    - Entrez (from Bio)
    Script usage: fillTableByGB <GB table>
    """
    for elem in GB_table:
        if (not elem.has_key("order")) or (not elem.has_key("family")):
            handle = Entrez.esearch(db="taxonomy", term=elem["Scientific_Name"])
            basicInfo = Entrez.read(handle)
            IdList = basicInfo["IdList"][0]
            handle = Entrez.efetch(db="Taxonomy", id=IdList, retmode="xml")
            record = Entrez.read(handle)
            record[0]["Lineage"] = record[0]["Lineage"].split(";")
            for taxon in record[0]["Lineage"]:
                if (not elem.has_key("order")) and taxon.endswith("formes"):
                    elem["order"] = taxon.strip()
                if (not elem.has_key("family")) and (taxon.endswith("idae") or taxon.endswith("aceae")):
                    elem["family"] = taxon.strip()
    return GB_table

def sortGB_table(GB_table, taxa='order'):
    """
    This function sort the table alphabetically by order or family 
    Script usage: sortGB_table <GB table> <sort by (dafaul='order')>
    """
    GB_table = sorted(GB_table, key=lambda k: k[taxa]) 
    return GB_table
        
def makeGB_LaTeX_table(GB_table, page, tableTitle, taxa, titlePage=False, section=False, comments=False):
    """
    This function returns a string objects which could by paste to LaTeX to build a table of taxonomy species.
    Usage: <table> <page number> < table title> <page title (default=False)> <section (default = False)> <comments (default=False)>
    """
    GB_table = sortGB_table(GB_table, taxa)
    LaTeX = "\\documentclass[12pt]{article}\n\\usepackage{caption}\n\\usepackage{longtable}\n\\usepackage{geometry}\n\
    \\geometry{legalpaper, total={21.59cm,27.947cm}, left=2cm, right=2cm, top=3cm, bottom=3cm}\n"
    if titlePage:
        LaTeX = LaTeX + "\\title{"+titlePage+"}\n"
    LaTeX = LaTeX + "\\begin{document}\n\\setcounter{page}{"+page+"}\n"
    if titlePage:
        LaTeX = LaTeX + "\\maketitle\n"
    if section:
        LaTeX = LaTeX + "\\section*{"+section+"}\n"
    if comments:
        LaTeX = LaTeX + comments + "\n"
    LaTeX = LaTeX + "\n\\centering\n\\begin{longtable}{ccccc}\n\\caption*{"+tableTitle+"}\\\\\n\
    GenBank code & GenBank ID & Order & Family & Scientific name \\\\ \\hline\n"
    for elem in GB_table:
        try:
            elem['order'] = elem['order']
        except:
            elem['order'] = 'No identified'
        LaTeX = LaTeX + elem["GB_accession"] + "&" + elem["id"] + "&" + elem["order"] + "&" + elem["family"] + "&" + elem["Scientific_Name"] + "\\\\\n"
    LaTeX = LaTeX + "\\end{longtable}\n\\end{document}"
    return LaTeX

def makeFE_table(inFile, constCov, constrain):
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = GBcode(openFasta(inFile))
    table = []
    for item in inFile:
        table.append([item.description, getFree_energy(str(item.seq)), getFree_energy(str(item.seq), constCov), getFree_energy(str(item.seq), constrain)])
    return table

def makeFE_LaTeX_table(FE_table, page, tableTitle, titlePage=False, section=False, comments=False):
    """
    This function returns a string objects which could by paste to LaTeX to build a table of free energy.
    Usage: <table> <page number> < table title> <page title (default=False)> <section (default = False)> <comments (default=False)>
    """
    LaTeX = "\\documentclass[12pt]{article}\n\\usepackage{caption}\n\\usepackage{longtable}\n\\usepackage{geometry}\n\
    \\geometry{legalpaper, total={21.59cm,27.947cm}, left=2cm, right=2cm, top=3cm, bottom=3cm}\n"
    if titlePage:
        LaTeX = LaTeX + "\\title{"+titlePage+"}\n"
    LaTeX = LaTeX + "\\begin{document}\n\\setcounter{page}{"+page+"}\n"
    if titlePage:
        LaTeX = LaTeX + "\\maketitle\n"
    if section:
        LaTeX = LaTeX + "\\section*{"+section+"}\n"
    if comments:
        LaTeX = LaTeX + comments + "\n"
    LaTeX = LaTeX + "\n\\centering\n\\begin{longtable}{cccc}\n\\caption*{"+tableTitle+"}\\\\\n\
    GenBank code & No constrained & constrained by covariance & Full constrained \\\\ \\hline\n"
    for elem in FE_table:
        LaTeX = LaTeX + elem[0] + "&" + elem[1] + "&" + "&" + elem[2] + "&" + elem[3] + "\\\\\n"
    LaTeX = LaTeX + "\\end{longtable}\n\\end{document}"
    return LaTeX

def signXij(dictCov, criticalValue):
    """
    Given a dictionary of covariance values for each base-pair, this
    function returns only those values greather than a specific value.
    Usage: signXij <dictionary> <critical value>
    """
    diction = {}
    for a in dictCov:
        if dictCov[a] > float(criticalValue):
            diction[a[0],a[1]] = dictCov[a]
    return diction

def graph_cova(dictResults, title, outFile, edge_labels=True):
    """
    Buil a graph that sumarize the covariation relationship of each positions
    Module required:
    - networkx (as nx)
    - matplotlib.pyplot (as plt)
    Usage: <dictResults> <title> <outFile> <edge_labels (default=True)>
    """
    for elem in dictResults:
        dictResults[elem] = round(dictResults[elem],2)
    G=nx.Graph()
    for elem in dictResults:
        G.add_edge(elem[0],elem[1], weight=1)
    pos = nx.circular_layout(G, scale=3)
    nx.draw_networkx_nodes(G,pos,node_size=400, node_color="white")
    nx.draw_networkx_edges(G,pos, width=1,alpha=1,edge_color='black')
    nx.draw_networkx_labels(G,pos,font_size=11,font_family='courier-new')
    if edge_labels:
        nx.draw_networkx_edge_labels(G, pos, dictResults, label_pos=0.5, font_size=8)
    plt.title(title)
    plt.axis('off')
    plt.figsize=(12,12)
    plt.savefig(outFile+'.png')
    plt.clf()

def _makeNGS_LateX_table(inFile, page, tableTitle, titlePage=False, section=False, comments=False):
    """
    This code is for personal use.
    """
    if (type(inFile) == str) and (str(type(inFile[0])) != "<class 'Bio.SeqRecord.SeqRecord'>"):
        inFile = openFasta(inFile)
    LaTeX = "\\documentclass[12pt]{article}\n\\usepackage{caption}\n\\usepackage{longtable}\n\\usepackage{geometry}\n\
    \\geometry{legalpaper, total={21.59cm,27.947cm}, left=2cm, right=2cm, top=3cm, bottom=3cm}\n"
    if titlePage:
        LaTeX = LaTeX + "\\title{"+titlePage+"}\n"
    LaTeX = LaTeX + "\\begin{document}\n\\setcounter{page}{"+page+"}\n"
    if titlePage:
        LaTeX = LaTeX + "\\maketitle\n"
    if section:
        LaTeX = LaTeX + "\\section*{"+section+"}\n"
    if comments:
        LaTeX = LaTeX + comments + "\n"
    LaTeX = LaTeX + "\n\\centering\n\\begin{longtable}{cccccc}\n\\caption*{"+tableTitle+"}\\\\\n\
    SRA code & conting ID & Order & Family & Scientific name & Strand direction \\\\ \\hline\n"
    for elem in inFile:
        conting_ID = re.findall(r'[A-Z]{1,3}[0-9]{3,9}.[0-9]{1,9}.[0-9]{1}?',elem.description)[0]
        orientation = re.findall('inverted', elem.description)
        if orientation:
            orientation = 'Positive'
        else:
            orientation = 'Negative'
        if conting_ID.startswith('SRR6847'):
            LaTeX = LaTeX + 'SRX228421' + "&" + conting_ID + "&" + 'Lamniformes' + "&" + 'Lamnidae' + "&" + 'Carcharodon carcharias' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR6845'):
            LaTeX = LaTeX + 'SRX228332' + "&" + conting_ID + "&" + 'Lamniformes' + "&" + 'Lamnidae' + "&" + 'Carcharodon carcharias' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('ERR375875'):
            LaTeX = LaTeX + 'ERX348252' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Hemiscylliidae' + "&" + 'Chiloscyllium griseum' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('ERR375876'):
            LaTeX = LaTeX + 'ERX348253' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Hemiscylliidae' + "&" + 'Chiloscyllium griseum' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR652972'):
            LaTeX = LaTeX + 'SRX219865' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Ginglymostomatidae' + "&" + 'Ginglymostoma cirratum' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR652971'):
            LaTeX = LaTeX + 'SRX219866' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Ginglymostomatidae' + "&" + 'Ginglymostoma cirratum' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521184') or conting_ID.startswith('SRR1521183'):
            LaTeX = LaTeX + 'SRX653790' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521190'):
            LaTeX = LaTeX + 'SRX657784' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521191'):
            LaTeX = LaTeX + 'SRX657785' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521192'):
            LaTeX = LaTeX + 'SRX657786' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521195'):
            LaTeX = LaTeX + 'SRX657787' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521197'):
            LaTeX = LaTeX + 'SRX657788' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521198'):
            LaTeX = LaTeX + 'SRX657789' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521199'):
            LaTeX = LaTeX + 'SRX657790' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521200'):
            LaTeX = LaTeX + 'SRX657791' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521201'):
            LaTeX = LaTeX + 'SRX657792' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1521204'):
            LaTeX = LaTeX + 'SRX657793' + "&" + conting_ID + "&" + 'Orectolobiformes' + "&" + 'Rhincodontidae' + "&" + 'Rhincodon typus' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1514129'):
            LaTeX = LaTeX + 'SRX651773' + "&" + conting_ID + "&" + 'Carcharhiniformes' + "&" + 'Scyliorhinidae' + "&" + 'Scyliorhinus canicula' + "&" + orientation + "\\\\\n"
        if conting_ID.startswith('SRR1514131'):
            LaTeX = LaTeX + 'SRX651775' + "&" + conting_ID + "&" + 'Carcharhiniformes' + "&" + 'Scyliorhinidae' + "&" + 'Scyliorhinus canicula' + "&" + orientation + "\\\\\n"
    LaTeX = LaTeX + "\\end{longtable}\n\\end{document}"
    return LaTeX
