#!/usr/bin/python
# coding: utf-8
import time
import table
import random
import logging, math, bisect
import matplotlib
import operator
import numpy as np
import objects
from scipy.stats.stats import pearsonr, spearmanr

#from numpy import random as nprand
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import xlwt
import re
from numpy import random as nprand
#from pylab import *
import src.constants as constants

MANY_TO_MANY = 0
ONE_GENE_TO_MANY_miRNA = 1
ONE_miRNA_TO_MANY_GENES = 2
ONE_TO_ONE = 3

def network_type_str(network_type):
    if network_type == MANY_TO_MANY:
        return "many_to_many"
    elif network_type == ONE_GENE_TO_MANY_miRNA:
        return "one_gene_to_many_miRs"
    elif network_type == ONE_miRNA_TO_MANY_GENES:
        return "one_miR_to_many_genes"
    elif network_type == ONE_TO_ONE:
        return "one_to_one"


def timer(timeList, restart, sec):
#Print the time elapsed between the first and second time this function is called.
    timeList.append(time.time())
    minInterval = 0
    timeElapsed = 0
    out_line = ''
    if len(timeList)%2 == 0:
        timeElapsed = round(timeList[-1] - timeList[-2],4)
        if(minInterval > 0):
            if(round(timeList[-1] - timeList[-2],4) > minInterval):
                logERROR("This is taking too long:")
                logERROR('Time elapsed: ' + str(round(timeList[-1] - timeList[-2],4)) + ' seconds.')            
        if(sec):
            out_line = 'Time elapsed: ' + str(round(timeList[-1] - timeList[-2],4)) + ' seconds.'
        else:
            out_line = 'Time elapsed: ' + str(round(timeList[-1] - timeList[-2],4)/60) + ' minutes.'
            out_line += time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        timeList.pop()
        timeList.pop()
    if(restart):
        timeList.append(time.time())
    return [timeList, timeElapsed, out_line]

def createProbabilities(data, dist):
    probability = {}
    if(dist == 'normalize'):
        normby = sum(data.values())
        for key in data:
            probability[key] = float(data[key])/normby
    return probability

def pwd():
    #pwd = '/cs/=/learning/shellym/imac/miRNA/'
    pwd = 'C:/miR/'
    return pwd

def initLogFile(LOG_FILENAME):
    logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)
    
def logINFO(msg):
    logging.info(msg)
    
def logERROR(msg):
    logging.error(msg)

def parseTable(filename, keyIdx, valIdx, header):
    f_in = open(pwd() + filename, 'r')
    i = 0
    data = {}
    for line in f_in:
        if(not header or i > 0):
            row = line.split('\t')
            data[row[keyIdx].strip()] = float(row[valIdx].strip())
        i += 1
    return data

# sheets is a map of sheet name to tuple:
# [0] header
# [1] data table
def write2excel(filename, sheets, mark_lines=None):
    book = xlwt.Workbook(encoding="utf-8")
    marked_style = xlwt.easyxf('pattern: pattern solid, fore_colour red;')
    for sheetname in sorted(sheets.keys()):
        sheet = book.add_sheet(sheetname)
        header = sheets[sheetname][0]
        data = sheets[sheetname][1]
        for i, title in enumerate(header):
            sheet.write(0, i, title)
        for i, data_tuple in enumerate(data):
            for j, curr_data in enumerate(data_tuple):
                #if(type(curr_data) is not str):
                #    print(curr_data)
                if(mark_lines and mark_lines in curr_data):
                    sheet.write(i+1,j, str(curr_data), marked_style)
                else:
                    sheet.write(i+1,j, str(curr_data))
    book.save(pwd() + filename)


def writeList(lst, out_file):
    f_out = open(pwd() + out_file, 'w')
    for item in sorted(lst):
        f_out.write(str(item) + '\n')
    f_out.close()

def writeMap(key2val, out_file):
    f_out = open(pwd() + out_file, 'w')
    for key in sorted(key2val):
        f_out.write(str(key) + '\t' + str(key2val[key]) + '\n')
    f_out.close()
    
def mergeIntDicts(dict1, dict2):
    for item in dict2:
        if(item in dict1):
            dict1[item] += dict2[item]
        else:
            dict1[item] = dict2[item]
    return dict1

def permute(lst):
    #permuteLst = list(itertools.permutations(lst))
    #idx = nprand.choice(range(len(lst)))
    #permuteLst = list(permuteLst[idx])
    random.shuffle(lst)
    #return permuteLst
    
def corr(list1, list2):
    return pearsonr(list1, list2)

def corr_spearman(list1, list2):
    return spearmanr(list1, list2)

def corr_columns(filename, ref):
    f_in = open(pwd() + filename, 'r')
    i = 0
    idx2key = {}
    key2col = {}
    for line in f_in:
        row = line.split('\t')
        for j in range(len(row)):
            if(j > 0 ):
                if(i == 0):
                    idx2key[j] = row[j].strip()
                elif(i == 1):
                    key2col[int(idx2key[j])] = [float(row[j].strip())]
                else:
                    key2col[int(idx2key[j])].append(float(row[j].strip()))
            j += 1
        i += 1
    f_in.close()
    corrMap = {}
    for key in key2col:
        corrMap[key] = corr(key2col[ref], key2col[key])
    return corrMap

def findIntersection(files, outfile):
    s = set()
    i = 0
    for f in files:
        lc = 0
        set1 = set()
        f_in = open(pwd() + 'output/HeLa/GSM1299957/' + f, 'r')
        for line in f_in:
            set1.add(line.strip())
            lc += 1
        f_in.close()
        print(f + " num of genes: " + str(lc))
        if(i == 0):
            s = set1
        else:
            s = s.intersection(set1)
        i += 1        
    f_out = open(pwd() + outfile, 'w')    
    for item in s:
        f_out.write(item + "\n")
    f_out.close()

def printException(dataMap):
    for ds in dataMap:
        logERROR(str(ds) + " : " + str(dataMap[ds]))


def load_data(cell_type, miR_exp_file, gene_exp_file, score_table, source="microrna.org",
              key="miR", name_format=True, network_type=MANY_TO_MANY, half_life_file=None, cell_type_gene=None):
    miR2exp = {}
    mRNA2exp = {}
    mRNAdist = {}
    clash_in_table = None
    miRmRNAexp_clash = None
    table, table_genes = scoreTable(score_table, key, 0, source, network_type)
    if cell_type_gene is None:
        cell_type_gene = cell_type
    if cell_type == 'HEK':
        # HEK293 cells data:
        # miR2exp = miRexpression('expression_data/GSM952037_Untreated_cells_annotated.txt')
        #miR2exp = miRexpression('expression_data/HEK/' + miR_exp_file)
        #[mRNA2exp, mRNAdist] = mRNAexpressionDetailed('mRNA_exp_hg19.tab')
        miR2exp = miRexpression('expression_data/HEK/' + miR_exp_file)
    if cell_type_gene == 'HEK':
        mRNAdist = mRNAexpression('HEK/' + gene_exp_file)
        # miRmRNAexp_clash, mRNA2chimera, clash_in_table = processCLASHdata('clashdataS1.txt', table=table, mir_exp=miR2exp, mRNA_exp=mRNAdist) # [0]
        #[table, table_genes] = scoreTable(score_table, 'miR', 1)
    if cell_type == 'RA':
        # right atrial data:
        miR2exp = miRexpression('expression_data/GSM792462_counts_s1_RA_miRNA.txt')
        mRNAdist = mRNAexpression('GSM792454_counts_s1_RA_mRNA.txt')
        #[table, table_genes] = scoreTable(score_table, 'miR', 0)
    if cell_type == 'HeLa':
        #miR2exp = miRexpression('expression_data/HeLa/miRnamexpGSM1299957.tab')
        miR2exp = miRexpression('expression_data/HeLa/' + miR_exp_file)
    if cell_type_gene == "HeLa":
        mRNAdist = mRNAexpression(gene_exp_file)
    if cell_type == "MCF7":
        miR2exp = miRexpression("expression_data/MCF7/" + miR_exp_file)
    if cell_type_gene == "MCF7":
        mRNAdist = mRNAexpression('MCF7/' + gene_exp_file)
    else:
        miR2exp = miRexpression("expression_data/HeLa/" + miR_exp_file)
        mRNAdist = mRNAexpression(gene_exp_file, gene_file="C:\\miR\\expression_data\\miRTarBase\\genes__common_targetscan.txt")
    half_life = {}
    if half_life_file is not None:
        f_half_life = open(half_life_file)
        for i, line in enumerate(f_half_life):
            if i > 0:
                row = line.split("\t")
                gene = row[0].strip()
                val = row[1].strip()
                if val == ">24":
                    half_life[gene] = 30
                elif val != "N.D.":
                    half_life[gene] = float(val)
    # creating target - miRNA table
    target_mir_table = {}
    for mir in table:
        for gene in table[mir]:
            if gene not in target_mir_table:
                target_mir_table[gene] = {}
            target_mir_table[gene][mir] = table[mir][gene]
    return miR2exp, mRNA2exp, mRNAdist, table, target_mir_table, half_life


def correct_dist(dist):
    """
    Corrects the given distribution so that the sum of probabilities would be 1
    :param dist:
    :return: list of keys and their corresponding probabilities
    """
    keys = list(dist.keys())
    probabilities = list(dist.values())
    total = sum(dist.values())
    if 1-total != 0:
        diff = 1-total
        if diff > 0:
            rand_idx = nprand.choice(range(len(keys)))
            probabilities[rand_idx] += diff
        else:
            fixed = False
            rand_idx = nprand.choice(range(len(probabilities)))
            while not fixed and rand_idx < len(probabilities):
                if probabilities[rand_idx] + diff >= 0:
                    probabilities[rand_idx] += diff
                    fixed = True
                else:
                    rand_idx += 1
                    if rand_idx == len(probabilities):
                        rand_idx = 0
    return [keys, probabilities]


def score2prob(score, source, exp_fold_change=-1, discrete_probabilities=False):
    """
    Converting table scores to probabilities.
    :param score:
    :param exp_fold_change: This is relevant to microrna.org table
    :param source:
    :return:
    """
    if source.startswith('both'):  # This is for the case that we are using intersection table of both sources
        from_source = source.split('_')[1]
        if from_source == 'm':
            source = "microrna.org"
            score = score[0]
        elif from_source == 't':
            source = "TargetScan"
            score = score[1]
        else:
            return random.random()  # return random probability
    if source == "microrna.org":
        if exp_fold_change >= -0.1:
            p1 = -0.3619
            p2 = -1.96
            p3 = -3.763
            p4 = -3.09
            p5 = -1.279
            p6 = 0.237        
        elif -0.5 <= exp_fold_change < -0.1:
            p1 = -0.2635
            p2 = -1.456
            p3 = -2.839
            p4 = -2.323
            p5 = -0.9859
            p6 = 0.1104
        elif -1 <= exp_fold_change < -0.5:
            p1 = -0.2274
            p2 = -1.24
            p3 = -2.324
            p4 = -1.733
            p5 = -0.6418
            p6 = 0.04136
        elif -1.5 <= exp_fold_change < -1:
            p1 = -0.1632
            p2 = -0.8814
            p3 = -1.624
            p4 = -1.151
            p5 = -0.3869
            p6 = 0.0131
        prob =p1*pow(score,5)+p2*pow(score,4)+p3*pow(score,3)+p4*pow(score,2)+p5*score+p6
    elif source == "TargetScan" or source == constants.RANDOM_SCORES:
        if score <= -20:
            prob = 1
        else:
            prob = 1 - (2**score)
    elif source == constants.MIR_TAR_BASE:
        if score >= 10:
            prob = 1
        else:
            prob = 0.01*score
    if discrete_probabilities:
        if prob > 0.75:
            prob = 1
        elif 0.25 <= prob <= 0.75:
            prob = 0.5
        else:
            prob = 0.2
    return prob


def scoreTable(table_filename, key, detailed, source, network_type=MANY_TO_MANY):
    if key == 'mRNA':
        return table.targetmiRScore(table_filename, source)
    if key == 'miR':
        if network_type == MANY_TO_MANY:
            if source == constants.MIR_TAR_BASE:
                return table.miRtargetScore_mirtarbase(table_filename, detailed, source)
            return table.miRtargetScore(table_filename, detailed, source)
        elif network_type == ONE_GENE_TO_MANY_miRNA:
            return table.one_gene_to_many_miRNA_table(table_filename, detailed, source)
        elif network_type == ONE_miRNA_TO_MANY_GENES:
            return table.one_miRNA_to_many_genes_table(table_filename, detailed, source)
        elif network_type == ONE_TO_ONE:
            return table.one_to_one_table(table_filename, detailed, source)
        else:
            logERROR("Error in Methods.scoreTable \n network type is not defined, exiting...")
        exit(0)
    else:
        logERROR("Error in Methods.scoreTable \n key type is not mRNA or miR, exiting...")
        exit(0)

def plotTimeStats(timeMap, factor):
    colors = {1:'red', 2:'blue', 3:'green'}
    i = 1
    for item in timeMap:
        x = range(len(timeMap[item]))
        y = timeMap[item]
        for idx in x:
            y[idx] /= factor            
        plt.plot(x, y, label=item, color=colors[i])
        i += 1
    plt.legend(loc='lower right') 
    plt.savefig(pwd() + "/images/plotTimeStats.png")
    plt.cla()
    #plt.show()
    
def colors(idx, num_colors=15):
    if(idx == 0):
        return 'red'
    if(idx == 1):
        return 'blue'
    if(idx == 2):
        return 'green'
    if(idx == 3):
        return 'cyan'
    if(idx == 4):
        return 'magenta'
    if(idx == 5):
        return 'yellow'
    if(idx == 6):
        return 'gold'
    cm = get_cmap('rainbow')
    color = cm(1.*idx/num_colors)
    return color


def plotDict(miRs, mRNAs,  xticksLables, figname):
    x = sorted(miRs.keys())
    ymiRs = []
    ymRNAs = []
    ticks = [0 for j in range(len(x))]
    tc = 0
    prev_count = 0
    for i in range(len(x)):
        if(i == 1 and isinstance(xticksLables[1], str)):
            count = int(xticksLables[1].split('_')[1])
        elif(i > 0):
            prev_count = count
            count = int(xticksLables[i])            
        else:
            count = 0
        if(prev_count <= count):
            tc += count - prev_count
        else:
            tc += count 
        ticks[i] = tc
        ymiRs.append(miRs[x[i]])
        ymRNAs.append(mRNAs[x[i]])
    #[fig, ax] = plt.subplots()
    #for t in range(len(ticks)):
    #    print(str(ticks[t]) + '\t' + str(ymiRs[t]) + '\t' + str(ymRNAs[t]))
    plt.plot(ticks, ymiRs, label='miRs', color='red')
    plt.plot(ticks, ymRNAs, label='mRNAs', color='blue')
    tick1 = []
    ticklable1 = []
    for i in range(len(x)):
        if((i > 2 and i%10 == 0) or i <= 2):
            tick1.append(ticks[i])
            ticklable1.append(xticksLables[i])
    plt.xticks(tick1, ticklable1, rotation='vertical')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.title(figname)
    #ax.set_xticklabels(xticksLables)
    plt.savefig(pwd() + "/images/" + figname + ".png")
    plt.cla()
    
def plotDict1(dataLsts, xticksLables, figname):
    ticks = [0 for j in range(len(xticksLables))]
    tc = 0
    prev_count = 0
    for i in range(len(xticksLables)):
        if(i == 1 and isinstance(xticksLables[1], str)):
            count = int(xticksLables[1].split('_')[1])
        elif(i > 0):
            prev_count = count
            count = int(xticksLables[i])            
        else:
            count = 0
        if(prev_count <= count):
            tc += count - prev_count
        else:
            tc += count 
        ticks[i] = tc        
    i = 0
    for key in dataLsts:
        if(not(len(ticks) == len(dataLsts[key]))):
            logERROR("in plotDict1 len(ticks): " + str(len(ticks)) + " , len(dataLsts[key]): " + str(len(dataLsts[key])))
        else:
            plt.plot(ticks, dataLsts[key], label=key, color=colors(i))
        i += 1
    tick1 = []
    ticklable1 = []
    for i in range(len(xticksLables)):
        if((i > 2 and i%10 == 0) or i <= 2):
            tick1.append(ticks[i])
            ticklable1.append(xticksLables[i])
    plt.xticks(tick1, ticklable1, rotation='vertical')
    plt.legend(loc='upper right')
    plt.grid(True)
    plt.title(figname)
    #ax.set_xticklabels(xticksLables)
    plt.savefig(pwd() + "/images/" + figname + ".png")
    plt.cla()
    
def plotHist(dist, bins, title, x_label='', y_label=''):
    x = {}
    count248 = 0
    for item in dist:
        x[item] = math.floor(dist[item])
        if(x[item] == 248):
            count248 += 1
    print("count248: " + str(count248))
    n, bins, patches = plt.hist(list(x.values()), bins, facecolor='g', alpha=0.75)
    plt.grid(True)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.savefig(pwd() + "/images/" + title + ".png")
    plt.cla()
    return [n, bins, patches]

def plotBar(vals, labels=None, title='', subplots=1):
    fig = plt.figure()
    if(type(vals[0]) is not list):
        vals = [vals]
        labels = [labels]
    for p in range(subplots):
        curr_vals = vals[p]
        if(labels):
            curr_labels = labels[p]
        curr_subplot = 100*subplots + 10 + p+1
        ax = fig.add_subplot(curr_subplot)
        ind = np.arange(len(curr_vals)) 
        width = 0.35
        ax.bar(ind,curr_vals , width)
        ax.set_xticks(ind+width)
        if(labels):
            xtickNames = ax.set_xticklabels(curr_labels)
        plt.setp(xtickNames, rotation=45, fontsize=10)
        ax.set_xlim(-width,len(ind)+width)
        #ax.set_title
    plt.title(title)
    plt.savefig(pwd() + 'images/' + title + ".png")
    plt.cla()
    
    
def dist2absolute(dist, total, threshold):
    absDist = {}
    updatedTotal = 0
    for item in dist:
        count = math.floor(dist[item]*total)
        if(count > 0):
            absDist[item] = count
            updatedTotal += count
    return [absDist, updatedTotal]


def absolute2dist(abs_dist, total):
    dist = {}
    for item in abs_dist:
        dist[item] = abs_dist[item] / total
    return dist
    
def fixDist(dist):
    distsum = 0
    for item in dist:
        distsum += dist[item]
    for item in dist:
        dist[item] /= distsum
    return dist

def affimetrix2geneSymbol(aff_filename):
    aff2symbol = {}
    f_affgenes = open(pwd() + aff_filename, 'r')
    i = 0
    for line in f_affgenes:
        if(i > 25):
            row = line.split('\t')
            probeID = row[0].strip()
            genesymbol = row[14].strip()
            aff2symbol[probeID] = genesymbol
        i += 1
    f_affgenes.close()
    return aff2symbol

def affID2miRname():
    aff2miR = {}
    f_namemap = open(pwd() + "expression_data/HeLa/miRnameID.tab", 'r')
    i = 0
    for line in f_namemap:
        if(i > 0):
            try:
                row = line.split('\t')
                affID = row[0].strip()
                mirID = row[2].strip()
                aff2miR[affID] = miRnameFormat(mirID)
            except:
                print(line)
        i += 1
    f_namemap.close()
    return aff2miR

def geneslist(filename):
    # assume that the format is: Ensembl Gene ID    Gene Start (bp)    Gene End (bp)    Chromosome Name    Strand    Associated Gene Name
    f_genes = open(pwd() + filename, 'r')
    i = 0
    col2idx = {}
    chr2genes = {}
    for line in f_genes:
        row = line.split('\t')
        if(i == 0):
            for idx in range(len(row)):
                col2idx[row[idx].strip()] = idx
        else:
            curr_chr = row[col2idx['Chromosome Name']]
            strand = row[col2idx['Strand']]
            genename = row[col2idx['Associated Gene Name']]
            start = row[col2idx['Gene Start (bp)']]
            end = row[col2idx['Gene End (bp)']]
            if(curr_chr in chr2genes):
                strand2genes = chr2genes[curr_chr]                
            else:
                strand2genes = {1:{}, -1:{}}
            gene2pos = strand2genes[strand]
            gene2pos[genename] = start + '_' + end
            strand2genes[strand] = gene2pos
            chr2genes[curr_chr] = strand2genes
    f_genes.close()
    return chr2genes

def genesMap(ensembl_filename):
    # Ensembl Gene ID Chromosome Name Gene Start (bp) Gene End (bp)   Strand  Associated Gene Name    EntrezGene ID   HGNC symbol     HGNC ID(s)
    strandLabel = {'1' : '+', '-1': '-'}
    ensembl2genename = {}
    genes = {}
    f_genemap = open(pwd() + ensembl_filename, 'r')
    i = 0
    for line in f_genemap:
        if(i > 0):
            row = line.split('\t')
            ensembl2genename[row[0].strip()] = row[5].strip()
            curr_chr = row[1].strip()
            gene_start = int(row[2].strip())
            gene_end = int(row[3].strip())
            gene_strand = strandLabel[row[4].strip()]
            genes[row[5].strip()] = objects.Gene(row[5].strip(), curr_chr, gene_strand, gene_start, gene_end)
            #genes[row[0].strip()] = curr_chr + '_' + gene_start + '_' + gene_end + '_' + gene_strand
        i += 1
    f_genemap.close()
    return [ensembl2genename, genes]


def miRexpression(miRexp_file, name_format=True):
    """
    :param miRexp_file:
    :return: {miR:rpm} map
    """
    miR2exp = {}
    f_miRexp = open(pwd() + miRexp_file, 'r')
    i = 0
    total = 0
    for i, line in enumerate(f_miRexp):
        if i > 0:
            row = line.split('\t')
            miR = row[0].strip()
            rpm = float(row[1].strip())
            if name_format:
                miR = miRnameFormat(miR)
            if miR in miR2exp:
                miR2exp[miR] += rpm
            else:
                miR2exp[miR] = rpm
            total += rpm
    f_miRexp.close()
    for miR in miR2exp:
        miR2exp[miR] /= total
    return miR2exp

def removeUnderThreshold(dist, total, threshold, shouldLOG, miRormRNA):
    totalPercentage = 0
    removedCount = 0
    for item in dist.keys():
        if(math.floor(dist[item]*total) < threshold):
            del dist[item]
            removedCount += 1
        else:
            totalPercentage += dist[item]
    errorCount = 0
    #updatedTotal = 0
    for item in dist.keys():
        dist[item] /= totalPercentage
        #updatedTotal += math.floor(dist[item]*total)
        if(math.floor(dist[item]*total) < threshold):
            errorCount += 1
    if(shouldLOG):
        logINFO(miRormRNA + " removed below threshold " + str(threshold) + ": " + str(removedCount) + ", total: " + str(total))
        if(errorCount):
            logERROR("error count: " + str(errorCount))
    return dist

# Returns a descending list of sorted variance of the given input lists
def getVarianceSorted(key2list):
    varsAbsMap = {}
    varsRelativeMap = {}
    #threshold = 10
    #varsRelativeThreshold = {}
    for key in key2list:        
        lst = key2list[key]
        varsAbsMap[key] = np.var(lst)
        total = lst[0]
        relative_lst = lst[:]
        for i in range(len(lst)):
            relative_lst[i] = 100*float(lst[i])/float(total)
        varsRelativeMap[key] = np.var(relative_lst)
        # Adding the total value at the beginning to the relative lists
        #relative_lst.insert(0,total)
        #if(total >= threshold):
        #    varsRelativeThreshold[key] = varsRelativeMap[key]
    sortedAbsKeys = sorted(varsAbsMap.iteritems(), key=operator.itemgetter(1), reverse=True)
    sortedRelativeKeys = sorted(varsRelativeMap.iteritems(), key=operator.itemgetter(1), reverse=True)
    #sortedRelativeThreshold = sorted(varsRelativeThreshold.iteritems(), key=operator.itemgetter(1), reverse=True)
    return [sortedAbsKeys, sortedRelativeKeys]

# Return a descending list of sorted retention of the the given input list
def sort_by_retention(key2list, th):
    retention = {}
    for key in key2list:
        lst = key2list[key]
        # check that the initial amount is above th
        if(float(lst[1]) >= th):
            curr_ret = float(lst[-1])/float(lst[0])*100
            retention[key] = curr_ret
    retention_sorted = sorted(retention.items(), key=operator.itemgetter(1))
    return retention_sorted
        

# Returns:
# 1. mRNA2exp: {gene:{'[chr:strand]':rpm/mRNAcount}}
# 2. mRNAdist: {mRNA+pos:rpm}
# 3. miRdist: {miR+pos:rpm}
def mRNAexpressionDetailed(mRNAexp_file):
    # {gene: {chr:strand : rpm}}
    mRNA2exp = {}
    miR2exp = {}
    f_mRNAexp = open(pwd()  + 'expression_data/'  + mRNAexp_file, 'r')
    mRNAcount = 0
    miRcount = 0
    for line in f_mRNAexp:
        row = line.split('\t')
        gene = row[3].strip()
        curr_chr = row[0].strip()
        if(curr_chr.startswith('chr')):
            curr_chr = curr_chr[3:]
        strand = row[4].strip()
        rpm = float(row[6].strip())
        if(gene.startswith('MIR')):
            miRcount += rpm
            gene2exp = miR2exp
        else:
            mRNAcount += rpm
            gene2exp = mRNA2exp
        if(gene in gene2exp):
            pos2rpm = gene2exp[gene]
        else:
            pos2rpm = {}
        pos2rpm[curr_chr + ':' + strand] = rpm
        gene2exp[gene] = pos2rpm
    f_mRNAexp.close()
    mRNAdist = {}
    miRdist = {}
    for mRNA in mRNA2exp:
        pos2rpm = mRNA2exp[mRNA]
        for pos in pos2rpm:
            pos2rpm[pos] = pos2rpm[pos]/mRNAcount
            mRNAdist[mRNA+'_'+pos] = pos2rpm[pos]
    for miR in miR2exp:
        pos2rpm = miR2exp[miR]
        for pos in pos2rpm:
            pos2rpm[pos] = pos2rpm[pos]/miRcount
            miRdist[miR+'_'+pos] = pos2rpm[pos]
    return [mRNA2exp, mRNAdist]

def mRNAexpression(mRNAexp_file, exp_col=1, gene_file=None):
    relevant_genes = []
    if gene_file:
        with open(gene_file) as f:
            for line in f:
                relevant_genes.append(line.strip())
    mRNA2exp = {}
    miR2exp = {}
    f_mRNAexp = open(pwd() + 'expression_data/'  + mRNAexp_file, 'r')
    mRNAcount = 0
    miRcount = 0
    for i, line in enumerate(f_mRNAexp):
        if i > 0:
            row = line.split('\t')
            gene = row[0].strip().split(',')[0].strip()
            gene = gene.split('///')[0].strip()
            rpm = float(row[exp_col].strip())
            if gene.startswith('MIR'):
                miRcount += rpm
                gene2exp = miR2exp
            else:
                mRNAcount += rpm
                gene2exp = mRNA2exp
            if gene_file is not None and gene not in relevant_genes:
                rpm = 0
            gene2exp[gene] = rpm
    f_mRNAexp.close()
    for mRNA in mRNA2exp:
        mRNA2exp[mRNA] /= mRNAcount
    for miR in miR2exp:
        miR2exp[miR] /= miRcount
    return mRNA2exp


def miRnameFormat(miR, only_pre=False):
    orig_miR = miR
    if only_pre:
        miR = re.sub(r'-.p.*$', "", miR)
    if(miR.endswith('*')):
        miR = miR[:-1]
    if miR.startswith('hsa') or miR.startswith('mmu'):
        miR = miR[4:]
    miR = miR.lower()
    miRarr = miR.split('-')
    #miRarr[1] = re.sub("\D", "", miRarr[1])
    #miR = '-'.join(miRarr)
    if(len(miRarr) >= 2):
        #miRarr[1] = re.sub("\D", "", miRarr[1])
        if len(miRarr) > 2:
            if "p" not in miRarr[2]:
                miRarr.remove(miRarr[2])
            if "." in miRarr[-1]:
                miRarr[2] = miRarr[2].split(".")[0]
        if only_pre:
            miR = ('-').join([miRarr[0], miRarr[1]])
        else:
            miR = ('-').join(miRarr)
    else:
        logERROR("error in methods.miRnameFormat: " + miR + ", orig miR: " + orig_miR)
    return miR


def processCLASHdata(clash_file, name_format=True, table=None, mir_exp=None, mRNA_exp=None):
    [ensembl2genename, geneMap] = genesMap('expression_data/ensembl_extend.tab')
    miRmRNAexp = {}
    #{miR : {mRNA : {bs : chimera_reads+non_chimera/total}}}
    mRNA2chimera = {}
    #{mRNA : [chimera_count , non_chimera_count]}
    # total = chimer + non_chimera
    total = 0
    i = 0
    f_in = open(pwd() + 'expression_data/' + clash_file)
    nonchimera_missingdata = 0
    mRNAcount = {}
    genecount = {}
    #ensembl2genename = {}
    for line in f_in:
        if i > 30:
            gene = ''
            row = line.split('\t')
            miR = row[1].strip().split('_')[2].lower()
            mRNA = row[5]
            mRNAstart = int(row[6].strip())
            mRNAend = int(row[7].strip())
            chimera_reads = float(row[9].strip())
            nonchimeraRNA = row[24].strip()
            genename = mRNA.split('_')[2]
            ensembl = mRNA.split('_')[0]
            if name_format:
                miR = miRnameFormat(miR)
            #===================================================================
            # if(ensembl in ensembl2genename):
            #    if not (ensembl2genename[ensembl] == genename):
            #        print mRNA
            # else:
            #    ensembl2genename[ensembl] = genename
            #===================================================================
            #genename = mRNA
            #mRNA = ensembl + '_' + genename
            if(genename in geneMap):
                gene = geneMap[genename]
            elif(ensembl in ensembl2genename and ensembl2genename[ensembl] in geneMap):
                genename = ensembl2genename[ensembl]
                gene = geneMap[genename]
            else:
                gene = objects.Gene(genename, '', '', 0,0)
            bs = objects.BS(gene, str(mRNAstart) + '_' + str(mRNAend) , 'relativeBS')
            #bs = mRNAstart + '_' + mRNAend
            if miR not in miRmRNAexp:
                miRmRNAexp[miR] = {genename : {bs : chimera_reads}}
            elif(genename not in miRmRNAexp[miR]):
                miRmRNAexp[miR][genename] = {bs : chimera_reads}
            else:
                miRmRNAexp[miR][genename][bs] = chimera_reads
            try:
                chimera_reads = float(chimera_reads)
            except:
                logERROR(miR + '\t' + genename + '\t' + str(chimera_reads) + '\t' + str(nonchimeraRNA))
            try:
                nonchimeraRNA = float(nonchimeraRNA)
                miRmRNAexp[miR][genename][bs] += nonchimeraRNA
                total += chimera_reads + nonchimeraRNA
            except:
                nonchimera_missingdata += 1
            if(genename not in mRNA2chimera):
                mRNA2chimera[genename] = [0,0]
            if(type(chimera_reads) is float):
                mRNA2chimera[genename][0] += chimera_reads
            if(type(nonchimeraRNA) is float):
                mRNA2chimera[genename][1] += nonchimeraRNA
            if(genename not in mRNAcount):
                mRNAcount[genename] = 1
            if(genename not in genecount):
                genecount[genename] = 1
        i += 1
    f_in.close()
#    print "nonchimera_missingdata: " + str(nonchimera_missingdata)
#    print "total: " + str(total)
#    print 'mRNA count: ' + str(len(mRNAcount.keys()))
#    print 'gene count: ' + str(len(genecount.keys()))
    clash_in_table = {}
    clash_in_exp = {}
    total_pairs = 0
    f_out = open("C:/miR/output/HEK/clash_pairs_expressed.tab", "w")
    f_out.write("miRNA\tGene\tmiRNA exp.\tgene exp.\n")
    for miR in miRmRNAexp:
        for mRNA in miRmRNAexp[miR]:
            sum_bs = 0
            for bs in miRmRNAexp[miR][mRNA]:
                sum_bs += miRmRNAexp[miR][mRNA][bs]
                miRmRNAexp[miR][mRNA][bs] = miRmRNAexp[miR][mRNA][bs]/total
                total_pairs += 1
            #print(miR, "\t", mRNA, "\t", sum_bs)
            if table is not None and miR in table and mRNA in table[miR]:
                clash_in_table[miR + '_' + mRNA] = sum_bs
                if mir_exp is not None and mRNA_exp is not None:
                    if miR in mir_exp and mRNA in mRNA_exp:
                        clash_in_exp[miR + '_' + mRNA] = (mir_exp[miR], mRNA_exp[mRNA])
                        f_out.write(miR + '\t' + mRNA + '\t' + str(mir_exp[miR]) + "\t" + str(mRNA_exp[mRNA]) + "\n")
            else:
                pass
                # f_out.write(miR + '\t' + mRNA + '\t' + str(sum_bs) + "\n")
    f_out.close()
    logINFO("total pairs in CLASH: " + str(total_pairs))
    logINFO("total CLASH pairs in table: " + str(len(clash_in_table)))
    logINFO("total CLASH pairs expressed: " + str(len(clash_in_exp)))
    return miRmRNAexp, clash_in_table, clash_in_exp

# Given miR, mRNA, and the score table it returns a BS2score map, with all possible BSs and their scores
def BSscore(gene, miR, table):
    BSobj2score = {}
    if(miR in table):
        if(gene in table[miR]):
            bs2score = table[miR][gene]
            for bs in bs2score:
                BSobj = objects.BS(gene, bs, '[hg19:chr:start-end:-]')
                BSobj2score[BSobj] = bs2score[bs]
        else:
            logERROR("gene is not found in the table")
    else:
        logERROR("miR is not found in the table")
    if(not BSobj2score):
        logERROR("BS2score is empty")
    return BSobj2score

def compareCLASH2miRBASE(clash_file, miRBase_file):
    miR2exp = {}
    # {miR : [miRBase, clash]}
    mature2primary = mat2primary('hsa.gff3')
    miRinCLASH = {}
    miRinmiRBase = {}
    miRinboth = 0    
    miRfound = {}
    f_clash = open(pwd() + 'expression_data/' + clash_file, 'r')
    i = 0
    for line in f_clash:
        if(i > 30):
            row = line.split('\t')
            miR = row[1].strip()
            chimera_reads = float(row[9].strip())
            if(miR not in miR2exp):
                miR2exp[miR] = [0,0]
            miR2exp[miR][0] += chimera_reads
            miRinCLASH[miR] = 1
            miRfound[miR] = 0
        i += 1
    f_clash.close()
    f_miRBase = open(pwd() + 'expression_data/' + miRBase_file, 'r')
    i = 0
    for line in f_miRBase:
        if(i > 0):
            row = line.split('\t')
            miRID = row[0].strip()[4:]
            miRacc = row[1].strip()
            rpm = row[2].strip()
            try:
                rpm = float(rpm)
            except:
                if(rpm == '-'):
                    rpm = 0
                else:
                    logERROR(miRID + '\t' + rpm)
#            found = 0
            for miR in miR2exp:
                miR_lst = miR.split('_')
                if(len(miR_lst) > 2):
                    curr_miRacc = miR_lst[0].strip()
                    if(curr_miRacc in mature2primary):
                        curr_miRacc = mature2primary[miR_lst[0].strip()]
                    curr_miRID = miR_lst[2].strip()
                    if(miRID == curr_miRID or miRacc == curr_miRacc):
                        miR2exp[miR][1] += rpm
                        miRinboth += 1
                        miRfound[miR] = 1
#                        found = 1
#            if not found:
#                if(miRID + '_' + miRacc in miR2exp):
#                    print miRID + '_' + miRacc
#                miR2exp[miRID + '_' + miRacc] = [0, rpm]
            miRinmiRBase[miRacc] = 1
        i += 1
    f_miRBase.close()
#    for miR in miRfound:
#        if(miRfound[miR] == 0):
#            print miR
#    print 'miRinCLASH: ' + str(len(miRinCLASH.keys()))
#    print 'miRinmiRBase: ' + str(len(miRinmiRBase.keys()))
#    print 'miRinboth: ' + str(miRinboth)
    return miR2exp

def mat2primary(miRBase_file):
    mature2primary = {}
    i = 0
    f_miRBase = open(pwd() + 'expression_data/' + miRBase_file, 'r')
    primaryID = ''
    for line in f_miRBase:
        if(i > 12):
            row = line.split('\t')
            miRtype = row[2].strip()
            id_data = row[8].strip()
            if(miRtype == 'miRNA_primary_transcript'):
                primaryID = id_data.split(';')[0].split('=')[1]
            else:
                matID = id_data.split(';')[0].split('=')[1]
                mature2primary[matID] = primaryID
            if not (id_data.split(';')[0].startswith('ID')):
                logERROR(line)
        i += 1
    f_miRBase.close()
    return mature2primary

def id2name(miRBase_file):
    #miRBase_file = 'hsa.gff3'
    id2name = {}
    ID2matID = {}
    matID2name = {}
    i = 0
    f_miRBase = open(pwd() + 'expression_data/' + miRBase_file, 'r')
#    primaryID = ''
    for line in f_miRBase:
        if(i > 12):
            row = line.split('\t')
            miRtype = row[2].strip()
            id_data = row[8].strip()
            data_row = id_data.split(';')
            data_map = {}
            for data in data_row:
                curr_info = data.split('=')
                data_map[curr_info[0]] = curr_info[1]
            if(miRtype == 'miRNA_primary_transcript'):                
                primaryID = data_map['ID']
                name = data_map['Name']
                id2name[primaryID] = name
                if(primaryID in ID2matID):
                    logERROR(primaryID + " already exists")
                else:
                    ID2matID[primaryID] = []
            else:
                matID = data_map['ID']
                name = data_map['Name']
                primaryID = data_map['Derives_from']
                matID2name[matID] = name
                if(primaryID in ID2matID):
                    ID2matID[primaryID].append(matID)
                else:
                    logERROR(primaryID + " does not exist")
                    ID2matID[primaryID] = [matID]                    
        i += 1
    f_miRBase.close()
    return [id2name, ID2matID, matID2name]

def compareCLASH2miRexp(clash_file, miRexp_file):
    miR2exp = {}
    # {miR : [miRBase, clash]}
    norm_factor_clash = 11.721864
    ID2name = id2name('hsa.gff3.txt')[0]
    miRinCLASH = {}
    miRinexp = {}
    miRinboth = 0    
    miRfound = {}
    f_clash = open(pwd() + 'expression_data/' + clash_file, 'r')
    i = 0
    for line in f_clash:
        if(i > 30):
            row = line.split('\t')
            miR = row[1].strip()
            chimera_reads = float(float(row[9].strip())/norm_factor_clash)
            if(miR not in miR2exp):
                miR2exp[miR] = [0,0]
            miR2exp[miR][0] += chimera_reads
            miRinCLASH[miR] = 1
            miRfound[miR] = 0
        i += 1
    f_clash.close()
    f_miRexp = open(pwd() + 'expression_data/' + miRexp_file, 'r')
    i = 0
    for line in f_miRexp:
        if(i > 0):
            row = line.split('\t')
            miRname = row[0].strip()
            rpm = row[1].strip()
            try:
                rpm = float(rpm)
            except:
                if(rpm == '-'):
                    rpm = 0
                else:
                    logERROR(miRname + '\t' + rpm)
#            found = 0
            for miR in miR2exp:
                miR_lst = miR.split('_')
                if(len(miR_lst) > 2):
                    curr_miRid = miR_lst[0].strip()
                    if(curr_miRid in ID2name):
                        curr_miRname = ID2name[miR_lst[0].strip()]
                    else:
                        curr_miRname = miR_lst[2].strip().lower()
                    if(curr_miRname.endswith('*')):
                        curr_miRname = curr_miRname[:-1]
                    if(curr_miRname.endswith('-3p') or curr_miRname.endswith('-5p')):
                        curr_miRname = curr_miRname[:-3]
                    if(miRname == curr_miRname):
                        miR2exp[miR][1] += rpm
                        miRinboth += 1
                        miRfound[miR] = 1
#                        found = 1
#            if not found:
#                if(miRname + '_' + miRid in miR2exp):
#                    print miRname + '_' + miRid
#                miR2exp[miRname + '_' + miRid] = [0, rpm]
            miRinexp[miRname] = 1
        i += 1
    f_miRexp.close()
#    not_found_count = 0
#    for miR in miRfound:
#        if(miRfound[miR] == 0):
#            print miR
#            not_found_count += 1
#    print 'miRinCLASH: ' + str(len(miRinCLASH.keys()))
#    print 'miRinexp: ' + str(len(miRinexp.keys()))
#    print 'miRinboth: ' + str(miRinboth)
    return miR2exp


def compareCLASH2mRNAexp(clash_file, exp_file):
    mRNA2exp = {}
    # {mRNA : [exp, clash]}
    mRNA = {}
    mRNAinExp = {}
    mRNAinboth = 0    
    mRNAfound = {}
    f_clash = open(pwd() + 'expression_data/' + clash_file, 'r')
    i = 0
    for line in f_clash:
        if(i > 30):
            row = line.split('\t')
            genename = row[5].strip().split('_')[2]
            chimera_reads = float(row[9].strip())
            if(genename not in mRNA2exp):
                mRNA2exp[genename] = [0,0]
            mRNA2exp[genename][0] += chimera_reads
            mRNA[genename] = 1
            mRNAfound[genename] = 0
        i += 1
    f_clash.close()
    exp_file = open(pwd() + 'expression_data/' + exp_file, 'r')
    i = 0
    for line in exp_file:
        row = line.split('\t')
        genename = row[3].strip()
        rpm = row[6].strip()
        try:
            rpm = float(rpm)
        except:
            if(rpm == '-'):
                rpm = 0
            else:
                logERROR(genename + '\t' + rpm)
#            found = 0
        if(genename in mRNA2exp.keys()):
            mRNA2exp[genename][1] += rpm
            mRNAinboth += 1
            mRNAfound[genename] = 1
        mRNAinExp[genename] = 1
        i += 1
    exp_file.close()
#    for target in mRNAfound:
#        if(mRNAfound[target] == 0):
#            print target
#    print 'mRNA: ' + str(len(mRNA.keys()))
#    print 'mRNAinExp: ' + str(len(mRNAinExp.keys()))
#    print 'mRNAinboth: ' + str(mRNAinboth)
    return mRNA2exp

def comparemiRexp(data1, data2):
    # the format in this data is MIR4316[chr17:-]
    miRrpm = {}
    familiesrpm = {}
    f_data1 = open(pwd() + 'expression_data/' + data1, 'r')
    i = 0
    for line in f_data1:
        if(i > 0):
            row = line.split('\t')
            miR = row[0].strip().split('[')[0]
            rpm = float(row[1].strip())
            c = 0
            if(miR[0:3] == 'MIR'):
                if(miR[3:4].isdigit()):
                    curr_miR = "mir-"
                    c = 3
                elif(miR[3:6] == 'LET'):
                    curr_miR = "let-"
                    c = 6
                while(miR[c:c+1].isdigit()):
                    curr_miR += miR[c:c+1]
                    c += 1
                miRfamily = curr_miR
                if(len(miR) > c):
                    if(miR[c:c+1] == '-'):
                        curr_miR += miR[c:]
                    else:
                        curr_miR += miR[c:c+1].lower()
                        c += 1
                        if(len(miR) > c and miR[c:c+1].isdigit()):
                            curr_miR += '-' + miR[c:]
            else:
                logERROR(miR)
            miRrpm[curr_miR] = [rpm, 0]
            if(miRfamily not in familiesrpm):
                    familiesrpm[miRfamily] = [0, 0]
            familiesrpm[miRfamily][0] += rpm
        i += 1
    f_data1.close()
    # the format in this data is mir-103-1
    f_data2 = open(pwd() + 'expression_data/' + data2, 'r')
    i = 0
    for line in f_data2:
        if(i > 0):
            row = line.split('\t')
            miR = row[0].strip()
            rpm = float(row[1].strip())
            if(miR in miRrpm):
                miRrpm[miR][1] = rpm
            else:
                miRrpm[miR] = [0, rpm]
            c = 4
            miRfamily = miR[0:4]
            while(miR[c:c+1].isdigit()):
                miRfamily += miR[c:c+1]
                c += 1
            if(miRfamily in familiesrpm):
                familiesrpm[miRfamily][1] = rpm
            else:
                familiesrpm[miRfamily] = [0, rpm]
        i += 1
    f_data2.close()
    return [miRrpm, familiesrpm]