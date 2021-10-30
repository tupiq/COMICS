#!/usr/bin/python
# coding: utf-8
import csv
import time
import random
import methods
import sys, getopt
from numpy import random as nprand
from copy import deepcopy
import src.constants as constants
import pandas as pd
 
out_file = ''
in_file = ''


# python miRtargetTable.py -i human_predictions_S_C_aug2010.txt -o targetmiRtable_human.tab
def main(argv):   
    try:
        opts, args = getopt.getopt(argv,"o:i:",["in=","out="])
    except getopt.GetoptError:
        print('miRtargetTable.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('miRtargetTable.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-o", "--out"):
            global out_file
            out_file = arg
        elif opt in ("-i", "--in"):
            global in_file
            in_file = arg

if __name__ == "__main__":
    main(sys.argv[1:])
    
timeList = []
miRs = {}
targets = {}
gene2coordinate = {}


def target_mir_score_table_old(table_filename, source):
    timeList.append(time.time())
    #{genename: {transcript: {bs: {miR : score}}}}
    target2trans2bs2miR = {}
    miR2mRNAlst = {}
    col2idx = {}
    genes = set()
    miRslst = set()
    f_table = open(methods.pwd() + "expression_data/" + table_filename, 'r')
    i = 0
    for line in f_table:
        row = line.strip().split('\t')
        if i == 0:
            for idx in range(len(row)):
                col2idx[row[idx].strip()] = idx
        else:
            if source == constants.MICRORNA_ORG:
                miRname = row[col2idx['mirna_name']].strip()
                miRacc = row[col2idx['#mirbase_acc']].strip()
                genename = row[col2idx['gene_symbol']].strip()
                score = float(row[col2idx['mirsvr_score']].strip())
                transcript = row[col2idx['ext_transcript_id']].strip()
                coordinate = row[col2idx['genome_coordinates']].strip()
            elif source == constants.TARGET_SCAN:
                miRname = row[col2idx['miR']].strip()
                genename = row[col2idx['gene']].strip()
                score = float(row[col2idx['score']].strip())
                coordinate = row[col2idx['coordinates']].strip()
                miRacc = miRname
                transcript = genename
            elif source == constants.RANDOM_SCORES:
                miRname = row[col2idx['miRNA']].strip()
                genename = row[col2idx['target']].strip()
                score = float(row[col2idx['score']].strip())
                coordinate = row[col2idx['MBS']].strip()
                miRacc = miRname
                transcript = genename
            elif source.startswith("both"):
                miRname = row[col2idx['miR']].strip()
                genename = row[col2idx['gene']].strip()
                score1 = float(row[col2idx['score1']].strip())
                score2 = float(row[col2idx['score2']].strip())
                score = (score1,score2)
                coordinate = row[col2idx['coordinates']].strip()
                miRacc = miRname
                transcript = genename
            # removing 'hsa-' prefix:
            #miRname = miRname[4:].lower()
            miRname = methods.miRnameFormat(miRname)
            if genename in gene2coordinate:
                coo_lst = gene2coordinate[genename]
                if coordinate not in coo_lst:
                    coo_lst.append(coordinate)
                    gene2coordinate[genename] = coo_lst
            else:
                gene2coordinate[genename] = [coordinate]
            genes.add(genename)
            miRslst.add(miRname)
            if miRname in miR2mRNAlst:
                miR2mRNAlst[miRname].add(genename)
            else:
                miR2mRNAlst[miRname] = {genename}
            if miRname not in miRs:
                miRs[miRname] = miRacc
            #elif(not (miRs[miRname] == miRacc)):
            #    print(miRname + ': ' + miRs[miRname] + miRacc)
            if genename not in targets:
                targets[genename] = {transcript: [coordinate]}
            else:
                if transcript not in targets[genename]:
                    targets[genename][transcript] = [coordinate]
                else:
                    if coordinate not in targets[genename][transcript]:
                        targets[genename][transcript].append(coordinate)
            if genename in target2trans2bs2miR:
                trans2bsmap = target2trans2bs2miR[genename]
                if transcript in trans2bsmap:
                    bs_map = trans2bsmap[transcript]
                else:
                    bs_map = {}
                if coordinate in bs_map:
                    miRscoreMap = bs_map[coordinate]
                else:
                    miRscoreMap = {}
            else:
                trans2bsmap = {}
                bs_map = {}
                miRscoreMap = {}
            miRname = methods.miRnameFormat(miRname)
            if miRname in miRscoreMap:
                #print "different scores for the same miR: " + miRname + " in transcript: " + refseq
                miRscoreMap[miRname] = min([score, miRscoreMap[miRname]])
            else:
                miRscoreMap[miRname] = score
            bs_map[coordinate] = miRscoreMap
            trans2bsmap[transcript] = bs_map
            target2trans2bs2miR[genename] = trans2bsmap           
        i += 1
    f_table.close()
    methods.timer(timeList, 0, 0)
    for target in target2trans2bs2miR:
        transcripts = target2trans2bs2miR[target].keys()
        curr_transcript = nprand.choice(list(transcripts))
        target2trans2bs2miR[target] = target2trans2bs2miR[target][curr_transcript]
    return [target2trans2bs2miR, miRslst]

def target_mir_score_table(table_filename, source):
    timeList.append(time.time())
    #{genename: {transcript: {bs: {miR : score}}}}
    target2mir2bs2score= {}
    col2idx = {}
    f_table = open(methods.pwd() + "expression_data/" + table_filename, 'r')
    for i, line in enumerate(f_table):
        row = line.strip().split('\t')
        if i == 0:
            for idx in range(len(row)):
                col2idx[row[idx].strip()] = idx
        else:
            if source == constants.MICRORNA_ORG:
                miRname = row[col2idx['mirna_name']].strip()
                miRacc = row[col2idx['#mirbase_acc']].strip()
                genename = row[col2idx['gene_symbol']].strip()
                score = float(row[col2idx['mirsvr_score']].strip())
                transcript = row[col2idx['ext_transcript_id']].strip()
                coordinate = row[col2idx['genome_coordinates']].strip()
            elif source == constants.TARGET_SCAN:
                miRname = row[col2idx['miR']].strip()
                genename = row[col2idx['gene']].strip()
                score = float(row[col2idx['score']].strip())
                coordinate = row[col2idx['coordinates']].strip()
                miRacc = miRname
                transcript = genename
            elif source == constants.RANDOM_SCORES:
                miRname = row[col2idx['miRNA']].strip()
                genename = row[col2idx['target']].strip()
                score = float(row[col2idx['score']].strip())
                coordinate = row[col2idx['MBS']].strip()
                miRacc = miRname
                transcript = genename
            elif source.startswith("both"):
                miRname = row[col2idx['miR']].strip()
                genename = row[col2idx['gene']].strip()
                score1 = float(row[col2idx['score1']].strip())
                score2 = float(row[col2idx['score2']].strip())
                score = (score1,score2)
                coordinate = row[col2idx['coordinates']].strip()
                miRacc = miRname
                transcript = genename
            mir = methods.miRnameFormat(miRname)
            if genename not in target2mir2bs2score:
                target2mir2bs2score[genename] = {}
            if mir not in target2mir2bs2score[genename]:
                    target2mir2bs2score[genename][mir] = {}
            target2mir2bs2score[genename][mir][coordinate] = score
    f_table.close()
    methods.timer(timeList, 0, 0)
    return target2mir2bs2score

def miRtargetScore(table_filename, detailed, source):   
    timeList = [time.time()]
    #{miR: {gene : transcript: {bs: score}}}}
    miR2gene2trans2bs2score = {}
    miRs = {}
    targets = {}
    col2idx = {}
    bs2id = {}
    counterID = 1
    f_table = open(methods.pwd() + "expression_data/" + table_filename, 'r')
    i = 0
    orig_mirs = set()
    miRs_set = set()
    for line in f_table:
        if i%1000000 == 0:
            print("curr i: " + str(i))       
        row = line.strip().split('\t')
        if i == 0:
            for idx in range(len(row)):
                col2idx[row[idx].strip()] = idx
        else:
            if(source == "microrna.org"):
                miRname = row[col2idx['mirna_name']].strip()
                miRacc = row[col2idx['#mirbase_acc']].strip()
                genename = row[col2idx['gene_symbol']].strip()
                score = float(row[col2idx['mirsvr_score']].strip())
                transcript = row[col2idx['ext_transcript_id']].strip()
                coordinate = row[col2idx['genome_coordinates']].strip()
            elif(source == "TargetScan"):
                miRname = row[col2idx['miR']].strip()
                genename = row[col2idx['gene']].strip()
                score = float(row[col2idx['score']].strip())
                coordinate = row[col2idx['coordinates']].strip()
                miRacc = miRname
                transcript = genename
            elif source == constants.RANDOM_SCORES:
                miRname = row[col2idx['miRNA']].strip()
                genename = row[col2idx['target']].strip()
                score = float(row[col2idx['score']].strip())
                coordinate = row[col2idx['MBS']].strip()
                miRacc = miRname
                transcript = genename
            elif(source.startswith("both")):
                miRname = row[col2idx['miR']].strip()
                genename = row[col2idx['gene']].strip()
                score1 = float(row[col2idx['score1']].strip())
                score2 = float(row[col2idx['score2']].strip())
                score = (score1,score2)
                coordinate = row[col2idx['coordinates']].strip()
                miRacc = miRname
                transcript = genename
            # [hg19:chr:start-end:-]
            curr_chr = coordinate.split(':')[1]
            strand = coordinate.split(':')[3][0]
            #if('hg19' not in coordinate):
            #    methods.logERROR("hg19 not in " + coordinate)
            if(coordinate not in bs2id):
                currID = counterID
                counterID += 1
                bs2id[coordinate] = currID
            else:
                currID = bs2id[coordinate]          
            coordinate = coordinate.replace('hg19', str(currID))            
            if(detailed):
                genename += '_' + curr_chr + ':' + strand
            # removing 'hsa-' prefix:
            #miRname = miRname[4:].lower()
            if(source == "microrna.org"):
                if miRname not in miRs:
                    miRs[miRname] = miRacc
                elif miRs[miRname] != miRacc:
                    print(miRname + ': ' + miRs[miRname] + miRacc)            
            orig_mirs.add(miRname)
            miRname = methods.miRnameFormat(miRname)
            miRs_set.add(miRname)
            if genename not in targets:
                targets[genename] = {transcript: [coordinate]}                
            else:
                if transcript not in targets[genename]:
                    targets[genename][transcript] = [coordinate]
                elif(coordinate not in targets[genename][transcript]):
                    targets[genename][transcript].append(coordinate)                
            if(miRname in miR2gene2trans2bs2score):
                gene2trans2bs2score = miR2gene2trans2bs2score[miRname]
            else:
                gene2trans2bs2score = {}
            if(genename in gene2trans2bs2score):
                trans2bsmap = gene2trans2bs2score[genename]
            else:
                trans2bsmap = {}
            if(transcript in trans2bsmap):
                bs_map = trans2bsmap[transcript]
            else:
                bs_map = {}
            if(coordinate in bs_map):
                #print("found coordinate: " + coordinate + " miR: " + miRname + " gene: " + genename + " transcript: " + transcript)
                bs_map[coordinate] = min(bs_map[coordinate], score)
            else:
                bs_map[coordinate] = score
            trans2bsmap[transcript] = bs_map
            gene2trans2bs2score[genename] = trans2bsmap
            miR2gene2trans2bs2score[miRname] = gene2trans2bs2score
        i += 1
    f_table.close()
    print("number of miRs: " + str(len(miRs_set)))
    print("number of orig miRs: " + str(len(orig_mirs)))
    for miR in miR2gene2trans2bs2score:
        #mir_arr = miR.split("-")
        #for item in mir_arr:
        #    print(item, end="\t")
        #print()
        for gene in miR2gene2trans2bs2score[miR]:
            transcripts = miR2gene2trans2bs2score[miR][gene].keys()
            curr_transcript = nprand.choice(list(transcripts))
            miR2gene2trans2bs2score[miR][gene] = miR2gene2trans2bs2score[miR][gene][curr_transcript]
    [timeList, timeElapsed, out_line] = methods.timer(timeList, 0, 0)
    methods.logINFO("uploading score table: " + out_line)
    return [miR2gene2trans2bs2score, targets.keys()]

def miRtargetScore_mirtarbase(table_filename, detailed, source):
    timeList = [time.time()]
    #{miR: {gene : transcript: {bs: score}}}}
    miR2gene2trans2bs2score = {}
    miRs = {}
    targets = {}
    col2idx = {}
    bs2id = {}
    counterID = 1
    #f_table = open(methods.pwd() + "expression_data/" + table_filename, 'r')
    with open(methods.pwd() + "expression_data/" + table_filename) as csvfile:
        f_table = csv.reader(csvfile)
        i = 0
        orig_mirs = set()
        miRs_set = set()
        for row in f_table:
            if i%1000000 == 0:
                print("curr i: " + str(i))
            # row = line.strip().split('\t')
            if i == 0:
                for idx in range(len(row)):
                    col2idx[row[idx].strip()] = idx
            else:
                if(source == "microrna.org"):
                    miRname = row[col2idx['mirna_name']].strip()
                    miRacc = row[col2idx['#mirbase_acc']].strip()
                    genename = row[col2idx['gene_symbol']].strip()
                    score = float(row[col2idx['mirsvr_score']].strip())
                    transcript = row[col2idx['ext_transcript_id']].strip()
                    coordinate = row[col2idx['genome_coordinates']].strip()
                elif(source == "TargetScan"):
                    miRname = row[col2idx['miR']].strip()
                    genename = row[col2idx['gene']].strip()
                    score = float(row[col2idx['score']].strip())
                    coordinate = row[col2idx['coordinates']].strip()
                    miRacc = miRname
                    transcript = genename
                elif source == constants.MIR_TAR_BASE:
                    miRname = row[col2idx['miRNA']].strip()
                    genename = row[col2idx['Target Gene']].strip()
                    score = float(row[col2idx['count']].strip())
                    coordinate = row[col2idx['coordinates']].strip()
                    miRacc = miRname
                    transcript = genename
                elif source == constants.RANDOM_SCORES:
                    miRname = row[col2idx['miRNA']].strip()
                    genename = row[col2idx['target']].strip()
                    score = float(row[col2idx['score']].strip())
                    coordinate = row[col2idx['MBS']].strip()
                    miRacc = miRname
                    transcript = genename
                elif(source.startswith("both")):
                    miRname = row[col2idx['miR']].strip()
                    genename = row[col2idx['gene']].strip()
                    score1 = float(row[col2idx['score1']].strip())
                    score2 = float(row[col2idx['score2']].strip())
                    score = (score1,score2)
                    coordinate = row[col2idx['coordinates']].strip()
                    miRacc = miRname
                    transcript = genename
                # [hg19:chr:start-end:-]
                curr_chr = coordinate.split(':')[1]
                strand = coordinate.split(':')[3][0]
                #if('hg19' not in coordinate):
                #    methods.logERROR("hg19 not in " + coordinate)
                if(coordinate not in bs2id):
                    currID = counterID
                    counterID += 1
                    bs2id[coordinate] = currID
                else:
                    currID = bs2id[coordinate]
                coordinate = coordinate.replace('hg19', str(currID))
                if(detailed):
                    genename += '_' + curr_chr + ':' + strand
                # removing 'hsa-' prefix:
                #miRname = miRname[4:].lower()
                if(source == "microrna.org"):
                    if miRname not in miRs:
                        miRs[miRname] = miRacc
                    elif miRs[miRname] != miRacc:
                        print(miRname + ': ' + miRs[miRname] + miRacc)
                orig_mirs.add(miRname)
                miRname = methods.miRnameFormat(miRname)
                miRs_set.add(miRname)
                if genename not in targets:
                    targets[genename] = {transcript: [coordinate]}
                else:
                    if transcript not in targets[genename]:
                        targets[genename][transcript] = [coordinate]
                    elif(coordinate not in targets[genename][transcript]):
                        targets[genename][transcript].append(coordinate)
                if(miRname in miR2gene2trans2bs2score):
                    gene2trans2bs2score = miR2gene2trans2bs2score[miRname]
                else:
                    gene2trans2bs2score = {}
                if(genename in gene2trans2bs2score):
                    trans2bsmap = gene2trans2bs2score[genename]
                else:
                    trans2bsmap = {}
                if(transcript in trans2bsmap):
                    bs_map = trans2bsmap[transcript]
                else:
                    bs_map = {}
                if(coordinate in bs_map):
                    #print("found coordinate: " + coordinate + " miR: " + miRname + " gene: " + genename + " transcript: " + transcript)
                    bs_map[coordinate] = min(bs_map[coordinate], score)
                else:
                    bs_map[coordinate] = score
                trans2bsmap[transcript] = bs_map
                gene2trans2bs2score[genename] = trans2bsmap
                miR2gene2trans2bs2score[miRname] = gene2trans2bs2score
            i += 1
        #f_table.close()
        print("number of miRs: " + str(len(miRs_set)))
        print("number of orig miRs: " + str(len(orig_mirs)))
        for miR in miR2gene2trans2bs2score:
            #mir_arr = miR.split("-")
            #for item in mir_arr:
            #    print(item, end="\t")
            #print()
            for gene in miR2gene2trans2bs2score[miR]:
                transcripts = miR2gene2trans2bs2score[miR][gene].keys()
                curr_transcript = nprand.choice(list(transcripts))
                miR2gene2trans2bs2score[miR][gene] = miR2gene2trans2bs2score[miR][gene][curr_transcript]
        [timeList, timeElapsed, out_line] = methods.timer(timeList, 0, 0)
        methods.logINFO("uploading score table: " + out_line)
    return [miR2gene2trans2bs2score, targets.keys()]

def one_gene_to_many_miRNA_table(table_filename, detailed, source):
    orig_table, targets = miRtargetScore(table_filename, detailed, source)
    new_targets = []
    for miRNA in orig_table:
        rand_target = nprand.choice(list(orig_table[miRNA].keys()))
        orig_table[miRNA] = {rand_target:orig_table[miRNA][rand_target]}
        new_targets.append(rand_target)
    for miRNA in orig_table:
        methods.logINFO(miRNA + "\t"+ str(orig_table[miRNA]))
    return orig_table, new_targets

def one_miRNA_to_many_genes_table(table_filename, detailed, source):
    orig_table, targets = miRtargetScore(table_filename, detailed, source)
    new_targets = set()
    new_table = {}
    for miRNA in orig_table:
        if len(orig_table[miRNA]) == 1:
            new_table[miRNA] = deepcopy(orig_table[miRNA])
            for target in new_table[miRNA]:
                new_targets.add(target)
    for miRNA in orig_table:
        if len(orig_table[miRNA]) > 1:
            for target in orig_table[miRNA]:
                if target not in new_targets:
                    new_targets.add(target)
                    if miRNA in new_table:
                        new_table[miRNA][target] = deepcopy(orig_table[miRNA][target])
                    else:
                        new_table[miRNA] = {target:deepcopy(orig_table[miRNA][target])}
        if len(orig_table[miRNA]) == 0:
            methods.logERROR("no unique target for miRNA: "+ miRNA)
    for miRNA in new_table:
        methods.logINFO(miRNA + "\t"+ str(new_table[miRNA]))
    return new_table, new_targets

def one_to_one_table(table_filename, detailed, source):
    orig_table, targets = miRtargetScore(table_filename, detailed, source)
    new_targets = set()
    new_table = {}
    for miRNA in orig_table:
        if len(orig_table[miRNA]) == 1:
            new_table[miRNA] = deepcopy(orig_table[miRNA])
            for target in new_table[miRNA]:
                new_targets.add(target)
    for miRNA in orig_table:
        if len(orig_table[miRNA]) > 1:
            found = False
            while len(orig_table[miRNA].keys()) > 0 and not found:
                rand_target = nprand.choice(list(orig_table[miRNA].keys()))
                if rand_target in new_targets:
                    del orig_table[miRNA][rand_target]
                else:
                    new_table[miRNA] = {rand_target:orig_table[miRNA][rand_target]}
                    found = True
                    new_targets.add(rand_target)
            if not found:
                print("no unique target for miRNA: ", miRNA)
    for miRNA in new_table:
        methods.logINFO(miRNA + "\t"+ str(new_table[miRNA]))
    return new_table, new_targets

if in_file:
    timeList = [time.time()]
    target2trans2bs2miR = target_mir_score_table(in_file)
    f_out = open(methods.pwd() + out_file, 'w')
    f_targetcount = open(methods.pwd() + 'target_count.tab', 'w')
    f_tarnscriptcount = open(methods.pwd() + 'transcript_count.tab', 'w')
    miRs_sorted = sorted(miRs.keys())
    targets_sorted = sorted(targets.keys())
    target_count = 0
    header_line = 'targets'
    for miR in miRs_sorted:
        header_line += '\t' + miR
    f_out.write(header_line + '\n')
    for target in targets_sorted:
        transcripts = target2trans2bs2miR[target].keys()
        f_tarnscriptcount.write(target + '\t' + str(len(transcripts)) + '\n')
        rand_idx = random.randint(0,(len(transcripts)-1))
        curr_transcript = transcripts[rand_idx]
        bs_lst = targets[target][curr_transcript]
        bs_map = target2trans2bs2miR[target][curr_transcript]
        #bs_map = target2trans2bs2miR[target]
        target_count += 1
        bs_count = 0
        for bs in bs_lst:
            bs_count += 1
            out_line = target + ':' + curr_transcript + ':' + bs + ':[' + ('|').join(gene2coordinate[target]) + ']'
            if bs in bs_map.keys():
                miRscoreMap = bs_map[bs]
                for miR in miRs_sorted:
                    miRinlst = 0
                    for curr_miR in miRscoreMap:
                        miR_score = miRscoreMap[curr_miR]
                        if(miR == curr_miR):
                            if(miRinlst):
                                print('miR ' + miR +' appeared more than once ' + target + ' bs: ' + bs)
                                for miR in miRscoreMap:
                                    print(miR)
                                exit(8)
                            out_line += '\t' + miR_score
                            miRinlst = 1
                    if(not miRinlst):
                        out_line += '\t 0'
            else:
                for miR in miRs_sorted:
                    out_line += '\t 0'
            f_out.write(out_line + '\n')
        f_targetcount.write(target + '\t ' + str(bs_count) + '\n')
    print('number of targets:' + str(target_count))
    print('number of miRs:' + str(len(miRs)))
    f_targetcount.close()
    f_out.close()
    f_tarnscriptcount.close()
    
