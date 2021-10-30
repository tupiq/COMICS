#!/usr/bin/python
# coding: utf-8
import methods
import objects
import sys
import getopt
import math
from numpy import random as nprand
import logging
import src.constants as constants


# Setting here all relevant parameter
# we assume that there is about 50000 miRs molecules in the cell
# so there is 0.5*50000 = 6000 mRNA in the cell
total_miR = 50000
mRNA2miR_ratio = 0.5
update_quanta = 1  # The number of mRNA participate in each iteration
num_of_iter = 100000  # Number of total iterations
remove_occupied_mRNA_interval = 1000  # Number of iterations, after which an occupied mRNA should be removed
over_exp_miR_factor = 1  # The factor by which the over/under expressed miRNA should raise
over_exp_miRNA = ''  # The name of the miRNA that is over/under expression
report_threshold = 5  # The threshold of number of mRNA molecules that should be reported at the end of simulation

# Input and output files and data:
cell_type = 'HeLa'
cell_type_gene = "HeLa"
# cell_type = 'mmu_ES'
# cell_type_gene = "mmu_ES"
# Source of score table:
# source = 'TargetScan'
# source = 'microrna.org'
# score_table = 'human_predictions_S_C_aug2010.txt'  # microrna.org
# in case we are using microrna.org:
exp_fold_change = -1  # The fold change of over expression, lower number will result lower probabilities
# score_table = 'targetScan_scoretable.tab'  # TargetScan
# score_table = 'score_table_mirnaorg_targetscan.tab'  # intersection of microran.org and TargetScan
# source = constants.TARGET_SCAN
#score_table = constants.TARGET_SCAN_TABLE
# score_table = constants.MMU_TARGET_SCAN_TABLE
source = constants.MIR_TAR_BASE
score_table = constants.MIR_TAR_BASE_TABLE
discrete_probabilities = False
compare_to_clash = False
clash_file = 'clashdataS1.txt'
# miRNA and mRNA expression files:
if cell_type == "HeLa":
    # miR_exp_file = 'miRnamexpGSM1299957.tab'
    miR_exp_file = 'our_miRNA/HeLa_no_mir/0h/miRNA_mature_exp.tab'# 'our_miRNA/miR_exp_total_no_mir.tab'
    # miR_exp_file = 'our_miRNA/alpha_amanitin/HeLa_no_mir/0h/miRNA_mature_exp.tab'
    # gene_exp_file = 'genenameControlGSM546991.tab'
if cell_type_gene == "HeLa":
    gene_exp_file = 'HeLa/our_mRNA/HeLa_no_mir/gene_exp_total.tab'
    # gene_exp_file = 'our_mRNA/alpha_amanitin/HeLa_no_mir/gene_exp_total.tab'
if cell_type == "HEK":
    miR_exp_file = 'HEK_no_mir/0h/miRNA_mature_exp.tab'# 'our_miRNA/miR_exp_total_no_mir.tab'
if cell_type_gene == "HEK":
    gene_exp_file = 'HEK_no_mir/gene_exp_total.tab'
    compare_to_clash = True
if cell_type == "MCF7":
    miR_exp_file = "mcf7_miRNA_expression.tab"
if cell_type_gene == "MCF7":
    gene_exp_file = "GSE103520_mcf7_geneexp.tab"#'mcf7_nci60_ave.tab'
if cell_type_gene == "mmu_ES":
    gene_exp_file = "mouse_ES_cells/ES_gene_expression.txt"
if cell_type == "mmu_ES":
    miR_exp_file = "mouse_ES_cells/ES_CCE_mir_expression.txt"

network_type = methods.MANY_TO_MANY
half_life_file = "C:/miR/expression_data/HeLa/gene_half_life_time.tab"
# Output:
out_dir = cell_type + '/' #'/our_data/'  # The output directory, should be located under miR/output/cell_type/
# out_dir = 'HeLa/GSM1299957/allmiRs/control/'
# out_file = "not_removing_mRNA_ts"
out_file = over_exp_miRNA + '_' + str(over_exp_miR_factor) + '_th' + str(report_threshold) + source + \
           "mir_" + cell_type + "_gene_" + cell_type_gene + "_shuffled"
# output files:
# 1. log file: methods.pwd() + 'output/' + out_dir + out_file + '.log'
# 2. miRNA expression at time 0:
# methods.pwd() + "output/" + out_dir + over_exp_miRNA + '_' + str(over_exp_miR_factor) + '.tab'

# More features:
remove_free_mRNA = False
iter_add_mRNA = -1  # The number of iterations, after which a random new mRNA should be added
iter_remove_free_mRNA = -1  # The number of iterations, after which a free mRNA should be removed
lam = 1000  # The value of lambda for the poisson probability of removing mRNA

# Initializing parameters:
use_CLASH = False  # percentage of interaction using CLASH data, only in case of HEK cells
random_init = 0  # The percentage of random initiation of interactions

#mRNAs2plot = set()
debug = False
mir_simulation = True
random_pairing = False

def main(argv):
    # Getting the parameters from the user:
    try:
        opts, args = getopt.getopt(argv,"o:g:m:t:u:f:r:z:q:n:k:curr_iter:c:",["out=", "gene_exp_file=", "miR_exp_file=", "th="])
        command_help = 'main.py -o <out_file,out_dir> -g <gene_exp_file> -m <miR_exp_file> -t <total_miR> ' \
                       '-u <mRNA2miR_ratio> -f <over_exp_miR_factor -r <over_exp_miRNA> -z <report_threshold>' \
                       '-q <update_quanta> -n <num_of_iter> -k <remove_occupied_mRNA_interval>' \
                       '-curr_iter <random_init> -c <use_CLASH>'
    except getopt.GetoptError:
        print(command_help)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(command_help)
            sys.exit()
        elif opt in ("-o", "--out"):
            global out_file
            global out_dir
            out = arg.split(',')
            out_file = out[0]
            out_dir = out[1]
        elif opt in ("-g", "--gene_exp_file"):
            global gene_exp_file
            gene_exp_file = arg
        elif opt in ("-m", "--miR_exp_file"):
            global miR_exp_file
            miR_exp_file = arg
        elif opt in ("-t", "--total_miR"):
            global total_miR
            total_miR = int(arg)
        elif opt in ("-u", "--mRNA2miR_ratio"):
            global mRNA2miR_ratio
            mRNA2miR_ratio = float(arg)
        elif opt in ("-f", "--over_exp_miR_factor"):
            global over_exp_miR_factor
            over_exp_miR_factor = int(arg)
        elif opt in ("-r", "--over_exp_miRNA"):
            global over_exp_miRNA
            over_exp_miRNA = arg
        elif opt in ("-d", "--directory"):
            global directory
            directory = arg
        elif opt in ("-z", "--report_threshold"):
            global report_threshold
            report_threshold = int(arg)
        elif opt in ("-q", "--update_quanta"):
            global update_quanta
            update_quanta = int(arg)
        elif opt in ("-n", "--num_of_iter"):
            global num_of_iter
            num_of_iter = int(arg)
        elif opt in ("-curr_iter", "--random_init"):
            global random_init
            random_init = float(arg)
        elif opt in ("-c", "--use_CLASH"):
            global use_CLASH
            use_CLASH = int(arg)
        elif opt in ("-k", "--remove_occupied_mRNA_interval"):
            global remove_occupied_mRNA_interval
            remove_occupied_mRNA_interval = int(arg)

    try:
        print("starting...")
        log_filename = methods.pwd() + 'output/' + out_dir + out_file + '.log'
        print(log_filename)
        methods.initLogFile(log_filename)
        methods.logINFO("=================================================")
        methods.logINFO("out_file: " + out_file + "over/under expressed miRNA: " + over_exp_miRNA +
                        " over_exp_miR_factor: " + str(over_exp_miR_factor))
        # setting variables for tracing total miRNA and mRNA counts
        total_miRs_map = {}
        total_mRNAs_map = {}
        total_labels = []
        record_idx_gen = generate_record_idx(num_of_iter)  # index of recording cell state (default recording after 1000 iterations)
        time_list = methods.timer([], 0, 0)[0]
        # Loading the data:
        print("loading..")
        [miR_expression, mRNA2exp, mRNAdist, miRtargetTable, table_genes, half_life] = \
             methods.load_data(cell_type, miR_exp_file, gene_exp_file, score_table,
                               source=source, half_life_file=half_life_file, cell_type_gene=cell_type_gene)
        print("done loading")
        cell_state = objects.CellState(miR_expression, mRNAdist, total_miR, total_miR * mRNA2miR_ratio,
                                       remove_occupied_mRNA_interval, update_quanta, exp_fold_change,
                                       miRtargetTable, source,
                                       half_life=half_life,
                                       discrete_probabilities=discrete_probabilities)
        methods.logINFO("miRNAinTable: " + str(len(miRtargetTable)))
        methods.logINFO("mRNAinTable: " + str(len(table_genes)))
        [in_table_miR, zero_count_miR, not_in_table_miR] = cell_state.remove_not_in_table('miR', miRtargetTable.keys())
        [in_table_mRNA, zero_count_mRNA, not_in_table_mRNA] = cell_state.remove_not_in_table('mRNA', table_genes)
        cell_state.change_dist_to_abs_count(table=miRtargetTable)
        if cell_type_gene == "HEK":
            clash, clash_in_table, clash_in_exp = \
                methods.processCLASHdata(clash_file, table=miRtargetTable, mir_exp=cell_state.miR2dist, mRNA_exp=cell_state.mRNA2dist)
        methods.logINFO("total_mRNA: " + str(cell_state.totalmRNA) + " number of mRNAs: " + str(len(in_table_mRNA)) +
                        " zero count mRNA: " + str(zero_count_mRNA) + " mRNA not in table: " + str(len(not_in_table_mRNA)))
        methods.logINFO("total_miR: " + str(cell_state.totalmiR) + " number of miRs: " + str(len(in_table_miR)) +
                        " zero count miR: " + str(zero_count_miR) + " miRNAs not in table: " + str(len(not_in_table_miR)))
        if remove_free_mRNA:
            for mRNA in cell_state.mRNA2dist:
                if mRNA not in half_life:
                    methods.logINFO(mRNA + " not in half life")
        #for mir in not_in_table_miR:
        #    print(mir)

        if over_exp_miR_factor != 1:  # change miRNAs counts according to the over expressed miRNA
            over_exp_miRNAs = cell_state.change_miR_counts(over_exp_miR_factor, over_exp_miRNA, miRtargetTable)
            over_exp_miRNA = "_".join(over_exp_miRNAs)
        cell_state.create_pairs_probabilities(miRtargetTable, source)
        # Reporting the initial count of miRNA:
        f_miRs = open(methods.pwd() + "output/" + out_dir + "mir_exp_" + over_exp_miRNA + '_' + str(over_exp_miR_factor) + '.tab', 'w')
        f_miRs.write("miR \t count \n")
        for mirt in cell_state.miR2dist:
            f_miRs.write(mirt + " \t " + str(cell_state.miR2dist[mirt]) + '\n')
        f_miRs.close()
        record_idx = next(record_idx_gen)
        total_miRs_map[record_idx] = cell_state.totalmiR
        total_mRNAs_map[record_idx] = cell_state.totalmRNA
        cell_state.recordmRNAcount(record_idx)
        total_labels.append("0")
        pass_quanta = 0
        # Report how many miRNA and mRNA are relevant:
        for mRNA in cell_state.mRNA2dist:
            if cell_state.mRNA2dist[mRNA]*cell_state.totalmRNA >= update_quanta:
                pass_quanta += 1
        methods.logINFO("mRNA pass quanta: " + str(pass_quanta))
        pass_quanta = 0
        for miR in cell_state.miR2dist:
            if cell_state.miR2dist[miR]*cell_state.totalmiR >= update_quanta:
                pass_quanta += 1
        methods.logINFO("miR pass quanta: " + str(pass_quanta))
        # updating CLASH data:
        curr_iter = 0
        if use_CLASH:
            record_idx = next(record_idx_gen)
            curr_iter = use_clash_data(mRNA2exp,clash, cell_state, miRtargetTable, table_genes)
            total_miRs_map[record_idx] = cell_state.totalmiR
            total_mRNAs_map[record_idx] = cell_state.totalmRNA
            cell_state.recordmRNAcount(record_idx)
            total_labels.append("CLASHinit_" + str(curr_iter))
            [time_list, timeElapsed, out_line] = methods.timer(time_list, 1,0)
            methods.logINFO(out_line)
        # Check if initiation is needed:
        updated_count = 0
        miRnotinTable = set()
        mRNAnotinTable = set()
        no_matching_mRNA = 0
        # Check if we should make some random interactions before the simulation
        if random_init > 0 and not use_CLASH:
            num_miR_init = int(math.floor(random_init * cell_state.totalmiR))
            for curr_iter in range(num_miR_init):
                rand_miR = cell_state.get_random('miR')
                if rand_miR in miRtargetTable:
                    [rand_mRNA, bs, score,  BSfree, currgeneNotFound] = cell_state.get_random_mRNA(rand_miR, miRtargetTable)
                    if rand_mRNA and rand_miR and bs:
                        cell_state.force_update(rand_mRNA, rand_miR, bs, BSfree, update_quanta, 0)
                        updated_count += 1
                    elif rand_mRNA:
                        gene = rand_mRNA.split('_')[0]
                        methods.logINFO(rand_mRNA + ' : ' + rand_miR )
                        bs2score = methods.BSscore(gene, rand_miR, miRtargetTable)
                    else:
                        no_matching_mRNA += 1
            record_idx = next(record_idx_gen)
            total_miRs_map[record_idx] = cell_state.totalmiR
            total_mRNAs_map[record_idx] = cell_state.totalmRNA
            cell_state.recordmRNAcount(record_idx)
            total_labels.append("randInit_" + str(curr_iter))
            methods.logINFO("updated_count: " + str(updated_count))
            methods.logINFO("no_matching_mRNA: " + str(no_matching_mRNA))
            [time_list, timeElapsed, out_line] = methods.timer(time_list, 1,0)
            methods.logINFO(out_line)
        # Starting the simulation
        updated_count = 0
        no_matching_mRNA = 0
        clash_pairs_cnt = {}
        pair_cnt_bu = {}
        pair_cnt = {}
        curr_iter = 0
        tl = []
        iter_time_lst = []
        iter_time_lst = methods.timer(iter_time_lst,0,1)[0]
        methods_time_count = {"cell_state.getRandommRNA":[0], "cell_state.update":[0], "cell_state.removeOccupiedmRNA":[0]}
        methods_idx = 0
        totalmiRremoved = 0
        totalmRNAremoved = 0
        genesNotFound = set()
        for curr_iter in range(num_of_iter):
            iter_timer = methods.timer([],0,0)[0]
            # Adding random mRNA (if this feature is on)
            if curr_iter > 0 and iter_add_mRNA > 0 and curr_iter%iter_add_mRNA == 0:
                cell_state.addRandommRNA(update_quanta)
            # Remove random free mRNA
            if curr_iter > 0 and iter_remove_free_mRNA > 0 and curr_iter%iter_remove_free_mRNA == 0:
                cell_state.removeRandomFreemRNA(update_quanta)
            if mir_simulation:
                if random_pairing:
                    cell_state.random_interaction(update_quanta, curr_iter)
                else:
                    rand_miR = cell_state.get_random('miR')
                    if rand_miR in miRtargetTable:
                        tl = methods.timer(tl,0,0)[0]
                        [rand_mRNA, bs, score,  BSfree, currgeneNotFound] = cell_state.get_random_mRNA(rand_miR, miRtargetTable)
                        genesNotFound.union(currgeneNotFound)
                        [tl, methodTime1, out_line] = methods.timer(tl,0,1)
                        methods_time_count["cell_state.getRandommRNA"][methods_idx] += methodTime1
                        if rand_mRNA and rand_miR:
                            if bs and score:
                                tl = methods.timer(tl,0,0)[0]
                                updated = cell_state.update(rand_miR, rand_mRNA, {bs:score}, BSfree, update_quanta, curr_iter, source)
                                [tl, methodTime2, out_line] = methods.timer(tl,0,1)
                                methods_time_count["cell_state.update"][methods_idx] += methodTime2
                            else:
                                gene = rand_mRNA.split('_')[0]
                                if debug:
                                    methods.logINFO(rand_mRNA + ' : ' + rand_miR)
                                bs2score = methods.BSscore(gene, rand_miR, miRtargetTable)
                                updated = cell_state.update(rand_miR, rand_mRNA, bs2score, 0, update_quanta, curr_iter, source)
                            pair = rand_miR + '_' + rand_mRNA
                            if updated:
                                updated_count += 1
                                if pair in pair_cnt:
                                    pair_cnt[pair] += 1
                                else:
                                    pair_cnt[pair] = 1
                            if compare_to_clash and rand_miR + '_' + rand_mRNA in clash_in_table:
                                if pair in clash_pairs_cnt:
                                    clash_pairs_cnt[pair] += 1
                                else:
                                    clash_pairs_cnt[pair] = 1
                            if pair in pair_cnt_bu:
                                pair_cnt_bu[pair] += 1
                            else:
                                pair_cnt_bu[pair] = 1
                        else:
                            if not rand_mRNA:
                                no_matching_mRNA += 1
                    else:
                        miRnotinTable.add(rand_miR)
                tl = methods.timer(tl,0,0)[0]
                [currmiRremoved, currmRNAremoved] = cell_state.removeOccupiedmRNA(curr_iter)
                totalmiRremoved += currmiRremoved
                totalmRNAremoved += currmRNAremoved
                [tl, methodTime3, out_line] = methods.timer(tl,0,1)
                methods_time_count["cell_state.removeOccupiedmRNA"][methods_idx] += methodTime3
                [iter_timer, iterTimeElapsed, out_line] = methods.timer(iter_timer,0,1)
                if iterTimeElapsed > 30:
                    methods.logERROR("Slow running time: " + str(iterTimeElapsed))
                    methods.logERROR("cell_state.getRandommRNA: " + str(methodTime1))
                    methods.logERROR("cell_state.update: " + str(methodTime2))
                    methods.logERROR("cell_state.removeOccupiedmRNA: " + str(methodTime3))
                    methods.logERROR("exiting..")
                # exit()
            if cell_state.totalmRNA <= 1000 or cell_state.totalmiR <= 100:
                methods.logINFO("exiting:")
                methods.logINFO("current iteration: " + str(curr_iter))
                methods.logINFO("total mRNA: " + str(cell_state.totalmRNA))
                methods.logINFO("total miR: " + str(cell_state.totalmiR))
                break
            if curr_iter%1000 == 0 and curr_iter > 0:
                record_idx = next(record_idx_gen)
                methods.logINFO("current iteration: " + str(curr_iter))
                methods.logINFO("total mRNA: " + str(cell_state.totalmRNA))
                methods.logINFO("total miR: " + str(cell_state.totalmiR))
                methods.logINFO("totalmiRremoved: " + str(totalmiRremoved))
                methods.logINFO("totalmRNAremoved: " + str(totalmRNAremoved))
                methods.logINFO("total pairs: " + str(len(pair_cnt)))
                methods.logINFO("total pairs before update: " + str(len(pair_cnt_bu)))
                total_miRs_map[record_idx] = cell_state.totalmiR
                total_mRNAs_map[record_idx] = cell_state.totalmRNA
                cell_state.recordmRNAcount(record_idx)
                total_labels.append(curr_iter)
                for method in methods_time_count:
                    runningTime = methods_time_count[method][methods_idx]
                    if debug:
                        methods.logINFO(method + " : " + str(runningTime/1000))
                    if runningTime > 0.02*1000:
                        methods.plotTimeStats(methods_time_count, 1000)
                        methods.logERROR("Slow running time " + str(runningTime/1000) )
                        #exit()
                    methods_time_count[method].append(0)
                methods_idx += 1
                if remove_free_mRNA:
                    cell_state.remove_mRNA()
                if compare_to_clash:
                    pearson, spearman , pearson_before_update, spearman_before_update = \
                        cell_state.compare_to_CLASH(clash_in_table, clash_pairs_cnt)
                    methods.logINFO("Correlation with CLASH, Pearson: " + str(pearson) + ", current iteration: " + str(curr_iter))
                    methods.logINFO("Correlation with CLASH, Spearman: " + str(spearman) + ", current iteration: " + str(curr_iter))
                    methods.logINFO("Correlation with CLASH before update, Pearson: " + str(pearson_before_update) + ", current iteration: " + str(curr_iter))
                    methods.logINFO("Correlation with CLASH before update, Spearman: " + str(spearman_before_update) + ", current iteration: " + str(curr_iter))
        methods.logINFO("updated_count: " + str(updated_count))
        methods.logINFO("miRnotintable: " + str(len(miRnotinTable)))
        methods.logINFO("mRNAnotintable: " + str(len(mRNAnotinTable)))
        methods.logINFO("no_matching_mRNA: " + str(no_matching_mRNA))
        methods.logINFO("total mRNA: " + str(cell_state.totalmRNA))
        methods.logINFO("total miR: " + str(cell_state.totalmiR))
        methods.logINFO("genes not found in the exp file: " + str(len(genesNotFound)))
        methods.logINFO(str(set))
        if over_exp_miRNA:
            title = over_exp_miRNA + '_' + str(over_exp_miR_factor)
        else:
            title = out_file
        methods.logINFO("output title: " + out_dir + title)
        sorted_genes = methods.sort_by_retention(cell_state.mRNA_counts, report_threshold)
        cell_state.outputCellState('output/' + out_dir, title, sorted_genes)
        [time_list, timeElapsed, out_line] = methods.timer(time_list, 0,0)
        methods.logINFO(out_line)
    except:
        logging.exception('Got exception on main handler')
        print("FAILED:\n" + "over expressed miRNA: " + over_exp_miRNA + " by factor " + str(over_exp_miR_factor) +
              "\nupdate_quanta: " + str(update_quanta) + " removeOccupiedmRNAinterval: " + str(remove_occupied_mRNA_interval))


def generate_record_idx(max_idx):
    idx = 0
    while idx <= max_idx:
        yield idx
        idx += 1


def use_clash_data(mRNA2exp, clash, cell_state, miRtargetTable, tableGenes):
    i = 0
    CLASHmRNAnotFound = {}
    CLASHmiRnotFound = {}
    CLASHmiRfound = set()
    CLASHmRNAfound = set()
    CLASHmiRtotal = set()
    CLASHmRNAtotal = set()
    f_outCLASH = open(methods.pwd() + 'expression_data/pairNotinCLASH.tab', 'w')
    f_outCLASH.write("miR \t mRNA \t miRexistinTable \t mRNAexistinTable \t BS \t count \n")
    for miR in clash:
        mRNA2bs = clash[miR]
        miR = methods.miRnameFormat(miR)
        CLASHmiRtotal.add(miR)
        for mRNA in mRNA2bs:
            CLASHmRNAtotal.add(mRNA)
            if(i%1000 == 0):
                methods.logINFO("i: " + str(i) )
            i += 1
            if(mRNA in mRNA2exp.keys()):
                pos_lst = mRNA2exp[mRNA]
                pos = nprand.choice(list(pos_lst.keys()))
                curr_chr = pos.split(':')[0]
                curr_strand = pos.split(':')[1]
            else:
                curr_chr = '*'
                curr_strand = '*'
            for bs in mRNA2bs[mRNA]:
                quantity = math.floor(mRNA2bs[mRNA][bs] * total_miR * random_init)
                if quantity > 0:
                    cell_state.force_update(mRNA + '_' + bs.chr + ':' + bs.strand, miR, bs, -1, quantity, 0)
                # This is only for checking if the CLASH data appears in the prediction table
                bs2score = methods.BSscore(mRNA, miR, miRtargetTable)
                miR_table = 0
                mRNA_table = 0
                if(miR in miRtargetTable):
                    CLASHmiRfound.add(miR)
                    miR_table = 1
                else:
                    if(miR in CLASHmiRnotFound):
                        CLASHmiRnotFound[miR] += quantity
                    else:
                        CLASHmiRnotFound[miR] = quantity
                if mRNA not in tableGenes:
                    if mRNA in CLASHmRNAnotFound:
                        CLASHmRNAnotFound[mRNA] += quantity
                    else:
                        CLASHmRNAnotFound[mRNA] = quantity
                if len(bs2score) == 0:
                    if not(miR_table and mRNA_table):
                        f_outCLASH.write(miR + '\t' + mRNA + '\t' + str(miR_table) + '\t' + str(mRNA_table) + '\t' + bs.toString() + str(quantity) + "\n")
    f_outCLASH.close()
    methods.logINFO("CLASHmiRfound: " + str(len(CLASHmiRfound)) )
    methods.logINFO("CLASHmiRnotFound: " + str(len(CLASHmiRnotFound)) )
    methods.logINFO("CLASHmRNAnotFound: " + str(len(CLASHmRNAnotFound)) )
    methods.logINFO("CLASHmiRtotal: " + str(len(CLASHmiRtotal)) )
    methods.logINFO("CLASHmRNAtotal: " + str(len(CLASHmRNAtotal)) )
    methods.writeMap(CLASHmiRnotFound, 'expression_data/CLASHmiRnotFoundinTable.tab')
    methods.writeMap(CLASHmRNAnotFound, 'expression_data/CLASHmRNAnotFoundinTable.tab')
    methods.logINFO("i: " + str(i))
    methods.logINFO("total miR: " + str(cell_state.totalmiR) + " total mRNA: " + str(cell_state.totalmRNA))
    return i

if __name__ == "__main__":
    main(sys.argv[1:])