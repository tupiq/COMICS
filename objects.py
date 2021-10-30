from numpy import random as nprand
import math
import methods
import operator
import numpy as np
import traceback


class CellState:
        
    def __init__(self, miRexpression, mRNAdist, totalmiR, totalmRNA, del_interval, update_quanta,
                 exp_fold_change, score_table, source, half_life=None, discrete_probabilities=False,
                 lam=0):
        self.miR2dist = miRexpression
        self.mRNA2dist = mRNAdist
        self.totalmiR = totalmiR
        self.totalmRNA = totalmRNA
        self.del_interval = del_interval
        self.update_quanta = update_quanta # how many interaction should be upon one interaction event
        self.exp_fold_change = exp_fold_change
        self.lam = lam
        self.updated_miRs = {}  # counts the number of updates per miRNA
        self.updated_mRNAs = {}  # counts the number of updates per mRNA
        self.mRNA_occupied = {}  # map of mRNA name to mRNA object, of occupied mRNAs
        self.mRNA_counts = {}
        self.mRNA_counts_idx = []
        self.mRNAoccupBScount = {}  # mRNA to list of tow counts:
                                    # number of occupied, and number of newly bounds (to free mRNAs)
        self.half_life = half_life
        self.dead_by_half_life = 0
        self.discrete_probabilities = discrete_probabilities
        self.pairs = {}
        self.table_pairs = None
        self.probabilities = None
        self.score_table = score_table
        self.source = source

    def remove_not_in_table(self, miRormRNA, table_genes):
        """
        Removing genes that are not in the miRNA-mRNA interaction table,
        and those that have total counter that is lower than 1
        :param miRormRNA: type of gene - mRNA or miRNA
        :param table_genes: the list of relevant genes that are in the score table
        :return: in_table_dist, zero_count, not_in_table
        """
        if miRormRNA == 'mRNA':
            exp_dist = self.mRNA2dist
            total = self.totalmRNA
        elif miRormRNA == 'miR':
            exp_dist = self.miR2dist
            total = self.totalmiR
            table_pre_mir = {}
            for mir in table_genes:
                if mir.endswith("p"):
                    pre_mir = mir[:-3]
                    if pre_mir in table_pre_mir:
                        table_pre_mir[pre_mir].add(mir)
                    else:
                        table_pre_mir[pre_mir] = {mir}
            methods.logINFO("len miRDist: " + str(len(exp_dist)) + " total: " + str(total))
        in_table_dist = {}
        not_in_table = {}
        total_prob = 0
        not_in_table_total_prob = 0
        for gene in exp_dist:
            pre_mir = gene
            if miRormRNA == "miR" and gene.endswith("p"):
                pre_mir = gene[:-3]
            if gene in table_genes:
                in_table_dist[gene] = exp_dist[gene]
                total_prob += exp_dist[gene]
            elif pre_mir in table_genes:
                if pre_mir in in_table_dist:
                    in_table_dist[pre_mir] += exp_dist[gene]
                    print(pre_mir, "already found")
                else:
                    in_table_dist[pre_mir] = exp_dist[gene]
                total_prob += exp_dist[gene]
            elif miRormRNA == "miR" and gene in table_pre_mir:
                print(gene, "\t", table_pre_mir[gene], "pre mir in expression, mature in table")
            elif miRormRNA == "miR" and pre_mir in table_pre_mir:
                print(gene, "\t", pre_mir, table_pre_mir[pre_mir], "\tpre mir in expression, mature in table")
            else:
                not_in_table[gene.split('_')[0]] = exp_dist[gene]
                not_in_table_total_prob += exp_dist[gene]
        methods.logINFO(miRormRNA + ": not_in_table_total_prob: " + str(not_in_table_total_prob))
        if miRormRNA == "miR":
            print("-------miRNA not in table---------------")
            for mir in not_in_table:
                print(mir)
            print("----------------------------------------")
        zero_count = 0  # count the number of genes that have counter lower than1
        total_prob_corrected = 0
        genes2delete = set()
        for gene in in_table_dist.keys():
            # Removing genes that their counter is lower than 1:
            if in_table_dist[gene]*total/total_prob < 1:
                genes2delete.add(gene)
                zero_count += 1
            else:
                in_table_dist[gene] /= total_prob
                total_prob_corrected += in_table_dist[gene]
        for gene in genes2delete:
            del in_table_dist[gene]
        for gene in in_table_dist:
            in_table_dist[gene] /= total_prob_corrected
        total = math.floor(total*total_prob_corrected)
        methods.logINFO(miRormRNA + ": number of genes below threshold: " + str(zero_count))
        methods.logINFO(miRormRNA + ": number of passed threshold: " + str(len(in_table_dist)))
        # updating the current object:
        if miRormRNA == 'mRNA':
            self.mRNA2dist = in_table_dist
            self.totalmRNA = total
        elif miRormRNA == 'miR':
            self.miR2dist = in_table_dist
            self.totalmiR = total
        return [in_table_dist, zero_count, not_in_table]
    
    def change_dist_to_abs_count(self, table=None):
        """
        Changes the miRNA and mRNA distribution to absolute values
        :return:
        """
        [self.mRNA2dist, self.totalmRNA] = methods.dist2absolute(self.mRNA2dist, self.totalmRNA, 1)
        [self.miR2dist, self.totalmiR] = methods.dist2absolute(self.miR2dist, self.totalmiR, 1)
        if table is not None:
            pairs_cnt = 0
            for mir in self.miR2dist:
                for gene in table[mir]:
                    if gene in self.mRNA2dist:
                        pairs_cnt += 1
            methods.logINFO("total table pairs expressed: " + str(pairs_cnt))

    def update(self, miRNA, mRNA, bs2score, BS_free_in_occupied, quantity, iter_num, source):
        """
        Updates the interaction of the given miRNA and mRNA, if it occurs according to the binding probability
        :param miRNA:
        :param mRNA:
        :param bs2score:
        :param BS_free_in_occupied: number of free BS in occupied mRNA
        :param quantity:
        :param iter_num:
        :param source:
        :return: True if there was a change and False otherwise
        """
        if not BS_free_in_occupied:
            currBS = nprand.choice(list(bs2score.keys()))
        else:  # there are free BS in the occupied mRNA
            currBS = list(bs2score.keys())[0]
        binding_prob = methods.score2prob(bs2score[currBS], source, self.exp_fold_change, self.discrete_probabilities)
        try:
            should_update = nprand.choice([0,1], p=[1-binding_prob, binding_prob])
            if should_update:
                # need to update
                self.force_update(mRNA, miRNA, currBS, BS_free_in_occupied, quantity, iter_num)
                # update the number of interaction per the miRNA and the mRNA:
                if miRNA in self.updated_miRs:
                    self.updated_miRs[miRNA] += 1
                else:
                    self.updated_miRs[miRNA] = 1
                if mRNA in self.updated_mRNAs:
                    self.updated_mRNAs[mRNA] += 1
                else:
                    self.updated_mRNAs[mRNA] = 1
                return True
        except:
            methods.logERROR("Could not force update " + miRNA + " and " + mRNA + ", binding_prob: " + str(binding_prob))
            return False
        return False
    
    def force_update(self, mRNA, miRNA, bs, BS_free_in_occupied, quantity, iter_num):
        """
        Create interaction between the given miRNA and mRNA.
        There is preference for updating occupied mRNA, but it is because occupied mRNA was randomly chosen.
        :param mRNA:
        :param miRNA:
        :param bs:
        :param BS_free_in_occupied: number of free BS in occupied mRNA
        :param quantity:
        :param iter_num:
        :return: None
        """
        if miRNA not in self.pairs:
            self.pairs[miRNA] = {}
        if mRNA not in self.pairs[miRNA]:
            self.pairs[miRNA][mRNA] = 0
        self.pairs[miRNA][mRNA] += 1
        BSobj = None
        if BS_free_in_occupied%1 !=0:
            print("BS_free_in_occupied is not an integer!")
        occupied_quanta = 0
        free_quanta = 0
        if miRNA in self.miR2dist:
            freemiRs = int(self.miR2dist[miRNA])  # get the number of available miRNA
        else:
            methods.logERROR("miRNA: " + miRNA + " not in miR2dist. mRNA: " + mRNA)
            freemiRs = quantity  # initiate it to the given quantity if it is not exist
        if mRNA in self.mRNA2dist:
            freemRNAs = int(self.mRNA2dist[mRNA])
        else:
            methods.logERROR("mRNA: " + mRNA + " not in mRNA2dist. miRNA: " + miRNA)
            freemRNAs = quantity
        if freemiRs > 0:
            quantity = min(quantity, freemiRs)  # the  minimum of the update quanta and available miRNAs
            mRNAobj = ''
            if mRNA in self.mRNA_occupied:
                mRNAobj = self.mRNA_occupied[mRNA]
            if BS_free_in_occupied < 0 and mRNAobj:
                # this means we need to check if there are occupied mRNAs with free BS
                if(isinstance(bs, BS)):
                    BSobj = bs
                else:
                    BSobj = BS(mRNA.split('_')[0], bs, '[hg19:chr:start-end:-]')
                bs = str(BSobj)
                BS_free_in_occupied = mRNAobj.isBSavailable(BSobj)
            if BS_free_in_occupied >= quantity:
                occupied_quanta = quantity
            elif BS_free_in_occupied > 0:
                occupied_quanta = BS_free_in_occupied
                free_quanta = min((quantity - occupied_quanta), freemRNAs)
            else:
                free_quanta = min(quantity, freemRNAs)
            quantity = occupied_quanta + free_quanta
            if not mRNAobj:
                mRNAobj = mRNAObj(mRNA, self.del_interval)
            if free_quanta > 0:
                if mRNA not in self.mRNA2dist:
                    self.mRNA2dist[mRNA] = 0
                else:
                    self.mRNA2dist[mRNA] -= free_quanta
                    self.totalmRNA -= free_quanta
                if self.mRNA2dist[mRNA] < 0:
                    methods.logERROR("count value is less than 0 for mRNA: " + mRNA)
                BSobj = mRNAobj.updateBSfree(bs, miRNA, iter_num, free_quanta)  # bind the miRNA to the free mRNA on the relevant BS
            if occupied_quanta > 0:
                BSobj = mRNAobj.updateBSoccupied(bs, miRNA, iter_num, occupied_quanta)  # bind the miRNA to the occupied mRNA on the relevant BS
            if not BSobj:
                methods.logERROR("in force_update: occupied_quanta = " + str(occupied_quanta) + ", free_quanta = " + str(free_quanta))
            elif BSobj.free_bounded + BSobj.occupied_bounded > mRNAobj.occupied:
                methods.logERROR("in forceUpdate: bs.bounded > mRNAobj.occupied, free_quanta: " + str(free_quanta) + ", occupied_quanta: " + str(occupied_quanta))
            self.mRNA_occupied[mRNA] = mRNAobj
            # count BS on each mRNA:
            if mRNA in self.mRNAoccupBScount:
                self.mRNAoccupBScount[mRNA][0] += occupied_quanta + free_quanta
                self.mRNAoccupBScount[mRNA][1] += free_quanta
            else:
                self.mRNAoccupBScount[mRNA] = [occupied_quanta + free_quanta, free_quanta]
            # update the free miR distribution:
            if miRNA not in self.miR2dist:
                self.miR2dist[miRNA] = 0
            else:
                self.miR2dist[miRNA] -= quantity
                self.totalmiR -= quantity
            if self.miR2dist[miRNA] < 0:
                methods.logERROR("count value is less than 0 for miR: " + miRNA)
            if not(quantity%1==0):
                print("quantity for update is not an integer")
            
   
    # return true or false if input BS are overlapping
    @staticmethod
    def BSoverlap(bs1, bs2):
        if bs1.chr == bs2.chr and bs1.strand == bs2.strand:
            if abs(int(bs1.start) - int(bs2.start)) < 50 or abs(int(bs1.end) - int(bs2.end)) < 50:
                return True
        return False

    def create_pairs_probabilities(self, table, source):
        pairs = []
        probabilities = []
        total = 0
        target2idx = {}
        mir2idx = {}
        i = 0
        for mir in self.miR2dist:
            mir2idx[mir] = set()
            mir_count = self.miR2dist[mir]
            for gene in table[mir]:
                if gene not in self.mRNA2dist:
                    continue
                if gene not in target2idx:
                    target2idx[gene] = set()
                gene_count = self.mRNA2dist[gene]
                for bs in table[mir][gene]:
                    target2idx[gene].add(i)
                    mir2idx[mir].add(i)
                    score = table[mir][gene][bs]
                    p = methods.score2prob(score, source)
                    value = p*mir_count*gene_count
                    pairs.append((mir, gene, bs))
                    probabilities.append(value)
                    total += value
                    i += 1
        self.table_pairs = pairs
        self.probabilities = np.array(probabilities)
        self.target2idx = target2idx
        self.mir2idx = mir2idx
        return pairs, probabilities

    def random_interaction(self, quantity, iter_num):
        remove = False
        if len(self.table_pairs) != len(self.probabilities):
            print(len(self.table_pairs),len(self.probabilities), iter_num)
        idx = nprand.choice(list(range(len(self.table_pairs))), p=self.probabilities/sum(self.probabilities))
        if iter_num%1000 == 0:
            print(self.score_table)
            print(self.probabilities)
            print(self.table_pairs)
        if self.probabilities[idx] == 0:
            #print("index is zero",idx, self.table_pairs[idx])
            #print("sum probabilities:", sum(self.probabilities))
            print(self.score_table)
            print(self.probabilities)
            print(self.table_pairs)
            exit()
            return
        mir, mRNA, bs = self.table_pairs[idx]
        BSobj = BS(mRNA.split('_')[0], bs, '[hg19:chr:start-end:-]')
        curr_free = 0
        free_in_occupied = 0
        # see if to get if from occupied or free
        if mRNA in self.mRNA2dist:
            curr_free = self.mRNA2dist[mRNA]  # get number of free mRNA
        # we need to know if the update should be on the occupied or the free mRNA:
        # check if this mRNA is occupied and if it is in full capacity:
        # the occupied BS map:
        if mRNA in self.mRNA_occupied:
            mRNAobj = self.mRNA_occupied[mRNA]
            free_in_occupied = mRNAobj.isBSavailable(BSobj)
        total_gene = curr_free + free_in_occupied
        total_mir = self.miR2dist[mir]
        target_idx = self.target2idx[mRNA]
        mir_idx = self.mir2idx[mir]
        idx_to_remove = set()
        for idx in target_idx:
            if total_gene - 1 == 0:
                idx_to_remove.add(idx)
                self.probabilities[idx] = 0
            elif total_gene > 0:
                self.probabilities[idx] *= (total_gene-1)/total_gene
            else:
                self.probabilities[idx] = 0
        for idx in mir_idx:
            if total_mir -1 ==0:
                idx_to_remove.add(idx)
                self.probabilities[idx] = 0
            elif total_mir > 0 :
                self.probabilities[idx] *= (total_mir-1)/total_mir
            else:
                self.probabilities[idx] = 0
        if remove and len(idx_to_remove) > 0:
            print("length before",len(self.probabilities))
            print("length before",len(self.table_pairs))
            removed_idx = len(idx_to_remove)
            items = map(self.table_pairs.__getitem__, idx_to_remove)
            for i, curr_idx in enumerate(idx_to_remove):
                if curr_idx > removed_idx:
                    curr_idx -= 1*i
                removed_idx = curr_idx
                if self.table_pairs[curr_idx] not in self.table_pairs[curr_idx]:
                    print("error in remove indexing")
                del self.table_pairs[curr_idx]
            self.probabilities = np.delete(self.probabilities, list(idx_to_remove))
            print(len(self.probabilities))
            print(len(self.table_pairs))
            print(idx_to_remove)
            print(iter_num)
        self.force_update(mRNA, mir, bs, free_in_occupied, quantity, iter_num)

    def get_random(self, miRormRNA):
        """
        Returns a random miR/mRNA according to their distribution
        :param miRormRNA:
        :return:
        """
        if miRormRNA == 'miR':
            abs_dist = self.miR2dist
            dist = methods.absolute2dist(self.miR2dist, self.totalmiR)            
        elif miRormRNA == 'mRNA':
            abs_dist = self.mRNA2dist
            dist = methods.absolute2dist(self.mRNA2dist, self.totalmRNA)            
        else:
            methods.logERROR("Error in molecule type, exiting...")
            exit()
        genes = list(dist.keys())
        probabilities = list(dist.values())
        total = sum(dist.values())
        # fixing the probabilities so it will sum to 1
        if 1-total != 0:
            diff = 1-total
            if diff > 0:
                genes.append("")
                probabilities.append(diff)
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
        try:
            rand_mol = nprand.choice(genes, p=probabilities)
            if abs_dist[rand_mol] < 1:
                methods.logERROR("molecule count < 1: " + rand_mol + ", total count: " + str(abs_dist[rand_mol]))
            return rand_mol
        except:
            if not total == 1:
                methods.logERROR("probabilities are not equal to 1 for: " + miRormRNA + " total: " + str(total))
                return ""

    def get_random_mRNA(self, miR, miRtargetTable):
        """
        Chooses one random free mRNA according to their distribution
        :param miR:
        :param miRtargetTable:
        :return: list: [mRNA, bs, score,
                        free_in_occupied (if an occupied mRNA was chosen how many free bs in occupied, else it os 0),
                        genes_not_found]
        """
        mRNAs = {}  # get the number of relevant free mRNAs
        genes_not_found = set()
        sum_mRNA = 0  # sum of relevant free mRNAs
        table_genes = miRtargetTable[miR]
        mRNA2BS2score = {}
        mRNA2bs2free_in_occupied = {}
        for mRNA in table_genes:
            if mRNA in self.mRNA2dist:
                bs2score = miRtargetTable[miR][mRNA]
                curr_free = self.mRNA2dist[mRNA]  # get number of free mRNA
                if curr_free > 0:
                    mRNAs[mRNA] = curr_free
                    sum_mRNA += curr_free
                    mRNA2BS2score[mRNA] = bs2score
                # we need to know if the update should be on the occupied or the free mRNA:
                # check if this mRNA is occupied and if it is in full capacity:
                # the occupied BS map:
                if mRNA in self.mRNA_occupied:
                    mRNAobj = self.mRNA_occupied[mRNA]
                    BScopy = list(bs2score.keys())
                    for origBS in BScopy:
                        BSobj = BS(mRNA.split('_')[0], origBS, '[hg19:chr:start-end:-]')
                        bs = str(BSobj)
                        if bs in mRNAobj.id2BS:  # check if the current BS is occupied
                            BSobj = mRNAobj.id2BS[bs]
                        free_in_occupied = mRNAobj.isBSavailable(BSobj)
                        if free_in_occupied > 0:
                            if mRNA in mRNA2BS2score:
                                mRNA2BS2score[mRNA][bs] = bs2score[origBS]
                            else:
                                mRNA2BS2score[mRNA] = {bs:bs2score[origBS]}
                            if mRNA in mRNA2bs2free_in_occupied:
                                mRNA2bs2free_in_occupied[mRNA][bs] = free_in_occupied
                            else:
                                mRNA2bs2free_in_occupied[mRNA] = {bs:free_in_occupied}
                            mRNAs[mRNA + '|' + bs + '|occupied'] = free_in_occupied
                            sum_mRNA += free_in_occupied
            else:
                genes_not_found.add(mRNA)
        if sum_mRNA > 0:
            for mRNA in mRNAs:
                mRNAs[mRNA] /= sum_mRNA
            [mRNAsfixed, probabilities] = methods.correct_dist(mRNAs)
            randmRNA = nprand.choice(mRNAsfixed, p=probabilities)
            if randmRNA.endswith('occupied'):
                mRNA_arr = randmRNA.split('|')
                mRNA = mRNA_arr[0]
                bs = mRNA_arr[1]
                free_in_occupied = mRNA2bs2free_in_occupied[mRNA][bs]
                return [mRNA, bs, mRNA2BS2score[mRNA][bs], free_in_occupied, genes_not_found]
            bs2score = mRNA2BS2score[randmRNA]
            if bs2score:
                bs = nprand.choice(list(bs2score.keys()))
                return [randmRNA, bs, bs2score[bs], 0, genes_not_found]
        return['',0,0,0,genes_not_found]


    
    def outputCellState(self, outdir, title, sortedGenes, outmRNA=False, outmiR=False, outRetention=True):
        if outmRNA:
            mRNAoutfile = outdir + 'mRNA' + title + ".tab"
            mRNAupdatedcount = 0
            miRupdatedcount = 0
            f_out = open(methods.pwd() + mRNAoutfile, 'w')
            f_out.write("gene \t chr \t strand \t count \t updated \t aveOccupBS \n")
            for mRNA in self.mRNA2dist:
                if('_' in mRNA):
                    gene = mRNA.split('_')[0]
                    if(len(mRNA.split('_')[1].split(':')) == 2):
                        curr_chr = mRNA.split('_')[1].split(':')[0]
                        curr_strand = mRNA.split('_')[1].split(':')[1]
                    else:
                        curr_chr = '*'
                        curr_strand = '*'
                else:
                    gene = mRNA
                    curr_chr = '*'
                    curr_strand = '*'
                count = self.mRNA2dist[mRNA]
                if(mRNA in self.updated_mRNAs):
                    updated = self.updated_mRNAs[mRNA]
                    mRNAupdatedcount += 1
                else:
                    updated = 0
                if(mRNA in self.mRNAoccupBScount):
                    totalBS = self.mRNAoccupBScount[mRNA][0]
                    totalmRNA = self.mRNAoccupBScount[mRNA][1]
                    aveOccupBS = float(totalBS/totalmRNA)
                else:
                    aveOccupBS = 0
                f_out.write(gene + "\t" + curr_chr + "\t" + curr_strand + "\t" + str(count) + "\t" + str(updated) + "\t" + str(aveOccupBS) + "\n")
            f_out.close()
            methods.logINFO("total mRNA: " + str(self.totalmRNA))
            methods.logINFO("updated mRNA count: " + str(mRNAupdatedcount))
            methods.logINFO("Removed due to half life time: " + str(self.dead_by_half_life))
        if outmiR:
            miRoutfile = outdir + 'miR' + title + ".tab"
            f_out = open(methods.pwd() + miRoutfile, 'w')
            f_out.write("miR \t count \t updated \n")
            for miR in self.miR2dist:
                count = self.miR2dist[miR]
                if(miR in self.updated_miRs):
                    updated = self.updated_miRs[miR]
                    miRupdatedcount += 1
                else:
                    updated = 0            
                f_out.write(miR + "\t" + str(count) + "\t" + str(updated) + "\n")
            f_out.close()
            methods.logINFO("total miR: " + str(self.totalmiR))
            methods.logINFO("updated miR count: " + str(miRupdatedcount))
        if outRetention:
            mRNA_retentionfile = outdir + 'retention' + title + '.tab'
            # This is for output the genes according to their retention rate.
            # I am checking that the threshold is 10
            f_out = open(methods.pwd()  + mRNA_retentionfile, 'w')
            header  = "Gene "
            for i in range(len(self.mRNA_counts_idx)):
                header += "\t " + str(self.mRNA_counts_idx[i])
            header += "\t retention \n"
            f_out.write(header)
            for i in range(len(sortedGenes)):
                gene = sortedGenes[i][0]
                line = gene 
                counts = self.mRNA_counts[gene]
                for j in range(len(counts)):
                    line += "\t " + str(self.mRNA_counts[gene][j])
                retention = 100*counts[-1]/counts[0]
                line += "\t " + str(retention) + " \n"
                f_out.write(line)
            f_out.close()
            methods.logINFO("Removed due to half life time: " + str(self.dead_by_half_life))
        #mRNAcountsfile = outdir + 'counts' + title + ".tab"        
        #=======================================================================
        # This is for genes that are sorted by their variance 
        # f_out = open(methods.pwd() + mRNAcountsfile, 'w')
        # header = "Gene "
        # for i in range(len(self.mRNAcountsIdx)):
        #    header += "\t " + str(self.mRNAcountsIdx[i])
        # header += " \n"
        # f_out.write(header)
        # for i in range(len(sortedGenes)):
        #    gene = sortedGenes[i][0]
        #    line = gene 
        #    for j in range(len(self.mRNAcounts[gene])):
        #        line += "\t " + str(self.mRNAcounts[gene][j])
        #    line += " \n"
        #    f_out.write(line)
        # f_out.close()
        #=======================================================================
        
    def remove_mRNA(self, r=100, T=24):
        for mRNA in self.mRNA2dist:
            if mRNA in self.half_life:
                for i in range(self.mRNA2dist[mRNA]):
                    t = self.half_life[mRNA]
                    epsilon = 1 - 0.5**(T/(r*t))
                    should_remove = nprand.choice([0,1], p=[1-epsilon, epsilon])
                    if should_remove:
                        self.mRNA2dist[mRNA] -= 1
                        #methods.logINFO(mRNA + " was removed due to life time")
                        self.dead_by_half_life += 1
        
    def compare_to_CLASH(self, clash_pairs, clash_pairs_cnt):
        pair_counts_sim = []
        pair_counts_clash = []
        for miRNA in self.pairs:
            for mRNA in self.pairs[miRNA]:
                curr_pair = miRNA + '_' + mRNA
                if curr_pair in clash_pairs:
                    pair_counts_sim.append(self.pairs[miRNA][mRNA])
                    pair_counts_clash.append(clash_pairs[curr_pair])
        corr_spearman = methods.corr_spearman(pair_counts_sim, pair_counts_clash)
        corr = methods.corr(pair_counts_sim, pair_counts_clash)
        methods.logINFO("numbers of pairs to compare with CLASH: " + str(len(pair_counts_sim)))
        pair_counts_clash = []
        pair_counts_sim_before_update = []
        for curr_pair in clash_pairs_cnt:
            if curr_pair in clash_pairs:
                pair_counts_clash.append(clash_pairs[curr_pair])
                pair_counts_sim_before_update.append(clash_pairs_cnt[curr_pair])
            else:
                methods.logERROR("Could not find pair: " + curr_pair +" in CLASH")
        corr_spearman_before_update = methods.corr_spearman(pair_counts_sim_before_update, pair_counts_clash)
        corr_before_update = methods.corr(pair_counts_sim_before_update, pair_counts_clash)
        methods.logINFO("numbers of pairs to compare with CLASH before update: " + str(len(pair_counts_sim_before_update)))
        return corr, corr_spearman, corr_before_update, corr_spearman_before_update

    def correlateBS2length(self):
        gene2len = {}
        f_lendata = open(methods.pwd() + "expression_data/genesLengthInfo.tab", 'r')
        i = 0
        for line in f_lendata:
            if(i > 0):
                row = line.split('\t')
                gene2len[row[0].strip()] = [int(row[4].strip()), int(row[5].strip())]
            i += 1
        f_lendata.close()
        aveBS = []
        UTRlen = []
        totalLen = []
        notFound = 0
        for mRNAfull in self.mRNAoccupBScount:
            mRNA = mRNAfull.split('_')[0]
            if(mRNA in gene2len):
                BScount = self.mRNAoccupBScount[mRNAfull][0]
                mRNAcount = self.mRNAoccupBScount[mRNAfull][1]
                aveBS.append(float(BScount/mRNAcount))
                UTRlen.append(gene2len[mRNA][0])
                totalLen.append(gene2len[mRNA][1])
            else:
                notFound += 1
        methods.logINFO("Correlation between average BS and 3' UTR length: " + str(methods.corr(aveBS, UTRlen)))
        methods.logINFO("Correlation between average BS and total length: " + str(methods.corr(aveBS, totalLen)))
        methods.logINFO("genes not found: " + str(notFound))
                
        
    def plotUpdateCount(self, title):
        types = ['miRs', 'mRNAs']
        for mol_type in types:   
            genes2plot = {}
            if(mol_type == 'miRs'):
                updatedCounts = self.updated_miRs
            else:
                updatedCounts = self.updated_mRNAs
            for item in updatedCounts:
                if(updatedCounts[item] > 1):
                    genes2plot[item] = updatedCounts[item]
            sortedmiRs = sorted(genes2plot.iteritems(), key=operator.itemgetter(1), reverse=True)
            genes = []
            counts = []
            i = 0            
            for item in sortedmiRs:
                if(i > 20):
                    break
                genes.append(item[0])
                counts.append(item[1])
                i += 1
            methods.plotBar(counts, genes, title + '_' + mol_type)
    # Choosing uniformly an mRNA and removing quantity of it
    def removeRandomFreemRNA(self, quantity):
        mRNAs = set()
        for mRNA in self.mRNA2dist:
            if(self.mRNA2dist[mRNA]> 0):
                mRNAs.add(mRNA)
        rand_mRNA = nprand.choice(list(mRNAs))
        #if(self.mRNA2dist[rand_mRNA] < quantity):
        #    quantity = math.floor(self.mRNA2dist[rand_mRNA]*self.totalmRNA)
        #[self.mRNA2dist, self.totalmRNA] = self.updateDist(rand_mRNA, self.mRNA2dist, self.totalmRNA, quantity, -1)
        quantity = min(quantity, self.mRNA2dist[rand_mRNA])
        self.mRNA2dist[rand_mRNA] -= quantity
        self.totalmRNA -= quantity
        return [rand_mRNA, quantity]
    
    def removeOccupiedmRNA(self, iter_num):
        mRNAs= list(self.mRNA_occupied.keys())
        miR2release = {}
        total_miR_removed = 0
        total_mRNA_removed = 0
        for mRNA in mRNAs:
            # remove mRNA that are not occupied:
            if self.mRNA_occupied[mRNA].occupied <= 0:
                del self.mRNA_occupied[mRNA]
            else:
                mRNAobj = self.mRNA_occupied[mRNA]
                [currmiR2release, totalmiR, totalmRNA] = mRNAobj.releasemRNAs(iter_num)
                miR2release = methods.mergeIntDicts(miR2release, currmiR2release)
                total_miR_removed += totalmiR
                total_mRNA_removed += totalmRNA
                if mRNAobj.occupied == 0:
                # Release all mRNAobj from self.mRNAoccupied
                    miR2release = methods.mergeIntDicts(miR2release,mRNAobj.totalmiRcount())
                    del self.mRNA_occupied[mRNA]
        if mRNAs:
            # Add the removed miRs back to the pool:
            for miR in miR2release:
                add_mir_cnt = miR2release[miR]
                prev_mir_cnt = self.miR2dist[miR]
                for idx in self.mir2idx[miR]:
                    if prev_mir_cnt > 0:
                        self.probabilities[idx] *= (prev_mir_cnt+add_mir_cnt)/prev_mir_cnt
                    else:
                        mir, target, bs = self.table_pairs[idx]
                        if mir != miR:
                            methods.logERROR("error in table index: "+mir + ", " + miR)
                        else:
                            if target in self.mRNA2dist:
                                total_gene = self.mRNA2dist[target]
                                score = self.score_table[mir][target][bs]
                                p = methods.score2prob(score, self.source)
                                self.probabilities[idx] = p*total_gene
                self.miR2dist[miR] += miR2release[miR]
                self.totalmiR += miR2release[miR]

                #[self.miR2dist, self.totalmiR] = self.updateDist(miR, self.miR2dist, self.totalmiR, miR2release[miR], 1)
        return [total_miR_removed, total_mRNA_removed]
    
    def addRandommRNA(self, quantity):
        rand_mRNA = nprand.choice(list(self.mRNA2dist.keys()))
        self.mRNA2dist[rand_mRNA] += quantity
        self.totalmRNA += quantity
        #self.updateDist(rand_mRNA, self.mRNA2dist, self.totalmRNA, quantity, 1)
               
    #===========================================================================
    # def removemiR(self, percentage):
    #    miRsorted = sorted(self.miR2dist.iteritems(), key=operator.itemgetter(1), reverse=True)
    #    miRs = []
    #    count = 0
    #    i = 0
    #    while(count <= percentage*self.totalmiR):
    #        miRs.append(miRsorted[i][0])
    #        count += miRsorted[i][1]
    #        i += 1
    #    methods.logINFO("miRs for random remove: " + ",".join(miRs))
    #    toRemoveIdx = nprand.choice(range(len(miRs)))
    #    miR2remove = miRs[toRemoveIdx]
    #    quantity = self.miR2dist[miR2remove]
    #    #self.miR2dist[miR2remove] = 0
    #    # I do not change the total miRs
    #    dist = methods.absolute2dist(self.miR2dist, self.totalmiR)
    #    dist[miR2remove] = 0
    #    dist = methods.fixDist(dist)
    #    [self.miR2dist, self.totalmiR] = methods.dist2absolute(dist, self.totalmiR, 1)
    #    methods.logINFO("Removed miR: " + miR2remove + ", total: " + str(quantity))
    #    return miR2remove
    #===========================================================================
    
    def change_miR_counts(self, factor, miR2change, table, percentage=0.5):
        """
        Changes miRNA count according to the given factor,
        or randomly chooses one miRNA that its count is lower than the given
        percentage of total miRNAs
        :param percentage:
        :param factor:
        :param miR2change:
        :return:
        """
        if not miR2change:  # randomly choose one miRNA that has relatively low expression
            miR_sorted = sorted(self.miR2dist.iteritems(), key=operator.itemgetter(1), reverse=True)
            miRs = []
            count = 0
            i = 0
            while count <= percentage*self.totalmiR:
                miRs.append(miR_sorted[i][0])
                count += miR_sorted[i][1]
                i += 1
            methods.logINFO("miRs for random change: " + ",".join(miRs))
            tochangeidx = nprand.choice(range(len(miRs)))
            miR2change = miRs[tochangeidx]
        input_mirs2change = set(miR2change.split(","))
        miRs2change = []
        for miR2change in input_mirs2change:
            if miR2change not in table:
                for miR in table:
                    if miR.startswith(miR2change):
                        miR_arr = miR.split("-")
                        if "-".join(miR_arr[:2]) == miR2change:
                            miRs2change.append(miR)
            else:
                miRs2change.append(miR2change)
            for miR in miRs2change:
                if miR not in self.miR2dist:  # if the miRNA is not expressed we initiate its counter 0.1% of total miRNA
                    miR_quanta = math.floor(self.totalmiR*0.0001)
                    self.miR2dist[miR] = miR_quanta
                    self.totalmiR += miR_quanta
                    methods.logINFO("miR not found in dist: " + miR)
                quantity = self.miR2dist[miR]
                methods.logINFO("quantity of miRNA to change: " + str(quantity))
                self.miR2dist[miR] *= factor
                self.totalmiR = sum(self.miR2dist.values())
                dist = methods.absolute2dist(self.miR2dist, self.totalmiR)
                #dist[miR] *= factor
                #dist = methods.fixDist(dist)
                [self.miR2dist, self.totalmiR] = methods.dist2absolute(dist, self.totalmiR, 1)
                if miR not in self.miR2dist:
                    quanta = 0
                else:
                    quanta = self.miR2dist[miR]
        methods.logINFO("Changed miRs: " + str(miRs2change) + ", factor: " + str(factor) + ", prev quanta: " + str(quantity) + ", curr quanta: " + str(quanta) +", totalmiRs: " + str(self.totalmiR))
        return miRs2change
    
    def recordmRNAcount(self, idx):
        self.mRNA_counts_idx.append(idx)
        for mRNA in self.mRNA2dist:
            if(mRNA in self.mRNA_counts):
                self.mRNA_counts[mRNA].append(self.mRNA2dist[mRNA])
            else:
                self.mRNA_counts[mRNA] = [self.mRNA2dist[mRNA]]
                
    
class BS:
    
    BSnotinGene = set()
    
    def __init__(self, gene, input_BS, form):
        if form == '[hg19:chr:start-end:-]':
            bs_arr = input_BS.split(':')
            self.chr = bs_arr[1]
            self.strand = bs_arr[3][0]
            self.id = int(input_BS.split(':')[0][1:])
            try:
                self.start = int(bs_arr[2].split(',')[0].split('-')[0])
                self.end = int(bs_arr[2].split(',')[0].split('-')[1])
            except:
                methods.logERROR(gene + " : " + input_BS)
                self.start = 0
                self.end = 22
            self.genename = gene
        elif form == 'relativeBS':
            #methods.logINFO("relativeBS form")
            self.genename = gene.name
            self.chr = gene.chr
            self.strand = gene.strand
            bs_arr = input_BS.split('_')
            self.start = int(gene.start) + int(bs_arr[0])
            self.end = int(gene.start) + int(bs_arr[1])
            if self.end > gene.end:
                BS.addBSnew(self.genename)
        self.free_bounded = 0  # how many miRs are bounded at this BS to free mRNA
        self.occupied_bounded = 0  # how many miRs are bounded at this BS to occupied mRNA
        # This should follow the number of times it was bounded on a free mRNA: {miR:{iterNum:quantity}}
        self.freemRNA = {}
        # This should follow the number of times it was bounded on an occupied mRNA: {miR:{iterNum:quantity}}
        self.occupiedmRNA = {}
        self.free2occupied = {}
        
    
    @staticmethod
    def addBSnew(genename):
        BS.BSnotinGene.add(genename)
        
    @staticmethod
    def printBSnew():
        methods.logINFO(BS.BSnotinGene)
       
    def compareBS(self, BSobj):
        if BSobj.genename == self.genename and \
                        BSobj.chr == self.chr and BSobj.strand == self.strand and \
                        BSobj.start == self.start and BSobj.end == self.end:
            return True
        return False
    
    def overlap(self, BSobj):
        if abs(int(self.start) - int(BSobj.start)) < 50 or abs(int(self.end) - int(BSobj.end)) < 50:
                return True
        return False
    
    def copy(self):
        copyBS = BS(self.genename, str(self) , '[hg19:chr:start-end:-]')
        copyBS.free_bounded = self.free_bounded
        copyBS.occupied_bounded = self.occupied_bounded
        return copyBS        
    
    def __str__(self):
        # return "[hg19:" + self.chr + ":" + str(self.start) + "-" + str(self.end) + ':' + self.strand + "]"
        return "[" + str(self.id) + ":" + self.chr + ":" + str(self.start) + "-" + str(self.end) + ':' + self.strand + "]"
    
    def getID(self):
        return self.id
    
    def setID(self, bsID):
        self.id = bsID
    
class Gene:
    
    def __init__(self, genename, in_chr, strand, start, end):
        self.name = genename
        self.chr = in_chr
        self.strand = strand
        self.start = start
        self.end = end
        self.lifetime = -1

        
    #===========================================================================
    # @property
    # def name(self):
    #     return self.name
    # 
    # @property
    # def chr(self):
    #     return self.chr
    # 
    # @property
    # def strand(self):
    #     return self.strand
    # 
    # @property
    # def start(self):
    #     return self.start
    # 
    # @property
    # def end(self):
    #     return self.end    
    #===========================================================================
    
class mRNAObj:
    
    # BSoccupied: {bs:miR:{iterNum of insertion:quanta}}}
    # test: BSoccupied: {regionBS:{bs:{miR:{iter:quanta}}}}
    
    def __init__(self, name, del_interval):
        self.name = name
        self.occupied = 0
        self.del_interval = del_interval
        self.id2BS = {}
        self.id2rBS = {}
        #self.freeBS2occupiedBS = {}
        self.miRscount = {}
        self.rBS2BS = {}
        self.BS2rBS = {}
        self.rBScounter = 1

        
    #===========================================================================
    # @property
    # def name(self):
    #     return self.name
    # @property
    # def occupied(self):
    #     return self.occupied
    # @property
    # def free(self):
    #     return self.free
    # @property
    # def BSoccupied(self):
    #     return self.BSoccupied
    # @property
    # def updatequanta(self):
    #     return self.updatequanta
    #===========================================================================
    
    def copyBSmap(self, BSmap):
        copyBSmap = {}
        for bs in BSmap:
            copyBS = bs.copy()
            copyBSmap[copyBS] = BSmap[bs].copy()
        return copyBSmap
    
    def copyrBSmap(self, rBSmap):
        copyrBSmap = {}
        for rBS in rBSmap:
            copyrBS = rBS.copy()
            copyBSmap = self.copyBSmap(rBSmap[rBS])
            copyrBSmap[copyrBS] = copyBSmap
        return copyrBSmap
                
    
    def updateOccupiedCount(self, quantity):
        self.occupied += quantity
        return self.occupied
    
    def updateBSfree(self, bs, miRNA, iter_num, free_quanta):
        if isinstance(bs, BS):
            BSobj = bs
            #bs = BSobj.toString()
            bsID = BSobj.getID()
        else:
            bsID = int(bs.split(':')[0][1:])
            if bsID in self.id2BS:  # check if the BS is already found
                BSobj = self.id2BS[bsID]
                #bs = BSobj.toString()
            else:
                BSobj = BS(self.name, bs, '[hg19:chr:start-end:-]')
        if iter_num in BSobj.freemRNA:
            if miRNA in BSobj.freemRNA[iter_num]:
                BSobj.freemRNA[iter_num][miRNA] += free_quanta
            else:
                BSobj.freemRNA[iter_num][miRNA] = free_quanta
        else:
            BSobj.freemRNA[iter_num] = {miRNA : free_quanta}
        BSobj.free_bounded += free_quanta
        #self.updateBSoccupiedField(BSobj, miR, iterNum, free_quanta)
        self.id2BS[bsID] = BSobj
        self.updaterBS2BS(BSobj, free_quanta,1)
        #self.BS2rBS[bs] = rBS
        if(miRNA in self.miRscount):
            self.miRscount[miRNA] += free_quanta
        else:
            self.miRscount[miRNA] = free_quanta
        self.occupied += free_quanta
        return BSobj
    
    def updaterBS2BS(self, BSobj, quantity, freebound):
        foundRegion = 0
        foundBS = 0
        #bs = BSobj.toString()
        bsID = BSobj.getID()
        for rBSid in self.rBS2BS:
            rBSobj = self.id2rBS[rBSid]
            if(BSobj.overlap(rBSobj)):
                foundRegion = 1
                if(freebound):
                    rBSobj.free_bounded += quantity
                else:
                    rBSobj.occupied_bounded += quantity
                for currBSid in self.rBS2BS[rBSid]:
                    currBSobj = self.id2BS[currBSid]
                    if(BSobj.compareBS(currBSobj)):
                        if(not(currBSobj.free_bounded == BSobj.free_bounded or currBSobj.occupied_bounded == BSobj.occupied_bounded)):
                            methods.logERROR("not(currBSobj.free_bounded == BSobj.free_bounded or currBSobj.occupied_bounded == BSobj.occupied_bounded)")
                            currBSobj.free_bounded = BSobj.free_bounded
                            currBSobj.occupied_bounded = BSobj.occupied_bounded
                        if(not(bsID == currBSid)):
                            methods.logERROR("bsID != currBSid: " + str(BSobj) + ", " + str(currBSobj))
                        BSobj = currBSobj
                        foundBS = 1
                        break
                if(not(foundBS)):
                    oldstart = rBSobj.start
                    oldend = rBSobj.end
                    rBSobj.start = min(rBSobj.start, BSobj.start)
                    rBSobj.end = max(rBSobj.end, BSobj.end)
                    #updatedrBS = rBSobj.toString()
                    self.rBS2BS[rBSid].add(bsID)
                    #if(not(updatedrBS == rBSid)):
                    if(not(rBSobj.start == oldstart and rBSobj.end == oldend)):
                        BSset = self.rBS2BS[rBSid].copy()
                        updatedrBSid = self.rBScounter
                        rBSobj.setID(updatedrBSid)
                        self.rBScounter += 1
                        for irBSid in list(self.rBS2BS.keys()):
                            irBSobj = self.id2rBS[irBSid]
                            if(not(irBSid == rBSid) and rBSobj.overlap(irBSobj)):
                                rBSobj.start = min(rBSobj.start, irBSobj.start)
                                rBSobj.end = max(rBSobj.end, irBSobj.end)
                                #updatedrBS = rBSobj.toString()                                
                                BSset.update(self.rBS2BS[irBSid].copy())
                                rBSobj.free_bounded += irBSobj.free_bounded
                                rBSobj.occupied_bounded += irBSobj.occupied_bounded
                                del self.rBS2BS[irBSid]
                                del self.id2rBS[irBSid]
                        del self.rBS2BS[rBSid]
                        del self.id2rBS[rBSid]
                        self.rBS2BS[updatedrBSid] = BSset
                        self.id2rBS[updatedrBSid] = rBSobj
                        rBSobj.setID(updatedrBSid)
                        #self.id2BS[bs] = BSobj
                        for currBSid in BSset:
                            self.BS2rBS[currBSid] = updatedrBSid
                    else:
                        self.BS2rBS[bsID] = rBSid
                break
        if not foundRegion:
            rBSobj = BSobj.copy()
            #rBSid = rBSobj.toString()
            rBSid = self.rBScounter
            self.rBScounter += 1
            rBSobj.setID(rBSid)
            self.rBS2BS[rBSid] = {bsID}
            self.id2rBS[rBSid] = rBSobj
            self.BS2rBS[bsID] = rBSid
            #self.id2BS[bs] = BSobj
        if bsID not in self.BS2rBS:
            methods.logERROR("in updaterBS2BS bsID not in self.BS2rBS")
        return rBSobj   

    def updateBSoccupied(self, bs, miR, iter_num, occupied_quanta):
        """
        Binds the given miRNA  to the already occupied mRNA
        :param bs:
        :param miR:
        :param iter_num:
        :param occupied_quanta:
        :return:
        """
        if isinstance(bs, BS):
            BSobj = bs
            bsID = BSobj.getID()
        else:
            bsID = int(bs.split(':')[0][1:])
            if bsID in self.id2BS:
                BSobj = self.id2BS[bsID]
                #bs = BSobj.toString()
            else:
                BSobj = BS(self.name, bs, '[hg19:chr:start-end:-]')
        quanta_updated = 0
        curr_total = 0
        # we count how many total miRs are bounded to this BS
        for curriter in BSobj.occupiedmRNA:
            for curr_miR in BSobj.occupiedmRNA[curriter]:
                curr_total += BSobj.occupiedmRNA[curriter][curr_miR]
        if curr_total != BSobj.occupied_bounded:  # sanity check
            methods.logERROR("Begin: curr_total != BSobj.occupied_bounded: " + str(BSobj.toString()) + ", iterNUm: " + str(iter_num))
        # In the case of initiation we do not need to update free2occupied field:
        if iter_num == 0:
            quanta_updated = occupied_quanta
            if iter_num in BSobj.occupiedmRNA:
                if miR in BSobj.occupiedmRNA[0]:
                    BSobj.occupiedmRNA[0][miR] += quanta_updated
                else:
                    BSobj.occupiedmRNA[0][miR] = quanta_updated
            else:
                BSobj.occupiedmRNA[0] = {miR: quanta_updated}
            BSobj.occupied_bounded += quanta_updated
        else:
            rBSlst = []
            total_available = 0
            for rBSid in list(self.rBS2BS.keys()):
                rBSobj = self.id2rBS[rBSid]
                if not rBSobj.overlap(BSobj):
                    rBSlst.append(rBSid)
            if BSobj.free_bounded + BSobj.occupied_bounded >= self.occupied:
                methods.logERROR("in updateBSoccupied: bs.bounded > mRNAobj.occupied")
            if rBSlst:  # choose now mRNA that was bounded on free mRNA
                methods.permute(rBSlst)
                #permuterBSs = rBSlst
                for rBSid in rBSlst:
                    if self.rBS2BS[rBSid]:
                        rBSobj = self.id2rBS[rBSid]
                        BSlist = list(self.rBS2BS[rBSid])
                        methods.permute(BSlist)
                        #permuteBSs = list(self.rBS2BS[rBSid])
                        for currBSid in BSlist:
                            curr_freeBS_obj = self.id2BS[currBSid]
                            if curr_freeBS_obj.freemRNA:
                                freelst = list(curr_freeBS_obj.freemRNA.keys())
                                methods.permute(freelst) # for choose randomly
                                #permuteFreelst = list(curr_freeBS_obj.freemRNA.keys())
                                for freeiter in freelst:
                                    available = 0
                                    for currmiR in curr_freeBS_obj.freemRNA[freeiter]:
                                        available += curr_freeBS_obj.freemRNA[freeiter][currmiR]
                                    if not available:
                                        methods.logERROR("available == 0")
                                    total_available += available
                                    foundBSobj = 0
                                    if freeiter in curr_freeBS_obj.free2occupied:
                                        #methods.logINFO("updateBSoccupied in curr_freeBS_obj.free2occupied")
                                        for BSoccupiedID in curr_freeBS_obj.free2occupied[freeiter]:  # remove the occupied BS counter from the total available
                                            BSoccobj = self.id2BS[BSoccupiedID]
                                            if BSoccobj.overlap(BSobj):
                                                if BSoccupiedID == bsID:
                                                    foundBSobj = BSoccobj
                                                for occupiter in curr_freeBS_obj.free2occupied[freeiter][BSoccupiedID]:
                                                    # when occupiter == 0 we do not update curr_freeBS_obj.free2occupied so it is not relevant
                                                    if occupiter > 0:
                                                        for currmiR in BSoccobj.occupiedmRNA[occupiter]:
                                                            available -= BSoccobj.occupiedmRNA[occupiter][currmiR]
                                                            total_available -= BSoccobj.occupiedmRNA[occupiter][currmiR]
                                                            if available < 0:
                                                                methods.logERROR("available < 0")
                                    q2update = min(max(available,0), (occupied_quanta - quanta_updated))
                                    #methods.logINFO("need to update: " + str(occupied_quanta) + ", q2update: " + str(q2update))
                                    if q2update < 0:
                                        methods.logERROR("q2update < 0")
                                    quanta_updated += q2update
                                    if q2update:
                                        if not foundBSobj:
                                            foundBSobj = BSobj
                                        if iter_num in foundBSobj.occupiedmRNA:
                                            if miR in foundBSobj.occupiedmRNA[iter_num]:
                                                foundBSobj.occupiedmRNA[iter_num][miR] += q2update
                                            else:
                                                foundBSobj.occupiedmRNA[iter_num][miR] = q2update
                                        else:
                                            foundBSobj.occupiedmRNA[iter_num] = {miR : q2update}
                                        if freeiter in curr_freeBS_obj.free2occupied:
                                            if bsID in curr_freeBS_obj.free2occupied[freeiter]:
                                                #curr_freeBS_obj.free2occupied[freeiter][bs].add(iterNum)
                                                curr_freeBS_obj.free2occupied[freeiter][bsID][iter_num] = miR + '_' + str(q2update)
                                            else:
                                                curr_freeBS_obj.free2occupied[freeiter][bsID] = {iter_num: miR + '_' + str(q2update)}
                                        else:
                                            curr_freeBS_obj.free2occupied[freeiter] = {bsID: {iter_num: miR + '_' + str(q2update)}}
                                        foundBSobj.occupied_bounded += q2update
                                        #rBSobj.occupiedbounded += q2update                                        
                                    if quanta_updated == occupied_quanta:
                                        break                                                       
                                if quanta_updated == occupied_quanta:
                                    break                                                
                        if quanta_updated == occupied_quanta:
                            break
            if quanta_updated < occupied_quanta:
                methods.logERROR("mRNA: " + self.name + ", miR: " + miR + ", quanta_updated: " + str(quanta_updated) + ", occupied_quanta: " + str(occupied_quanta) + ", self.isBSavailable(BSobj): " + str(self.isBSavailable(BSobj)) + ", total_available: " + str(total_available) )
        self.id2BS[bsID] = BSobj
        rBSobj = self.updaterBS2BS(BSobj, quanta_updated, 0)
        self.BS2rBS[bsID] = rBSobj.getID()
        if miR in self.miRscount:
            self.miRscount[miR] += quanta_updated
        else:
            self.miRscount[miR] = quanta_updated
        curr_total = 0
        # Another sanity check:
        for curriter in BSobj.occupiedmRNA:
            for miR in BSobj.occupiedmRNA[curriter]:
                curr_total += BSobj.occupiedmRNA[curriter][miR]
        if curr_total != BSobj.occupied_bounded:
            methods.logERROR("curr_total != BSobj.occupied_bounded: " + str(BSobj.toString()) + ", iterNUm: " + str(iter_num) + ", curr_total: " + str(curr_total) + ", BSobj.occupied_bounded: " + str(BSobj.occupied_bounded) + ", q2update: " + str(q2update) + ", occupied_quanta: " + str(occupied_quanta) + ", available: " + str(total_available) + ", quanta_updated: " + str(quanta_updated))
        return BSobj

    def isBSavailable(self, BSobj):
        """
        Checks if the given BS is available, by checking also overlapping occupied BS
        :param BSobj:
        :return: number of available bs???
        """
        available = self.occupied
        for rBS in self.rBS2BS:
            rBSobj = self.id2rBS[rBS]
            if rBSobj.overlap(BSobj):
                available -= rBSobj.free_bounded + rBSobj.occupied_bounded
                if(not(available%1 ==0)):
                    methods.logERROR("isBSavailable: available%1 != 0")
                return self.occupied - (rBSobj.free_bounded + rBSobj.occupied_bounded)
        return available
    
    def totalmiRcount(self):
        return self.miRscount
            
    def releasemRNAs(self, iter_num):
        miR2release = {}
        totalmRNAremoved = 0
        totalmiRremoved = 0
        BS2remove = set()
        for rBS in list(self.rBS2BS.keys()):
            rBSobj = self.id2rBS[rBS]
            BSsetcopy = self.rBS2BS[rBS].copy()
            for BS in self.rBS2BS[rBS]:
                BSobj = self.id2BS[BS]
                # remove all occupied that were bounded at initiation (0 is the number of iteration):
                if (iter_num == self.del_interval/2) and (0 in BSobj.occupiedmRNA.keys()):
                    curr_total = 0
                    try:
                        for miR in BSobj.occupiedmRNA[0]:
                            quanta = BSobj.occupiedmRNA[0][miR]
                            curr_total += quanta
                            if miR in miR2release:
                                miR2release[miR] += quanta
                            else:
                                miR2release[miR] = quanta
                            self.miRscount[miR] -= quanta
                        totalmiRremoved += curr_total
                        rBSobj.occupied_bounded -= curr_total
                        BSobj.occupied_bounded -= curr_total
                        del BSobj.occupiedmRNA[0]
                        # find occupied bs that can appear on the same mRNA:
                        nonoverlaprBSs = set(self.rBS2BS.keys()).copy()
                        nonoverlaprBSs.remove(rBS)
                        methods.permute(list(nonoverlaprBSs))
                        free_need2remove = curr_total
                        for nonoverlaprBS in nonoverlaprBSs:
                            BSs = self.rBS2BS[nonoverlaprBS]
                            for currBS in BSs:
                                currBSobj = self.id2BS[currBS]
                                if 0 in currBSobj.freemRNA:
                                    for miR in list(currBSobj.freemRNA[0].keys()):
                                        currRemove = min(free_need2remove, currBSobj.freemRNA[0][miR])
                                        currBSobj.freemRNA[0][miR] -= currRemove
                                        if currBSobj.freemRNA[0][miR] == 0:
                                            del currBSobj.freemRNA[0][miR]
                                        if miR in miR2release:
                                            miR2release[miR] += currRemove
                                        else:
                                            miR2release[miR] = currRemove
                                        free_need2remove -= currRemove
                                        totalmiRremoved += currRemove
                                        totalmRNAremoved += currRemove
                                        currBSobj.free_bounded -= currRemove
                                        self.id2rBS[nonoverlaprBS].free_bounded -= currRemove
                                        if not(free_need2remove):
                                            break
                                if not free_need2remove:
                                    break
                            if not free_need2remove:
                                break                                                                          
                    except:
                        methods.logERROR("failed in removing initiation: miR: " + miR + ", BSobj.occupiedmRNA: " + str(BSobj.occupiedmRNA))
                        exit()                
                # Now for bounded mRNA after the simulation started
                for free_iter in list(BSobj.freemRNA.keys()):  # need to get by the first that was bounded
                    is_occupied = False
                    if iter_num - free_iter >= self.del_interval/2:
                        visitedBS = set()
                        if free_iter in BSobj.free2occupied:
                            # if there is more than one BS we would remove after self.delInterval/2:
                            if(len(BSobj.free2occupied[free_iter]) > 0):
                                is_occupied = 1
                                for BSoccupiedID in list(BSobj.free2occupied[free_iter].keys()):
                                    if(BSoccupiedID in visitedBS):
                                        methods.logERROR("already been in this BS")
                                    else:
                                        visitedBS.add(BSoccupiedID)
                                    BSoccobj = self.id2BS[BSoccupiedID]
                                    for occupiter in BSobj.free2occupied[free_iter][BSoccupiedID]:
                                        curr_total = 0
                                        if(occupiter not in BSoccobj.occupiedmRNA):
                                            methods.logERROR("occupiter not in BSoccobj.occupiedmRNA: " + BSoccupiedID  + ", BS: " + BS + ", occupiter: " + str(occupiter) + ", free_iter: " + str(free_iter))
                                            continue
                                        miR_quanta = BSobj.free2occupied[free_iter][BSoccupiedID][occupiter]
                                        miR = miR_quanta.split('_')[0]
                                        quanta = int(float(miR_quanta.split('_')[1]))
                                        if(miR in BSoccobj.occupiedmRNA[occupiter]):
                                            BSoccobj.occupiedmRNA[occupiter][miR] -= quanta
                                            curr_total += quanta
                                            if(miR in miR2release):
                                                miR2release[miR] += quanta
                                            else:
                                                miR2release[miR] = quanta
                                            self.miRscount[miR] -= quanta
                                        else:
                                            methods.logERROR("miR was not found in BSoccobj.occupiedmRNA[occupiter]: " + str(BSoccobj.occupiedmRNA[occupiter]))
                                        if(BSoccobj.occupiedmRNA[occupiter][miR] == 0):                                        
                                            del BSoccobj.occupiedmRNA[occupiter]
                                        elif(BSoccobj.occupiedmRNA[occupiter][miR] < 0):
                                            methods.logERROR("BSoccobj.occupiedmRNA[occupiter][miR] < 0: " + str(BSoccobj.occupiedmRNA[occupiter][miR] < 0))
                                        #methods.logINFO("removing from occupied: " + BSoccupiedID + ", BS: " + BS + ", occupiter: " + str(occupiter) + ", free_iter: " + str(free_iter) + ", BSoccobj.occupiedbounded: " + str(BSoccobj.occupiedbounded) + " - " + str(curr_total))
                                        BSoccobj.occupied_bounded -= curr_total
                                        totalmiRremoved += curr_total
                                        rBSID = self.BS2rBS[BSoccupiedID]
                                        occrBSobj = self.id2rBS[rBSID]
                                        occrBSobj.occupied_bounded -= curr_total
                                    if(BSoccobj.occupied_bounded + BSoccobj.free_bounded == 0):
                                        if(BSoccupiedID in BS2remove):
                                            methods.logERROR("BS is already exist in BS2remove: " + BSoccupiedID)
                                        else:
                                            BS2remove.add(BSoccupiedID)
                                        if(BSoccupiedID not in self.BS2rBS):
                                            methods.logERROR("BSoccupiedID not in BS2rBS: " + BSoccupiedID)    
                                    copyset = BSobj.free2occupied[free_iter][BSoccupiedID].copy()
                                    for occupiter in copyset:
                                        del BSobj.free2occupied[free_iter][BSoccupiedID][occupiter]
                                    if(not(BSobj.free2occupied[free_iter][BSoccupiedID])):
                                        del BSobj.free2occupied[free_iter][BSoccupiedID]
                            if(not(BSobj.free2occupied[free_iter])):
                                del BSobj.free2occupied[free_iter]
                        if(is_occupied or iter_num - free_iter >= self.del_interval):
                            curr_total = 0
                            for miR in list(BSobj.freemRNA[free_iter].keys()):
                                quanta = BSobj.freemRNA[free_iter][miR]
                                curr_total += quanta
                                if(miR in miR2release):
                                    miR2release[miR] += quanta
                                else:
                                    miR2release[miR] = quanta
                                self.miRscount[miR] -= quanta
                                del BSobj.freemRNA[free_iter][miR]
                            if(not(BSobj.freemRNA[free_iter])):
                                del BSobj.freemRNA[free_iter]
                            BSobj.free_bounded -= curr_total
                            rBSobj.free_bounded -= curr_total
                            totalmiRremoved += curr_total
                            totalmRNAremoved += curr_total
                        if(BSobj.occupied_bounded + BSobj.free_bounded == 0):
                            if(BS in BS2remove):
                                methods.logERROR("BS is already exist in BS2remove: " + str(BS))
                            else:
                                BS2remove.add(BS)
                                if(BS not in self.BS2rBS):
                                    methods.logERROR("BS not in BS2rBS: " + BS + ", iterNum: " + str(iter_num))
                            BSsetcopy.remove(BS)                            
            self.rBS2BS[rBS] = BSsetcopy
        for BS in BS2remove:
            rBS = self.BS2rBS[BS]            
            del self.BS2rBS[BS]
            del self.id2BS[BS]
            if(rBS in self.rBS2BS and BS in self.rBS2BS[rBS]):
                self.rBS2BS[rBS].remove(BS)
                if(not(self.rBS2BS[rBS])):
                    del self.rBS2BS[rBS]
                    del self.id2rBS[rBS]         
        self.occupied -= totalmRNAremoved
        #shouldexit = 0
        if(totalmRNAremoved > 0):
            for rBS in self.rBS2BS:
                rBSobj = self.id2rBS[rBS]
                #methods.logINFO("rBS: " + rBS + ", rBSobj.freebounded: " + str(rBSobj.freebounded) + ", rBSobj.occupiedbounded: " + str(rBSobj.occupiedbounded) + ", self.occupied: " + str(self.occupied))
                if(rBSobj.free_bounded + rBSobj.occupied_bounded > self.occupied):
                    methods.logERROR("rBS.free_bounded + rBS.occupied_bounded > self.occupied: " + str(rBS) +", iterNum: " + str(iter_num) and abs(rBSobj.end - rBSobj.start) <= 70)
                    #shouldexit = 1
                for bs in self.rBS2BS[rBS]:
                    BSobj = self.id2BS[bs]
                    #methods.logINFO("BS: " + bs + ", BSobj.freebounded: " + str(BSobj.freebounded) + ", BSobj.occupiedbounded: " + str(BSobj.occupiedbounded) + ", self.occupied: " + str(self.occupied))
                    if(BSobj.free_bounded + BSobj.occupied_bounded > self.occupied):
                        methods.logERROR("BS.free_bounded + BS.occupied_bounded > self.occupied: " + bs + ", iterNum: " + str(iter_num))
                        methods.logERROR("BS.freemRNA length: " + str(len(BSobj.freemRNA)))
                        for curriter in BSobj.freemRNA:
                            methods.logERROR("freemRNA iter: " + str(curriter) + ":")
                            for miR in BSobj.freemRNA[curriter]:
                                methods.logERROR(miR + " : " + str(BSobj.freemRNA[curriter][miR]))
                        methods.logERROR("BS.occupiedmRNA length: " + str(len(BSobj.occupiedmRNA)))
                        for curriter in BSobj.occupiedmRNA:
                            methods.logERROR("occupiedmRNA iter: " + str(curriter) + ":")
                            for miR in BSobj.occupiedmRNA[curriter]:
                                methods.logERROR(miR + " : " + str(BSobj.occupiedmRNA[curriter][miR]))
                        #shouldexit = 1
        #if(shouldexit):
        #    exit()
        return [miR2release, totalmiRremoved , totalmRNAremoved]
            
                                
                            
                            
    
    
