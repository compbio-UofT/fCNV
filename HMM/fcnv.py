#!/usr/bin/pypy

import math
import argparse
import numpypy as np
import sys
import fcnvHMM
import cvrgHMM
from datetime import datetime
import os


def readInput(in_file_name):
    """Read the input from files.
    Read maternal alleles, paternal alleles, and samples per position.
    
    Arguments:
    file name -- string
    returns lists of parsed data    
    
    """
    in_file = open(in_file_name, 'r')
    positions = []
    samples = []
    M = []; P = [];
    MC = []; PC = [];
    while True:
        line = in_file.readline()
        if not line: break
        if line[0] == '#': continue  #skip comment
        line = line.rstrip('\n').split('\t')
        
        #genomic positions and allele support in plasma samples
        positions.append(int(line[0]))
        samples.append(tuple(map(int, line[1:5])))
        
        #maternal and paternal alleles
        M.append(tuple(line[5:7]))
        MC.append(tuple(map(float, line[7:9])))
        
        P.append(tuple(line[9:11]))
        PC.append(tuple(map(float, line[11:13])))     
    
    in_file.close()
    return positions, samples, M, P, MC, PC

def computeEval(reference, prediction, sensitivity, num_patt):
    """Compute evaluation of referenced anotation versus predicted one.
    
    Prediction for one continuous region R of an inheritance pattern is considered
    correct if there is a continuous part with length at least '1/sensitivity'*len(R)
    that is predicted correctly.
    
    """
    num_real = 0
    num_called = 0
    ok_called = []
    not_called = []
    
    ind = 0
    while ind < len(reference):
        ctype = reference[ind]
        length = 0
        begin_pos = ind
        while ind < len(reference) and reference[ind] == ctype:
            length += 1
            ind += 1
        end_pos = ind
        
        max_ok_length = 0
        tmp_length = 0
        stats = [0 for x in range(num_patt)]
        for i in range(begin_pos, end_pos):
            stats[prediction[i]] += 1
            if prediction[i] == ctype:
                tmp_length += 1
            else:
                #print i, prediction[i],  ctype, type(prediction[i]), type(ctype)
                max_ok_length = max(max_ok_length, tmp_length)
                tmp_length = 0
        max_ok_length = max(max_ok_length, tmp_length)
        
        num_real += 1
        if max_ok_length >= length/float(sensitivity):
            num_called += 1
            ok_called.append((length, ctype))
        else:
            max_predic = -1
            max_type = -1
            for i in range(len(stats)):
                if stats[i] > max_predic:
                    max_predic = stats[i]
                    max_type = i
            not_called.append((length, ctype, max_type))
        #else: print  max_ok_length, length, reference[ind], prediction[ind], ctype
    return num_called, num_real, ok_called, not_called

def categorize(data, by_length, by_type):
    lengths = [100, 500, 1000, 5000, 10000, 20000, 50000, 100000]
    for item in data:
        length = item[0]
        ptype = item[1]
        gen_len = sorted([(abs(length-l), l) for l in lengths])[0][1]
        #print length, gen_len
        
        #by length
        if not gen_len in by_length: by_length[gen_len] = 0
        if ptype != 3: by_length[gen_len] += 1 #don't count normal IP into CNV length stats
        #by type
        if not ptype in by_type: by_type[ptype] = 0
        by_type[ptype] += 1

def computeStats(ok, wrong, pref, num_patt):
    lengths = [100, 500, 1000, 5000, 10000, 20000, 50000, 100000]
    ok_by_length = dict([(i,0) for i in lengths])
    ok_by_type = dict([(i,0) for i in range(num_patt)])
    wr_by_length = dict([(i,0) for i in lengths])
    wr_by_type = dict([(i,0) for i in range(num_patt)])
    categorize(ok, ok_by_length, ok_by_type)
    categorize(wrong, wr_by_length, wr_by_type)
    for l in lengths:
        o = ok_by_length[l]
        w = wr_by_length[l]
        ratio = 0
        if o+w > 0: ratio = o/float(o+w) *100
        print pref, l, ": ", o, w, ratio, '%'
    for t in range(num_patt):
        o = ok_by_type[t]
        w = wr_by_type[t]
        ratio = 0
        if o+w > 0: ratio = o/float(o+w) *100
        print pref, t, ": ", o, w, ratio, '%'

def test(fcnv, snp_positions, samples, M, P, MSC, PSC, mixture, ground_truth, file_name_prefix):
    #el = fcnv.extendedLabeling(samples, M, P, mixture)
    #el = fcnv.mixedDecoding(samples, M, P, mixture) 
    vp, v_state_path = fcnv.viterbiPath(samples, M, P, MSC, PSC, mixture)
#    posterior = fcnv.posteriorDecoding(samples, M, P, MSC, PSC, mixture)
#    byLL = fcnv.likelihoodDecoding(samples, M, P, MSC, PSC, mixture)    
    
    date = datetime.now().strftime('%m-%d-%H-%M')
    #fout = file(file_name_prefix + ".prediction" + date + ".txt", 'w')
    annot_out = file(file_name_prefix + ".annotation" + date + ".txt", 'w')
    
    
    print fcnv.inheritance_patterns
    num_patt = fcnv.getNumIP()
    
    #labeling_correct = 0
    viterbi_correct = 0
#    max_posterior_correct = 0
#    ll_correct = 0
    #states distributions for posterior and pure likelihood
#    posterior_dist = [0 for x in range(num_patt)]
#    ll_dist = [0 for x in range(num_patt)]
    #other stats
#    state_stats = [0 for x in range(num_patt)]
    #state_correct_el = [0 for x in range(num_patt)]
    state_correct_vp = [0 for x in range(num_patt)]
#    state_correct_pp = [0 for x in range(num_patt)]
#    state_correct_ll = [0 for x in range(num_patt)]
#    pp = []
#    tableLH = []
#    avg_distance = 0.
#    distance_set = {}
#    ll_misclass = [[0]*num_patt for x in range(num_patt)]
#    state3_stats = [{ (0,0):0, (0,1):0, (1,0):0, (1,1):0 } for x in range(num_patt)]
    for i in xrange(len(vp)):
#        state_stats[ground_truth[i]] += 1
#        post = []
#        for x in posterior[i]:   
#            post.append(x[1])
#        pp.append(post[0])
#        
#        ll_state = []
#        ll_value = []
#        for x in byLL[i]:   
#            ll_state.append(x[1])
#            ll_value.append(x[0])
#        tableLH.append(ll_value)
        
        #print the results
        #print >>fout, i+1, samples[i], M[i], list(MSC[i]), P[i], list(PSC[i]), vp[i], post[0], [ (ll_state[x], int(ll_value[x]*1e5)/1.e5) for x in range(len(ll_state))]
        print >>annot_out, snp_positions[i], ground_truth[i], vp[i] #, v_state_path[i]
            
        #print ground_truth[i], vp[i], pp[i], '|', post
        #labeling_correct += int(ground_truth[i] == el[i])
        #state_correct_el[ground_truth[i]] += int(ground_truth[i] == el[i])
        viterbi_correct += int(ground_truth[i] == vp[i])
        state_correct_vp[ground_truth[i]] += int(ground_truth[i] == vp[i])
#        max_posterior_correct += int(ground_truth[i] == post[0])
#        state_correct_pp[ground_truth[i]] += int(ground_truth[i] == post[0])
#        x = post.index(ground_truth[i])
#        posterior_dist[x] += 1
        
#        #for pure likelihood
#        threshold = .1
#        x = ll_state.index(ground_truth[i])
#        ll_diff = abs(ll_value[x]-ll_value[0])
#        if ll_diff < threshold:
#            ll_dist[0] += 1
#            #ll_correct += int(ground_truth[i] == ll_state[0])
#            ll_correct += 1
#            #state_correct_ll[ground_truth[i]] += int(ground_truth[i] == ll_state[0])
#            state_correct_ll[ground_truth[i]] += 1
#        else:            
#            ll_dist[x] += 1
#            for j in range(num_patt):
#                if ll_value[x] - ll_value[j] < 0.01:
#                #if abs(ll_value[0] - ll_value[j]) < 1.:
#                    ll_misclass[ground_truth[i]][ll_state[j]] += 1
#                    if ground_truth[i] == 3:
#                        m_poly = M[i][0] == M[i][1]
#                        p_poly = P[i][0] == P[i][1]
#                        state3_stats[ll_state[j]][(m_poly, p_poly)] += 1
#        
#        avg_distance += ll_diff
#        int_distance = int(round(ll_diff))
#        if int_distance not in distance_set:
#            distance_set[int_distance] = 1
#        else:
#            distance_set[int_distance] += 1
#        #print ground_truth[i], vp[i], post, '|', post.index(int(ground_truth[i]))
    
    
#    posterior_dist = np.array(posterior_dist)
#    posterior_dist = (posterior_dist*100.)/len(vp)    
#    print posterior_dist  
#    ll_dist = np.array(ll_dist)
#    ll_dist = (ll_dist*100.)/len(vp)    
#    print ll_dist
#    print "stats of misclassified by LL:"
#    ll_misclass = np.array(ll_misclass, float)
#    for i in range(len(ll_misclass)):
#        ll_misclass[i] = 100.*ll_misclass[i]/float(state_stats[i])
#        for j in range(len(ll_misclass[i])): ll_misclass[i][j] = round(ll_misclass[i][j]*100)/100.
#    print ll_misclass
    #print "mistakes for state3: "
    #for x in range(7):
    #    print x,":", state3_stats[x]
#    print "stats  : ", state_stats
    #print "mplabel: ", state_correct_el
    print "viterbi: ", state_correct_vp
#    print "mposter: ", state_correct_pp
#    print "byLHood: ", state_correct_ll
#    print "avgDist: ", avg_distance / len(vp) 
#    print "dist_set:", distance_set
    
    #print 'Labeling : ',(labeling_correct*100.)/len(vp), '% OK'
    print 'Viterbi  : ',(viterbi_correct*100.)/len(vp), '% OK'
#    print 'Posterior: ',(max_posterior_correct*100.)/len(vp), '% OK'
#    print 'LikeliH. : ',(ll_correct*100.)/len(vp), '% OK'
    
    if len(vp) != len(ground_truth):
        print "UNEXPECTED ERROR: different prediction lengths:", len(vp), len(ground_truth)
        return
    
    for i in [2]: #, 2, 4, 8, 16]:
        #precision and recall
        print "sensitivity: 1 /", i
        
        #recall Labeling
        #called_l, real, ok_called, not_called = computeEval(ground_truth, el, i, num_patt)
        #print "mplabel recall   : ", called_l, '/',  real, " => ", 
        #print "%0.3f" % (called_l*100./real), "%"
        #print ok_called
        #print not_called
        #computeStats(ok_called, not_called, "LR", num_patt)
        
        #recall Viterbi
        called_v, real, ok_called, not_called = computeEval(ground_truth, vp, i, num_patt)
        print "viterbi recall   : ", called_v, '/',  real, " => ", 
        print "%0.3f" % (called_v*100./real), "%"
        print ok_called
        print not_called
        computeStats(ok_called, not_called, "VR", num_patt)
        
#        #recall max posterior
#        called_p, real, ok_called, not_called = computeEval(ground_truth, pp, i, num_patt)       
#        print "mposter recall   : ", called_p, '/',  real, " => ",
#        print "%0.3f" % (called_p*100./real), "%"
#        print ok_called
#        print not_called
#        computeStats(ok_called, not_called, "PR", num_patt)
        
        #precision Labeling 
        #correct_l, claimed_l, ok_prediction, wr_prediction = computeEval(el, ground_truth, i, num_patt)
        #print "mplabel precision: ", correct_l , '/',  claimed_l, " => ", 
        #print "%0.3f" % (correct_l*100./claimed_l), "%"
        #print ok_prediction
        #print wr_prediction
        #computeStats(ok_prediction, wr_prediction, "LP", num_patt)
        
        #precision Viterbi 
        correct_v, claimed_v, ok_prediction, wr_prediction = computeEval(vp, ground_truth, i, num_patt)
        print "viterbi precision: ", correct_v , '/',  claimed_v, " => ", 
        print "%0.3f" % (correct_v*100./claimed_v), "%"
        print ok_prediction
        print wr_prediction
        computeStats(ok_prediction, wr_prediction, "VP", num_patt)
        
#        #precision max posterior
#        correct_p, claimed_p, ok_prediction, wr_prediction = computeEval(pp, ground_truth, i, num_patt)
#        print "mposter precision: ", correct_p , '/',  claimed_p, " => ", 
#        print "%0.3f" % (correct_p*100./claimed_p), "%"
#        print ok_prediction
#        print wr_prediction
#        computeStats(ok_prediction, wr_prediction, "PP", num_patt)
        
        #format for LaTeX tabular
        #print "%0.3f" % ((viterbi_correct*100.)/len(vp)), "\%"
        #print "%0.3f" % ((max_posterior_correct*100.)/len(vp)), "\%"
        #print "%0.3f" % (called_l*100./real), "\%"
        print "%0.3f" % (called_v*100./real), "\%"
#        print "%0.3f" % (called_p*100./real), "\%"
        #print "%0.3f" % (correct_l*100./claimed_l), "\%"
        print "%0.3f" % (correct_v*100./claimed_v), "\%"
#        print "%0.3f" % (correct_p*100./claimed_p), "\%"
        
        
def main():
    #parse command line arguments
    parser = argparse.ArgumentParser(description='Performs fetal CNV analysis from maternal plasma and phased parental data.')
    parser.add_argument('input', type=str, nargs=1, help='path to input file with allele counts in plasma and parental haplotypes')
    parser.add_argument('target', type=str, nargs=1, help='path to file with background truth - "target"')
    parser.add_argument('plasma', type=str, nargs=1, help='path to file with plasma sequencing DOC for all chromosomal positions')
    parser.add_argument('ref', type=str, nargs=1, help='path to file with reference plasma sequencing DOC for all chromosomal positions')
    parser.add_argument('seq', type=str, nargs=1, help='path to ref. genomic sequence in fasta format')
    parser.add_argument('--ff', type=float, help='fetal mixture ratio', default=-1.)
    parser.add_argument('--useCvrg', help='use coverage flag', action="store_true")
    args = parser.parse_args()
    
    in_file_name = args.input[0]
    target_file_name = args.target[0]
    plasma_doc_file = open(args.plasma[0], "r")
    ref_doc_file = open(args.ref[0], "r")
    seq_file = open(args.seq[0], "r")
    if args.ff > 0: mix = args.ff
    
    #print input info
    print "------------------------------------------"
    print "Running fCNV, input parameters:"
    print "input:", in_file_name
    print "target:", target_file_name
    print "plasma:", plasma_doc_file
    print "refDOC:", ref_doc_file
    print "seq:", seq_file
    print "--ff:", args.ff
    print "--useCvrg:", args.useCvrg
    print "------------------------------------------"
    os.system("hostname")
    
    #read the pre-processed input
    snp_positions, samples, M, P, MSC, PSC = readInput(in_file_name)
    
    #get genomic positions on the last lines of the pileup files to estimate the length of the chromosome
    with open(args.plasma[0], 'rb') as fh:
        fh.seek(-256, 2)
        last_pos_plasma = int(fh.readlines()[-1].decode().split(' ')[0])
        fh.close()
#    with open(args.ref[0], 'rb') as fh:
#        fh.seek(-256, 2)
#        last_pos_ref = int(fh.readlines()[-1].decode().split(' ')[0])
#        fh.close()
    chr_length = last_pos_plasma + 4742
    
    gc_sum = [0] * chr_length
    prefix_sum_plasma = [0] * chr_length
    prefix_count_plasma = [0] * chr_length
    prefix_sum_ref = [0] * chr_length
    prefix_count_ref = [0] * chr_length
    
    #get GC content prefix sums from the reference
    gen_pos = 0
    keep_reading = True
    while keep_reading:
        line = seq_file.readline().strip().upper()
        if len(line) == 0: break
        if line[0] == '>': continue
        for i in range(len(line)):
            gc_sum[gen_pos] = gc_sum[max(gen_pos - 1, 0)]
            if line[i] in 'GC': gc_sum[gen_pos] += 1
            gen_pos += 1
            if gen_pos >= chr_length:
                keep_reading = False
                break
    seq_file.close()
    
    last = 0
    while True: 
        line = plasma_doc_file.readline()
        if not line: break
        row = map(int, line.split(' '))
        if row[0] >= chr_length: break
        prefix_sum_plasma[row[0]] = prefix_sum_plasma[last] + row[1]
        prefix_count_plasma[row[0]] = prefix_count_plasma[last] + 1
        last = row[0]
    plasma_doc_file.close()
    
    last = 0
    while True: 
        line = ref_doc_file.readline()
        if not line: break
        row = map(int, line.split(' '))
        if row[0] >= chr_length: break
        prefix_sum_ref[row[0]] = prefix_sum_ref[last] + row[1]
        prefix_count_ref[row[0]] = prefix_count_ref[last] + 1
        last = row[0]
    ref_doc_file.close()
    
    #snp_positions = []
    ground_truth = []
    target_file = open(target_file_name, 'r')
    while True: 
        line = target_file.readline()
        if not line: break
        line = line.rstrip("\n").split("\t")
        #snp_positions.append(int(line[0]))
        ground_truth.append(int(line[-1]))
    target_file.close()
    
    fcnv = fcnvHMM.FCNV(None, None, False)
    mix, mix_median, ct = fcnv.estimateMixture(samples, M, P)
    print "Est. Mixture: ", mix, mix_median, '(', ct ,')',
    #mix = 0.13 #proportion of fetal genome in plasma
    if args.ff > 0: mix = args.ff
    print "used:", mix
    
    cnv_prior = None
    if args.useCvrg:
        cvrg = cvrgHMM.coverageFCNV(snp_positions, prefix_sum_plasma, prefix_count_plasma, prefix_sum_ref, prefix_count_ref, gc_sum)
        cvrg_posterior = cvrg.posteriorDecoding(mix)
#       byLL = cvrg.likelihoodDecoding(mix)
        del prefix_sum_plasma, prefix_count_plasma, prefix_sum_ref, prefix_count_ref, gc_sum
        
        cnv_prior = [ [0., 0., 0.] for x in range(len(snp_positions)) ]
        for pos in range(len(cvrg_posterior)):
            for cp_num_posterior in cvrg_posterior[pos]:
                cnv_prior[pos][cp_num_posterior[1]] = cp_num_posterior[0]
        del cvrg, cvrg_posterior
    
#        ll_state = []
#        ll_value = []
#        for x in byLL[pos]:   
#            ll_state.append(x[1])
#            ll_value.append(x[0])
#        
#        ll_str = ''
#        for j in range(len(ll_state)):
#            cvrg.logNormalize(ll_value)
#            ll_str += "%.8f"%math.exp(ll_value[j])+' '+str(ll_state[j])+' | '
#            
#        posterior_str = ''
#        for x in cvrg_posterior[pos]:
#            posterior_str += "%.8f"%math.exp(x[0])+' '+str(x[1])+' | '
#        posterior_str+='\t'

#        print snp_positions[pos], 'PP:', posterior_str, 'LL:', ll_str
#    del cvrg, cvrg_posterior, byLL
    
    del fcnv
    fcnv = fcnvHMM.FCNV(snp_positions, cnv_prior, args.useCvrg)
    
    ground_truth = []
    target_file = open(target_file_name, 'r')
    for line in target_file.readlines():
        line = line.rstrip("\n").split("\t")
        ground_truth.append(int(line[-1]))
    target_file.close()
    
    #res_file = open(out_file_name, 'w')
    file_name_prefix = target_file_name.split('/')[-1].split('.')[0].replace(':', '-')
    print "------------------ w/o TRAINING -------------------"
    test(fcnv, snp_positions, samples, M, P, MSC, PSC, mix, ground_truth, file_name_prefix)
    #test(fcnv, snp_positions[:1000], samples[:1000], M[:1000], P[:1000], MSC[:1000], PSC[:1000], mix, ground_truth[:1000], file_name_prefix)
    
    '''
    print "------------------ BEFORE TRAINING -------------------"
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TRAINING DATA: >>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples_t, M_t, P_t, mix_t, ground_truth_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TESTING DATA: >>>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples, M, P, mix, ground_truth)
    print "=================================================================================="
    
    
    def_trans = copy.deepcopy(fcnv.transitions)
    print "----------------- AFTER up to 5 rounds of B-W -------------------"
    fcnv.baumWelshForTransitions(samples_t, M_t, P_t, mix_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TRAINING DATA: >>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples_t, M_t, P_t, mix_t, ground_truth_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TESTING DATA: >>>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples, M, P, mix, ground_truth)
    print "=================================================================================="
    
    fcnv.transitions = copy.deepcopy(def_trans)
    print "----------- AFTER up to 5 rounds of Viterbi Training ------------"
    fcnv.viterbiTrainingForTransitions(samples_t, M_t, P_t, mix_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TRAINING DATA: >>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples_t, M_t, P_t, mix_t, ground_truth_t)
    print ">>>>>>>>>>>>>>>>>>>>>>>>> ON TESTING DATA: >>>>>>>>>>>>>>>>>>>>>>"
    test(fcnv, samples, M, P, mix, ground_truth)
    print "=================================================================================="
    '''
    
    #fcnv.baumWelshForTransitions(samples, M, P, mixture)
    #fcnv.viterbiTrainingForTransitions(samples, M, P, mixture)
    #print fcnv.computeForward(samples, M, P, mixture)
    #print fcnv.computeBackward(samples, M, P, mixture)
    #pp = fcnv.maxPosteriorDecoding(samples, M, P, mixture)
    
if __name__ == "__main__":
    #import doctest
    #doctest.testmod()
    
    main()
    



