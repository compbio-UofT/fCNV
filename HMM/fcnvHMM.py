import random
import math
#import numpypy as np
import itertools
import copy
import functools
from bisect import bisect_left

class memoized(object):
       """Decorator that caches a function's return value each time it is called.
       If called later with the same arguments, the cached value is returned, and
       not re-evaluated.
       """
       def __init__(self, func):
          self.func = func
          self.cache = {}
          
       def __call__(self, *args):
          try:
             return self.cache[args]
          except KeyError:
             value = self.func(*args)
             self.cache[args] = value
             return value
          except TypeError:
             # uncachable -- for instance, passing a list as an argument.
             # Better to not cache than to blow up entirely.
             print 'memo arg TypeError', args
             return self.func(*args)
             
       def __repr__(self):
          """Return the function's docstring."""
          return self.func.__doc__
          
       def __get__(self, obj, objtype):
          """Support instance methods."""
          fn = functools.partial(self.__call__, obj)
          fn.reset = self.reset
          return fn
          
       def reset(self):
          self.cache = {}

class HMMState(object):
    """Class to represent a particular phased inheritance pattern."""
    
    def __init__(self, ip, pp):
        """Initialize new HMMState object
        
        Arguments:
        ip -- tuple(int, int) -- general inheritance pattern
        pp -- tuple(tuple(ints), tuple(ints)) -- phased inheritance pattern
        
        """
        super(HMMState, self).__init__()
        
        self.inheritance_pattern = ip
        self.phased_pattern = pp
        
    def __str__(self):
        return "IP: {0}, PP: {1}".format(self.inheritance_pattern, self.phased_pattern)
        '''
        dic = ["1", "2"]
        nms = ["M", "P"]
        result = ""
        for h in range(2):
            for x in self.phased_pattern[h]:
                result += " " + nms[h] + dic[x]
            if h == 0: result += " |"
        return result
        '''
        
    def __eq__(self, other):
        return self.inheritance_pattern == other.inheritance_pattern and \
            self.phased_pattern == other.phased_pattern

class FCNV(object):
    """Hidden Markov Model for fetal CNV calling"""
    
    nucleotides = ['A', 'C', 'G', 'T']
    
    def __init__(self, positions, cnv_prior, use_prior):
        """Initialize new FCNV object"""
        super(FCNV, self).__init__()
        
        #store genomic positions of the SNPs
        self.positions = positions
        
        #CNV per position prior
        self.use_prior = use_prior
        self.cnv_prior = cnv_prior
        print "Use priors: ", self.use_prior
        
        #run intern tests
        self.neg_inf = float('-inf')
        self.inf = float('inf')
        self.nan = float('nan')
        self.testLogSum()
        
        #cache for log-likelihood values
        self.logLikelihoodCache = {}
        #self.distributionCache = {}
           
        '''#precompute lookup table for logSum function: 
        self.logTable1 = []
        for i in range(501):
            self.logTable1.append(math.log(1 + math.exp(-i/100.)) )
        self.logTable2 = []
        for i in range(51):
            self.logTable2.append(math.log(1 + math.exp(-i)) )
        '''
        
        #generate inheritance patterns
        self.inheritance_patterns = []
        for mats in range(3):
            for pats in range(3):
                if mats+pats in [0,4]: continue
                self.inheritance_patterns.append((mats, pats))
        
               
        #generate phased variants of the inheritance patterns and list of HMM States
        self.states = []
        self.phased_patterns = []
        self.max_component_size = 0
        for ip in self.inheritance_patterns:
            num_phased_states = 0
            for mPhased in itertools.combinations_with_replacement([0,1], ip[0]):
                for pPhased in itertools.combinations_with_replacement([0,1], ip[1]):
                    phased_pattern = (tuple(mPhased), tuple(pPhased))
                    self.phased_patterns.append(phased_pattern)
                    new_state = HMMState(ip, phased_pattern)
                    self.states.append(new_state)
                    num_phased_states += 1
            self.max_component_size = max(self.max_component_size, num_phased_states)
        
        #generate silent states
        self.states.append( HMMState("s", "s") ) #start state
        
        normalIP = (1,1)
        for ip in self.inheritance_patterns:
            if ip == normalIP: continue
            new_state = HMMState(ip, "out")
            self.states.append(new_state)
        
        self.states.append( HMMState(normalIP, "in") )
        self.states.append( HMMState(normalIP, "out") )
        
        for ip in self.inheritance_patterns:
            if ip == normalIP: continue
            new_state = HMMState(ip, "in")
            self.states.append(new_state)
        
        self.states.append( HMMState("t", "t") ) #end state
        
        #DEPRECATED
        self.main_state = 3
        self.main_states = []
        for i, s in enumerate(self.states):
            if s.inheritance_pattern == (1,1): self.main_states.append(i)
        
        #generate the table of transition probaboilities
        self.transitions = self.generateTransitionProb()

        #print the table
        num_states = self.getNumStates()
        #for k in xrange(num_states):
        #    for l in xrange(num_states):
        #        print "%.4f" % (math.exp(self.transitions[k][l])),
        #    print " "
            
        #create list of neighbors (adjacency list) from the transitions (matrix)
        self.adjacency_out = [ [] for x in range(num_states)]
        self.adjacency_in = [ [] for x in range(num_states)]
        for i in range(num_states):
            for j in range(num_states):
                if self.transitions[i][j] != self.neg_inf:
                    self.adjacency_out[i].append(j)
                    self.adjacency_in[j].append(i)
    
    @memoized     
    def getNumIP(self):
        """Return number of inheritance patterns"""
        return len(self.inheritance_patterns)
    
    @memoized     
    def getNumPP(self):
        """Return number of phased inheritance patterns"""
        return len(self.phased_patterns)
    
    @memoized    
    def getNumStates(self):
        """Return number of HMM states per one SNP position"""
        return len(self.states)   
    
    @memoized 
    def isNormal(self, state):
        return state.inheritance_pattern == (1,1)
    
    @memoized   
    def isReal(self, state):
        return isinstance(state.inheritance_pattern, tuple) and isinstance(state.phased_pattern, tuple)
    
    def areRecombination(self, state1, state2):
        if state1.inheritance_pattern != state2.inheritance_pattern:
            return False
        eq = 0
        for hap in range(2):
            eq += (state1.phased_pattern[hap] == state2.phased_pattern[hap])
        return (eq > 0)
    
    @memoized
    def getCorrespondingInState(self, state):
        num_real_states = self.getNumPP()
        for st_id, st in enumerate(self.states[num_real_states:]):
            if st.inheritance_pattern == state.inheritance_pattern and \
               st.phased_pattern == "in":
                return num_real_states + st_id, st
        return None
    
    @memoized    
    def getCorrespondingOutState(self, state):
        num_real_states = self.getNumPP()
        for st_id, st in enumerate(self.states[num_real_states:]):
            if st.inheritance_pattern == state.inheritance_pattern and \
               st.phased_pattern == "out":
                return num_real_states + st_id, st
        return None
    
    @memoized    
    def getStartState(self):
        for i, state in enumerate(self.states):
            if state.inheritance_pattern == "s": return i, state
        return None
    
    @memoized   
    def getExitState(self):
        for i, state in enumerate(self.states):
            if state.inheritance_pattern == "t": return i, state
        return None
    
    def areAdjacent(self, state1, state2):
        """Find out whether the two states are adjacent.
        Adjacency is defined as hamming distance <= 1.
        
        """
        #convert from tuples to lists
        s = list(state1.phased_pattern)
        t = list(state2.phased_pattern)
        for i in range(2):
            s[i] = list(s[i])
            t[i] = list(t[i])
          
        #enumerate the  neighborhood of s with hamming distance 1
        #and look for discovery of t
        match = (s == t)
        for h in range(2): #over both haplotypes
            for i, x in enumerate(s[h]): #over all alleles in a haplotype
                temp = copy.deepcopy(s)
                temp[h][i] = (x + 1) % 2 #flip this allele
                match += (temp == t)
                
                del temp[h][i] #delete this allele
                match += (temp == t)
            for x in range(2):
                temp = copy.deepcopy(s)
                temp[h].append(x) #add a new allele
                match += (temp == t)
                
        return (match >= 1)
    
    def generateTransitionProb(self):
        """Generate transition probabilities between HMM states"""
        #number of possible states per one SNP position
        num_states = self.getNumStates()
        num_real_states = self.getNumPP() 
        
        #precompute statistics
        self.num_recombs = [0 for x in range(self.getNumIP())]
        self.max_recombs = 0
        for ip_id, ip in enumerate(self.inheritance_patterns):
            for i, state1 in enumerate(self.states[:num_real_states]):
                if state1.inheritance_pattern != ip: continue
                #states "inside the IP component"
                self.num_recombs[ip_id] = 0
                for j, state2 in enumerate(self.states[:num_real_states]):
                    if i == j: #the exactly same state
                        self.num_recombs[ip_id] += 1
                    elif state2.inheritance_pattern == ip: #to recombination of the same IP
                        if state1.inheritance_pattern == (2, 1):
                            if state1.phased_pattern[0][0] ^ state1.phased_pattern[0][1] != state2.phased_pattern[0][0] ^ state2.phased_pattern[0][1]:
                                continue
                        if state1.inheritance_pattern == (1, 2):
                            if state1.phased_pattern[1][0] ^ state1.phased_pattern[1][1] != state2.phased_pattern[1][0] ^ state2.phased_pattern[1][1]:
                                continue
                        self.num_recombs[ip_id] += 1
                self.max_recombs = max(self.max_recombs, self.num_recombs[ip_id])
        
        
        trans = [[0. for x in range(num_states)] for y in range(num_states)]
        #first generate transitions from the real states
        for ip_id, ip in enumerate(self.inheritance_patterns):
            
            if ip == (1, 1): #generate Normal states transitions\
                pstay = 0.9799
                #precomb = 0.5 / (num_recombs[self.inheritance_patterns.index(ip)] - 1)
                #pstay = precomb = 0.999 / 6.
                precomb = 0.02 / (self.max_recombs - 1)
                pgo = 0.0001
                for i, state1 in enumerate(self.states[:num_real_states]):
                    if state1.inheritance_pattern != ip: continue
                    #states "inside the IP component"
                    for j, state2 in enumerate(self.states[:num_real_states]):
                        if i == j: #stay in the state
                            trans[i][i] = pstay
                        elif state2.inheritance_pattern == ip: #to recombination of the same IP
                            trans[i][j] = precomb
                    #to the silent exit node    
                    outState_id= self.getCorrespondingOutState(state1)[0]
                    trans[i][outState_id] = pgo
                    
            else: #generate CNV states transitions
                pstay = 0.979
                #precomb = 0.5 / (num_recombs[self.inheritance_patterns.index(ip)] - 1)
                #pstay = precomb = 0.999 / 6.
                precomb = 0.02 / (self.max_recombs - 1)
                pgo = 0.001
                for i, state1 in enumerate(self.states[:num_real_states]):
                    if state1.inheritance_pattern != ip: continue
                    #states "inside the IP component"
                    for j, state2 in enumerate(self.states[:num_real_states]):
                        if i == j: #stay in the state
                            trans[i][i] = pstay
                        elif state2.inheritance_pattern == ip: #to recombination of the same IP
                            if state1.inheritance_pattern == (2, 1):
                                if state1.phased_pattern[0][0] ^ state1.phased_pattern[0][1] != state2.phased_pattern[0][0] ^ state2.phased_pattern[0][1]:
                                    trans[i][j] = 0.
                                    continue
                            if state1.inheritance_pattern == (1, 2):
                                if state1.phased_pattern[1][0] ^ state1.phased_pattern[1][1] != state2.phased_pattern[1][0] ^ state2.phased_pattern[1][1]:
                                    trans[i][j] = 0.
                                    continue
                            trans[i][j] = precomb
                    
                    #to the silent exit node    
                    outState_id = self.getCorrespondingOutState(state1)[0]
                    trans[i][outState_id] = pgo
        
        #now generate transitions from the silent states
        inNormal_id = self.getCorrespondingInState( HMMState((1,1), ()) )[0]
        for i, state1 in enumerate(self.states[num_real_states:]):
            i += num_real_states
            #if it is the start node
            if state1.phased_pattern == "s":
                for j, state2 in enumerate(self.states[num_real_states:]):
                    if state2.phased_pattern == "in":
                        trans[i][num_real_states + j] = 1. / self.getNumIP()
            
            #if it is a silent exit node
            elif state1.phased_pattern == "out":
                prob = 1. / self.getNumIP()
                j = self.getExitState()[0]
                trans[i][j] = prob
                        
                if self.isNormal(state1):
                    #prob = 1. / self.getNumIP()
                    for j, state2 in enumerate(self.states[num_real_states:]):
                        if state2.phased_pattern == "in" and \
                         state2.inheritance_pattern != state1.inheritance_pattern \
                         and state2.inheritance_pattern != (0,2) \
                         and state2.inheritance_pattern != (2,0):
                            trans[i][num_real_states + j] = prob
                else:
                    trans[i][inNormal_id] = prob
                    
            #if it is a silent starting node
            elif state1.phased_pattern == "in":
                prob = 1. / self.max_component_size
                for j, state2 in enumerate(self.states[:num_real_states]):
                    if state2.inheritance_pattern == state1.inheritance_pattern:
                        trans[i][j] = prob
        
#        for i in range(num_states):
#            print i, self.states[i]
#        print '\n'
#        for i in range(num_states):
#            print "%2d"%i,
#            for j in range(26,42):
#                print "%.3f "%trans[i][j],
#            print '\n',
#        print '\n',
        
        #normalize and take logarithm
        for i in range(num_states):
            for j in range(num_states):
                if trans[i][j] < 10e-10: trans[i][j] = self.neg_inf
                else: trans[i][j] = math.log(trans[i][j])
            self.logPseudoNormalizeTrans(i, trans[i])
#            self.logNormalize(trans[i])
        
#        for i in range(num_states):
#            print "%2d"%i,
#            for j in range(26,42):
#                print "%.3f "%math.exp(trans[i][j]),
#            print '\n',
#        print 0/0
        
        return trans
    
    def logPseudoNormalizeTrans(self, state_id, l):
        """Normalize the given list of transition log-probabilities for state 'state_id'
        
        Normalizes the list in place and returns the used scale factor.
        
        """
        num_real_states = self.getNumPP()
        state1 = self.states[state_id]
        if self.isReal(state1):
            # 3types of transitions: stay, recomb, leave
            # special case: normalize to max number of recombs
            
            sum_ = self.neg_inf
            #recombs
            recomb_num = 0
            for j, state2 in enumerate(self.states[:num_real_states]):
                if state_id != j and state2.inheritance_pattern == state1.inheritance_pattern:
                    #recombination of the same IP
                    if state1.inheritance_pattern == (2, 1):
                        if state1.phased_pattern[0][0] ^ state1.phased_pattern[0][1] != state2.phased_pattern[0][0] ^ state2.phased_pattern[0][1]:
                            continue
                    if state1.inheritance_pattern == (1, 2):
                        if state1.phased_pattern[1][0] ^ state1.phased_pattern[1][1] != state2.phased_pattern[1][0] ^ state2.phased_pattern[1][1]:
                            continue
                    sum_ = self.logSum(sum_, l[j])
                    recomb_num += 1
            sum_ = sum_ + math.log(self.max_recombs-1) - math.log(recomb_num)
            
            #stay
            sum_ = self.logSum(sum_, l[state_id])
            #leave
            sum_ = self.logSum(sum_, l[self.getCorrespondingOutState(state1)[0]])
            #print state1, sum_, math.exp(sum_), recomb_num
            
            #normalize
            for i in range(len(l)): 
                if l[i] != self.neg_inf: l[i] -= sum_
            return sum_
        elif state1.phased_pattern == "in":
            # special case: normalize to max number of phased 'substates'
            sum_ = self.neg_inf
            num = 0
            for x in l:
                if x != self.neg_inf:
                    num += 1
                    sum_ = self.logSum(sum_, x)
            sum_ = sum_ + math.log(self.max_component_size) - math.log(num)
            for i in range(len(l)): 
                if l[i] != self.neg_inf: l[i] -= sum_
            return sum_
        else:
            return self.logNormalize(l)
        
    
    def logNormalize(self, l):
        """Normalize the given list of log-probabilities.
        
        Normalizes the list in place and returns the used scale factor.
        
        """
        sum_ = self.neg_inf
        for x in l: sum_ = self.logSum(sum_, x)
        for i in range(len(l)): 
            if l[i] != self.neg_inf: l[i] -= sum_
        return sum_
        
    def logSum(self, x, y):
        """Return sum of two numbers in log space
        i.e. equivalent to log(exp(x)+exp(y))
        
        """
        '''
        def precise(x):
            return math.log(1 + math.exp(x) )

        def lookup(x):
            #return math.log(1 + math.exp(x) )
            x = -x
                    
            if x < 5:
                x *= 100
                fx = int(math.floor(x))
                return (x-fx)*(self.logTable1[fx+1] - self.logTable1[fx]) + self.logTable1[fx]
            elif x < 50:
                fx = int(math.floor(x))
                return (x-fx)*(self.logTable2[fx+1] - self.logTable2[fx]) + self.logTable2[fx]
            else: return 0.
        '''
        
        if x == self.neg_inf and y == self.neg_inf: return self.neg_inf
        elif x == self.inf and y == self.inf: return self.inf
        elif math.isnan(x) or math.isnan(y): return self.nan
        
        a = max(x, y)
        b = min(x, y)
        return a + math.log(1 + math.exp(b-a) )
    
    def testLogSum(self):
        # From the numpy documentation
        prob1 = math.log(1e-50)
        prob2 = math.log(2.5e-50)
        prob12 = self.logSum(prob1, prob2)
        assert prob12 == -113.87649168120691
        assert math.exp(prob12) == 3.5000000000000057e-50
        assert self.logSum(0, 0) == math.log(2)
        assert self.logSum(float('-inf'), 0) == 0
        #assert self.logSum(12345678, 12345678) == float('inf')

        assert math.isnan(self.logSum(float('nan'), 1))
        assert math.isnan(self.logSum(1, float('nan')))
        assert math.isnan(self.logSum(float('nan'), float('inf')))
        assert math.isnan(self.logSum(float('inf'), float('nan')))
        assert self.logSum(float('-inf'), float('-inf')) == float('-inf')
        assert self.logSum(float('-inf'), float('inf')) == float('inf')
        assert self.logSum(float('inf'), float('-inf')) == float('inf')
        assert self.logSum(float('inf'), float('inf')) == float('inf')
    
    #memoized
    def allelesDistribution(self, maternal, fetal, mix):
        """Compute nucleotides distribution in maternal plasma for a position with
        given maternal and fetal alleles, assuming the given fetal admixture.
        
        Arguments:
        maternal -- maternal alleles
        paternal -- paternal alleles
        mix -- fetal admixture ratio
        returns list(floats) -- list of nucleotides probabilities
        
        >>> f = FCNV()
        >>> f.allelesDistribution(['A', 'A'], ['T', 'A'], 0.1)
        [0.9230769230769231, 0.009615384615384616, 0.009615384615384616, 0.057692307692307696]
        """
        '''#Caching:
        code = (tuple(maternal), tuple(fetal), mix)
        val = self.distributionCache.get(code)
        if val != None: return val
        '''

        adjusted_fetal_admix = mix/2. * len(fetal)
        adjusted_maternal_admix = (1.-mix)/ 2. * len(maternal)
        cmix = adjusted_fetal_admix / (adjusted_maternal_admix + adjusted_fetal_admix)
        dist = []
        for nuc in self.nucleotides:
            p = maternal.count(nuc) / float(len(maternal)) * (1.-cmix)
            p += fetal.count(nuc) / float(len(fetal)) * (cmix)
            #if p < 0.01: p = 0.01
            dist.append(p)
            
        #normalize
        #summ = sum(dist)
        #dist = [dist[i] / summ for i in range(len(dist)) ]
        
        '''self.distributionCache[code] = dist'''
        return dist
    
    def allelesDistributionExtended(self, Malleles, Mcounts, Palleles, Pcounts, pattern, mix):
        """Compute nucleotides distribution in maternal plasma for a position with
        given maternal and fetal alleles, assuming the given fetal admixture.
        
        Arguments:
        maternal -- maternal alleles
        paternal -- paternal alleles
        mix -- fetal admixture ratio
        returns list(floats) -- list of nucleotides probabilities
        """

#        Falleles = []
#        #maternaly inherited
#        for mHpatt in pattern[0]:
#            Falleles.append(Malleles[mHpatt])
#        #paternaly inherited
#        for pHpatt in pattern[1]:
#            Falleles.append(Palleles[pHpatt])
            
        numFalleles = float(len(pattern[0]) + len(pattern[1]))
        numMinherited = len(pattern[0])
        numPinherited = len(pattern[1])
        numMalleles = float(len(Malleles))
        sumMcount = float(sum(Mcounts))
        sumPcount = float(sum(Pcounts))
        
        #adjust mixture to currently considered event (deletion, duplication, normal)
        adjusted_fetal_admix = mix/2. * numFalleles
        adjusted_maternal_admix = (1.-mix)/2. * numMalleles
        cmix = adjusted_fetal_admix / (adjusted_maternal_admix + adjusted_fetal_admix)
        
        dist = {}
        for nuc in self.nucleotides: dist[nuc] = 0.
        alphaM = 32/2. + sum(Mcounts)/2.
        alphaP = 39/2. + sum(Pcounts)/2.
        
        #fraction that is trully from mother DNA
        for i, nuc in enumerate(Malleles):
            dist[nuc] += (alphaM + Mcounts[i]) / float(2*alphaM + sumMcount) * (1.-cmix)
            
        #fraction that is fetal but maternaly inherited
        for mHpatt in pattern[0]:
            nuc = Malleles[mHpatt]
            if numMinherited == 2 and pattern[0][0] != pattern[0][1]:
                #if both haplotypes are inherited, use the maternal seq. ratio
                dist[nuc] += (alphaM + Mcounts[mHpatt]) / float(2*alphaM + sumMcount) * (cmix*2./numFalleles)
            else:
                dist[nuc] += cmix / numFalleles
            
        #fraction that is fetal but paternaly inherited
        for pHpatt in pattern[1]:
            nuc = Palleles[pHpatt]
            # dist[nuc] += cmix / numFalleles
            # use paternal ratio correction
            if numPinherited == 2 and pattern[1][0] != pattern[1][1]:
                #if both haplotypes are inherited, use the paternal seq. ratio
                dist[nuc] += (alphaP + Pcounts[pHpatt]) / float(2*alphaP + sumPcount) * (cmix*2./numFalleles)
            else:
                dist[nuc] += cmix / numFalleles
        
        dist_list = [dist[nuc] for nuc in self.nucleotides]
        #print pattern, Malleles, Mcounts, Palleles, Pcounts,':', dist_list, sum(dist_list)
        
        #normalize
        #summ = float(sum(dist_list))
        #dist_list = [dist_list[i] / summ for i in range(len(dist_list)) ]

        return dist_list
    
    def logGaussian(self, x, mus, cov_diagonal):
        '''
        log probability of x ~ N(mu, cov)
        '''
        D = len(x)
        N = sum(x)
        xm = [ float(x[i] - mus[i]) for i in range(D) ]
        
        det = 1.
        for i in range(D): det *= cov_diagonal[i]
            
        px = - D/2. * math.log(2.*math.pi)
        px += - math.log(det) /2.
        
        #compute the 'exponent' part #px += -(np.dot(np.dot(xm.T, inv), xm)) /2.
        part1 = [ xm[i] / cov_diagonal[i] for i in range(D) ]
        part2 = sum([ part1[i] * xm[i] for i in range(D) ])
        px += - part2 /2.
        
        return px
    
    def logLHGivenStateWCoverage(self, pos_ind, nuc_counts, maternal_alleles, paternal_alleles,\
           maternal_sq_counts, paternal_sq_counts, mix, state):
        
        if state.inheritance_pattern == (2, 0) or state.inheritance_pattern == (0, 2): return -50
        
        #pattern = state.phased_pattern
        #mus = self.allelesMeans(nuc_counts, maternal_alleles, paternal_alleles, maternal_sq_counts, paternal_sq_counts, pattern, mix)
        
        pattern = state.phased_pattern
#        fetal_alleles = []
#        fetal_pseudocounts = []
#        for mHpatt in pattern[0]:
#            fetal_alleles.append(maternal_alleles[mHpatt])
#            fetal_pseudocounts.append(maternal_sq_counts[mHpatt])
#        for pHpatt in pattern[1]:
#            fetal_alleles.append(paternal_alleles[pHpatt])
#            fetal_pseudocounts.append(paternal_sq_counts[pHpatt])
#        ps = self.allelesDistribution(maternal_alleles, tuple(fetal_alleles), mix)
        ps = self.allelesDistributionExtended(maternal_alleles, maternal_sq_counts, paternal_alleles, paternal_sq_counts, pattern, mix)

        N = sum(nuc_counts)
        avg_doc = 70. #TODO: make this a parameter
        norm_coef = avg_doc / N
        
        mus = [ avg_doc * ps[i] for i in range(4) ]
        cov_diagonal = [ max(0.8, (N * ps[x]) * norm_coef**2) for x in range(4)]
        nuc_counts = [ x * norm_coef for x in nuc_counts ]
        
        ratios_prob = self.logGaussian(nuc_counts, mus, cov_diagonal)
        result = ratios_prob

        #print sum(nuc_counts), ratios_prob
        if result < -15: result = -15
        
        return result
    
    def estimateMixture(self, samples, M, P):
        """Estimate mixture of fetal genome in the plasma samples.
        
        Arguments:
        samples -- number of nucleotides sampled per individual SNP positions
        M -- maternal alleles
        P -- parental alleles
        returns float
        
        """
        mix = 0.
        num = 0
        estimates = []
        for i in xrange(len(samples)):
            #if mother is homozygous A and father has at least one alternative B
            if M[i][0] == M[i][1] and (P[i][0] != M[i][0] and P[i][1] != M[i][0]):
                #print M[i], P[i], samples[i]            
                type1 = M[i][0]
                type2 = P[i][0]
                if type1 == type2: type2 = P[i][1]
                try:
                    num_type2 = samples[i][self.nucleotides.index(type2)]
                except: 
                    print "ERROR unexpected allele: ", type2
                    
                #only if it looks like the child inherited the allele B
                if num_type2 != 0: 
                    sum_all = sum(samples[i])
                    num += 1
                    mix += (2.*num_type2)/sum_all
                    estimates.append((2.*num_type2)/sum_all)
                else:
                    num += 1.
                    estimates.append(0.)

        estimates.sort()
        return mix/num, estimates[len(estimates)/2], len(estimates)
        
    def estimateCoverage(self, samples):
        """Estimate coverage of the SNP loci in the plasma samples as empirical mean.
        """
        result = sum(map(sum, samples)) / float(len(samples))
        return result
    
    '''
    def adjustInTransitionProbForPos(self, pos, trans, inplace=False):
        """Multiply the CNV prior into transition probabilities of incoming 
        edges of the 'in' silent states at the specified position
        """
        return
        #if not inplace: trans = copy.deepcopy(self.transitions)
        num_states = self.getNumStates()
        pos = min(pos, len(self.cnv_prior) - 1)
        for state_id, state in enumerate(self.states):
            if state.phased_pattern == "in":
                for prev_id in range(num_states):
                    if trans[prev_id][state_id] == self.neg_inf: continue
                    trans[prev_id][state_id] += self.cnv_prior[pos][sum(state.inheritance_pattern) - 1]
        for st_id in range(num_states):
            self.logNormalizeTrans(st_id, trans[st_id])
    
    def adjustOutTransitionProbForPos(self, pos, trans, inplace=False):
        """Multiply the negation of CNV prior into transition probabilities
        of incoming edges of the 'out' silent states at the specified position
        """
        return 
        #if not inplace: trans = copy.deepcopy(self.transitions)
        num_states = self.getNumStates()
        pos = min(pos, len(self.cnv_prior) - 1)
        for state_id, state in enumerate(self.states):
            if state.phased_pattern == "out":
                for prev_id in range(num_states):
                    if trans[prev_id][state_id] == self.neg_inf: continue
                    trans[prev_id][state_id] -= math.log(1. - math.exp(self.cnv_prior[pos][sum(state.inheritance_pattern) - 1]))
                    if trans[prev_id][state_id] == self.inf:
                        trans[prev_id][state_id] = self.neg_inf
#            if self.isReal(state):
#                outState_id = self.getCorrespondingOutState(state)[0]
#                for next_id in range(num_states):
#                    if trans[state_id][next_id] == self.neg_inf or next_id == outState_id: continue
#                    trans[state_id][next_id] += math.log(1. - math.exp(self.cnv_prior[pos][sum(state.inheritance_pattern) - 1]))
#                    if trans[state_id][next_id] == self.inf:
#                        trans[state_id][next_id] = self.neg_inf
        for st_id in range(num_states):
            self.logNormalizeTrans(st_id, trans[st_id])
    '''
            
    def adjustTransitionProbForPos(self, pos, trans):
        """Multiply the CNV prior into transition probabilities at the specified position
        """
        if self.use_prior == False: return
        num_real_states = self.getNumPP()
        num_states = self.getNumStates()
        pos = min(pos, len(self.cnv_prior) - 1)
        pos = max(pos, 0)
        logP1n3 = self.logSum(self.cnv_prior[pos][0], self.cnv_prior[pos][2])
        for state_id, state in enumerate(self.states):
            if self.isReal(state):
                #inState_id = self.getCorrespondingInState(state)[0]
                for prev_id in range(num_real_states):
                    #if prev_id == inState_id: continue
                    if trans[prev_id][state_id] == self.neg_inf: continue
                    trans[prev_id][state_id] += self.cnv_prior[pos][sum(state.inheritance_pattern) - 1]
                    if trans[prev_id][state_id] == self.inf:
                        trans[prev_id][state_id] = self.neg_inf
            
            if state.phased_pattern == "out":
                for prev_id in range(num_real_states):
                    if trans[prev_id][state_id] == self.neg_inf: continue
                    if state.inheritance_pattern == (1, 1):
                        trans[prev_id][state_id] += logP1n3 #prior of going to a CNV
                    else:
                        trans[prev_id][state_id] += self.cnv_prior[pos][1] #prior of going to (1, 1)
                    if trans[prev_id][state_id] == self.inf:
                        trans[prev_id][state_id] = self.neg_inf
            
            if state.phased_pattern == "in" and state.inheritance_pattern != (1, 1):
                for prev_id in range(num_real_states, num_states):
                    if trans[prev_id][state_id] == self.neg_inf: continue
                    trans[prev_id][state_id] += self.cnv_prior[pos][sum(state.inheritance_pattern) - 1]
                    if trans[prev_id][state_id] == self.inf:
                        trans[prev_id][state_id] = self.neg_inf
       
        for st_id in range(num_states):
            self.logPseudoNormalizeTrans(st_id, trans[st_id])
#            self.logNormalize(trans[st_id])
            
            
    def viterbiPath(self, samples, M, P, MSC, PSC, mixture):
        """
        Viterbi decoding of the most probable path.
        """
        num_states = self.getNumStates()
        num_real_states = self.getNumPP() 
        transitions = copy.deepcopy(self.transitions)
        self.avgCoverage = self.estimateCoverage(samples)
        
        n = len(samples)
        predecessor = [[0 for i in range(num_states)] for j in xrange(n+1)] 
        #DP table dim: seq length+1 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)] 
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #multiply the CNV prior into transition probabilities
        self.adjustTransitionProbForPos(0, transitions)
        #the start state probability is 1 -> log(1)=0
        table[0][start_id] = 0. 
        #propagate over all silent states
        for state_id, state in enumerate(self.states[num_real_states:]):
            state_id += num_real_states
            if transitions[start_id][state_id] == self.neg_inf: continue #just for speed-up
            table[0][state_id] = table[0][start_id] + transitions[start_id][state_id]
            predecessor[0][state_id] = -1
         
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case
            
            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(self.states[:num_real_states]):
                #emission probability in the given state
                emis_p = self.logLHGivenStateWCoverage(pos-1, samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, state)
                
                #porobability of this state is `emission * max{previous state * transition}`
                max_prev = self.neg_inf
                arg_max = -1
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp = table[pos-1][prev_id] + transitions[prev_id][state_id]
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                table[pos][state_id] = max_prev + emis_p
                predecessor[pos][state_id] = arg_max
            
            #multiply the CNV prior into transition probabilities
            transitions = copy.deepcopy(self.transitions)
            self.adjustTransitionProbForPos(pos, transitions)
            
            #(ii, iii) transitions from 'real' states and silent states with lower id to silent states
            #note: the states in self.states are already ordered such that a transistion from a silent 
            #   state to another silent state always points to a state with heigher index
            for state_id, state in enumerate(self.states[num_real_states:]):
                state_id += num_real_states
                #porobability of this state is `max{(actual real state or silent state with lower id) * transition}`
                max_prev = self.neg_inf
                arg_max = -1
                for prev_id in range(state_id + 1):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp = table[pos][prev_id] + transitions[prev_id][state_id]
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                table[pos][state_id] = max_prev
                predecessor[pos][state_id] = arg_max
            
            #normalize to sum to 1
            self.logNormalize(table[pos])
        
        '''RESULT'''
        path = [-1 for i in xrange(n+1)]
        path[n] = exit_id 
        while path[n] >= num_real_states: path[n] = predecessor[n][path[n]]
        for i in reversed(xrange(n)):
            path[i] = predecessor[i+1][path[i+1]]
            while path[i] >= num_real_states: path[i] = predecessor[i][path[i]]
        
        state_path = [[] for i in xrange(n+1)]
        for i in xrange(1, n+1):
            state = self.states[ path[i] ]
            path[i] = self.inheritance_patterns.index(state.inheritance_pattern)
            state_path[i] = (state.inheritance_pattern, state.phased_pattern)
        
        print "Viterbi path probability: ", table[n][exit_id]    
        return path[1:], state_path[1:]
        
    '''    
    def preCompute(self, v, u, j, i, samples, M, P, mixture):
        code = (v, u, j)
        if i == j:
            self.PCache[code] = float(v == u)
        else:
            num_real_states = self.getNumPP()
            transitions = self.transitions
            states = self.states
            logSum = self.logSum
            
            #for silent *v*
            if v >= num_real_states:
                summ = self.neg_inf
                for w in range(num_real_states):
                    if transitions[v][w] == self.neg_inf: continue
                    if states[v].inheritance_pattern != \
                       states[w].inheritance_pattern: continue
                    tmp = transitions[v][w] + self.PCache.get((w, u, j))
                    summ = logSum(summ, tmp)
                result = summ
                
            #for real *v*
            else:
                emis_p = self.logLHGivenState(samples[j-1], M[j-1], P[j-1], mixture, states[v])
                summ = self.neg_inf
                for w in range(num_real_states):
                    if transitions[v][w] == self.neg_inf: continue
                    #if states[v].inheritance_pattern != \
                    #   states[w].inheritance_pattern: continue
                    tmp = transitions[v][w] + self.PCache.get((w, u, j+1))
                    summ = logSum(summ, tmp)
                result = summ + emis_p
            self.PCache[code] = result
    '''
    @memoized
    def getPRec(self, v, u, j, i):
        transitions = self.transitions
        logSum = self.logSum
        adjacency_out = self.adjacency_out
        num_real_states = self.getNumPP()
        states = self.states
        
        if states[v].inheritance_pattern != states[u].inheritance_pattern: return self.neg_inf
        #base case
        if j == i:
            if u == v and v < num_real_states:
                return self.logLHGivenState(self.samples[j-1], self.M[j-1], self.P[j-1], self.mixture, states[u])
            elif v >= num_real_states:
                return transitions[v][u] + self.logLHGivenState(self.samples[j-1], self.M[j-1], self.P[j-1], self.mixture, states[u])
            else:
                return self.neg_inf
        
        #if a silent state
        if v >= num_real_states:
            summ = self.neg_inf
            for w in adjacency_out[v]:
                if states[v].inheritance_pattern != \
                   states[w].inheritance_pattern: continue
                tmp = transitions[v][w] + self.getPRec(w, u, j, i)
                summ = logSum(summ, tmp)
            return summ
                
        #real state
        emis_p = self.logLHGivenState(self.samples[j-1], self.M[j-1], self.P[j-1], self.mixture, states[v])
        summ = self.neg_inf
        for w in adjacency_out[v]:
            if states[v].inheritance_pattern != \
               states[w].inheritance_pattern: continue
            tmp = emis_p + transitions[v][w] + self.getPRec(w, u, j+1, i)
            summ = logSum(summ, tmp)
        return summ
        
    
    def computeValueFor(self, v, u, j, i):
        #caching:
        code = (v, u, j)
        PCache = self.PCache
        try: return PCache[code]
        except KeyError: pass
        
        transitions = self.transitions
        logSum = self.logSum
        adjacency_out = self.adjacency_out
        num_real_states = self.getNumPP()
        
        emis_p = self.logLHGivenState(self.samples[j-1], self.M[j-1], self.P[j-1], self.mixture, self.states[v])
        summ = self.neg_inf
        for w in adjacency_out[v]:
            if w >= num_real_states: continue
            #if transitions[v][w] == self.neg_inf: continue
            #if states[v].inheritance_pattern != \
            #   states[w].inheritance_pattern: continue
            if j+1 == i:
                if u == w and u < num_real_states:
                    probP = self.logLHGivenState(self.samples[j], self.M[j], self.P[j], self.mixture, self.states[u])
                else:
                    probP = self.neg_inf
            else:
                probP = PCache[(w, u, j+1)]
            tmp = transitions[v][w] + probP #self.getP(w, u, j+1, i, samples, M, P, mixture)
            summ = logSum(summ, tmp)
        result = summ + emis_p
        
        PCache[code] = result
        return result
        
    
    def getP(self, v, u, j, i):
        #caching:
        code = (v, u, j)
        PCache = self.PCache
        try: return PCache[code]
        except KeyError: pass
        
        num_real_states = self.getNumPP()
        transitions = self.transitions
        adjacency_out = self.adjacency_out
        
        #base case
        if j == i:
            if u == v and v < num_real_states:
                result = self.logLHGivenState(self.samples[j-1], self.M[j-1], self.P[j-1], self.mixture, self.states[u])
            elif v >= num_real_states:
                result = transitions[v][u] + self.logLHGivenState(self.samples[j-1], self.M[j-1], self.P[j-1], self.mixture, self.states[u])
            else:
                result = self.neg_inf
            PCache[code] = result
            return result
        
        #for silent *v*
        if v >= num_real_states:
            '''print "X", 
            summ = self.neg_inf
            for w in adjacency_out[v]:
                #if transitions[v][w] == self.neg_inf: continue
                #if states[v].inheritance_pattern != \
                #   states[w].inheritance_pattern: continue
                if w < num_real_states: tmp = transitions[v][w] + self.computeValueFor(w, u, j, i)
                else: tmp = transitions[v][w] + PCache[(w, u, j)]
                summ = self.logSum(summ, tmp)
            result = summ
            '''
            result = self.neg_inf            
            
        #for real *v*
        else:
            result = self.computeValueFor(v, u, j, i)

        PCache[code] = result
        return result
    
    def extendedLabeling(self, samples, M, P, mixture):
        """
        Most probable extended labeling.
        """
        #import sys
        #print >> sys.stderr, "BEGIN"
        self.avgCoverage = self.estimateCoverage(samples)
        self.samples = samples
        self.M = M
        self.P = P
        self.mixture = mixture
        num_states = self.getNumStates()
        num_real_states = self.getNumPP() 
        transitions = self.transitions
        adjacency_out = self.adjacency_out
        adjacency_in = self.adjacency_in
        states = self.states
        getP = self.getP
        
        
        n = len(samples)
        predecessor = [[0 for i in range(num_states)] for j in xrange(n+1)]
        jump = [[0 for i in range(num_states)] for j in xrange(n+1)]  
        #DP table dim: seq length+1 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)] 
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #the start state probability is 1 -> log(1)=0
        table[0][start_id] = 0. 
        #propagate over all silent states
        for state_id, state in enumerate(self.states[num_real_states:]):
            state_id += num_real_states
            if transitions[start_id][state_id] == self.neg_inf: continue #just for speed-up
            table[0][state_id] = table[0][start_id] + transitions[start_id][state_id]
            predecessor[0][state_id] = start_id
            jump[0][state_id] = 0
         
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case
            #self.PCache = {}
            self.getPRec.reset()
            #print >> sys.stderr, pos, 
            
            '''
            preCompute = self.preCompute
            for j in reversed(xrange(1,  pos+1)):
                i = pos
                for v in xrange(num_states):
                    for u in xrange(num_real_states):
                        preCompute(v, u, j, i, samples, M, P, mixture)
            '''
            
            
            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(states[:num_real_states]):
                max_prev = self.neg_inf
                arg_max = -1
                jump_from = -1
                
                '''
                #if the pravious transition is the critical edge:
                #emission probability in current state
                emis_p = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], mixture, state)
                #porobability of this state is `emission * max{previous state * transition}`
                for prev_id in adjacency_in[state_id]: #range(num_states)
                    #if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp = table[pos-1][prev_id] + transitions[prev_id][state_id] + emis_p
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                        jump_from = pos-1             
                '''
                
                #where the previous critical edge is:
                for j in reversed(xrange(1, pos+1)):
                    #between which states the critical edge is:
                    #if abs(pos-j)>200 : break
                    
                    #the 'jumped' region after critical edge must have one lable
                    to_id = self.getCorrespondingInState(states[state_id])[0]
                           
                    paths_prob = self.getPRec(to_id, state_id, j, pos) #getP(to_id, state_id, j, pos)
                    #backup_prob = self.getPRec(to_id, state_id, j, pos)
                    #if abs(paths_prob - backup_prob)>1e-4: print "PRUSEEEERRRR", to_id, state_id, j, pos, " :", paths_prob, backup_prob
                    for from_id in adjacency_in[to_id]:
                        #to be a critical edge, it must change the label
                        if states[from_id].inheritance_pattern == \
                           states[to_id].inheritance_pattern: continue
                        
                        #tmp = table[j-1][from_id] + transitions[from_id][to_id] + self.PCache.get((to_id, state_id, j))
                        tmp = table[j-1][from_id] + transitions[from_id][to_id] + paths_prob
                        if tmp > max_prev:
                            max_prev = tmp
                            arg_max = from_id
                            jump_from = j-1
                
                #record the max value and traceback info
                table[pos][state_id] = max_prev
                predecessor[pos][state_id] = arg_max
                jump[pos][state_id] = jump_from
            
            #(ii, iii) transitions from 'real' states and silent states with lower id to silent states
            #note: the states in self.states are already ordered such that a transistion from a silent 
            #   state to another silent state always points to a state with heigher index
            for state_id, state in enumerate(states[num_real_states:]):
                state_id += num_real_states
                #porobability of this state is `max{(actual real state or silent state with lower id) * transition}`
                max_prev = self.neg_inf
                arg_max = -1
                jump_from = -1
                for prev_id in adjacency_in[state_id]: #prev_id in range(state_id + 1):
                    #if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    if prev_id > state_id: continue #update silent states in increasing ID order
                    tmp = table[pos][prev_id] + transitions[prev_id][state_id]
                    if tmp > max_prev:
                        max_prev = tmp
                        arg_max = prev_id
                        jump_from = pos
                table[pos][state_id] = max_prev
                predecessor[pos][state_id] = arg_max
                jump[pos][state_id] = jump_from
            
            #normalize to sum to 1
            #self.logNormalize(table[pos])
        
        '''RESULT'''
        #print >> sys.stderr, "COMPUTING BACKTRACE"
        path = [-1 for i in xrange(n+1)]
        
        path[n] = exit_id 
        pos = n
        while pos > 0:
            #print >> sys.stderr, "pos:", pos, "jump:", jump[pos][path[pos]]
            new_pos = jump[pos][path[pos]]
            path[new_pos] = predecessor[pos][path[pos]]
            pos = new_pos
            while path[pos] >= num_real_states and path[pos]!=start_id:
                #print >> sys.stderr, "pos:", pos, ": ", path[pos], "->", predecessor[pos][path[pos]]
                path[pos] = predecessor[pos][path[pos]]
            #print >> sys.stderr, path[pos]
            for x in range(jump[pos][path[pos]]+1, pos):
                path[x] = path[pos]
                #print >> sys.stderr, "path", x, ": ", path[x]
        
        #print path[1:]
        for i in xrange(1, n+1):
            state = states[ path[i] ]
            path[i] = self.inheritance_patterns.index(state.inheritance_pattern)
        
        #print path[1:]
        print "labeling probability: ", table[n][exit_id]
        return path[1:]
    
    def computeForward(self, samples, M, P, MSC, PSC, mixture):
        '''
        Posterior decoding: forward algorithm
        '''
        num_states = self.getNumStates()
        num_real_states = self.getNumPP()
        patterns = self.inheritance_patterns
        transitions = copy.deepcopy(self.transitions)
        self.avgCoverage = self.estimateCoverage(samples)
        
        n = len(samples)
        #scale factors
        scale = [0. for i in xrange(n+1)]
        #DP table dim: seq length+1 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)] 
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #multiply the CNV prior into transition probabilities
        self.adjustTransitionProbForPos(0, transitions)
        #the start state probability is 1 -> log(1)=0
        table[0][start_id] = 0. 
        #propagate over all silent states
        for state_id, state in enumerate(self.states[num_real_states:]):
            state_id += num_real_states
            if transitions[start_id][state_id] == self.neg_inf: continue #just for speed-up
            table[0][state_id] = table[0][start_id] + transitions[start_id][state_id]
        
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case
            
            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(self.states[:num_real_states]):
                #emission probability in the given state
                emis_p = self.logLHGivenStateWCoverage(pos-1, samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, state)
                
                #probability of this state is `emission * \sum {previous state * transition}`
                summ = self.neg_inf
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp =  table[pos-1][prev_id] + transitions[prev_id][state_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = emis_p + summ
            
            #multiply the CNV prior into transition probabilities
            transitions = copy.deepcopy(self.transitions)
            self.adjustTransitionProbForPos(pos, transitions)
            
            #(ii, iii) transitions from 'real' states and silent states with *lower* id to silent states
            #note: the states in self.states are already ordered such that a transistion from a silent 
            #   state to another silent state always points to a state with heigher index
            for state_id, state in enumerate(self.states[num_real_states:]):
                state_id += num_real_states
                
                #probability of this state is `\sum {(actual real state or silent state with lower id) * transition}`
                summ = self.neg_inf
                for prev_id in range(state_id):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp =  table[pos][prev_id] + transitions[prev_id][state_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = summ
                
            #normalize to sum to 1
            scale[pos] = self.logNormalize(table[pos])
        
        '''RESULT'''
        #p(X) - probability of the observed sequence (samples)
        pX = table[n][exit_id]
        #print "fwd: ", math.exp(pX), pX

        return table, pX, scale
        
        
    def computeBackward(self, samples, M, P, MSC, PSC, mixture, scale):
        '''
        Posterior decoding: backward algorithm
        '''
        num_states = self.getNumStates()
        num_real_states = self.getNumPP()
        patterns = self.inheritance_patterns
        transitions = copy.deepcopy(self.transitions)
        self.avgCoverage = self.estimateCoverage(samples)
        
        n = len(samples)
        #DP table dim: seq length+2 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)]
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
        #multiply the CNV prior into transition probabilities
        self.adjustTransitionProbForPos(n-1, transitions)
        #the exit state probability is 1 -> log(1)=0
        table[n][exit_id] = 0. 
        #(iii) add transitions from silent to silent states with *higher* id (*no emission*)
        for state_id in reversed(range(num_real_states, num_states)):
            #\sum {'next'_state * transition_from_here} (no emission)
            summ = self.neg_inf
            for next_id in reversed(range(state_id, num_states)):
                if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                tmp = transitions[state_id][next_id] + table[n][next_id]
                summ = self.logSum(summ, tmp)
            table[n][state_id] = self.logSum(summ, table[n][state_id])
        
        #(ii) add transitions from 'real' states to silent states (*no emission*)
        for state_id in range(num_real_states):
            #\sum {'next'_state * transition_from_here} (no emission)
            summ = self.neg_inf
            for next_id in range(num_real_states, num_states):
                if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                tmp = transitions[state_id][next_id] + table[n][next_id]
                summ = self.logSum(summ, tmp)
            table[n][state_id] = self.logSum(summ, table[n][state_id])
        
        '''DYNAMIC PROGRAMMING'''    
        #for all SNP positions do:
        for pos in reversed(xrange(0, n)):
            #real positions are <1..n+1), the 0 and n+1 are sentinels
            
            #precompute emission probabilities
            emis_p = []
            for state in self.states[:num_real_states]:
                emis_p.append(self.logLHGivenStateWCoverage(pos, samples[pos], M[pos], P[pos], MSC[pos], PSC[pos], mixture, state))
                
            #multiply the CNV prior into transition probabilities
            #transitions = copy.deepcopy(self.transitions)
            #self.adjustTransitionProbForPos(pos, transitions)
            
            #(i) transitions from all states to 'real' states of the HMM:
            for state_id in range(num_states):
                # \sum {emission_in_'next' * 'next'_state * transition_from_here}
                summ = self.neg_inf
                for next_id in range(num_real_states):
                    if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                    tmp = transitions[state_id][next_id] + table[pos+1][next_id] + emis_p[next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = summ
            
            #(iii) add transitions from silent to silent states with *higher* id (*no emission*)
            for state_id in reversed(range(num_real_states, num_states)):
                #\sum {'next'_state * transition_from_here} (no emission)
                summ = self.neg_inf
                for next_id in reversed(range(state_id, num_states)):
                    if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                    tmp = transitions[state_id][next_id] + table[pos][next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = self.logSum(summ, table[pos][state_id])
            
            #(ii) add transitions from 'real' states to silent states (*no emission*)
            for state_id in range(num_real_states):
                #\sum {'next'_state * transition_from_here} (no emission)
                summ = self.neg_inf
                for next_id in range(num_real_states, num_states):
                    if transitions[state_id][next_id] == self.neg_inf: continue #just for speed-up
                    tmp = transitions[state_id][next_id] + table[pos][next_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = self.logSum(summ, table[pos][state_id])
            
            #multiply the CNV prior into transition probabilities
            transitions = copy.deepcopy(self.transitions)
            self.adjustTransitionProbForPos(pos-1, transitions)
            
            #pseudonormalize by scaling factor from fwd
            for state_id in range(num_states):
                table[pos][state_id] -= scale[pos+1] 
        
        '''RESULT'''
        #p(X) - probability of the observed sequence (samples)
        pX = table[0][start_id]
        #print "bck: ", math.exp(pX), pX
        
        return table, pX
        
        
    def maxPosteriorDecoding(self, samples, M, P, MSC, PSC, mixture):
        """Maximum posterior probability decoding. 
        
        Returns list of states -- for each position the one with highest posterior
        probability.
        
        """
        n = len(samples)
        
        fwd, pX1, scale = self.computeForward(samples, M, P, MSC, PSC, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, MSC, PSC, mixture, scale)
        
        if abs(pX1-pX2) > 10e-8:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        path = []
        for pos in xrange(1, n+1):
            max_posterior = self.neg_inf
            arg_max = -1
            for s in range(self.getNumPP()):
                # fwd, bck are in log space (therefore + instead *); p(X) is constant
                tmp = fwd[pos][s] + bck[pos][s]
                if tmp > max_posterior:
                    max_posterior = tmp
                    arg_max = s
            #print max_posterior
            max_IP = self.states[ arg_max ].inheritance_pattern
            path.append( self.inheritance_patterns.index(max_IP) )
        
        return path
    
    
    def posteriorDecoding(self, samples, M, P, MSC, PSC, mixture):
        """Posterior probability decoding.
        
        For each position returns a list of states ordered by their posterior
        porobability. Returns list(list(tuple(probability, state)))
        
        """
        n = len(samples)
        num_real_states = self.getNumPP()
        
        fwd, pX1, scale = self.computeForward(samples, M, P, MSC, PSC, mixture)
        bck, pX2 = self.computeBackward(samples, M, P, MSC, PSC, mixture, scale)
        
        #print pX1, pX2, sum(scale)
        
        if abs(pX1-pX2) > 10e-8:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        table = []
        over90 = 0
        numall = 0
        for pos in xrange(1, n+1):
            tmp = []
            tmp2 = []
            for ip_id, ip in enumerate(self.inheritance_patterns):
                #max over corresponding phased patterns
                max_ = self.neg_inf
                for s_id, s in enumerate(self.states[:num_real_states]):
                    if s.inheritance_pattern == ip:
                        #fwd, bck are in log space (therefore + instead *); p(X) is constant
                        #max_ = max(max_, fwd[pos][s_id] + bck[pos][s_id])
                        max_ = max(max_, fwd[pos][s_id] + bck[pos][s_id] - pX1)
                
                '''#sum over corresponding phased patterns
                sum_ = self.neg_inf
                for s_id, s in enumerate(self.states):
                    if s.inheritance_pattern == ip:
                        #fwd, bck are in log space (therefore + instead *); p(X) is constant
                        sum_ = self.logSum(sum_, fwd[pos][s_id] + bck[pos][s_id])
                max_ = sum_
                '''
                    
                tmp.append((max_, ip_id))
                tmp2.append([max_, ip_id])
            
            sum_ = self.neg_inf
            for x in tmp2:
                sum_ = self.logSum(sum_, x[0])
            for i in range(len(tmp2)):
                tmp2[i][0] -= sum_
                tmp2[i][0] = int(math.exp(tmp2[i][0])*1000)
            tmp2.sort(reverse=True)
            #print tmp2
            if tmp2[0][0] >= 900: over90 += 1
            numall += 1
            
            tmp.sort(reverse=True)
            table.append(tmp)
        
        print over90, numall, "that is:", over90/float(numall)
        return table
    
    def mixedDecoding(self, samples, M, P, MSC, PSC, mixture):
        posterior = self.likelihoodDecoding(samples, M, P, MSC, PSC, mixture)
        n = len(M)
        prediction = [-1 for i in xrange(n)]
        
        last = -1
        for i, pp in enumerate(posterior):
            column = pp
            
            #compute sum
            sum_ = self.neg_inf
            for x in column:
                sum_ = self.logSum(sum_, x[0])
            
            #normalize
            for j in range(len(column)):
                column[j] = list(column[j])
                column[j][0] -= sum_
                column[j][0] = math.exp(column[j][0])
            #column.sort(reverse=True)
            
            #if high confidency, then take that prediction
            if column[0][0] >= 0.5 or i == n-1:
                if last + 1 == i:
                    #print column
                    prediction[i] = column[0][1]
                    last = i
                else: #else run max prob labeling
                    last += 1
                    #a = max(last-5, 0)
                    #b = min(n, i+6)
                    #path = self.extendedLabeling(samples[a:b], M[a:b], P[a:b], mixture)
                    print i - last + 1, 
                    path = self.viterbiPath(samples[last:i+1], M[last:i+1], P[last:i+1], MSC[last:i+1], PSC[last:i+1], mixture)
                    for pos in range(last, i+1):
                        prediction[pos] = path[pos-last]
                    last = i
                    
        print prediction
        return prediction    
            
        
        
    #TODO: NEEDS REWRITING !!!!!!!!!!!!!!!!!    
    def baumWelshForTransitions(self, samples, M, P, mixture):
        """Baum-Welsh EM algorithm for HMM parameter estimation.
        
        Training implemented only for transition probabilities.
        
        Doesn't return anything, directly modifies self.transitions
        
        """
        num_states = self.getNumStates()
        n = len(samples)
        A = [[self.neg_inf for x in range(num_states)] for y in range(num_states)]
        num_iter = 0
        while True and num_iter < 5:
            num_iter += 1
            #expectation
            fwd, pX1, scale = self.computeForward(samples, M, P, mixture)
            bck, pX2 = self.computeBackward(samples, M, P, mixture, scale)
            if abs(pX1-pX2) > 10e-8:
                print "p(X) from forward and backward DP doesn't match", pX1, pX2
                return  
            
            
            for k in range(num_states):            
                for l in range(num_states):
                    if self.transitions[k][l] == self.neg_inf:
                        A[k][l] = self.neg_inf
                        continue
                        
                    val = self.neg_inf
                    for pos in xrange(1, n):
                        emis_p = self.logLHGivenState(samples[pos], M[pos], P[pos], mixture, self.states[l])
                        tmp = fwd[pos][k] + self.transitions[k][l] + emis_p + bck[pos+1][l]
                        if tmp == self.neg_inf: continue
                        val = self.logSum(val, tmp)
                    
                    if val == self.neg_inf:
                        A[k][l] = self.neg_inf
                    else:
                        val -= pX1
                        A[k][l] = val
                    #if math.exp(A[k][l]) < 0.1:
                    #    A[k][l] = math.log(0.1)
                         
            #maximization
            old_trans = copy.deepcopy(self.transitions)
            for k in range(num_states):
                summ = self.neg_inf
                for m in range(num_states):
                    if A[k][m] == self.neg_inf: continue
                    summ = self.logSum(summ, A[k][m])
                #print summ, math.exp(summ)
                for l in range(num_states):
                    self.transitions[k][l] = A[k][l] - summ
                    #if self.transitions[k][l] < math.log(0.0005):
                    #    self.transitions[k][l] = math.log(0.0005)
            
            change = 0.
            for k in range(num_states):
                summ = self.neg_inf
                for m in range(num_states):
                    if A[k][m] == self.neg_inf: continue
                    summ = self.logSum(summ, self.transitions[k][m])
                #print summ
                for l in range(num_states):
                    #self.transitions[k][l] -= summ
                    if not (self.transitions[k][l] == self.neg_inf or old_trans[k][l] == self.neg_inf):
                        change += (self.transitions[k][l] - old_trans[k][l])**2
                    
            print "change: ", change 
            if change < 1: break
                    
        #print the table    
        for k in xrange(num_states):
            for l in xrange(num_states):
                print "%0.4f" % (math.exp(self.transitions[k][l])), 
            print " "
    
    
    #TODO: NEEDS REWRITING !!!!!!!!!!!!!!!!!         
    def viterbiTrainingForTransitions(self, samples, M, P, mixture):
        """Viterbi training algorithm for HMM parameter estimation.
        
        Training implemented only for transition probabilities.
        
        Doesn't return anything, directly modifies self.transitions
        
        """
        num_states = self.getNumStates()
        n = len(samples)
        A = [[0. for x in range(num_states)] for y in range(num_states)]
        num_iter = 0
        while True and num_iter < 5:
            num_iter += 1
            #expectation
            path, table = self.viterbiPath(samples, M, P, mixture)
            p_path = max(table[len(samples)])
            print "P of the max probable path: ", p_path
            
            for k in range(num_states):
                for l in range(num_states):
                    val = 0
                    for pos in xrange(1, n-1):
                        if path[pos] == k and path[pos+1] == l:
                            val += 1
                    #if val < 1: val = 1
                    A[k][l] = val + int(k == self.main_state)
                A[k][k] += 1
                A[k][self.main_state] += 1 
                
                    
                         
            #maximization
            old_trans = copy.deepcopy(self.transitions)
            for k in range(num_states):
                summ = sum(A[k])
                if summ == 0:
                    summ = float("inf")
                else: 
                    summ = math.log(sum(A[k]))
                #for m in range(num_states):
                #    summ = self.logSum(summ, A[k][m])
                #print summ, math.exp(summ)
                for l in range(num_states):
                    if A[k][l] == 0:
                        self.transitions[k][l] = self.neg_inf
                    else:
                        self.transitions[k][l] = math.log(A[k][l]) - summ
                    #if self.transitions[k][l] < math.log(0.0005):
                    #    self.transitions[k][l] = math.log(0.0005)
                        
            change = 0.
            for k in range(num_states):
                summ = self.neg_inf
                for m in range(num_states):
                    if A[k][m] == self.neg_inf: continue
                    summ = self.logSum(summ, self.transitions[k][m])
                for l in range(num_states):
                    #self.transitions[k][l] -= summ
                    if not (self.transitions[k][l] == self.neg_inf or old_trans[k][l] == self.neg_inf):
                        change += (self.transitions[k][l] - old_trans[k][l])**2
                    
            print "change: ", change 
            if change < 1: break
                   
        #print the table    
        for k in xrange(num_states):
            for l in xrange(num_states):
                print "%0.4f" % (math.exp(self.transitions[k][l])), 
            print " "
        
            
    def likelihoodDecoding(self, samples, M, P, MSC, PSC, mixture):
        '''
        for each position returns a list of states ordered by the likelihood
        of the data in corresponding state, i.e. by emission porobabilities
        '''
        n = len(samples)
        num_real_states = self.getNumPP()
        
        table = []
        for pos in xrange(1, n+1):
            tmp = []
            for ip_id, ip in enumerate(self.inheritance_patterns):
                #max over corresponding phased patterns
                max_ = self.neg_inf
                for s in self.states[:num_real_states]:
                    if s.inheritance_pattern == ip:
                        ll = self.logLHGivenStateWCoverage(pos-1, samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, s)
                        max_ = max(max_, ll)
                
                '''#sum over corresponding phased patterns
                sum_ = self.neg_inf
                for s in self.states:
                    if s.inheritance_pattern == ip:
                        ll = self.logLHGivenState(samples[pos-1], M[pos-1], P[pos-1], MSC[pos-1], PSC[pos-1], mixture, s)
                        sum_ = self.logSum(sum_, ll)
                '''        

                tmp.append((max_ , ip_id))
                
            tmp.sort(reverse=True)    
            table.append(tmp)
        return table
    
    def doTheShit(self):
        daMem = []
        alleles = [ 'A', 'T' ] 
        for M1 in alleles:
            for M2 in alleles:
                for F1 in alleles:
                    for F2 in alleles:
                        for P1 in alleles:
                            for P2 in alleles:
                                M = [M1, M2]
                                P = [P1, P2]
                                F = [F1, F2]
                                MSC = PSC = (10, 10)
                                counts = [90*(M1=='A')+90*(M2=='A')+10*(F1=='A')+10*(F2=='A'), 0, 0, 90*(M1=='T')+90*(M2=='T')+10*(F1=='T')+10*(F2=='T')]
                                tmp = self.maxLikelyState(counts, M, P, MSC, PSC, 0.1);
                                daMem.append((M, P, F, counts, tmp))
        return daMem
                        
    
    def maxLikelyState(self, nuc_counts, M, P, MSC, PSC, mixture):
        '''
        blah blah
        import fcnvHMM; f = fcnvHMM.FCNV()
        tmp = f.maxLikelyState((100,0,0,100), ('A','T'), ('A','T'), 0.1)

        '''
        num_real_states = self.getNumPP()
        
        tmp = []
        for ip_id, ip in enumerate(self.inheritance_patterns):
            #max over corresponding phased patterns
            max_ = self.neg_inf
            for s in self.states[:num_real_states]:
                if s.inheritance_pattern == ip:
                    ll = self.logLHGivenState(nuc_counts, M, P, MSC, PSC, mixture, s)
                    ll2 =self.logLHGivenState(nuc_counts, M, P, MSC, PSC, mixture, s, True)
                    max_ = max(max_, ll)

                    tmp.append((ll, ll2, s.phased_pattern, s.inheritance_pattern))
            
        tmp.sort(reverse=True)    

        return tmp

if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    
    f = FCNV()
    #print "\n".join(map(str, f.doTheShit()))
    mem = f.doTheShit()
    for i, x in enumerate(mem):
        print i, ">>  ", " ".join(map(str,x[:4]))
        for a in x[4]:
            print "         ", a[0], a[1], str(a[2])+"\t"+str(a[3]), '**'*(a[3] == (1,1))
    
    
    
