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
    
    def __init__(self, ip):
        """Initialize new HMMState object
        
        Arguments:
        ip -- tuple(int, int) -- general inheritance pattern
        pp -- tuple(tuple(ints), tuple(ints)) -- phased inheritance pattern
        
        """
        super(HMMState, self).__init__()
        
        self.inheritance_pattern = ip
        
    def __str__(self):
        return "IP: {0}".format(self.inheritance_pattern)
        
    def __eq__(self, other):
        return self.inheritance_pattern == other.inheritance_pattern

class coverageFCNV(object):
    """Hidden Markov Model for fetal CNV calling"""
    
    nucleotides = ['A', 'C', 'G', 'T']
    normal_scale_factor = 1.4826022185
    
    win_size = 1000
    num_neighbors = 200.
    tag_size = 200. # equal to 2 * read length  
    magic_scale_factor = win_size * 10
    # learned line of correlation between sample and reference WRVs
    wrv_coef_c0 = 0. # intercept
    wrv_coef_c1 = 1. # slope
    #----- I1 plasma vs G1 plasma WRV, IRLS
    #wrv_coef_c0 = 0.1082 # intercept
    #wrv_coef_c1 = 0.7834 # slope
    #---------------
    #----- I1 plasma vs I1 mother WRV, IRLS
    #wrv_coef_c0 = 0.3719 # intercept
    #wrv_coef_c1 = 0.2528 # slope
    #---------------
    #----- I1 plasma vs I1 mother WRV, pure LS
    #wrv_coef_c0 = 0.1553 # intercept
    #wrv_coef_c1 = 0.6927 # slope
    #---------------
    
    
    def __init__(self, positions, prefix_sum_plasma, prefix_count_plasma, prefix_sum_ref, prefix_count_ref, gc_sum):
        """Initialize new FCNV object"""
        super(coverageFCNV, self).__init__()
        
        #store the DOC prefix data
        self.positions = positions
        self.prefix_sum_plasma = prefix_sum_plasma
        self.prefix_count_plasma = prefix_count_plasma
        self.prefix_sum_ref = prefix_sum_ref
        self.prefix_count_ref = prefix_count_ref
        self.gc_sum = gc_sum
        
        #precompute WRV related stats
        self.plasma_wins = self.getGCWindows(self.win_size, gc_sum, prefix_sum_plasma, prefix_count_plasma)
        self.ref_wins = self.getGCWindows(self.win_size, gc_sum, prefix_sum_ref, prefix_count_ref)
        self.brv_diff_mean, self.brv_diff_var = self.computeBRVDiffEstimate( \
            self.win_size, gc_sum, \
            prefix_sum_plasma, prefix_count_plasma, self.plasma_wins, \
            prefix_sum_ref, prefix_count_ref, self.ref_wins )
        self.brv_diff_mean = 0.
        
        #compute list of reference and plasma WRV values for all bins
        #self.computeWRVlist( \
        #    self.win_size, gc_sum, \
        #    prefix_sum_plasma, prefix_count_plasma, self.plasma_wins, \
        #    prefix_sum_ref, prefix_count_ref, self.ref_wins )
        
        #run intern tests
        self.neg_inf = float('-inf')
        self.inf = float('inf')
        self.nan = float('nan')
        self.testLogSum()
        
        #cache for log-likelihood values
        self.logLikelihoodCache = {}
        #self.distributionCache = {}
        
        
        #generate inheritance patterns
        self.inheritance_patterns = [1, 2, 3]
     
        #generate HMM States
        self.states = []
        for ip in self.inheritance_patterns:
            new_state = HMMState(ip)
            self.states.append(new_state)
        
        #generate silent states
        self.states.append( HMMState("s") ) #start state
        self.states.append( HMMState("t") ) #end state
        
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
    
    def getGCWindows(self, win_size, gc_sum, prefix_sum, prefix_count):
        """
        Create a sorted list of (GCratio, fragments) for bins of size win_size
        """
        windows = []
        #get the DOC means for GC-ratio bins
        pos = 0
        while pos < len(prefix_sum):
            if prefix_sum[pos]==0:
                pos += 1
                continue
            
            end = pos + win_size
            if end >= len(prefix_sum): break
            
            while end > pos and prefix_sum[end]==0: end -= 1
            if end - pos < win_size * 0.9:
                pos += 1
                continue 
            
            arrivals = (prefix_sum[end]-prefix_sum[pos]) / (prefix_count[end]-prefix_count[pos]) * win_size / self.tag_size
            gc_ratio = (gc_sum[end] - gc_sum[pos]) / float(end-pos)
            windows.append((gc_ratio, arrivals))
            pos = end + 1
            
        windows.sort()
        return windows

    def computeWRVlist(self, win_size, gc_sum, pl_sum, pl_count, pl_wins, ref_sum, ref_count, ref_wins):
        """
        Compute WRV values in the same bins of plasma and reference sequencing and output them to a file
        """
        wrv_out_file = open("wrv_out_file" + str(random.randint(1, 99999)) + ".txt", "w")
        print wrv_out_file
        diffs = []
        pos = 0
        while pos < len(pl_sum):
            if pl_sum[pos]==0 or ref_sum[pos]==0:
                pos += 1
                continue
            
            end = pos + win_size
            if end >= len(pl_sum): break
            
            while end > pos and (pl_sum[end]==0 or ref_sum[end]==0): end -= 1
            if end - pos < win_size * 0.9:
                pos += 1
                continue 
            
            pl_arrivals = (pl_sum[end] - pl_sum[pos]) / (pl_count[end] - pl_count[pos]) * win_size / self.tag_size
            pl_gc_ratio = (gc_sum[end] - gc_sum[pos]) / float(end-pos)
            pl_close_arrivals, pl_var = self.getCloseGCArrivalsSum(pl_gc_ratio, pl_wins)
            pl_brv = pl_arrivals / pl_close_arrivals
            
            try:
                ref_arrivals = (ref_sum[end] - ref_sum[pos]) / (ref_count[end] - ref_count[pos]) * win_size / self.tag_size
                ref_gc_ratio = (gc_sum[end] - gc_sum[pos]) / float(end-pos)
                ref_close_arrivals, ref_var = self.getCloseGCArrivalsSum(ref_gc_ratio, ref_wins)
                ref_brv = ref_arrivals / ref_close_arrivals
            except ZeroDivisionError:
                pos += 1
                continue
            
            print >>wrv_out_file, pl_brv, ref_brv, pos, end, end-pos
            
            pos = end + 1
            
        wrv_out_file.close()
        return
        
    def computeBRVDiffEstimate(self, win_size, gc_sum, pl_sum, pl_count, pl_wins, ref_sum, ref_count, ref_wins):
        """
        Estimate mean and variance of BRV differences between bins of plasma and reference sequencing
        """
        diffs = []
        pos = 0
        while pos < len(pl_sum):
            if pl_sum[pos]==0 or ref_sum[pos]==0:
                pos += 1
                continue
            
            end = pos + win_size
            if end >= len(pl_sum): break
            
            while end > pos and (pl_sum[end]==0 or ref_sum[end]==0): end -= 1
            if end - pos < win_size * 0.9:
                pos += 1
                continue 
            
            pl_arrivals = (pl_sum[end] - pl_sum[pos]) / (pl_count[end] - pl_count[pos]) * win_size / self.tag_size
            pl_gc_ratio = (gc_sum[end] - gc_sum[pos]) / float(end-pos)
            pl_close_arrivals, pl_var = self.getCloseGCArrivalsSum(pl_gc_ratio, pl_wins)
            pl_brv = pl_arrivals / pl_close_arrivals
            #pl_brv = self.wrv_coef_c0 + self.wrv_coef_c1*pl_brv
            
            try:
                ref_arrivals = (ref_sum[end] - ref_sum[pos]) / (ref_count[end] - ref_count[pos]) * win_size / self.tag_size
                ref_gc_ratio = (gc_sum[end] - gc_sum[pos]) / float(end-pos)
                ref_close_arrivals, ref_var = self.getCloseGCArrivalsSum(ref_gc_ratio, ref_wins)
                ref_brv = ref_arrivals / ref_close_arrivals
            except ZeroDivisionError:
                pos += 1
                continue
            
            diffs.append((pl_brv - ref_brv) * self.magic_scale_factor)
            
            pos = end + 1
            
        mean = sum(diffs) / float(len(diffs))
        var = sum([(x-mean)**2 for x in diffs]) / float(len(diffs) - 1)
        print mean, var, len(diffs)
        
        diffs.sort()
        print diffs[len(diffs)/2], pl_brv*self.magic_scale_factor, ref_brv*self.magic_scale_factor
        median = diffs[len(diffs)/2]
        mad = sorted([abs(x - median) for x in diffs])[len(diffs)/2] #median absolute deviation (MAD)
        mad_estim_var = (self.normal_scale_factor * mad)**2
        print median, mad, mad_estim_var
        
        #return median, mad_estim_var
        return mean, var
    
    @memoized     
    def getNumIP(self):
        """Return number of inheritance patterns"""
        return len(self.inheritance_patterns)
    
    @memoized    
    def getNumStates(self):
        """Return number of HMM states per one SNP position"""
        return len(self.states)   
    
    @memoized 
    def isNormal(self, state):
        return state.inheritance_pattern == 2
    
    @memoized   
    def isReal(self, state):
        return isinstance(state.inheritance_pattern, int)
    
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
    
    def generateTransitionProb(self):
        """Generate transition probabilities between HMM states"""
        #number of possible states per one SNP position
        num_states = self.getNumStates()
        num_real_states = self.getNumIP()
        
        trans = [[0. for x in range(num_states)] for y in range(num_states)]
        #first generate transitions from the real states
        for ip in self.inheritance_patterns:
            
            if ip == 2: #generate Normal states transitions
                pstay = 1./6. #0.4
                pgo = 1./6. #0.6 / (num_real_states - 1)
                for i, state1 in enumerate(self.states[:num_real_states]):
                    if state1.inheritance_pattern != ip: continue
                    #states "inside the IP component"
                    for j, state2 in enumerate(self.states[:num_real_states]):
                        if i == j: #stay in the state
                            trans[i][i] = pstay
                        else:
                            trans[i][j] = pgo
                    #to the silent exit node    
                    outState_id = self.getExitState()[0]
                    trans[i][outState_id] = pgo
                    
            else: #generate CNV states transitions
                pstay = 1./6. #0.4
                pgo = 1./6. #0.6 / (num_real_states - 1)
                for i, state1 in enumerate(self.states[:num_real_states]):
                    if state1.inheritance_pattern != ip: continue
                    #states "inside the IP component"
                    for j, state2 in enumerate(self.states[:num_real_states]):
                        if i == j: #stay in the state
                            trans[i][i] = pstay
                        elif state1.inheritance_pattern + state2.inheritance_pattern != 4:
                            trans[i][j] = pgo
                        else:
                            trans[i][j] = 0.
                    #to the silent exit node    
                    outState_id = self.getExitState()[0]
                    trans[i][outState_id] = pgo
        
        #now generate transitions from the silent states
        for i, state1 in enumerate(self.states[num_real_states:]):
            i += num_real_states
            #if it is the start node
            if state1.inheritance_pattern == "s":
                for j, state2 in enumerate(self.states[:num_real_states]):
                    trans[i][j] = 1. / self.getNumIP()
        
        
        #normalize and take logarithm
        for i in range(num_states):
            for j in range(num_states):
                if trans[i][j] < 10e-10: trans[i][j] = self.neg_inf
                else: trans[i][j] = math.log(trans[i][j])
            self.logNormalize(trans[i])
        
        return trans
    
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
    
    def getCloseGCArrivalsSum(self, gc_ratio, wins):
        """
        Computes sum of sequencing fragments in bins with closest GC ratio
        """
        num_neighbors = self.num_neighbors
        
        #binary search the position
        pos_l = bisect_left(wins, (gc_ratio, 0.)) 
        pos_r = bisect_left(wins, (gc_ratio, float("inf")))
        pos = (pos_l + pos_r) / 2
        if pos == len(wins): pos -= 1
        
        close_arrivals = wins[pos][1]
        have = 1
        dis = 0
        left = pos
        right = pos+1
        while have < num_neighbors:
            dis += 1
            if pos+dis < len(wins):
                close_arrivals += wins[pos+dis][1]
                right += 1
                have += 1
            if pos-dis > 0:
                close_arrivals += wins[pos-dis][1]
                left -= 1
                have += 1
        
        var = 0.
        mean = close_arrivals / num_neighbors
        for x in range(left, right):
            var += (wins[x][1] - mean)**2
        var /= (num_neighbors - 1)
            
        return close_arrivals, var
    
    def logLHGivenStateWCoverage(self, pos_ind, mix, state):
        """
        Compute the likelihood of the observed coverage conditional on CNV event
        i.e. emission probability given the state
        """
        win_size = self.win_size
        begin_ind = max(0, pos_ind - 1)
        end_ind = min(pos_ind + 1, len(self.positions) - 1)
        
        #mid_left = max(int((self.positions[pos_ind] + self.positions[begin_ind]) / 2.), self.positions[pos_ind] - win_size/2)
        #mid_right = min(int((self.positions[pos_ind] + self.positions[end_ind]) / 2.), self.positions[pos_ind] + win_size/2)
        mid_left = min(int((self.positions[pos_ind] + self.positions[begin_ind]) / 2.), self.positions[pos_ind] - win_size/10)
        mid_right = max(int((self.positions[pos_ind] + self.positions[end_ind]) / 2.), self.positions[pos_ind] + win_size/10)
        
        #leftmost = self.positions[pos_ind] - win_size/2
        #rightmost = self.positions[pos_ind] + win_size/2
        
        #compute stats for REFERENCE
        b = mid_left   #leftmost
        e = mid_right  #rightmost
        while b < 0 or self.prefix_sum_ref[b] == 0: b += 1
        while e >= len(self.prefix_sum_ref) or self.prefix_sum_ref[e] == 0: e -= 1
        ref_win_size = e - b
        if ref_win_size <= 0: return 0.
        mu_doc = self.prefix_sum_ref[e] - self.prefix_sum_ref[b]
        mu_doc /= float(self.prefix_count_ref[e] - self.prefix_count_ref[b])
        ref_gc_ratio = (self.gc_sum[e] - self.gc_sum[b]) / float(ref_win_size)
        
        #adjust DOC conditional on inheritance pattern
        ref_doc = mu_doc
        ref_doc += (state.inheritance_pattern - 2) * (mu_doc * mix/2.)
        #get arrivals rate
        mu_arrivals = (mu_doc * win_size) / self.tag_size
        ref_arrivals = (ref_doc * win_size) / self.tag_size
        ref_close_arrivals, ref_var = self.getCloseGCArrivalsSum(ref_gc_ratio, self.ref_wins)
        
        #compute bin ratio values
        mu_brv  = mu_arrivals / ref_close_arrivals
        ref_brv = ref_arrivals / ref_close_arrivals
        
        #compute stats for PLASMA
        b = mid_left   #leftmost
        e = mid_right  #rightmost
        while b < 0 or self.prefix_sum_plasma[b] == 0: b += 1
        while e >= len(self.prefix_sum_plasma) or self.prefix_sum_plasma[e] == 0: e -= 1
        pl_win_size = e - b
        if pl_win_size <= 0: return 0.
        pl_doc = self.prefix_sum_plasma[e] - self.prefix_sum_plasma[b]
        pl_doc /= float(self.prefix_count_plasma[e] - self.prefix_count_plasma[b])
        pl_arrivals = (pl_doc * win_size) / self.tag_size
        pl_gc_ratio = (self.gc_sum[e] - self.gc_sum[b]) / float(pl_win_size)
        pl_close_arrivals, pl_var = self.getCloseGCArrivalsSum(pl_gc_ratio, self.plasma_wins)
        pl_brv = pl_arrivals / pl_close_arrivals
        pl_brv = self.wrv_coef_c0 + self.wrv_coef_c1*pl_brv
        
        
        #noise_prob = self.logGaussian([pl_brv], [mu_brv], [mu_brv*20])
        coverage_prob = self.logGaussian([(pl_brv - ref_brv) * self.magic_scale_factor], [self.brv_diff_mean], [self.brv_diff_var])
        result = coverage_prob
        
        #print ref_doc, ref_arrivals, ref_win_size, ref_brv, ref_var, '|', pl_doc, pl_arrivals, pl_win_size, pl_brv, pl_var, '|', result
        
        if pl_win_size < win_size/2: result = 0.
        #if mid_right - mid_left < win_size/2: result = 0.
        if result < -15: result = -15
        
        return result
    
    def viterbiPath(self, mixture):
        """
        Viterbi decoding of the most probable path.
        """
        num_states = self.getNumStates()
        num_real_states = self.getNumIP() 
        transitions = self.transitions
        
        n = len(self.positions)
        predecessor = [[0 for i in range(num_states)] for j in xrange(n+1)] 
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
            predecessor[0][state_id] = -1
         
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case

            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(self.states[:num_real_states]):
                #emission probability in the given state
                emis_p = self.logLHGivenStateWCoverage(pos-1, mixture, state)
                
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
            state_path[i] = state.inheritance_pattern
        
        print "Viterbi path probability: ", table[n][exit_id]    
        return path[1:], state_path[1:]
        
    
    def computeForward(self, mixture):
        '''
        Posterior decoding: forward algorithm
        '''
        num_states = self.getNumStates()
        num_real_states = self.getNumIP()
        patterns = self.inheritance_patterns
        transitions = self.transitions
        
        n = len(self.positions)
        #scale factors
        scale = [0. for i in xrange(n+1)]
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
        
        '''DYNAMIC PROGRAMMING'''
        #for all SNP positions do:
        for pos in xrange(1, n+1):
            #real positions are <1..n+1), pos 1 is the base case
            
            #(i) compute new values for all phased patterns - "real" states of the HMM:
            for state_id, state in enumerate(self.states[:num_real_states]):
                #emission probability in the given state
                emis_p = self.logLHGivenStateWCoverage(pos-1, mixture, state)
                
                #probability of this state is `emission * \sum {previous state * transition}`
                summ = self.neg_inf
                for prev_id in range(num_states):
                    if transitions[prev_id][state_id] == self.neg_inf: continue #just for speed-up
                    tmp =  table[pos-1][prev_id] + transitions[prev_id][state_id]
                    summ = self.logSum(summ, tmp)
                table[pos][state_id] = emis_p + summ
            
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
        
        
    def computeBackward(self, mixture, scale):
        '''
        Posterior decoding: backward algorithm
        '''
        num_states = self.getNumStates()
        num_real_states = self.getNumIP()
        patterns = self.inheritance_patterns
        transitions = self.transitions
        
        n = len(self.positions)
        #DP table dim: seq length+2 x num_states
        table = [[self.neg_inf for i in range(num_states)] for j in xrange(n+1)]
        
        start_id, start_state = self.getStartState()
        exit_id, exit_state = self.getExitState()
        
        
        '''INITIALIZE'''
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
                emis_p.append(self.logLHGivenStateWCoverage(pos, mixture, state))
                
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
            
            #pseudonormalize by scaling factor from fwd
            for state_id in range(num_states):
                table[pos][state_id] -= scale[pos+1] 
        
        '''RESULT'''
        #p(X) - probability of the observed sequence (samples)
        pX = table[0][start_id]
        #print "bck: ", math.exp(pX), pX
        
        return table, pX
        
        
    def maxPosteriorDecoding(self, mixture):
        """Maximum posterior probability decoding. 
        
        Returns list of states -- for each position the one with highest posterior
        probability.
        
        """
        n = len(self.positions)
        
        fwd, pX1, scale = self.computeForward(mixture)
        bck, pX2 = self.computeBackward(mixture, scale)
        
        if abs(pX1-pX2) > 10e-8:
            print "p(X) from forward and backward DP doesn't match", pX1, pX2
            return
        
        path = []
        for pos in xrange(1, n+1):
            max_posterior = self.neg_inf
            arg_max = -1
            for s in range(self.getNumIP()):
                # fwd, bck are in log space (therefore + instead *); p(X) is constant
                tmp = fwd[pos][s] + bck[pos][s]
                if tmp > max_posterior:
                    max_posterior = tmp
                    arg_max = s
            #print max_posterior
            max_IP = self.states[ arg_max ].inheritance_pattern
            path.append( self.inheritance_patterns.index(max_IP) )
        
        return path
    
    
    def posteriorDecoding(self, mixture):
        """Posterior probability decoding.
        
        For each position returns a list of states ordered by their posterior
        porobability. Returns list(list(tuple(probability, state)))
        
        """
        n = len(self.positions)
        num_real_states = self.getNumIP()
        
        fwd, pX1, scale = self.computeForward(mixture)
        bck, pX2 = self.computeBackward(mixture, scale)
        
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
    
            
    def likelihoodDecoding(self, mixture):
        '''
        for each position returns a list of states ordered by the likelihood
        of the data in corresponding state, i.e. by emission porobabilities
        '''
        n = len(self.positions)
        num_real_states = self.getNumIP()
        
        table = []
        for pos in xrange(1, n+1):
            tmp = []
            for ip_id, ip in enumerate(self.inheritance_patterns):
                #max over corresponding phased patterns
                max_ = self.neg_inf
                for s in self.states[:num_real_states]:
                    if s.inheritance_pattern == ip:
                        ll = self.logLHGivenStateWCoverage(pos-1, mixture, s)
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


if __name__ == "__main__":
    pass


