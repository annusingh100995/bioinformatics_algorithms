# Hidden Markov Model
# REFERENCE: Bioinformatics Algorithms : Design and Implementation in Python
#            Miguel Rocha
#            University of Minho, Braga, Portugal
#            Pedro G. Ferreira
#            Ipatimup/i3S, Porto, Portugal

class HiddenMarkovModel:
    
    def __init__(self, initial_prob, emission_prob, transition_prob):
        
        self.initial_prob = initial_prob
        self.emission_prob = emission_prob
        self.transition_prob = transition_prob
        self.states = self.emission_prob.keys()
        self.symbols = self.emission_prob[list(self.emission_prob.keys())[0]].keys()
        #self.emission_prob[self.emission_prob.keys()[0]].keys()
     
    # Return initial prob of a given state
    def get_initial_prob(self, state):
        if state in self.states:
            return (self.initial_prob[state])
        else:
            return 0
    
    def set_intial_prob(self,state,prob):
        if state in self.states:
            self.initial_prob[state] = prob
     
    
    # Get emmission prob of a symbol given a state
    def get_emmision_prob(self,state,symbol):
        if state in self.states and symbol in self.symbols:
            return (self.emission_prob[state][symbol])
        else:
            return 0
        
    def set_emission_prob(self,state,symbol,prob):
        if state in self.states and symbol in self.symbols:
            self.emission_prob[state][symbol] = prob
        
    # Get transition prob from origin to destination state
    def get_transition_prob(self,orig_state, dest_state):
        if orig_state in self.states and dest_state in self.states:
            return (self.transition_prob[orig_state][dest_state])
        else:
            return 0
        
    def set_transition_prob(self,origin_state , dest_state, prob):
        if origin_state in self.states and dest_state in self.states:
            self.transition_prob[origin_state][dest_state] = prob
    
    # Calculate joint prob. given an obsevered sequence and a state path

    def joint_probibility(self,observed_seq, state_path):
        obs_seq_len = len(observed_seq)
        
        if obs_seq_len == 0:
            return None
        
        state_path_len = len(state_path)
        if obs_seq_len != state_path_len:
            print("Observed sequence and state path are of differernt length!!")

        prob = self.get_initial_prob(state_path[0])*(self.get_emmision_prob(state_path[0],observed_seq[0]))
        
        for i in range(1,len(state_path)):
            prob = prob*(self.get_transition_prob(state_path[i-1],state_path[i]))*(self.get_emmision_prob(state_path[i], observed_seq[i]))
        
        return prob
        
     
    def forward_probibilites(self,obsered_seq):
        obs_seq_len = len(observed_seq)
        if obs_seq_len == 0:
            return []
        
        prob_list = [{}]
                            
        # probilities for the first observation
        for state in self.states:
            prob_list[0][state] = self.get_initial_prob(state)*(self.get_emmision_prob(state, observed_seq[0]))
            
        for i in range(1,obs_seq_len):
            prob_list.append({})
            for dest_state in self.states:
                prob = 0
                for orig_state in self.states:
                    prob += prob_list[i-1][orig_state]*(self.get_transition_prob(orig_state,dest_state))
                
                prob_list[i][dest_state] = prob*self.get_emmision_prob(dest_state,observed_seq[i])
        
        return prob_list
    
    def backward_probilities(self,observed_seq):
        obs_seq_len = len(observed_seq)
        if obs_seq_len == 0:
            return []
        
        beta = [{}]
        # probabilty for last observations
        for state in self.states:
            beta[0][state] = 1
        
        for i in range(obs_seq_len-1,0,-1):
            beta.insert(0,{})
            
            for orig_state in self.states:
                prob = 0
                for dest_state in self.states:
                    prob += beta[1][dest_state]*(self.get_transition_prob(orig_state, dest_state))*(self.get_emmision_prob(dest_state,observed_seq[i]))
                    beta[0][orig_state] = prob
        return beta
            
    def viterbi_algorithm(self,observed_sequence):
        """Viterbi algorithm finds the most probable 
        state path for a given observed sequence"""
        obs_seq_len = len(observed_sequence)
        
        if obs_seq_len == 0:
            return []
        
        # Prob of observing state Si at time t given the sqeuence till ti-1
        delta = {}
        # origin state with maximum prob. for each Si at time t
        state_path = {}
        
        for state in self.states:
            delta[state] = self.get_initial_prob(state)* (self.get_emmision_prob(state,observed_sequence[0]))
            # yeah , because for first observation,state with max probility is state itself
            state_path[state] = [state]
            
        
        # Finding delta for t = 1 to T
        for t in range(1,obs_seq_len):
            new_state_path = {}
            new_path = {}
            delta_temp = {}
            for dest_state in self.states:
                intermediate_prob = []
                for origin_state in self.states:
                    """For a particular destination state dest_state, we can come from N differnt origina states
                    So here prob is the probility of coming from a origin state that iterates through 1 to N"""
                    #print("orgin stte:",origin_state,"destinaton state",dest_state)
                    #print("transition prob",self.get_transition_prob(origin_state,dest_state))
                    #print("delata of origin state",origin_state,delta[origin_state])
                    prob = (delta[origin_state])*(self.get_transition_prob(origin_state,dest_state))
                    #print("Prob: {}".format(prob))
                    """Intermiediate prob stores the prob of coming from a origin state and the state itself"""
                    
                    intermediate_prob.append((prob,origin_state))
                    #print("intermediate prob:" ,intermediate_prob)
                    
                    """So, here, for each destination state, we find the max probability of arriving to destination state
                    from a origin state and the respective state with the highest probability"""
                (max_prob , max_state) = max(intermediate_prob)
                
                #print("max_prob: ",max_prob, "max_state :" ,max_state)
                """This prob is the delta t of i delta[t][dest_state]"""
                prob = self.get_emmision_prob(dest_state, observed_sequence[t])*max_prob
                
                #print("prob: " , prob)
                
                delta_temp[dest_state] = prob
                #print("delta temp: od dest sate" ,delta_temp[dest_state],dest_state)
                """new_state_path is the state that is the origin state with maximum prob for a particular destination state.
                Or in other words, for this destination state, there are maximum prob that the sequenc is coming from this (max_state)
                which is the origin state"""
                new_state_path[dest_state] = max_state
                
                """Okay so , this is the actual path. rest are all just the states.
                    And this is the path of all the destination states and the max_prob origin state.
                    That is , the state path which has the max_state as the last state. which is equivalent to 
                    the last state being the origin state.
                    I think the length of each path will increase with t"""
                new_path[dest_state] = state_path[max_state]+[dest_state]
                #print("New path for dest_state:",dest_state,"PATH:",new_path[dest_state])
            delta = delta_temp
            
            state_path = new_path
            
        max_state = None
        max_prob = 0
        for state in self.states:
            if delta[state]> max_prob:
                max_prob = delta[state]
                max_state = state
        return (max_prob,state_path[max_state])
    
    def baum_welch(self,observed_seq):
        obs_seq_len = len(observed_seq)
        if obs_seq_len == 0:
            return []
        
        alpha = self.forward_probibilites(observed_seq)
        beta = self.backward_probilities(observed_seq)
        
        # gamma t of i
        gamma = [{} for t in range(obs_seq_len)]
        
        # Xi t for i to j, t= obs_seq_len -1 cuz, there is no transition from last state
        xi = [{} for t in range(obs_seq_len-1)]
        
        """Calculating gamma t of i, gamma(t)(i) = alpha(t)(i)*beta(t)(i)"""
        
        for t in range(obs_seq_len):
            sum_alpha_beta = 0
            for i in self.states:
                gamma[t][i] = (alpha[t][i])*(beta[t][i])
                sum_alpha_beta += gamma[t][i]
            
            """ Calculating the normalising factor of the gammas"""
            for i in self.states:
                gamma[t][i] = gamma[t][i]/sum_alpha_beta
                
            if t == 0:
                self.set_intial_prob(i,gamma[t][i])
            # for Xi ,no transition possible in the last state.
            # this was smart, to calculate the gamma before,then take care of the first and last states
            # and then calculate Xi for all the states
            if t == obs_seq_len-1:
                continue
            
            # sum_prob is the normalising factor.
            # It is the prob of all the transition possible for states.
            sum_prob = 0
            for origin_state in self.states:
                # for each origin state there are len(states) destination states
                xi[t][origin_state] = {}
                for dest_state in self.states:
                    p = (alpha[t][origin_state])*(self.get_transition_prob(origin_state,dest_state))*(self.get_emmision_prob(dest_state,observed_seq[t+1]))*(beta[t+1][dest_state])
                    xi[t][origin_state][dest_state] = p
                    sum_prob += p 
                
            for origin_state in self.states:
                for dest_states in self.states:
                    xi[t][origin_state][dest_states]/sum_prob
            
        # M-step
        
        """ This the prob that for each state i the expected number of times in appears in the t.
        Expected number of times i th state appears in the observed sequence"""
        for i in self.states:
            denominator = 0
            for t in range(obs_seq_len):
                denominator += gamma[t][i]
            
            """Expected number of times w was observed/emmited when the hidden state was i"""
            for w in self.symbols:
                numerator = 0
                for t in range(obs_seq_len):
                    if observed_seq[t] == w:
                        numerator += gamma[t][i]
                        
                        if denominator>0:
                            self.set_emission_prob(i,w,(numerator/denominator))
                        else:
                            self.set_emission_prob(i,w,0.0)
        
        #E- step
        for i in self.states:
            for j in self.states:
                denominator = 0.0
                for t in range(obs_seq_len-1):
                    denominator += gamma[t][i]
                    """Uptil here, this is the expected number of times state i was observed in the observed sequence
                    or number of times the origin state for a transition was i"""
                
                """Below, is the expected number transition that were made from state i to j"""
                numerator = 0 
                for t in range(obs_seq_len-1):
                    numerator += xi[t][i][j]
                    
                if denominator > 0:
                    self.set_transition_prob(i,j,(numerator/denominator))
                else:
                    self.set_transition_prob(i,j,0.0)
            
initial_prob = {'5': 0.8, 'M': 0.15, '3': 0.05}

emission_prob = {"5":{"A":0.20, "C": 0.30, "G":0.30, "T":0.20},
                "M":{"A":0.25, "C": 0.25, "G":0.25, "T":0.25},
                "3":{"A":0.35, "C": 0.15, "G":0.15, "T":0.35}}

transition_prob = {"5":{"5": 0.8, "M": 0.2,"3": 0.0},
                   "M":{"5": 0.0, "M": 0.9,"3": 0.1},
                    "3":{"5": 0.0, "M": 0.0,"3": 1.0}}

hmm = HiddenMarkovModel(initial_prob,emission_prob,transition_prob)

print(hmm.viterbi_algorithm(observed_sequence=['A','T','T','G']))