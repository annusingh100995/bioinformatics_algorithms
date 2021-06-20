from Graph import MyGraph

# Compostion of seq of size k
def composition(k,seq):
    res = []
    for i in range(len(seq)-k+1):
        res.append(seq[i:i+k])
        res.sort()
    return res

def test_composition():
    seq = "CAATCATGATG"
    k = 3
    print(composition(k,seq))

test_composition()

# Suffix and prefix
def suffix(seq):
    return seq[1:]

def prefix(seq):
    return seq[:-1]

class OverlapGraph(MyGraph):
    def __init__(self,frags, reps = True):
        MyGraph.__init__(self,{})
        if reps:
            self.create_overlap_graph_with_reps(frags)
        else:
            self.create_overlap_graph(frags)
        self.reps = reps
        
    def create_overlap_graph(self,frags):
        # Each fragement is a node 
        for seq in frags:
            self.add_vertex(seq)
        """If suffix of one node is equal to the prefix
        of the other, add an edge between them."""
        for seq in frags:
            suf = suffix(seq)
            for seq2 in frags:
                if prefix(seq2) == suf:
                    self.add_edge(seq,seq2)
    
    # When there are repeted fragments
    # Node is of the format: seq-1
    def create_overlap_graph_with_reps(self,frags):
        idnum = 1
        for seq in frags:
            self.add_vertex(seq+"-"+str(idnum))
            idnum += 1
        for seq in frags:
            suf = suffix(seq)
            for seq2 in frags:
                if prefix(seq2) == suf:
                    for x in self.get_instance(seq2):
                        self.add_edge(seq+"-"+str(idnum),x)
            idnum += 1
    
    def get_instance(self,seq):
        res = []
        for k in self.graph.keys():
            if seq in k:
                res.append(k)
        return res
    
    # Node with replicates: seq-1
    # node without replicates : seq
    def get_seq(self,node):
        # if node not in node list return none
        if node not in self.graph.keys():
            return None
        #if there are replicates, get the seq part of node
        if self.reps:
            return node.split("-")[0]
        else:
            # if no replicates, node is equal to seq
            # in this case just return the node
            return node
    
    """From a path set by the Hamiltonian Circuit,
    the original sequence can be retrived by taking
    the seq from the first node in the path and 
    concatenating the last character of the seq in the
    remaining nodes of teh path following its order"""
    # If the path is given, get the sequence back
    def seq_from_path(self,path):
        if not self.check_if_hamiltonian_path(path):
            return None
        seq = self.get_seq(path[0])
        for i in range(1,len(path)):
            nxt = self.get_seq(path[i])
            seq += nxt[-1]
        return seq
    
    """A path is hamiltonian, if it conatains all nodes
    and if there are no repeated edges.
    """
    # Checking if the nodes and edges are valid
    def check_if_valid_path(self,p):
        if p[0] not in self.graph.keys():
            return False
        for i in range(1,len(p)):
            if p[i] not in self.graph.keys() or p[i] not in self.graph[p[i-1]]:
                return False
        return True
    
    def check_if_hamiltonian_path(self,p):
        if not self.check_if_valid_path(p):
            return False
        to_visit = list(self.get_nodes())
        if len(p) != len(to_visit):
            return False
        """take a node from path and if that is in to_visit
        remove it from to_visit.
        if all nodes are romves,means all nodes are covered"""
        for i in range(len(p)):
            if p[i] in to_visit:
                to_visit.remove(p[i])
            else:
                return False
        # that is, if all the nodes are covered
        if not to_visit:
            return True
        else:
            return False
    
    def search_hamiltonian_path(self):
        for k in self.graph.keys():
            p = self.search_hamiltonian_path_from_node(k)
            if p != None:
                return p
        return None
    
    # IM TRYING TO GET WHAT THE HECK IS GOING ON 
    # AS I TYPE, IGNORE OBNOXIOUS OVER COMMENTING
    # it's a exhaustive dfs algo.
    def search_hamiltonian_path_from_node(self,start):
        current = start
        visited = {start:0}
        path = [start]
        
        #means there are still nodes remaining
        while len(path) < len(self.get_nodes()):
            nxt_index = visited[current]
            
            # this has soemthing to do with edges
            # and some how unvisited edges, anyways
            if len(self.graph[current]) > nxt_index:
                nxt_node = self.graph[current][nxt_index]
                # yeasss this is counting the number of 
                # outgoing unvisited egdes soemhow
                
                visited[current] += 1
                
                if nxt_node not in path:
                    path.append(nxt_node)
                    visited[nxt_node] = 0
                    current = nxt_node
                else:
                    if len(path) > 1:
                        removed_node = path.pop()
                        del visited[removed_node]
                        current = path [-1]
                    else:
                        return None
        return path

def test_overlap_graph():
    #frags = ["ACC","ATA","CAT","CCA","TAA"]
    #ovgr = OverlapGraph(frags)
    #ovgr.print_graph()
    #print(ovgr.get_nodes())
    #print(ovgr.get_edges())
    frags =  [ "ATA", "ACC", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA",
              "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    g1 = OverlapGraph(frags)
    #g1.print_graph()
    #print("Nodes \n",g1.get_nodes())
    #print("Edges \n",g1.get_edges())
    path = ["ACC-2", "CCA-8", "CAT-5", "ATG-3", "TGG-13", "GGC-10","GCA-9", 
            "CAT-6", "ATT-4", "TTT-15", "TTC-14", "TCA-12", "CAT-7","ATA-1", "TAA-11"]
    print(g1.seq_from_path(path))
    g1.check_if_valid_path(path)
    print(g1.check_if_hamiltonian_path(path))

test_overlap_graph()

def test_1():
    orig_seq = "CAATCATGAT"
    #GATGATC"
    frags = composition(3,orig_seq)
    print(frags)
    g = OverlapGraph(frags,True)
    g.print_graph()
    #path = g.search_hamiltonian_path()
    #print(path)
    #print(g.seq_from_path(path))

# Okay my pc cried , test_1 not possible here
