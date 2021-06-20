"""
DeBruijn Graph for genome assembly
The fragments are represented as edges.Here the solution 
will be paths over the graph containing all the edges exactly
once.

The node contain sequences that corresponds to either a 
suffix or a prefix of one of the fagments.
Each edge each edge corresponds to a fragment,connects a node
representing its orefis to a node representing its suffix

Repeated fragments are represented as multiple edges
connecting the same pair of nodes.
"""
from Graph import MyGraph

def suffix(seq):
    return seq[1:]
def prefix(seq):
    return seq[:-1]

class DeBruijnGraph(MyGraph):
    
    def __init__(self,frags):
        MyGraph.__init__(self, {})
        self.create_deBruijn_graph(frags)
        
    
    def add_edge(self,o,d):
        if o not in self.graph.keys():
            self.add_vertex(o)
        if d not in self.graph.keys():
            self.add_vertex(d)
        self.graph[o].append(d)
    
    def create_deBruijn_graph(self,frags):
        for seq in frags:
            suf = suffix(seq)
            self.add_vertex(suf)
            pref = prefix(seq)
            self.add_vertex(prefix)
            self.add_edge(pref,suf)
    
    def check_balanced_node(self,node):
        return self.in_degree(node) == self.out_degree(node)

    def check_balanced_graph(self):
        for n in self.graph.keys():
            if not self.check_balanced_node(n):
                return False
        
        return True
    
    def check_nearly_balanced_graph(self):
        res = None, None
        
        for n in self.graph.keys():
            indeg = self.in_degree(n)
            outdeg = self.out_degree(n)
            
            if indeg-outdeg == 1 and res[1] is None:
                res = res[0],n
            elif indeg-outdeg == -1 and res[0] is None:
                res = n,res[1]
            elif indeg == outdeg:
                pass
            else:
                return None, None
            
        return res
    # If it is none,none: balaned
    # If node(x), node(y) == means nearly balanced
    
    # Indegree is modified to account for
    # multiple edges from the same o and d.
    def in_degree(self,v):
        res = 0
        for k in self.graph.keys():
            if v in self.graph[k]:
                res += self.graph[k].count(v)
        return res
    
    def eulerian_cycle(self):
        if not self.is_connected() or not self.check_balanced_graph():
            return None
        edges_visit = list(self.get_edges())
        i = 1
        if res != []:
            # while there are edges to visit
            # or len(edge_list) > 0
            while edges_visit:
                pair = edges_visit[0]
                i = 1
                if res != []:
                    while pair[0] not in res:
                        pair = edges_visit[i]
                        i = i+1
                    
                edges_visit.remove(pair)
                # selecting a random vertex and its succesor
                start , nxt = pair
                cycle = [start,nxt]
                # until we reach the start node back
                while nxt != start:
                    # for successor of next 
                    for suc in self.graph[nxt]:
                        #if there is unvisited edge b/w nxt->suc
                        if (nxt,suc) in edges_visit:
                            # next pair = nxt,suc
                            pair = (nxt,suc)
                            # the nxt variable is the suc node
                            nxt = suc
                            # append the suc in the cylce
                            cycle.append(nxt)
                            # remove the current ege used from the edge list
                            edges_visit.remove(pair)
                
                # if res is empty
                if not res:
                    res = cycle
                # if res is not empty
                else:
                    #basically insert the new forund cycle
                    # at the correct position
                    pos = res.index(cycle[0])
                    for i in range(len(cycle)-1):
                        res.insert(pos+i+1, cycle[i+1])
        return res
        
def test_db_graph():
    frags = [ "ACC", "ATA", "ATG", "ATT", "CAT", "CAT", "CAT", "CCA",
             "GCA", "GGC", "TAA", "TCA", "TGG", "TTC", "TTT"]
    g = DeBruijnGraph(frags)
    g.print_graph()

test_db_graph()