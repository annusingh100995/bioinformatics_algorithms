# Tries are n-ary trees that allow to organize a set of
# string patterns

# {node: {character : destination node}}

class Trie:
    
    def __init__(self):
        # root node
        self.nodes = {0:{}}
        self.num = 0
        
    def print_trie(self):
        for k in self.nodes.keys():
            print(k, "-->", self.nodes[k])
    
    def add_node(self,origin,symbol):
        self.num += 1
        self.nodes[origin][symbol] = self.num
        self.nodes[self.num] = {}
        
    """From origin we are moving to node numbered self.num
    via the edge symbol.
    In tries edges are added when a node is created.Every
    node has to have a predecessor.
    There are no nodes hanging around alone in tries."""
    
    def add_pattern(self,p):
        pos = 0
        node = 0
        while pos < len(p):
            if p[pos] not in self.nodes[node].keys():
                #print("At pattern : ",p," and at position ",p[pos])
                self.add_node(node,p[pos])
            node = self.nodes[node][p[pos]]
            pos += 1
            
    def tries_from_pattern(self,patterns):
        for p in patterns:
            self.add_pattern(p)
    
    """We start at the root node.
    If the text[pos] is one of any edges going out of root
    then a this text[pos] to the match sequence.
    And also update the node (current node) to the
    destination node of text[pos]
    
    If there is not edge in the next node
    then return the match so far.
    else move one pos ahead in text
    """
    def prefix_trie_match(self,text):
        pos = 0
        match = ""
        node = 0
        
        while pos < len(text):
            if text[pos] in self.nodes[node].keys():
                node = self.nodes[node][text[pos]]
                match += text[pos]
                if self.nodes[node] == {}:
                    return match
                else:
                    pos += 1
            else:
                return None
        return None
    
    def trie_matches(self,text):
        res = []
        for i in range(len(text)):
            m = self.prefix_trie_match(text[i:])
            if m != None:
                res.append((i,m))
        return res

def test_tries():
    patterns = ["GAT","CCT","GAG"]
    t = Trie()
    t.tries_from_pattern(patterns)
    t.print_trie()

test_tries()

def test_prefix_find():
    patterns = ["AGAGAT", "AGC", "AGTCC", "CAGAT", "CCTA", "GAGAT", "GAT", "TC"]
    t = Trie()
    t.tries_from_pattern(patterns)
    #t.print_trie()
    print(t.prefix_trie_match("GAGATCCTA"))
    print(t.trie_matches("GAGATCCTA"))

test_prefix_find()
