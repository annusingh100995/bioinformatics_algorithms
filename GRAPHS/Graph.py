class MyGraph:
    
    def __init__(self,g ={}):
        self.graph = g
    
    def get_nodes(self):
        return list(self.graph.keys())
    
    def get_edges(self):
        edges = []
        for v in self.graph.keys():
            for d in self.graph[v]:
                edges.append((v,d))
        return edges
    
    def size(self):
        return len(self.get_nodes()) , len(self.get_edges())
    
    def print_graph(self):
        for v in self.graph.keys():
            print(v,"->", self.graph[v])
    
    def add_vertex(self,v):
        if v not in self.graph.keys():
            self.graph[v] = []
    
    #Adding edge, o: origin node, d:destination node
    def add_edge(self,o,d):
        if o not in self.graph.keys():
            self.add_vertex(o)
        if d not in self.graph.keys():
            self.add_vertex(d)
        # if the edge doesn't already exist, add the new edge
        if d not in self.graph[o]:
            self.graph[o].append(d)
    
    # for an ordered pair (s,v) in a graph
    # s is the predecessor of v, v is the successor
    def get_successors(self,v):
        return list(self.graph[v])
    
    def get_predecessors(self,v):
        result = []
        for k in self.graph.keys():
            if v in self.graph[k]:
                result.append(k)
        return result
    
    # Two nodes are adjacent if one is successor of another
    def get_adjacents(self,v):
        successors = self.get_successors(v)
        predessors = self.get_predecessors(v)
        result = predessors
        for s in successors:
            if s not in result:
                result.append(s)
        
        return result
    
    def out_degree(self,v):
        return len(self.graph[v])
    
    def in_degree(self,v):
        return len(self.get_predecessors(v))
    
    def degree(self,v):
        return (len(self.get_adjacents(v)))
    
    # Compute degree of all the nodes
    # deg_type can be "in", "out" or "inout"
    # Returns dictionary node ---> degree
    def all_degrees(self,deg_type= 'inout'):
        degrees = {}
        # for eavery node out degree is equal to number of 
        # connection it has to other node 3 -> [2, 4]
        # len(list) for each node 
        for v in self.graph.keys():
            if deg_type == "out" or deg_type == "inout":
                degrees[v] = len(self.graph[v])
            else:
                degrees[v] = 0
        if deg_type == "in" or deg_type == "inout":
            for v in self.graph.keys():
                for d in self.graph[v]:
                    if deg_type == "in" or v not in self.graph[d]:
                        degrees[d] = degrees[d]+1
        return degrees
    
    """BFS start by the source nodem then visit all it successors,
    folloed by thier sucessors , until all possible nodes are visited
    """
    
    def breadth_first_search(self,v):
        visited = [v]
        res = []
        while len(visited) > 0:
            node = visited.pop(0)
            if node != v:
                res.append(node)
            for element in self.graph[node]:
                if element not in res and element not in visited:
                    visited.append(element)
        return res
    
    """DFS: start at teh source node, explore first the successors
    tehen its first successor unitl no further exploration is possible
    then back trake to explore further alternatives"""
    
    def depth_first_search(self,v):
        visited = [v]
        res = []
        while len(visited)>0:
            node = visited.pop(0)
            if node != v:
                res.append(node)
            s = 0
            for elem in self.graph[node]:
                if elem not in res and elem not in visited:
                    visited.insert(s,elem)
                    s += 1
        return res
    
    # Finding the distance between the source and destination node
    def distance(self,s,d):
        if s == d:
            return 0
        l = [(s,0)]
        visited = [s]
        while len(l) > 0:
            node , dist = l.pop(0)
            for element in self.graph[node]:
                if element == d:
                    return dist+1
                elif element not in visited:
                    l.append((element, dist+1))
                    visited.append(element)
        return None
    
    
    # Find shortest path
    def shortest_path(self,s,d):
        if s == d:
            return 0
        l = [(s,[])]
        visited = [s]
        while len(l) > 0:
            node,preds = l.pop(0)
            for element in self.graph[node]:
                if element == d:
                    return preds+[node,element]
                elif element not in visited:
                    l.append((element,preds+[node]))
                    visited.append(element)
                    
        return None
    
    def is_in_tuple_list(self,tl,value):
        res = False
        for (x,y) in tl:
            if x == value:
                return True
        return res
    
    def reachable_with_distance(self,s):
        res = []
        l = [(s,0)]
        while len(l) > 0:
            node , dist = l.pop(0)
            if node != s:
                res.append((node,dist))
            for element in self.graph[node]:
                if not self.is_in_tuple_list(l,element) and not self.is_in_tuple_list(res,element):
                    l.append((element,dist+1))
        return res
    
    
    def node_has_cycle(self,v):
        l = [v]
        res = False
        visited = [v]
        while len(l)>0:
            node = l.pop()
            for element in self.graph[node]:
                if element == v:
                    return True
                elif element not in visited:
                    l.append(element)
                    visited.append(element)
        return res
    
    def has_cycle(self):
        res = False
        for v in self.graph.keys():
            if self.node_has_cycle(v):
                return True
        return res
    
    def is_connected(self):
        total = len(self.graph.keys())
        for v in self.graph.keys():
            reachable_node = self.breadth_first_search(v)
            if (len(reachable_node) < total):
                return False
        return True


def test_graph():
    gr = MyGraph()
    gr.add_vertex(1)
    gr.add_vertex(2)
    gr.add_vertex(3)
    gr.add_vertex(4)
    gr.add_edge(1,2)
    gr.add_edge(2,3)
    gr.add_edge(3,2)
    gr.add_edge(3,4)
    gr.add_edge(4,2)
    gr.print_graph()
    print("Successors of node 2: ",gr.get_successors(2))
    print("Predecessors of node 2: ",gr.get_predecessors(2))
    print("Adjacents of node 2: ",gr.get_adjacents(2))
    print("In degree of 2 : ", gr.in_degree(2))
    print("Out degree of 2 : ", gr.out_degree(2))
    print("Degree of 2: ", gr.degree(2))
    gr2 = MyGraph({1:[2,3,4],2:[5,6],3:[6,8],4:[8],5:[7],6:[],7:[],8:[]})
    print (gr2.breadth_first_search(1))
    print (gr2.depth_first_search(1))
    print (gr2.distance(1,7))
    print (gr2.shortest_path(1,7))
    print (gr2.distance(1,8))
    print (gr2.shortest_path(1,8))
    print (gr2.distance(6,1))
    print (gr2.shortest_path(6,1))
    gr2.print_graph()
    print(gr2.reachable_with_distance(2))
    print(gr.has_cycle())
    print(gr2.has_cycle())

test_graph()

def graph_distance_tests():
    gr2 = MyGraph({1:[2,3,4],2:[5,6],3:[6,8],4:[8],5:[7],6:[],7:[],8:[]})
    gr2.print_graph()
    print("BFS: ",gr2.breadth_first_search(1))
    print("DFS: ",gr2.depth_first_search(1))
    print("Distance b/w 1 and 7: ", gr2.distance(1,7))
    print("Shortest Path b/w 1 and 7", gr2.shortest_path(1,7))
    print("Reachable with distance fron 1: ", gr2.reachable_with_distance(1))
    print("Has cycle ?: ", gr2.has_cycle())
    print("Is connected ? :", gr2.is_connected())
    print("Adjecents of 2: ", gr2.get_adjacents(2))
    

graph_distance_tests()