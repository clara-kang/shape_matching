import numpy as np

class OrderCalc:
    def __init__(self, graph, group_closed):
        self.graph = graph
        self.group_closed = group_closed
        self.visited = np.zeros(len(self.group_closed), dtype=bool)
        self.ordered_group = []
        self.group_bndry_map = {}

    def findClosedNode(self):
        for idx in range(len(self.group_closed)):
            if self.group_closed[idx] and not (self.visited[idx]):
                self.visited[idx] = True
                return idx
        return -1

    def getOpenNodes(self):
        open_nodes = []
        for idx in range(len(self.group_closed)):
            if not (self.group_closed[idx]):
                open_nodes.append(idx)
        return open_nodes

    def getNbs(self, start_node):
        graph_entry = self.graph[start_node]
        nbs = []
        edges = []
        preferred_edges = [] # should it be used??

        for key in graph_entry:
            if not (self.visited[key]):
                edges.extend(graph_entry[key])
                # the open groups will go the last time slot
                if self.group_closed[key] :
                    self.visited[key] = True
                    nbs.append(key)
                else:
                    preferred_edges.extend(graph_entry[key])

        if len(preferred_edges) > 0:
            self.group_bndry_map[start_node] = preferred_edges[0]
        elif len(edges) > 0:
            self.group_bndry_map[start_node] = edges[0]
        return nbs

    def getOutgoingPairs(self, start_node):
        graph_entry = self.graph[start_node]
        edges = []

        for key in graph_entry:
            edges.extend(graph_entry[key])

        return edges

    def getNextLayer(self, this_layer):
        next_layer = set({})
        for v in this_layer:
            nbs = self.getNbs(v)
            for nb in nbs:
                next_layer.add(nb)
        return next_layer

    def add_to_ordered_group(self, nodes, slot):
        if len(self.ordered_group) < slot+1:
            for i in range(slot+1 - len(self.ordered_group)):
                self.ordered_group.append(set({}))
        for node in nodes:
            self.ordered_group[slot].add(node)

    def determine_order(self):
        while True:
            start_node_id = self.findClosedNode()
            if start_node_id < 0:
                break
            next_layer = [start_node_id]
            # print("next_layer: ", next_layer)
            slot = 0
            while len(next_layer) > 0:
                self.add_to_ordered_group(next_layer, slot)
                slot += 1
                this_layer = next_layer
                next_layer = self.getNextLayer(next_layer)
                # print("next_layer: ", next_layer)
        # add all open nodes in the end
        self.add_to_ordered_group(self.getOpenNodes(), len(self.ordered_group))
