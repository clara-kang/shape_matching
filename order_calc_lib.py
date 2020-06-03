import numpy as np

class OrderCalc:
    def __init__(self, graph, group_closed):
        self.graph = graph
        self.group_closed = group_closed
        self.visited = np.zeros(len(self.group_closed), dtype=bool)
        self.ordered_group = []
        self.group_bndry_map = {}

    def findClosedNode(self):
        self_closed_nodes = []
        closed_nodes = []
        for idx in range(len(self.group_closed)):
            if self.group_closed[idx] and not (self.visited[idx]):
                if len(self.graph[idx]) == 1:
                    self.visited[idx] = True
                    self_closed_nodes.append(idx)
                else:
                    closed_nodes.append(idx)
        if len(self_closed_nodes) > 0:
            return self_closed_nodes
        elif len(closed_nodes) > 0:
            idx = closed_nodes[0]
            self.visited[idx] = True
            return [idx]
        return []

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

        # print("start_node: ", start_node)
        for key in graph_entry:
            # print("key: ", key, ", visited: ", self.visited[key])
            if not (self.visited[key]):
                edges.extend(graph_entry[key])
                # the open groups will go the last time slot
                if self.group_closed[key] :
                    self.visited[key] = True
                    nbs.append(key)
                else:
                    preferred_edges.extend(graph_entry[key])

        # print("edges: ", edges)
        # if len(preferred_edges) > 0:
        #     self.group_bndry_map[start_node] = preferred_edges[0]
        # elif len(edges) > 0:
        #     self.group_bndry_map[start_node] = edges[0]
        return nbs

    def getOrder(self, node):
        for order in range(len(self.ordered_group)):
            if node in self.ordered_group[order]:
                return order
        return -1

    def move2NextOrder(self, node, order):
        self.ordered_group[order].remove(node)
        self.add_to_ordered_group([node], order+1)

    def setBrokenEdge(self, node, order, node_has_order):
        graph_entry = self.graph[node]
        # order = self.ordered_group
        nb_orders = {}
        for key in graph_entry:
            nb_ordr = self.getOrder(key)
            nb_orders[key]=nb_ordr
        for key in graph_entry:
            if nb_orders[key] > order:
                self.group_bndry_map[node] = graph_entry[key][0]
                return False
        # did not find nb of bigger order
        for key in graph_entry:
            # nb same time slot, but not processed yet
            if nb_orders[key] == order and node_has_order[key] == False:
                # move nb to next order
                self.move2NextOrder(key, order)
                self.group_bndry_map[node] = graph_entry[key][0]
                return True

    def setBrokenEdges(self):
        node_has_order = [False] * len(self.graph)
        for order in range(len(self.ordered_group)):
            while True:
                processed_num = 0
                for node in self.ordered_group[order]:
                    if not node_has_order[node]:
                        node_has_order[node] = True
                        processed_num += 1
                        order_changed = self.setBrokenEdge(node, order, node_has_order)
                        if order_changed:
                            break
                # finished with the set
                if processed_num == 0:
                    break

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
            start_nodes = self.findClosedNode()
            # print("start_nodes: ", start_nodes)
            if len(start_nodes) == 0:
                break
            next_layer = start_nodes
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
