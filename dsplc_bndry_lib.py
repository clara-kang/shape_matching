import shape_info_lib
import trnsfrm_func_lib
import order_calc_lib
import numpy as np
import math

def getDsplcmentArray(shape_info, curve_vid_1, curve_vid_2, trnsfrm_func) :
    curve_verts_num = len(curve_vid_1)
    dsplcmentArray = np.zeros([curve_verts_num, 2])

    for i in range(curve_verts_num):
        dsplcmentArray[i] = trnsfrm_func_lib.applyTrnsfrm(shape_info.bndry_verts_flat[curve_vid_1[i]], trnsfrm_func) - shape_info.bndry_verts_flat[curve_vid_2[i]]

    return dsplcmentArray

def getDsplcmntAdjSumSquare(shape_info, pair_info, trnsfrm_func_id):
    trnsfrm_func = trnsfrm_func_lib.trsfrm_funcs[trnsfrm_func_id]
    curve_vid_1, curve_vid_2 = shape_info.getPairCurveVids(pair_info)
    curve_verts_num = len(curve_vid_1)

    # print("curve_vid_1: ", curve_vid_1)
    # print("curve_vid_2: ", curve_vid_2)
    dsplcmentArray = getDsplcmentArray(shape_info, curve_vid_1, curve_vid_2, trnsfrm_func)
    # calculate average displacement
    avg_displcmt = np.sum(dsplcmentArray, axis=0) / curve_verts_num
    # displacment has to be integers
    avg_displcmt = np.round(avg_displcmt)

    # if trnsfrm_func_id == 1 or trnsfrm_func_id == 2:
    #     avg_displcmt[0] += 0.5
    dsplcmntAdjArray = dsplcmentArray - avg_displcmt
    dsplcmntSum = np.sum(np.square(dsplcmntAdjArray))
    return (avg_displcmt, dsplcmntSum)

def givenDsplcmntGetAdj(shape_info, curve_vid_1, curve_vid_2, trnsfrm_func, displcmt):

    dsplcmentArray = getDsplcmentArray(shape_info, curve_vid_1, curve_vid_2, trnsfrm_func)
    dsplcmntAdjArray = dsplcmentArray - displcmt
    dsplcmntSum = np.sum(np.square(dsplcmntAdjArray))
    return dsplcmntSum

def inv_rel_func(func):
    return [1 / func[0], -func[1] / func[0]]

def compose_rel_func(func1, func2):
    return [func1[0] * func2[0], func2[0] * func1[1] + func2[1]]

class DsplcCalc:

    def __init__(self, shape_info):
        self.pair_bndry_dsplcmnt = []
        self.shape_info = shape_info
        self.variable_ids = {} # key: vid, value: variable id
        self.var_groups = []
        self.vid_groups = []
        self.vid_group_rel_ptrs = []
        self.vid_pair_rels = {}
        self.group_closed = []
        self.var_values = {}

    def get_pair_bndry_dsplcmnt(self):
        # find matching curve adjustments
        for path_pair in self.shape_info.path_pairs:
            avg_displcmt = np.zeros(2)

            avg_displcmt_min = np.zeros(2)
            dsplcmntSum_min = np.Inf
            trnsfrm_func_id = -1

            # print("path_pair: ", path_pair)
            # try each transfer function
            for i in range(len(trnsfrm_func_lib.trsfrm_funcs)):
                avg_displcmt, dsplcmntSum = getDsplcmntAdjSumSquare(self.shape_info, path_pair, i)
                if dsplcmntSum < dsplcmntSum_min:
                    dsplcmntSum_min = dsplcmntSum
                    avg_displcmt_min = avg_displcmt
                    trnsfrm_func_id = i

            self.pair_bndry_dsplcmnt.append([avg_displcmt_min, trnsfrm_func_id])

    def getVarIds(self, vids):
        var_ids = [0] * len(vids)
        for i in range(len(vids)):
            var_ids[i] = self.variable_ids[vids[i]]
        return var_ids

    def getVidFromVar(self, var_id):
        x_var_id = var_id
        if x_var_id % 2 == 1:
            x_var_id -= 1
        for key in self.variable_ids:
            if self.variable_ids[key] == x_var_id:
                return key

    def get_var_ids(self):
        # x, y are separate vars
        cnt = 0
        for pair_id in range(len(self.shape_info.path_pairs)):
            path_pair = self.shape_info.path_pairs[pair_id]
            # for cid pairs
            for i in range(2):
                eps = self.shape_info.getEPs(path_pair[i])
                # for start and end of curve
                for j in range(2):
                    if not(eps[j] in self.variable_ids):
                        self.variable_ids[eps[j]] = cnt
                        cnt += 2

    def get_trnsfrm_id(self, pair_id):
        dsplcmnt_info = self.pair_bndry_dsplcmnt[pair_id]
        trnsfrm_id = dsplcmnt_info[1]
        return trnsfrm_id

    def get_var_from_vpair(self, vpair):
        x_1 = self.variable_ids[vpair[0]]
        x_2 = self.variable_ids[vpair[1]]
        rel_id = self.vid_pair_rels[vpair]
        # trnsfr func id
        if self.pair_bndry_dsplcmnt[rel_id][1] == 1 or self.pair_bndry_dsplcmnt[rel_id][1] == 2:
            var_pairs = [[x_1, x_2+1], [x_1+1, x_2]]
        else:
            var_pairs = [[x_1, x_2], [x_1+1, x_2+1]]
        return var_pairs

    # get variable pair with curve pair
    def get_var_pair(self, pair_id):
        eps, ep_pairs = self.shape_info.get_ep_pair(pair_id) #[start1, end1], [start2, end2]
        ep_var_ids = self.getVarIds(eps) #[start1_x, end1_x]
        eppair_var_ids = self.getVarIds(ep_pairs) #[start2_x, end2_x]
        trnsfrm_id = self.get_trnsfrm_id(pair_id)

        var_pairs = [None] * 2

        for i in range(2):
            x1_id = ep_var_ids[i] # start1_x, end1_x
            y1_id = x1_id + 1 # start1_y, end1_y
            x2_id = eppair_var_ids[i] # start2_x, end2_x
            y2_id = x2_id + 1 # start2_y, end2_y

            if trnsfrm_id == 0 or trnsfrm_id == 3:
                var_pairs[i] = [[x1_id, x2_id], [y1_id, y2_id]]
            else:
                var_pairs[i] = [[x1_id, y2_id], [y1_id, x2_id]]
        return var_pairs

    def get_group_id(self, var_id, group):
        for i in range(len(group)):
            if var_id in group[i]:
                return i
        return -1

    def spawn_var_group(self, vid_gid):
        vid_group_rel_ptr = self.vid_group_rel_ptrs[vid_gid]
        var_gs = []
        def get_gid(var_id, var_gs):
            for i in range(len(var_gs)):
                if var_id in var_gs[i]:
                    return i
            return -1

        for vpair in vid_group_rel_ptr:
            x_1 = self.variable_ids[vpair[0]]
            x_2 = self.variable_ids[vpair[1]]
            rel_id = self.vid_pair_rels[vpair]
            # trnsfr func id
            if self.pair_bndry_dsplcmnt[rel_id][1] == 1 or self.pair_bndry_dsplcmnt[rel_id][1] == 2:
                var_pairs = [[x_1, x_2+1], [x_1+1, x_2]]
            else:
                var_pairs = [[x_1, x_2], [x_1+1, x_2+1]]
            for i in range(2):
                var_pair = var_pairs[i]
                g1_id = get_gid(var_pair[0], var_gs)
                g2_id = get_gid(var_pair[1], var_gs)
                if g1_id >=0 and g2_id < 0:
                    var_gs[g1_id].append(var_pair[1])
                # merge var_id to g2
                elif g2_id >=0 and g1_id < 0:
                    var_gs[g2_id].append(var_pair[0])
                # both var are in groups
                elif g1_id >=0 and g2_id >= 0:
                    if g1_id != g2_id:
                        # in different group, merge them
                        var_gs[g1_id] = var_gs[g1_id] + var_gs[g2_id]
                        var_gs.pop(g2_id)
                # no relevant group exist
                else:
                    var_gs.append(var_pair)
        # print("")
        return var_gs

    def group_vids(self):
        for pair_id in range(len(self.shape_info.path_pairs)):
            eps, ep_pairs = self.shape_info.get_ep_pair(pair_id)
            # for start and end
            for i in range(2):
                vpair = (eps[i], ep_pairs[i])
                g1_id = self.get_group_id(vpair[0], self.vid_groups)
                g2_id = self.get_group_id(vpair[1], self.vid_groups)

                self.vid_pair_rels[vpair] = pair_id
                # merge var_pair_id to g1
                if g1_id >=0 and g2_id < 0:
                    self.vid_groups[g1_id].append(ep_pairs[i])
                    self.vid_group_rel_ptrs[g1_id].append(vpair)

                # merge var_id to g2
                elif g2_id >=0 and g1_id < 0:
                    self.vid_groups[g2_id].append(eps[i])
                    self.vid_group_rel_ptrs[g2_id].append(vpair)
                # both var are in groups
                elif g1_id >=0 and g2_id >= 0:
                    # if same group, this closes the loop
                    self.vid_group_rel_ptrs[g1_id].append(vpair)
                    if g1_id != g2_id:
                        # in different group, merge them
                        self.vid_groups[g1_id] = self.vid_groups[g1_id] + self.vid_groups[g2_id]
                        self.vid_groups.pop(g2_id)
                        self.vid_group_rel_ptrs[g1_id] = self.vid_group_rel_ptrs[g1_id] + self.vid_group_rel_ptrs[g2_id]
                        self.vid_group_rel_ptrs.pop(g2_id)

                # no relevant group exist
                else:
                    self.vid_groups.append([eps[i], ep_pairs[i]])
                    self.vid_group_rel_ptrs.append([vpair])

    def get_closed_groups(self):
        self.group_closed = [None] * len(self.vid_groups)
        for i in range(len(self.vid_groups)):
            v_num = len(self.vid_groups[i])
            if len(self.vid_group_rel_ptrs[i]) >= v_num:
                self.group_closed[i] = True
            else:
                self.group_closed[i] = False

    def group_vars(self):
        self.var_groups = [None] * len(self.vid_groups)
        for gid in range(len(self.vid_groups)):
            self.var_groups[gid] = self.spawn_var_group(gid)

    def buildGroupConnGraph(self):
        node_num = len(self.var_groups)
        self.group_graph = [None] * node_num

        def insert_from1to2(group_graph, node1, node2, edge):
            if node2 in self.group_graph[node1]:
                if not (edge in self.group_graph[node1][node2]):
                    self.group_graph[node1][node2].append(edge)
            else:
                self.group_graph[node1][node2] = [edge]

        def insert_to_graph(group_graph, node1, node2, edge):
            insert_from1to2(group_graph, node1, node2, edge)
            insert_from1to2(group_graph, node2, node1, edge)

        for i in range(node_num):
            self.group_graph[i] = {}

        for pair_id in range(len(self.shape_info.path_pairs)):
            var_pairs = self.get_var_pair(pair_id)
            eps, ep_pairs = self.shape_info.get_ep_pair(pair_id)


            start_g_id = self.get_group_id(eps[0], self.vid_groups)
            end_g_id = self.get_group_id(eps[1], self.vid_groups)

            insert_to_graph(self.group_graph, start_g_id, end_g_id, pair_id)

    def getOrder(self):
        self.order_calc = order_calc_lib.OrderCalc(self.group_graph, self.group_closed)
        self.order_calc.determine_order()
        self.order_calc.setBrokenEdges()
        # print("self.variable_ids: ", self.variable_ids)
        # print("self.order_calc.ordered_group: ", self.order_calc.ordered_group)
        # print("self.var_groups: ", self.var_groups)
        # print("self.group_closed: ", self.group_closed)
        print("self.order_calc.group_bndry_map: ", self.order_calc.group_bndry_map)
        # self.order_calc.group_bndry_map[0] = 9
        # self.order_calc.group_bndry_map[4] = 8
        # self.order_calc.group_bndry_map[2] = 11
        # self.order_calc.group_bndry_map[6] = 10

    def getRelationFunc(self, ep, ep_pair, dsplcmnt_info):
        rel_funcs = np.zeros([2,2])
        trnsfrm_id = dsplcmnt_info[1]

        # same xy
        if trnsfrm_id == 0:
            trslate_funcs = np.zeros(2)
            for dim in range(2):
                trslate_funcs[dim] = - ( dsplcmnt_info[0][dim] - self.shape_info.bndry_verts_flat[ep][dim] + self.shape_info.bndry_verts_flat[ep_pair][dim] )
            # for x
            rel_funcs[0] = [1, trslate_funcs[0]]
            # for y
            rel_funcs[1] = [1, trslate_funcs[1]]

        # 90 ccw
        elif trnsfrm_id == 1:
            trslate_funcs = np.zeros(2)
            trslate_funcs[0] = - ( dsplcmnt_info[0][1] - self.shape_info.bndry_verts_flat[ep][0] + self.shape_info.bndry_verts_flat[ep_pair][1])
            trslate_funcs[1] = - ( dsplcmnt_info[0][0] + self.shape_info.bndry_verts_flat[ep][1] + self.shape_info.bndry_verts_flat[ep_pair][0])
            # for x
            rel_funcs[0] = [1, trslate_funcs[0]]
            # for y
            rel_funcs[1] = [-1, trslate_funcs[1]]

        # 90 cw
        elif trnsfrm_id == 2:
            trslate_funcs = np.zeros(2)
            trslate_funcs[0] = - ( dsplcmnt_info[0][1] + self.shape_info.bndry_verts_flat[ep][0] + self.shape_info.bndry_verts_flat[ep_pair][1])
            trslate_funcs[1] = - ( dsplcmnt_info[0][0] - self.shape_info.bndry_verts_flat[ep][1] + self.shape_info.bndry_verts_flat[ep_pair][0])
            # for x
            rel_funcs[0] = [-1, trslate_funcs[0]]
            # for y
            rel_funcs[1] = [1, trslate_funcs[1]]

        # same xy, *-1
        elif trnsfrm_id == 3:
            trslate_funcs = np.zeros(2)
            for dim in range(2):
                trslate_funcs[dim] = - ( dsplcmnt_info[0][dim] + self.shape_info.bndry_verts_flat[ep][dim] + self.shape_info.bndry_verts_flat[ep_pair][dim] )
            # for x
            rel_funcs[0] = [-1, trslate_funcs[0]]
            # for y
            rel_funcs[1] = [-1, trslate_funcs[1]]
        return rel_funcs

    def calcVariables(self):
        for same_slot_nodes in self.order_calc.ordered_group:

            for node in same_slot_nodes:
                print("node: ", node)
                pair_rel_funcs = {}
                explcit_ep_ids = {}

                # # if node is odd, then the relationships already found
                # if node % 2 == 0:
                # find pairwise var relationship

                skip_a_rel = False
                for vpair in self.vid_group_rel_ptrs[node]:
                    pair_id = self.vid_pair_rels[vpair]

                    # check if this is broken pair
                    if self.group_closed[node] and node in self.order_calc.group_bndry_map \
                        and self.order_calc.group_bndry_map[node] == pair_id:
                            if len(self.var_groups[node]) == 1:
                                skip_a_rel = True
                            # if two chains, skip both rels
                            else:
                                continue
                    else:
                        skip_a_rel = False

                    var_pairs = self.get_var_from_vpair(vpair)
                    dsplcmnt_info = self.pair_bndry_dsplcmnt[pair_id]
                    rel_funcs = self.getRelationFunc(vpair[0], vpair[1], dsplcmnt_info)

                    # x and y
                    for i in range(2):
                        if i == 0 and skip_a_rel:
                            continue
                        pairs = var_pairs[i]
                        pair_rel_funcs[(pairs[0], pairs[1])] = rel_funcs[i]

                # TODO:: only skip 1 or 2 to break the cycle

                # for gid in range(len(self.var_groups[node])):
                #     vars = self.var_groups[node][gid]
                #     outgoing_pair_ids = self.order_calc.getOutgoingPairs(node)
                #     # print("node: ", node, ", outgoing_pair_ids: ", outgoing_pair_ids)
                #     # process all outgoing pairs
                #     # only break once for a group, in case both x, y are in same group
                #     broken_pair = False
                #     for pair_id in outgoing_pair_ids:
                #         # if this pair is to be skipped
                #         if self.group_closed[node] and node in self.order_calc.group_bndry_map \
                #             and self.order_calc.group_bndry_map[node] == pair_id:
                #             broken_pair = True
                #             # continue
                #         else:
                #             broken_pair = False
                #
                #         eps, ep_pairs = self.shape_info.get_ep_pair(pair_id)
                #         var_pairs = self.get_var_pair(pair_id)
                #         dsplcmnt_info = self.pair_bndry_dsplcmnt[pair_id]
                #
                #         # start and end v
                #         for i in range(2):
                #             pairs = var_pairs[i]
                #             # only one endpoint is in the group
                #             if pairs[0][0] in self.var_groups[node][gid] or pairs[1][0] in self.var_groups[node][gid]:
                #                 rel_funcs = self.getRelationFunc(eps[i], ep_pairs[i], dsplcmnt_info)
                #                 # x and y
                #                 for j in range(2):
                #                     # skip first dimension
                #                     if broken_pair and j == 0:
                #                         # print("skip pair: ", (pairs[j][0], pairs[j][1]))
                #                         continue
                #                     pair_rel_funcs[(pairs[j][0], pairs[j][1])] = rel_funcs[j]
                for key in pair_rel_funcs:
                    print(key, ", ", pair_rel_funcs[key])
            # print("over")
                # find relationship between the first of the group to the rest of a group, do search
                for gid in range(len(self.var_groups[node])):
                    # first var id in group
                    cur_var = self.var_groups[node][gid][0]
                    explct_var = cur_var

                    # create an entry in explcit_ep_ids
                    explcit_ep_ids[cur_var] = {}
                    next_vars = []
                    for i in range(1, len(self.var_groups[node][gid])):
                        rel_var = self.var_groups[node][gid][i]
                        # go through all pairs
                        for var_pair in pair_rel_funcs:
                            # directly related
                            if cur_var in var_pair and rel_var in var_pair:
                                # add to next vars to process
                                next_vars.append(rel_var)
                                # no need to invert rel_func
                                if var_pair.index(rel_var) == 1:
                                    explcit_ep_ids[explct_var][rel_var] = pair_rel_funcs[var_pair]
                                else:
                                    explcit_ep_ids[explct_var][rel_var] = inv_rel_func(pair_rel_funcs[var_pair])

                    done_vars = [cur_var]
                    while len(next_vars) > 0:
                        cur_var = next_vars.pop(0)
                        # go through all pairs
                        for var_pair in pair_rel_funcs:
                            if cur_var in var_pair:
                                if var_pair.index(cur_var) == 0:
                                    rel_var = var_pair[1]
                                else:
                                    rel_var = var_pair[0]

                                # already processed this pair
                                if rel_var in done_vars:
                                    continue
                                else:
                                    # process rel_var later
                                    next_vars.append(rel_var)
                                    if var_pair.index(cur_var) == 0:
                                        rel_func = pair_rel_funcs[var_pair]
                                    else:
                                        rel_func = inv_rel_func(pair_rel_funcs[var_pair])

                                    # compose function
                                    cmpsd_func = compose_rel_func(explcit_ep_ids[explct_var][cur_var], rel_func)
                                    explcit_ep_ids[explct_var][rel_var] = cmpsd_func
                        # done with cur_var
                        done_vars.append(cur_var)

                print("explcit_ep_ids: ", explcit_ep_ids)

                def roundPointFiveZero(num):
                    is_neg = num < 0
                    num_abs = abs(num)
                    frac, whole = math.modf(num_abs)
                    if abs(frac - 0.5) < frac and abs(frac - 0.5) < 1 - frac:
                        num_abs = whole + 0.5
                    else:
                        num_abs = round(num_abs)
                    if is_neg:
                        return -num_abs
                    return num_abs

                # def roundZero(num):
                #     is_neg = num < 0
                #     num_abs = abs(num)
                #     frac, whole = math.modf(num_abs)
                #     return whole
                def roundVarVal(shape_info, var_values, val, var_id):
                    vid = self.getVidFromVar(var_id)
                    if var_id % 2 == 1:
                        othr_var_id = var_id - 1
                    else:
                        othr_var_id = var_id + 1
                    vid_val = shape_info.bndry_verts_flat[vid][var_id % 2]
                    othr_vid_val = shape_info.bndry_verts_flat[vid][othr_var_id % 2]
                    adjusted_val = vid_val + val
                    # print("var_id: ", var_id)
                    # print("adjusted_val: ", adjusted_val)
                    # other dim already defined
                    if othr_var_id in var_values:
                        othr_adjsted_val = othr_vid_val + var_values[othr_var_id]
                        if math.modf(othr_adjsted_val)[0] == 0:
                            return round(adjusted_val) - vid_val
                        elif math.modf(othr_adjsted_val)[0] == 0.5:
                            return math.floor(adjusted_val) + 0.5 - vid_val
                    else:
                        rounded_adjusted_val = roundPointFiveZero(adjusted_val)
                        # print("rounded_adjusted_val: ", rounded_adjusted_val)
                        return rounded_adjusted_val - vid_val

                def getGroupId(var_groups, var_id):
                    for i in range(len(var_groups)):
                        for j in range(len(var_groups[i])):
                            if var_id in var_groups[i][j]:
                                return i
                    return -1

                def cal_vars(obj, var_values, explcit_ep_ids):
                    # calculate the explicit vars
                    for key in explcit_ep_ids:
                        coeff_1 = len(explcit_ep_ids[key]) + 1
                        coeff_2 = 0
                        for rel_var in explcit_ep_ids[key]:
                            coeff_2 += explcit_ep_ids[key][rel_var][1] * explcit_ep_ids[key][rel_var][0]
                        var_value = - coeff_2 / coeff_1
                        group_id = getGroupId(obj.var_groups, key)
                        if (obj.group_closed[group_id]):
                            # print("var_id: ", key, ", before var_value", var_value)
                            var_value = roundVarVal(obj.shape_info, var_values, var_value, key)
                            # print("var_id: ", key, ", after var_value", var_value)
                        # print("var_value: ", var_value)
                        var_values[key] = var_value

                    # calculate the implicit vars
                    for key in explcit_ep_ids:
                        for rel_var in explcit_ep_ids[key]:
                            var_value = var_values[key] *  explcit_ep_ids[key][rel_var][0] + explcit_ep_ids[key][rel_var][1]
                            var_values[rel_var] = var_value

                cal_vars(self, self.var_values, explcit_ep_ids)
                # print("self.var_values: ", self.var_values)

                # find relationships between broken pairs
                if self.group_closed[node] and node in self.order_calc.group_bndry_map:
                    broken_pair = self.order_calc.group_bndry_map[node]
                    print("node: ", node)
                    print("broken pair: ", broken_pair)
                    broken_pair_info = self.shape_info.path_pairs[broken_pair]
                    # var_pairs = self.get_var_pair(broken_pair)
                    eps, ep_pairs = self.shape_info.get_ep_pair(broken_pair)

                    # start or end of pair?
                    this_ep_id = 1
                    othr_ep_id = 0
                    if eps[0] in self.vid_groups[node]:
                    # if var_pairs[0][0][0] in self.var_groups[node][gid] or var_pairs[0][1][0] in self.var_groups[node][gid]:
                        this_ep_id = 0
                        othr_ep_id = 1
                    print("othr_ep_id: ", othr_ep_id)
                    print("this_ep_id: ", this_ep_id)

                    vid_1 = eps[this_ep_id]
                    vid_2 = ep_pairs[this_ep_id]
                    print("vid_1: ", vid_1)
                    print("vid_2: ", vid_2)
                    vid_1_dsplcmnt = np.array([self.var_values[self.variable_ids[vid_1]], self.var_values[self.variable_ids[vid_1]+1] ])
                    vid_2_dsplcmnt = np.array([self.var_values[self.variable_ids[vid_2]], self.var_values[self.variable_ids[vid_2]+1] ])

                    # print("vid_1_dsplcmnt: ", vid_1_dsplcmnt)
                    # print("vid_2_dsplcmnt: ", vid_2_dsplcmnt)
                    vid_1_loc = self.shape_info.bndry_verts_flat[vid_1] + vid_1_dsplcmnt
                    vid_2_loc = self.shape_info.bndry_verts_flat[vid_2] + vid_2_dsplcmnt
                    # print("diff: ", vid_1_loc - vid_2_loc)
                    curve_vid_1, curve_vid_2 = self.shape_info.getPairCurveVids(broken_pair_info)
                    # remove the endpoint
                    if this_ep_id == 0:
                        curve_vid_1 = curve_vid_1[1:]
                        curve_vid_2 = curve_vid_2[1:]
                    else:
                        curve_vid_1 = curve_vid_1[:len(curve_vid_1)-1]
                        curve_vid_2 = curve_vid_2[:len(curve_vid_2)-1]
                    # print("curve_vid_1: ", curve_vid_1)

                    # get best trnsfr function knowing the location of one pair of endpoints
                    best_trnsfr_func_id = -1
                    best_dsplcmnt = None
                    min_adj_sum = np.Inf

                    for trnsfr_func_id in range(4):
                        trsfrm = trnsfrm_func_lib.trsfrm_funcs[trnsfr_func_id]
                        dsplcmnt = trnsfrm_func_lib.applyTrnsfrm(vid_1_loc, trsfrm) - vid_2_loc
                        # print("dsplcmnt: ", dsplcmnt)
                        # if not (dsplcmnt[0].is_integer() and dsplcmnt[1].is_integer()):
                        #     continue

                        adj_sum = givenDsplcmntGetAdj(self.shape_info, curve_vid_1, curve_vid_2, trsfrm, dsplcmnt)
                        # print("dsplcmnt: ", dsplcmnt)
                        if adj_sum < min_adj_sum:
                            min_adj_sum = adj_sum
                            best_trnsfr_func_id = trnsfr_func_id
                            best_dsplcmnt = dsplcmnt
                    # print("before self.pair_bndry_dsplcmnt[broken_pair]: ", self.pair_bndry_dsplcmnt[broken_pair])
                    self.pair_bndry_dsplcmnt[broken_pair] = [best_dsplcmnt, best_trnsfr_func_id]
                    print("after self.pair_bndry_dsplcmnt[broken_pair]: ", self.pair_bndry_dsplcmnt[broken_pair])
                    # get the other group affected
                    othr_group_id = self.get_group_id(eps[othr_ep_id], self.vid_groups)
                    self.var_groups[othr_group_id] = self.spawn_var_group(othr_group_id)
                    print("new var groups: ")
                    for i in range(len(self.var_groups)):
                        print(i, ": ", self.var_groups[i])


    def intrp_curve(self, cid):
        start_vid, end_vid = self.shape_info.getEPs(cid)

        end_dsplcmnt = self.shape_info.bndry_vert_uv[end_vid] - self.shape_info.bndry_verts_flat[end_vid]
        start_dsplcmnt = self.shape_info.bndry_vert_uv[start_vid] - self.shape_info.bndry_verts_flat[start_vid]
        curve_v_num = len(self.shape_info.bndry_verts[cid])

        for i in range(1, curve_v_num):
            v_dsplcmnt = ( i / curve_v_num ) * end_dsplcmnt + ( 1 - i / curve_v_num ) * start_dsplcmnt
            self.shape_info.bndry_vert_uv[start_vid + i] = self.shape_info.bndry_verts_flat[start_vid + i] + v_dsplcmnt

    def itrpRestOfCurves(self):
        # move the endpoints first
        for vid in self.variable_ids:
            self.shape_info.bndry_vert_uv[vid][0] = self.shape_info.bndry_verts_flat[vid][0] + self.var_values[self.variable_ids[vid]]
            self.shape_info.bndry_vert_uv[vid][1] = self.shape_info.bndry_verts_flat[vid][1] + self.var_values[self.variable_ids[vid]+1]

        # move the curve pairs
        for pair_id in range(len(self.shape_info.path_pairs)):
            eps = self.shape_info.getEPs(self.shape_info.path_pairs[pair_id][0])
            curve_vid_1, curve_vid_2 = self.shape_info.getPairCurveVids(self.shape_info.path_pairs[pair_id])
            c1_len = len(curve_vid_1)
            dsplcmnt_info = self.pair_bndry_dsplcmnt[pair_id]

            self.intrp_curve(self.shape_info.path_pairs[pair_id][0])

            # use rel transform + translate to get second curve
            for t in range(1, c1_len-1):
                c1_vid = curve_vid_1[t]
                c2_vid = curve_vid_2[t]
                self.shape_info.bndry_vert_uv[c2_vid] = trnsfrm_func_lib.applyTrnsfrm(self.shape_info.bndry_vert_uv[c1_vid], trnsfrm_func_lib.trsfrm_funcs[dsplcmnt_info[1]]) - dsplcmnt_info[0]

        # gather indices of curves needed to be interpolated

        cid_intrp = set()
        cid_pair = set()

        for path_pair in self.shape_info.path_pairs:
            cid_pair.add(path_pair[0])
            cid_pair.add(path_pair[1])

        for path_pair in self.shape_info.path_pairs:
            for pid in range(2):
                next_cid = self.shape_info.getNextCurveInShape(path_pair[pid])
                prev_cid = self.shape_info.getPrevCurveInShape(path_pair[pid])
                if not (next_cid in cid_pair):
                    cid_intrp.add(next_cid)
                if not (prev_cid in cid_pair):
                    cid_intrp.add(prev_cid)

        # do the interpolation
        for cid in cid_intrp:
            self.intrp_curve(cid)
