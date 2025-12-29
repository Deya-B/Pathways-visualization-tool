from collections import defaultdict
import networkx as nx

from .model import ConversionInteraction
from .idgen import idgenerator

# Layout/Board Constants
BOARD_MARGIN = 200    # margin around everything
CENTERING_FACTOR = 0.8  # between 0.0 (no centering) and 1.0 (full centering)
LAYER_GAP = 120.0       # vertical separation
COL_GAP = 250.0         # approximate horizontal separation

ENZYME_OFFSET_X = 180.0     # horizontal distance from anchor to enzyme
ENZYME_OFFSET_Y = 40.0      # horizontal distance from anchor to enzyme
ENZYME_STACK_GAP = 22.0   # vertical separation when several enzymes share an anchor

class Layout:
    def __init__(self, nodes, interactions):
        self.nodes = nodes                # Dict of nodes, keyed by graph_id
        self.interactions = interactions  # List of Interactions
        self.boardwidth = None
        self.boardheight = None
        self._laid_out = False
        self._parent_of = {}
        self._node_depth = {}


    def compute_anchor_xy(self, interactions):
        """Compute anchor (ax, ay) for a given Interaction. Based on the
        source and target nodes positions."""
        source = interactions.source
        target = interactions.target

        # bottom border of source and top border of target
        x1 = source.x; y1 = source.y + source.height / 2.0
        x2 = target.x; y2 = target.y - target.height / 2.0

        ax = x1 + interactions.anchor_pos * (x2 - x1)
        ay = y1 + interactions.anchor_pos * (y2 - y1)

        if interactions.anchor_id is None:
            interactions.anchor_id = idgenerator.new("anchor")

        return (ax, ay)


    def layout_positions(self):
        """Pre-layout positions: Assign x,y to metabolites and pathway nodes only."""
        G = self._build_skeleton_graph()
        self._parent_of = self._compute_parent_map(G)
        in_deg, out_deg = self._precompute_degrees(G)

        sccs, node2scc, CG, topo_scc = self._build_scc_condensation(G)
        depth_scc, scc_layers = self._compute_scc_depths(CG, topo_scc)
        layers, node_depth = self._build_node_layers(sccs, scc_layers)

        placed_x, placed_scc_x = self._place_layers(
            G, CG, sccs, node2scc, layers, in_deg, out_deg
        )

        self._fix_vertical_single_parent_children(G, CG, node2scc)
        self._place_remaining_enzymes(placed_x, layers)
        self._node_depth = node_depth
        self._layout_conversions()
        self._shift_layout_into_view()
        self._laid_out = True


# Graph construction helpers
    def _build_skeleton_graph(self):
        """Skeleton graph with metabolite/pathway nodes and conversion edges only."""
        G = nx.DiGraph()
        for node in self.nodes.values():
            if node.node_type == "GeneProduct":
                continue
            G.add_node(node.graph_id)
        for inter in self.interactions:
            if isinstance(inter, ConversionInteraction):
                G.add_edge(inter.source.graph_id, inter.target.graph_id)
        return G


    def _compute_parent_map(self, G):
        """First parent per node (used later for vertical adjustments)."""
        parent_of = {}
        for u, v in G.edges:
            if v not in parent_of:
                parent_of[v] = u
        return parent_of


    def _precompute_degrees(self, G):
        """Precompute in/out degree maps for chain detection."""
        in_deg = {n: G.in_degree(n) for n in G.nodes}
        out_deg = {n: G.out_degree(n) for n in G.nodes}
        return in_deg, out_deg


    def _build_scc_condensation(self, G):
        """SCCs and condensation DAG."""
        sccs = list(nx.strongly_connected_components(G))
        node2scc = {}
        for idx, comp in enumerate(sccs):
            for n in comp:
                node2scc[n] = idx
        CG = nx.condensation(G, sccs)
        topo_scc = list(nx.topological_sort(CG))
        return sccs, node2scc, CG, topo_scc


    def _compute_scc_depths(self, CG, topo_scc):
        """Depth per SCC (longest path) and SCC layers."""
        depth_scc = {c: 0 for c in topo_scc}
        for u in topo_scc:
            for v in CG.successors(u):
                depth_scc[v] = max(depth_scc[v], depth_scc[u] + 1)

        max_depth = max(depth_scc.values(), default=0)
        scc_layers = [[] for _ in range(max_depth + 1)]
        for comp_id, d in depth_scc.items():
            scc_layers[d].append(comp_id)
        return depth_scc, scc_layers


    def _build_node_layers(self, sccs, scc_layers):
        """Expand SCC layers to node layers and assign base node depths."""
        layers = []
        node_depth = {}
        for d, comps in enumerate(scc_layers):
            layer_nodes = []
            for comp_id in comps:
                members = list(sccs[comp_id])
                layer_nodes.extend(members)
                if len(members) == 1:
                    node_depth[members[0]] = float(d)
                else:
                    for n in members:
                        node_depth[n] = float(d)
            layers.append(layer_nodes)
        return layers, node_depth


# Placement helpers
    def _place_layers(self, G, CG, sccs, node2scc, layers, in_deg, out_deg):
        """Horizontal spacing, SCC ordering, SCC and node placement."""
        k_max = max((len(L) for L in layers), default=1)
        span_max = max(0, (k_max - 1) * COL_GAP)

        placed_x = {}
        placed_scc_x = {}
        scc_child_count = defaultdict(int)

        scc_local_depth = self._compute_scc_local_depth(CG)

        for ly, layer_nodes in enumerate(layers):
            if not layer_nodes:
                continue
            self._place_single_layer(
                G, CG, sccs, node2scc, layer_nodes,
                ly, span_max, in_deg, out_deg,
                scc_local_depth, placed_x, placed_scc_x, scc_child_count
            )

        return placed_x, placed_scc_x


    def _compute_scc_local_depth(self, CG):
        """Local depth of each SCC within its parent’s subgraph."""
        scc_local_depth = defaultdict(int)
        for s in CG.nodes:
            if CG.out_degree(s) == 0:
                continue
            queue = [(s, 0)]
            seen = {s}
            while queue:
                u, d = queue.pop(0)
                for v in CG.successors(u):
                    if v not in seen:
                        seen.add(v)
                        scc_local_depth[v] = max(scc_local_depth[v], d + 1)
                        queue.append((v, d + 1))
        return scc_local_depth


    def _place_single_layer(self, G, CG, sccs, node2scc, layer_nodes, ly, span_max,
        in_deg, out_deg, scc_local_depth, placed_x, placed_scc_x, scc_child_count):
        """Place nodes in a single layer."""
        k = len(layer_nodes)
        span = max(0, (k - 1) * COL_GAP)
        start_x = BOARD_MARGIN + CENTERING_FACTOR * (span_max - span) / 2.0
        y = 80 + ly * LAYER_GAP

        # Group nodes in this layer by SCC
        by_scc = {}
        for n in layer_nodes:
            by_scc.setdefault(node2scc[n], []).append(n)

        def scc_parent_avg_x(scc_id):
            parent_sccs = list(CG.predecessors(scc_id))
            xs = []
            for ps in parent_sccs:
                if ps in placed_scc_x:
                    xs.append(placed_scc_x[ps])
                else:
                    for n in sccs[ps]:
                        preds = list(G.predecessors(n))
                        xs.extend(placed_x[p] for p in preds if p in placed_x)
            return sum(xs) / len(xs) if xs else float("inf")

        scc_ids_in_layer = list(by_scc.keys())
        scc_ids_in_layer.sort(key=scc_parent_avg_x)

        i = 0
        for comp_id in scc_ids_in_layer:
            members = by_scc[comp_id]
            size = len(members)
            col_x = start_x + i * COL_GAP

            center_x = self._choose_scc_center_x(
                CG, comp_id, col_x, scc_ids_in_layer,
                scc_local_depth, placed_scc_x, scc_child_count
            )

            if size == 1:
                self._place_singleton_in_scc(
                    G, members[0], center_x, y,
                    in_deg, out_deg, placed_x
                )
            else:
                self._layout_scc_as_tree(
                    members, G,
                    center_x=center_x,
                    base_y=y
                )
                for m in members:
                    placed_x[m] = self.nodes[m].x

            placed_scc_x[comp_id] = center_x
            i += 1


    def _choose_scc_center_x(
        self, CG, comp_id, col_x, scc_ids_in_layer,
        scc_local_depth, placed_scc_x, scc_child_count
    ):
        parent_sccs = list(CG.predecessors(comp_id))
        if len(parent_sccs) == 1 and parent_sccs[0] in placed_scc_x:
            p = parent_sccs[0]
            siblings = [cid for cid in scc_ids_in_layer if p in CG.predecessors(cid)]
            siblings.sort(key=lambda cid: (scc_local_depth[cid], cid))
            idx = siblings.index(comp_id)
            if idx == 0:
                center_x = placed_scc_x[p]
            else:
                offset_index = (idx + 1) // 2
                offset = offset_index * COL_GAP
                sign = -1 if idx % 2 == 1 else 1
                center_x = placed_scc_x[p] + sign * offset
            scc_child_count[p] += 1
            return center_x
        return col_x


    def _place_singleton_in_scc(self, G, graph_id, center_x, y, 
                                in_deg, out_deg, placed_x):
        node = self.nodes[graph_id]
        preds = list(G.predecessors(graph_id))
        if (preds and len(preds) == 1 and
            in_deg[graph_id] == 1 and out_deg[preds[0]] == 1 and
            preds[0] in placed_x):
            x = placed_x[preds[0]]
        else:
            x = center_x
        node.coords(x, y)
        placed_x[graph_id] = x


# Post‑placement helpers
    def _fix_vertical_single_parent_children(self, G, CG, node2scc):
        for child, parent in self._parent_of.items():
            if node2scc[child] == node2scc[parent]:
                continue
            parent_scc = node2scc[parent]
            if CG.in_degree(parent_scc) == 0:
                continue
            if G.in_degree(child) != 1:
                continue

            child_node = self.nodes[child]
            parent_node = self.nodes[parent]
            if parent_node.y is None:
                continue
            if child_node.y <= parent_node.y:
                child_node.coords(child_node.x, parent_node.y + LAYER_GAP)


    def _place_remaining_enzymes(self, placed_x, layers):
        rest = [n for graph_id, n in self.nodes.items() if graph_id not in placed_x]
        for j, node in enumerate(rest):
            node.coords(
                BOARD_MARGIN + j * COL_GAP,
                BOARD_MARGIN + (len(layers) + 1) * LAYER_GAP
            )


    def _shift_layout_into_view(self):
        xs = [node.x for node in self.nodes.values() if node.x is not None]
        if not xs:
            return
        min_x = min(xs)
        if min_x < BOARD_MARGIN:
            shift = BOARD_MARGIN - min_x
            for node in self.nodes.values():
                if node.x is not None:
                    node.coords(node.x + shift, node.y)


    def _layout_scc_as_tree(self, scc_nodes, G, center_x, base_y):
        """
        Layout SCC nodes as a vertical chain, expanding downwards.
        """
        # Subgraph restricted to SCC
        SG = G.subgraph(scc_nodes)

        def score(n):
            # Prefer nodes with FEWER incoming edges from outside the SCC
            return sum(1 for u in G.predecessors(n) if u not in scc_nodes)
        
        root = max(scc_nodes, key=score)

        # BFS tree
        tree = nx.bfs_tree(SG, root)

        levels = defaultdict(list)
        for n in tree.nodes:
            d = nx.shortest_path_length(tree, root, n)
            levels[d].append(n)
       
        # assign depths relative to root
        for depth, nodes in levels.items():
            for n in nodes:
                # 0.2 is arbitrary "within-SCC" depth increment
                self._node_depth[n] = self._node_depth.get(root, 0) + 0.2 * depth
                # WARNING:
                # This value is used only for layout decisions (arrow routing)
                # It must NOT be interpreted as biological or topological depth

            span = (len(nodes) - 1) * COL_GAP
            start_x = center_x - span / 2.0
            y = base_y + depth * LAYER_GAP

            for i, n in enumerate(nodes):
                x = start_x + i * COL_GAP
                self.nodes[n].coords(x, y)


    def layout_anchors(self):
        """Compute conversion anchor coordinates. To obtain anchor ID and xy"""
        for inter in self.interactions:
            if isinstance(inter, ConversionInteraction):
                xy = self.compute_anchor_xy(inter) # calculate anchor coords
                inter._anchor_xy = xy              # update atribute
        

    def _layout_conversions(self):
        """Decide attachment type and GPML points for each conversion."""
            # NOTE: Attachment points 
            # left edge: (-1, 0) 
            # right: (1, 0)
            # top: (0, -1) 
            # bottom: (0, 1)
        depth = getattr(self, "_node_depth", {})

        for inter in self.interactions:
            if not isinstance(inter, ConversionInteraction):
                continue

            src = inter.source
            tgt = inter.target

            # Decide vertical vs side based on y (source above target => vertical)
            d_src = depth.get(src.graph_id, 0.0)
            d_tgt = depth.get(tgt.graph_id, 0.0)

            if d_src < d_tgt:
                inter.attachment = "vertical"
            else:
                inter.attachment = "side"

            if inter.attachment == "vertical":
                # bottom of source -> top of target
                sx = src.x
                sy = src.y + src.height / 2.0
                s_relx, s_rely = "0.0", "1.0"

                tx = tgt.x
                ty = tgt.y - tgt.height / 2.0
                t_relx, t_rely = "0.0", "-1.0"

            else:
                # side-to-side for backward/return edges
                if src.x <= tgt.x:
                    # source left of target
                    sx = src.x + src.width / 2.0
                    s_relx, s_rely = "1.0", "0.0"

                    tx = tgt.x - tgt.width / 2.0
                    t_relx, t_rely = "-1.0", "0.0"
                else:
                    # source right of target
                    sx = src.x - src.width / 2.0
                    s_relx, s_rely = "-1.0", "0.0"

                    tx = tgt.x + tgt.width / 2.0
                    t_relx, t_rely = "1.0", "0.0"

                sy = src.y
                ty = tgt.y

            inter._src_point = {
                "X": str(sx),
                "Y": str(sy),
                "GraphRef": src.graph_id,
                "RelX": s_relx,
                "RelY": s_rely,
            }
            inter._tgt_point = {
                "X": str(tx),
                "Y": str(ty),
                "GraphRef": tgt.graph_id,
                "RelX": t_relx,
                "RelY": t_rely,
                "ArrowHead": "mim-conversion",
            }


    def layout_catalysis(self):
        """Place enzyme nodes next to their anchors once catalysis 
        interactions exist."""
        # Build mapping: anchor_id -> (ax, ay)
        anchor_xy_dict = {     
            inter.anchor_id: inter._anchor_xy
            for inter in self.interactions
            if isinstance(inter, ConversionInteraction) 
                and inter.anchor_id is not None
                and inter._anchor_xy is not None
        }

        # Group enzymes by anchor_id
        enzymes_by_anchor = defaultdict(list)
        for inter in self.interactions:
            if inter.type == "mim-catalysis":
            # Catalysis: source = enzyme/GeneProduct, target = anchor_id
                enz = inter.source       # GeneProduct node
                anchor_id = inter.target # anchor GraphId
                if enz is not None and anchor_id in anchor_xy_dict:
                    enzymes_by_anchor[anchor_id].append(enz)

        # Place enzymes next to their anchors, stacked vertically
        for anchor_id, enzymes in enzymes_by_anchor.items():
            ax, ay = anchor_xy_dict[anchor_id]

            # Put enzymes above the anchor, stacked
            base_x = ax + ENZYME_OFFSET_X
            base_y = ay - ENZYME_OFFSET_Y

            # Center the stack around the anchor in Y
            n = len(enzymes)
            if n == 1:
                offsets = [0.0]
            else:
                total_height = (n - 1) * ENZYME_STACK_GAP
                offsets = [(-total_height / 2.0) 
                           + i * ENZYME_STACK_GAP for i in range(n)]

            for enz, dy in zip(enzymes, offsets):
                enz.coords(base_x, base_y + dy)

        # Update catalysis interactions with the anchor_xy (for GPML Points)
        for inter in self.interactions:
            if inter.type == "mim-catalysis":
                anchor_id = inter.target
                if anchor_id in anchor_xy_dict:
                    inter.anchor_xy = anchor_xy_dict[anchor_id]
