# TODO: ADD Module docstring (what this file is about)

import datetime
import logging
import re

import pandas as pd
from collections import defaultdict
import xml.etree.ElementTree as ET  # create and manipulate XML structures
import networkx as nx
import matplotlib.pyplot as plt 


# TODO: Add PEP 257 style (https://peps.python.org/pep-0257/) Docstrings 
    #     + style guidelines (https://google.github.io/styleguide/pyguide.html)
    # * For classes: what it represents + key attrs
    # * For methods: what it does, args, returns, side-effects

# TODO: incorporate a yaml config file
# TODO: incorporate hypothetical/multi-step reaction


############################# CONFIGURATION ###################################
# Layout/Board
BOARD_MARGIN = 200    # margin around everything
CENTERING_FACTOR = 0.8  # between 0.0 (no centering) and 1.0 (full centering)
LAYER_GAP = 120.0       # vertical separation
COL_GAP = 250.0         # approximate horizontal separation

ENZYME_OFFSET_X = 180.0     # horizontal distance from anchor to enzyme
ENZYME_OFFSET_Y = 40.0      # horizontal distance from anchor to enzyme
ENZYME_STACK_GAP = 22.0   # vertical separation when several enzymes share an anchor

# Global variables with accepted variants
PATHWAY_DB_NAMES = {"wikipathways", "reactome", "kegg pathway"}

################################ NODE #########################################

class Node:

    def __init__(self, node_type, label, database, db_id, colour):
        self.node_type = node_type
        self.colour = colour
        self.font_weight = None
        self.label = label
        self.database = database
        self.db_id = db_id
        self.graph_id = idgenerator.new(self.node_type)
        self.x, self.y = None, None
        self.width = 80.0 + (len(self.label) * 4)
        self.height = 25.0
    
    
    def coords(self, x: float, y: float) -> None: # making sure to get floats
        """Set positions as center coordinates in pixels.

        Args:
            x: Center X in pixels.
            y: Center Y in pixels.
        """
        self.x, self.y = float(x), float(y)
        pass


    def to_gpml(self):
        "Serialize to XML"
        datanode = ET.Element("DataNode", {
            "TextLabel": self.label,
            "GraphId": self.graph_id,
            "Type": self.node_type
        })
        ET.SubElement(datanode, "Graphics", {
            "CenterX": str(self.x),
            "CenterY": str(self.y),
            "Width": str(self.width),
            "Height": str(self.height),
            "FontWeight": "Normal" if self.node_type == "Metabolite" 
                                   or self.node_type == "GeneProduct" 
                                   else "Bold",
            "FontSize": "12",
            "ShapeType": "RoundedRectangle" if self.node_type == "Pathway" 
                                else "Rectangle",
            "Valign": "Middle", 
            "Color": str(self.colour)
        })
        if self.database is not None:
            ET.SubElement(datanode, "Xref", {
                "Database": str(self.database),
                "ID": str(self.db_id)
            })
        
        return datanode


############################## INTERACTION ####################################

class Interaction:
    def __init__(self, source, target, interaction_type):
        self.graph_id = idgenerator.new("interaction")
        self.source = source                
        self.target = target               
        self.type = (interaction_type or "").lower()
        self.anchor_xy = None


    def bind_to_anchor(self, anchor_id, xy=None):
        self.target = anchor_id
        if xy is not None:
            self.anchor_xy = xy
            

    def to_gpml(self):
        """GPML for catalysis interactions."""
        inter_el = ET.Element("Interaction", {
            "GraphId": self.graph_id
        })
        graphics = ET.SubElement(inter_el, "Graphics", {
            "LineThickness": "1.0"
        })

        if self.type == "mim-catalysis":
            enz = self.source       # GeneProduct node
            anchor_id = self.target
            ex, ey = enz.x, enz.y   # Enzyme coordinates
        
            # Origin point: enzyme center
            ET.SubElement(graphics, "Point", {
                "X": str(ex),
                "Y": str(ey),
                "GraphRef": enz.graph_id,
                "RelX": "0.0",
                "RelY": "0.0"
            })
            # Destination point: anchor (use anchor_xy)
            if self.anchor_xy is None:
                # defensive: fall back to enzyme center
                ax, ay = ex, ey
            else:
                ax, ay = self.anchor_xy

            ET.SubElement(graphics, "Point", {
                "X": str(ax),
                "Y": str(ay),
                "GraphRef": anchor_id,
                "RelX": "0.0",
                "RelY": "0.0",
                "ArrowHead": "mim-catalysis"
            })

        ET.SubElement(inter_el, "Xref", {"Database": "", "ID":""})
        return inter_el

    
class ConversionInteraction(Interaction):
    def __init__(self, source, target, anchor_id=None):
        super().__init__(source, target, "mim-conversion")
        self.anchor_id = anchor_id
        self._anchor_xy = None # private: to be computed by Layout class
        self.anchor_pos = 0.8
        self.attachment = "vertical"  # attachment mode: "vertical" or "side"

        # GPML point info (filled by Layout)
        self._src_point = None  # dict with X,Y,RelX,RelY
        self._tgt_point = None  # dict with X,Y,RelX,RelY


    def to_gpml(self):
        """GPML for an anchored conversion."""
        inter_el = ET.Element("Interaction", {
            "GraphId": self.graph_id
        })
        graphics = ET.SubElement(inter_el, "Graphics", {
            "ConnectorType": "Elbow",
            "LineThickness": "1.0"
        })
        
        # Use points computed by Layout        
        ET.SubElement(graphics, "Point", self._src_point)
        ET.SubElement(graphics, "Point", self._tgt_point)

        # Anchor with GraphId: uses anchor_id and anchor_pos calculated in Layout                                        
        ET.SubElement(graphics, "Anchor", {
            "GraphId": self.anchor_id, 
            "Position": str(self.anchor_pos), 
            "Shape": "None"
        })

        ET.SubElement(inter_el, "Xref", {"Database": "", "ID":""})
        return inter_el


############################# LAYOUT CREATION #################################

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
        """Pre-layout positions: Assign x,y to metabolites and pathways 
        nodes only."""
       # Skeleton graph with source/target nodes and conversion interactions only
        G = nx.DiGraph()

        for node in self.nodes.values():
            if node.node_type == "GeneProduct": # Ignore enzymes for now
                continue
            G.add_node(node.graph_id)

        for inter in self.interactions:
            if isinstance(inter, ConversionInteraction):
                G.add_edge(inter.source.graph_id, inter.target.graph_id)

        parent_of = {}
        for u, v in G.edges:
            if v not in parent_of:
                parent_of[v] = u

        self._parent_of = parent_of

        # Precompute degrees for chain detection
        in_deg = {n: G.in_degree(n) for n in G.nodes}
        out_deg = {n: G.out_degree(n) for n in G.nodes}

    # Approach:   
    # Collapse each strongly connected component (SCC) (i.e. cycle) into a 
    #   single super‑node.
    # Run the topological‑depth layout on this condensation DAG.
    # Inside each SCC:
        # If it has just 1 node, place it exactly at its layer position.
        # If it has multiple nodes (cycle), place its members with a simple
        #   local layout (e.g. small horizontal/vertical cluster around that 
        #   layer position).

        # 1) Strongly connected components (SCC)
        sccs = list(nx.strongly_connected_components(G))  # set of node ids per SCC
        # Map node -> scc_id (index)
        node2scc = {}
        for idx, comp in enumerate(sccs):
            for n in comp:
                node2scc[n] = idx

        # 2) Condensation graph (DAG of SCCs)
        CG = nx.condensation(G, sccs)  # nodes are 0..len(sccs)-1, each is a SCC
        # Topological order on SCC DAG
        topo_scc = list(nx.topological_sort(CG))

        # 3) Depth per SCC (longest path)
        depth_scc = {c: 0 for c in topo_scc}
        for u in topo_scc:
            for v in CG.successors(u):
                depth_scc[v] = max(depth_scc[v], depth_scc[u] + 1)

        # 4) Layers of SCCs, then expand to node layers
        max_depth = max(depth_scc.values(), default=0)
        scc_layers = [[] for _ in range(max_depth + 1)]
        for comp_id, d in depth_scc.items():
            scc_layers[d].append(comp_id)

        # Turn into node-based layers: each SCC contributes its members
        layers = []
        node_depth = {}  # graph_id -> numeric depth (can be non-integer)

        for d, comps in enumerate(scc_layers):
            layer_nodes = []
            for comp_id in comps:
                members = list(sccs[comp_id])
                layer_nodes.extend(members)
                # inside the SCC, give each node a slight offset so edges
                # have a direction even within cycles
                if len(members) == 1:
                    node_depth[members[0]] = float(d)
                else:
                    # temporarily same base depth; refined after SCC layout
                    for n in members:
                        node_depth[n] = float(d)
            layers.append(layer_nodes)

        # 5) Horizontal spacing & placement
        k_max = max((len(L) for L in layers), default=1)
        span_max = max(0, (k_max - 1) * COL_GAP)

        placed_x = {}  # node_id -> x
        placed_scc_x = {}   # scc_id  -> x
        scc_child_count = defaultdict(int)  # parent_scc_id -> how many snapped children


        # Local depth of each SCC within its parent’s subgraph (for ordering)
        scc_local_depth = defaultdict(int)

        # Simple BFS from each SCC to assign increasing depth to descendants
        for s in CG.nodes:
            # skip if no outgoing edges
            if CG.out_degree(s) == 0:
                continue
            # BFS from s
            queue = [(s, 0)]
            seen = {s}
            while queue:
                u, d = queue.pop(0)
                for v in CG.successors(u):
                    if v not in seen:
                        seen.add(v)
                        # depth relative to s
                        scc_local_depth[v] = max(scc_local_depth[v], d + 1)
                        queue.append((v, d + 1))


        for ly, layer_nodes in enumerate(layers):
            if not layer_nodes:
                continue

            # Compute row geometry
            k = len(layer_nodes)
            span = max(0, (k - 1) * COL_GAP)
            start_x = BOARD_MARGIN + CENTERING_FACTOR * (span_max - span) / 2.0
            y = 80 + ly * LAYER_GAP

            # Group nodes in this layer by SCC
            by_scc = {}
            for n in layer_nodes:
                by_scc.setdefault(node2scc[n], []).append(n)

            # Helper: average parent x per SCC (using condensation graph)
            def scc_parent_avg_x(scc_id):
                parent_sccs = list(CG.predecessors(scc_id))
                xs = []
                for ps in parent_sccs:
                    if ps in placed_scc_x:
                        xs.append(placed_scc_x[ps])
                    else:
                        # fallback: use individual parent nodes if any are placed
                        for n in sccs[ps]:
                            preds = list(G.predecessors(n))
                            xs.extend(placed_x[p] for p in preds if p in placed_x)
                return sum(xs) / len(xs) if xs else float("inf")

            # Order SCC blocks by their parent barycenter
            scc_ids_in_layer = list(by_scc.keys())
            scc_ids_in_layer.sort(key=scc_parent_avg_x)

            i = 0
            for comp_id in scc_ids_in_layer:
                members = by_scc[comp_id]
                size = len(members)

                # Default column-based x for this SCC
                col_x = start_x + i * COL_GAP

                # Strong-chain rule at SCC level:
                # if SCC has exactly one incoming SCC and that parent SCC
                # has only this SCC as outgoing, snap under parent SCC
                parent_sccs = list(CG.predecessors(comp_id))

                if len(parent_sccs) == 1 and parent_sccs[0] in placed_scc_x:
                    p = parent_sccs[0]
                    
                    # Children of p in this layer, sorted by local depth then id
                    siblings = [
                        cid for cid in scc_ids_in_layer
                        if p in CG.predecessors(cid)
                    ]
                    siblings.sort(key=lambda cid: (scc_local_depth[cid], cid))

                    # Index of this comp_id among its siblings
                    idx = siblings.index(comp_id)

                    if idx == 0:
                        # first child: exactly under parent
                        center_x = placed_scc_x[p]
                    else:
                        # next children: alternate left/right, moving outward
                        offset_index = (idx + 1) // 2
                        offset = offset_index * COL_GAP
                        sign = -1 if idx % 2 == 1 else 1
                        center_x = placed_scc_x[p] + sign * offset

                    scc_child_count[p] += 1
                else:
                    center_x = col_x

                if size == 1:
                    graph_id = members[0]
                    node = self.nodes[graph_id]

                    # Optional: node-level strong-chain rule
                    # "strong chain" rule: single parent and parent has only this child
                    preds = list(G.predecessors(graph_id))
                    if (preds and len(preds) == 1 and
                        in_deg[graph_id] == 1 and out_deg[preds[0]] == 1 and
                        preds[0] in placed_x):
                        x = placed_x[preds[0]]
                    else:
                        x = center_x

                    node.coords(x, y)
                    placed_x[graph_id] = x

                else:
                    self._layout_scc_as_tree(
                        members, 
                        G,
                        center_x=center_x,
                        base_y=y,
                        scc_id=comp_id,
                        CG=CG,
                    )
                    for m in members:
                        placed_x[m] = self.nodes[m].x


                placed_scc_x[comp_id] = center_x
                i += 1  # one slot per SCC


        # 5.8) Fix vertical placement of non-SCC single-parent children
        for child, parent in self._parent_of.items():
            # Skip SCC-internal nodes
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

            # Force child below parent
            if child_node.y <= parent_node.y:
                child_node.coords(
                    child_node.x,
                    parent_node.y + LAYER_GAP
                )


        # 6) Temporarily place enzymes on extra row
        rest = [n for graph_id, n in self.nodes.items() if graph_id not in placed_x]
        for j, node in enumerate(rest):
            node.coords(
                BOARD_MARGIN + j * COL_GAP, 
                BOARD_MARGIN + (len(layers) + 1) * LAYER_GAP)

        # Save depths
        self._node_depth = node_depth

        # 7) Compute conversion attachment points
        self._layout_conversions()

        # 8) Shift whole layout into the visible canvas (fix negative/too-small X)
        xs = [node.x for node in self.nodes.values() if node.x is not None]
        if xs:
            min_x = min(xs)
            # Ensure the leftmost node is at least at BOARD_MARGIN
            if min_x < BOARD_MARGIN:
                shift = BOARD_MARGIN - min_x
                for node in self.nodes.values():
                    if node.x is not None:
                        node.coords(node.x + shift, node.y)

        self._laid_out = True


    def _layout_scc_as_tree(self, scc_nodes, G, center_x, base_y, scc_id, CG):
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


############################## BUILDER CLASSES ################################
################################## Parser #####################################


class Parser:
    def __init__(self, id_data_df, relations_df):
        self.id_data_df = id_data_df
        self.relations_df = relations_df
        self.db_index = {}
        self.nodes = {}
        self.interactions = []
        self.conversion_by_pair = {}  # (source_id, target_id) -> ConvInter 
                                      # to create only one ConvInter per pair
        # Get column names
        self.id_col, self.db_col, self.name_col = (
            self.id_data_df.columns[0:3])
        self.source_col, self.target_col, self.catalyser_col = (
            self.relations_df.columns[0:3])


    def _clean_field(self, value, *, empty_to_none=True):
        """Normalize string-like fields from TSV/CSV."""
        if isinstance(value, str):
            value = value.strip()
        if pd.isna(value):
            value = ""
        if empty_to_none and value == "":
            return None
        return value


    def _is_pathway(self, node_id):
        """
        Devuelve True si:
        - El ID tiene formato WikiPathways (WP5176, WP12345_r2, etc.), o
        - El ID tiene formato KEGG Pathway (hsa00120, map00121, etc.), o
        - En la columna db_col correspondiente a ese ID aparece la palabra
            'pathway'.
        """
        if pd.isnull(node_id):
            return False
        
        # WikiPathways: WP + 4-5 dígitos
        if re.fullmatch(r"WP\d{4,5}.*", str(node_id)):
            return True

        # KEGG Pathway típico: 3 letras (organismo) + 5 dígitos
        if re.fullmatch(r"[a-z]{3}\d{5}", str(node_id)):
            return True

        # Buscar el registro en el DataFrame
        match = self.id_data_df[self.id_data_df[self.id_col] == node_id]
        if match.empty:
            return False
        db_value = str(match.iloc[0][self.db_col]).lower()
        if any(name in db_value for name in PATHWAY_DB_NAMES):
            return True
        return False


    def _get_create_node (self, node_id):
        """For source and target nodes"""
        node_id = self._clean_field(node_id)  # normalize ID
        if pd.isnull(node_id):
            return None
        
        # If already created, return it
        if node_id in self.db_index:
            graph_id = self.db_index[node_id]
            return self.nodes[graph_id]
        
        # Find metadata row
        match = self.id_data_df[self.id_data_df[self.id_col] == node_id]
        if match.empty:
            return None
        
        info = match.iloc[0]
        is_pathway_flag = self._is_pathway(node_id)
        node_type = "Pathway" if is_pathway_flag else "Metabolite"
        colour = self.get_node_colour(node_type)
        
        label = self._clean_field(info[self.name_col], empty_to_none=False)
        database = self._clean_field(info[self.db_col], empty_to_none=False)

        # Extract info from tabular data and transform into a Node object
        node = Node(node_type, label, database, node_id, colour)
        self.nodes[node.graph_id] = node

        # Record that this metabolite is already mapped
        self.db_index[node_id] = node.graph_id     
        return node


    def build_conversions(self, row):
        # Obtain source and target ID
        raw_source = getattr(row, self.source_col)
        raw_target = getattr(row, self.target_col)

        source_id = self._clean_field(raw_source)
        target_id = self._clean_field(raw_target)

        # Build nodes
        source_node = self._get_create_node(source_id)
        target_node = self._get_create_node(target_id)

        if not (source_node and target_node):
            return None

        key = (source_id, target_id)

        # Reuse existing ConversionInteraction if present
        if key in self.conversion_by_pair:
            return self.conversion_by_pair[key]
        
        conversion = ConversionInteraction(source_node, target_node)
        self.interactions.append(conversion)
        self.conversion_by_pair[key] = conversion
        return conversion
        

    def _create_catalyser_node(self, catal_id):
        """For catalyser nodes"""
        catal_id = self._clean_field(catal_id)
        if pd.isnull(catal_id):
            return None
        
        # Find metadata row
        match = self.id_data_df[self.id_data_df[self.id_col] == catal_id]
        if match.empty:
            return None
        
        info = match.iloc[0]
        colour = self.get_node_colour("GeneProduct")
            
        label = self._clean_field(info[self.name_col], empty_to_none=False)
        database = self._clean_field(info[self.db_col], empty_to_none=False)

        node = Node("GeneProduct", label, database, catal_id, colour)
        self.nodes[node.graph_id] = node
        return node
    

    def build_catalysis(self, row, conv_interaction):
        # Obtain catalyser ID
        raw_catal = getattr(row, self.catalyser_col)
        catal_id = self._clean_field(raw_catal)
        if pd.isnull(catal_id):
            return None # no enzyme -> no interaction
        
        catal_node = self._create_catalyser_node(catal_id)
        if catal_node is None:
            return None
        
        anchor_id  = conv_interaction.anchor_id # shared anchor
        
        catal_inter = Interaction(
            source=catal_node, 
            target=anchor_id,
            interaction_type="mim-catalysis")

        self.interactions.append(catal_inter)
        return catal_inter


    def get_node_colour(self, node_type):
        """Get XML Graphics standards stablished by WikiPathways for the 
        different node types"""
        if node_type == "Metabolite": 
            return "0000ff"
        elif node_type == "Pathway": 
            return "14961e"
        else: # For Enzymes and others
            return "000000"

            
############################### XML build #####################################

class XMLBuilder:  
    def __init__(self, title, organism=None, nodes=None, interactions=None):
        self.title = title
        self.organism = organism
        self._laid_out = False 
        self.nodes = nodes or {}    # Node label
        self._nodes_by_id = {}      # Node graph_id
        self.interactions = interactions or []


    def to_etree(self):
        root = self.xml_header()

        # A) Nodes
        for node in self.nodes.values(): # iterar por los objetos Node
            # if node.x is None or node.y is None:
            #     raise ValueError(f"Nodo sin coordenadas: {node.label}")
            root.append(node.to_gpml())

        # B) Interactions
        for inter in self.interactions:
            root.append(inter.to_gpml())

        # C) InfoBox & Biopax xml finishing lines
        ET.SubElement(root, "InfoBox", {"CenterX": "0.0", "CenterY": "0.0"})
        ET.SubElement(root, "Biopax")
        return root
    

    # def xml_beginning(self, board_w, board_h):
    def xml_header(self):
        root = ET.Element("Pathway", {
            "xmlns": "http://pathvisio.org/GPML/2013a",
            "Name": self.title,
            "Version": datetime.date.today().isoformat(),
            "Organism": self.organism
        })
        ET.SubElement(root, "Graphics", {
                    # "BoardWidth": f"{board_w:.1f}", 
                    # "BoardHeight": f"{board_h:.1f}"})
                    "BoardWidth": "100", 
                    "BoardHeight": "100"
        })
        return root
    

################################ MAIN #########################################

def main(pathway_title, organism, 
         ID_data_file, relations_file, output_filename, 
         delimiter="\t"):
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    # DataFrames from Data and Relations files
    try:
        id_data_df = read_csv(ID_data_file, sep=delimiter)
        logging.info(f"Successfully read ID data file: {ID_data_file}")
    except Exception as e:
        logging.error(f"Failed to read ID data file {ID_data_file}: {e}")
        raise

    try:
        relations_df = read_csv(relations_file, sep=delimiter)
        logging.info(f"Successfully read relations file: {relations_file}")
    except Exception as e:
        logging.error(f"Failed to read relations file {relations_file}: {e}")
        raise

    # Parse DF and get node and interaction objects
    builder = Parser(id_data_df, relations_df)
    
    # Build source/target nodes + conversion interactions
    conversions_list = []
    for row in relations_df.itertuples():
        conversion = builder.build_conversions(row)
        if conversion is not None: 
            conversions_list.append(conversion)

    # Assign layout to compute x,y coordinates + anchors for catalysis
    layout = Layout(builder.nodes, builder.interactions)
    layout.layout_positions()
    layout.layout_anchors()  

    # Build catalytic nodes + interactions
    catalysis_list = []
    for row, conversion in zip(relations_df.itertuples(), conversions_list):
        if conversion is not None:
            catal = builder.build_catalysis(row, conversion)
            if catal is not None:     
                catalysis_list.append(catal)
    layout.layout_catalysis() # place enzymes

    # Build GPML
    xml_builder = XMLBuilder(
        title=pathway_title,
        organism=organism,
        nodes=builder.nodes,
        interactions=builder.interactions
    )
    root = xml_builder.to_etree()
    tree = ET.ElementTree(xml_builder.to_etree())
    
    try:
        ET.indent(tree, space="  ") 
    except AttributeError:
        pass
    tree.write(output_filename, encoding="utf-8", xml_declaration=True)


########################## Helper Functions ###################################

def read_csv(path, sep="\t", encodings=("utf-8", "utf-16", "cp1252")):
    last_err = None
    for enc in encodings:
        try:
            df = pd.read_csv(path, sep=sep, encoding=enc)
            logging.info(f"Successfully read {path} with encoding {enc}")
            return df
        except (UnicodeError, UnicodeDecodeError) as e:
            logging.warning(f"Failed to read {path} with encoding {enc}: {e}")
            last_err = e
        except Exception as e:
            logging.warning(f"Unexpected error reading {path} with encoding {enc}: {e}")
            last_err = e
    # If none of them work
    logging.error(f"Could not read {path} with any encoding")
    raise last_err


############################ ID GENERATOR #####################################

class IDGenerator: 
    """
    Generate unique IDs with custom prefixes.
    
    Attributes
    ----------
        counters (defaultdict): Maintains a count of IDs by prefix.

    Methods
    -------
        new(kind): Generates a new unique ID for the given kind.
    """
    def __init__(self):
        self.counters = defaultdict(int)        
 
    def _prefix(self, kind: str) -> str:         
        k = (kind or "").lower()
        if k in {"enzyme", "geneproduct", "gene", "protein", "e"}: return "e"
        if k in {"metabolite", "compound", "smallmolecule", "m"}:  return "m"
        if k in {"pathway"}:                                       return "p"
        if k in {"anchor", "a"}:                                   return "a"
        if k in {"interaction", "edge", "i"}:                      return "i"
        return "n"

    def new(self, kind: str) -> str:
        p = self._prefix(kind)
        self.counters[p] += 1
        return f"{p}{self.counters[p]:04d}"  # p=prefix + 4 digits with leading 0's 

# initiation of an instance for this class
idgenerator = IDGenerator()


############################# ENTRY POINT #####################################

if __name__ == "__main__":
    name = "AlternativeBA"
    # name = "7a-dehydroxylation"
    # name = "CA-DCA_UCA"
    # name = "CDCA-LCA_UDCA"
    # name = "MurineCDCA-MCA_MDCA"

    # ID_data_file = "c:/Users/dborrotoa/Desktop/TFM/src/examples/data/.tsv"
    # relations_file = "c:/Users/dborrotoa/Desktop/TFM/src/examples/data/_relationships.tsv"
    # output_filename = "c:/Users/dborrotoa/Desktop/TFM/src/examples/gpml/ruta.gpml"
    #home
    ID_data_file = f"C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/data/{name}.tsv"
    relations_file = f"C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/data/{name}_relationships.tsv"
    output_filename = f"C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/data/{name}_Test.gpml"
    
    pathway_title = f"{name}_Test"
    organism = "Homo sapiens" ## "Homo sapiens" "Mus-musculus"

    main(pathway_title, organism, 
         ID_data_file, relations_file, output_filename)




## NOTES:

# default delimiter="\t" 
    # can be changed as a 5th parameter. Example: 
        # main(pathway_title, organism, ID_data_file, relations_file, delimiter=";")

# Relations_file:
    # Must have the following columns in order:
    #   Source Node ID, Target Node ID, Catalyser Node ID
    # If the ID from the DataBase is not available (in cases where the  
    #   specific component is not found in a DataBase), then an ID such 
    #   as "CNIC-M001" (metabolite), "CNIC-E001" (enzyme)... must be provided
    #   and added to the ID metadata file.

# ID metadata file:
    # Must have the following columns in order:
    #   DataBase ID, DataBase Name, Label (Component Name)
    # IMPORTANT: Every ID from the relations file must be also present in 
    #   the ID_data_file provided. Otherwise this wont appear in the pathway.

# ET.indent(tree, space="  "), part of xml.etree.ElementTree requires 
#   Python ≥ 3.9 to work