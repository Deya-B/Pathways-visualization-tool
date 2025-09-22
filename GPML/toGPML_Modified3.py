import csv
import xml.etree.ElementTree as ET  # create and manipulate XML structures
import xml.dom.minidom as minidom
import networkx as nx
import datetime
import uuid
from collections import defaultdict


# TODO: Add PEP 257 style (https://peps.python.org/pep-0257/) Docstrings 
    #     + style guidelines (https://google.github.io/styleguide/pyguide.html)
    # * For classes: what it represents + key attrs
    # * For methods: what it does, args, returns, side-effects


# ----- layout / board -----
BOARD_MARGIN = 100.0    # margin around everything
LAYER_GAP = 100.0       # vertical separation between BFS layers
COL_GAP = 140.0         # approximate horizontal separation
ENZYME_OFFSET = 200.0     # distance from the anchor to the enzyme
ENZYME_STACK_GAP = 60.0 # separation between multiple enzymes on the same anchor

############################### ID GEN ########################################
class IDGenerator: 
    """Generate unique id"""
    def __init__(self):
        self.counters = defaultdict(int)        
    
    def _prefix(self, kind: str) -> str:         
        k = (kind or "").lower()
        if k in {"enzyme", "geneproduct", "gene", "protein", "e"}: return "e"
        if k in {"metabolite", "compound", "smallmolecule", "m"}:  return "m"
        if k in {"anchor", "a"}:                                   return "a"
        if k in {"interaction", "edge", "i"}:                      return "i"
        return "n"

    def new(self, kind: str) -> str:
        p = self._prefix(kind)
        self.counters[p] += 1
        return f"{p}{self.counters[p]:04d}"  # p=prefix + 4 digits with leading 0's 

idgenerator = IDGenerator()

################################ NODE #########################################
# Metabolite(Node) and enzyme(Node)

class Node:
    """GPML DataNode (metabolite or enzyme). Contains the Node properties.
    
    Attributes:
        node_type: Semantic type, e.g. "Metabolite" or "Enzyme".
        label: Display text.
        database: Xref DB name.
        db_id: Xref ID.
        graph_id: Unique GraphId in GPML.
        x, y: Center coordinates (PathVisio uses centers).
        width, height: Node box size in pixels.
    """ 

    def __init__(self, node_type, label, database, db_id):
        self.node_type = node_type      
        self.label = label
        self.database = database
        self.db_id = db_id
        self.graph_id = idgenerator.new(self.node_type)
        self.x, self.y = None, None
        self.width = 90.0 + len(label) * 2
        self.height = 25.0
    
    def coords(self, x: float, y: float) -> None: # makes sure we get floats
        """Set center coordinates in pixels.

        Args:
            x: Center X in pixels.
            y: Center Y in pixels.
        """
        self.x, self.y = float(x), float(y)

    def to_gpml(self):
        node_type = "GeneProduct" if self.node_type.lower() == "enzyme" else "Metabolite"
        datanode = ET.Element("DataNode", {
            "TextLabel": self.label,
            "GraphId": self.graph_id,
            "Type": node_type
        })
        ET.SubElement(datanode, "Graphics", {
            "CenterX": str(self.x),
            "CenterY": str(self.y),
            "Width": str(self.width),
            "Height": str(self.height),
            "FontSize": "12",
            "Valign": "Middle", 
            "Color": "0000ff" if node_type == "Metabolite" else "000000" 
        })
        ET.SubElement(datanode, "Xref", {
            "Database": self.database,
            "ID": self.db_id
        })
        return datanode

############################## INTERACTION ####################################
# Represents conversions or cathalysis from nodes

class Interaction: 
    default_anchor_pos = 0.4

    def __init__(self, source, target, interaction_type, anchor_pos=None, anchor_id=None):
                        # TODO: anchor_pos default 0.5 and then in parser changes to 0.4 
                            # The parser is overriding the default; pick a value to standardize on
        self.source = source                # Source Node
        self.target = target                # Target Node: if target=conversion | None: if target=anchor
        self.type = (interaction_type or "").lower()
        self.anchor_pos = float(anchor_pos) if anchor_pos is not None else self.default_anchor_pos
        self.anchor_id = anchor_id          # set for conversion > used as a destination in catalysis
        self._anchor_xy = None              # PRIVATE: to be computed ONLY when coordinates exist

    def _can_compute_anchor(self) -> bool:
        """Error handling"""
        return (
            self.type == "conversion"
            and self.source.x is not None and self.source.y is not None
            and self.target.x is not None and self.target.y is not None
        )

    def compute_anchor_xy(self):
        """Compute anchor coordinates on the source-target line."""
        if not self._can_compute_anchor():
            return None
        x1 = self.source.x; y1 = self.source.y + self.source.height/2.0                          
        x2 = self.target.x; y2 = self.target.y - self.target.height/2.0
        ax = x1 + self.anchor_pos * (x2 - x1)
        ay = y1 + self.anchor_pos * (y2 - y1)
        self._anchor_xy = (ax, ay)
        return self._anchor_xy
    
    @property
    def anchor_xy(self):
        """Safe, lazy read: computes after layout; otherwise returns None."""
        if self._anchor_xy is None and self._can_compute_anchor():
            self.compute_anchor_xy()
        return self._anchor_xy

    def bind_to_anchor(self, anchor_id, xy=None):
        """Public API to attach an interaction to an anchor."""
        self.anchor_id = anchor_id
        if xy is not None:
            self._anchor_xy = xy

    def to_gpml(self):
        interaction = ET.Element("Interaction")
        graphics = ET.SubElement(interaction, "Graphics", {"LineThickness": "1.0"})

        if self.type.lower() == "conversion": 
            # source point       
            ET.SubElement(graphics, "Point", {    
                "X": str(self.source.x),
                "Y": str(self.source.y + self.source.height/2.0), # lower edge of the origin node              
                "GraphRef": self.source.graph_id, "RelX": "0.0", "RelY": "1.0"
            })
                                                        # NOTE: Attachment point: 
                                                            # left edge: (-1, 0) 
                                                            # right: (1, 0)
                                                            # top: (0, -1) 
                                                            # bottom: (0, 1)
            # target point
            ET.SubElement(graphics, "Point", {
                "X": str(self.target.x),
                "Y": str(self.target.y- self.target.height/2.0), # top edge of the destination node               
                "GraphRef": self.target.graph_id, "RelX": "0.0", "RelY": "-1.0",
                "ArrowHead": "mim-conversion"
            })
            
            # anchor with GraphId (to reference from Catalysis)
            if self.anchor_id is None:
                self.anchor_id = idgenerator.new("anchor")
            if self._anchor_xy is None:
                self.compute_anchor_xy()
            ax, ay = self._anchor_xy    # anchor point coords                                          
            ET.SubElement(graphics, "Anchor", {
                "GraphId": self.anchor_id, 
                "Position": str(self.anchor_pos), "Shape": "None"
            })

        elif self.type.lower() == "catalysis": 
            # destination=ANCHOR of a conversion
            ax, ay = self._anchor_xy if self._anchor_xy else (self.source.x, self.source.y) 

            def side_relxy(node, tx, ty):
                """Compute nearest side of enzyme to select the source for its arrow"""
                dx, dy = tx - node.x, ty - node.y
                if abs(dx) > abs(dy):
                    # left/right
                    return (-1.0, 0.0) if dx < 0 else (1.0, 0.0)  # left or right edge
                else:
                    # top/bottom
                    return (0.0, -1.0) if dy < 0 else (0.0, 1.0)  # top or bottom edge
        
            rx, ry = side_relxy(self.source, ax, ay) # enzyme side

            # origin at chosen side of enzyme
            ET.SubElement(graphics, "Point", {    
                "X": str(self.source.x), "Y": str(self.source.y),          
                "GraphRef": self.source.graph_id, 
                "RelX": str(rx), "RelY": str(ry)
            })
            # end at the anchor (GraphRef=anchor_id)
            ET.SubElement(graphics, "Point", {
                "X": str(ax), "Y": str(ay),
                "GraphRef": self.anchor_id, "RelX": "0.0", "RelY": "0.0",
                "ArrowHead": "mim-catalysis"
            })

        ET.SubElement(interaction, "Xref", {"Database": "", "ID":""})
        return interaction

############################## PATHWAY ########################################
# Define the elements of "Pathway"
# Adds nodes and interactions
# Places enzymes close to anchors
# Creates and saves the pathway itself
# Computes canvas size

class Pathway:  
    def __init__(self, title, organism="Homo sapiens"):
        self.title = title
        self.organism = organism
        self.nodes = {}             # Node label
        self._nodes_by_id = {}      # Node graph_id
        self.interactions = []
        # conversion mappings -> pending anchor and catalysis
        self._conv_key_to_inter = {}     # (src_label, tgt_label) -> Interaction(conversion) 
        self._conv_key_to_catalysts = {} # (src_label, tgt_label) -> [enzyme_labels]


    def add_node(self, node: Node):
        self.nodes[node.label] = node
        self._nodes_by_id[node.graph_id] = node
        # These indices avoid O(n) scans and keep parsing, layout, and GPML writing simple and reliable

    def add_interaction(self, inter: Interaction):
        self.interactions.append(inter)

    def _compute_board_size(self):
        xs = []
        ys = []
        for n in self.nodes.values():                                
            if n.x is None: continue
            xs += [n.x - n.width/2.0, n.x + n.width/2.0]    # left & right box edges
            ys += [n.y - n.height/2.0, n.y + n.height/2.0]  # top & bottom box edges
        if not xs:
            return 1000.0, 1000.0
        w = (max(xs) - min(xs)) + 2*BOARD_MARGIN
        h = (max(ys) - min(ys)) + 2*BOARD_MARGIN
        return max(w, 300.0), max(h, 300.0)

    def xml_beginning(self, board_w, board_h):
        root = ET.Element("Pathway", {
            "xmlns": "http://pathvisio.org/GPML/2013a",
            "Name": self.title,
            "Version": datetime.date.today().isoformat(),
            "Organism": self.organism
        })
        ET.SubElement(root, "Graphics", {
                    "BoardWidth": f"{board_w:.1f}", 
                    "BoardHeight": f"{board_h:.1f}"})
        return root

    def _place_enzymes_near_anchors(self):
        """
        Place enzymes next to each conversion's anchor on a horizontal offset left/right.
        """          
        for (s_lbl, t_lbl), conv in self._conv_key_to_inter.items():
            axay = conv.anchor_xy
            if not axay: # do not place enzymes if there's no anchor
                continue
            ax, ay = axay

            enzymes = self._conv_key_to_catalysts.get((s_lbl, t_lbl), [])
            if not enzymes:
                continue

            # heuristic: place away from the graph centroid to reduce clutter
            xs = [n.x for n in self.nodes.values()
                if n.x is not None and n.node_type.lower() != "enzyme"]
            cx = (sum(xs) / len(xs)) if xs else ax
            side = 1.0 if ax < cx else - 1.0  # to the right if anchor is left of center
                                              # else to the left
                                            
            # base offset vector
            ux,uy = side, 0.0
            
            # stack enzymes VERTICALLY around the anchor y to avoid crossing the reaction line
            m = len(enzymes)
            for k, enz_lbl in enumerate(enzymes):
                n = self.nodes.get(enz_lbl)
                if not n:
                    continue
                # centered vertical stacking: ..., -1, 0, +1, ...
                vstack = (k - (m - 1) / 2.0) * ENZYME_STACK_GAP
                ex = ax + ux * ENZYME_OFFSET
                ey = ay + vstack
                n.coords(ex, ey)


            # # original perpendicular placement
            # x1 = conv.source.x; y1 = conv.source.y + conv.source.height/2.0
            # x2 = conv.target.x; y2 = conv.target.y - conv.target.height/2.0
            # dx, dy = (x2 - x1), (y2 - y1)
            # L = (dx*dx + dy*dy) ** 0.5 or 1.0 # normalize length 
            # px, py = (-dy/L, dx/L) # take a unit vector PERPENDICULAR to the reaction
            # enzymes = self._conv_key_to_catalysts.get((s_lbl, t_lbl), [])
            # for k, enz_lbl in enumerate(enzymes):
            #     n = self.nodes.get(enz_lbl)
            #     if not n: 
            #         continue
            #     off = ENZYME_OFFSET + k*ENZYME_STACK_GAP
            #     n.coords(ax + px*off, ay + py*off)


    def assign_layout(self):
        """Layout BFS for metabolites/products; 
        The enzymes are then placed next to the anchor."""

        board_w, board_h = self._compute_board_size()

        if not self.nodes:
            return
        # A) Build “skeleton” graph: only metabolites+conversion edges
        G = nx.Graph()
        for n in self.nodes.values():
            if n.node_type.lower() == "enzyme":
                continue
            G.add_node(n.graph_id)
        for e in self.interactions:
            if e.type == "conversion":
                G.add_edge(e.source.graph_id, e.target.graph_id)
        
        # B) BFS layering (topology > rows)                            
        layers = list(nx.bfs_layers(G, [next(iter(G.nodes))])) if G.number_of_nodes() else []
                     # returns lists of nodes by depth from a root

        # C) Place metabolites row-by-row (grid)
        # placed = set() # save graph_ids already positioned (metabolites)                                                    
        # for ly, layer_nodes in enumerate(layers):
        #     for i, gid in enumerate(layer_nodes):
        #         node = self._nodes_by_id[gid]
        #         x = BOARD_MARGIN + i * COL_GAP
        #         y = BOARD_MARGIN + ly * LAYER_GAP # LAYER_GAP places each layer in a row separated by COL_GAP
        #         node.coords(x, y)
        #         placed.add(gid)

        # C) Place metabolites row-by-row (grid)
        placed = set()
        for ly, layer_nodes in enumerate(layers):
            k = len(layer_nodes)
            for i, gid in enumerate(layer_nodes):
                node = self._nodes_by_id[gid]
                x = BOARD_MARGIN + (i + 1) * ((board_w - 2 * BOARD_MARGIN) / (k + 1))
                y = BOARD_MARGIN + ly * LAYER_GAP
                node.coords(x, y)
                placed.add(gid)
        

        # D) Temporarily place the rest (enzymes or others) on an extra row
        rest = [n for gid, n in self._nodes_by_id.items() if gid not in placed] # remaining graph_ids (enzymes) are left aside
        for j, node in enumerate(rest):
            node.coords(BOARD_MARGIN + j * COL_GAP, 
                        BOARD_MARGIN + (len(layers) + 1) * LAYER_GAP)          
        
        # E) Compute anchors on placed reactions
        for e in self.interactions:
            if e.type == "conversion":
                _ = e.anchor_xy

        # F) Synchronize catalysis interactions with those anchor coordinates
        anchor_xy = {
            e.anchor_id: e.anchor_xy
            for e in self.interactions
            if e.type == "conversion" and e.anchor_id and e.anchor_xy
        }
        for e in self.interactions:
            if e.type == "catalysis"and e.anchor_id in anchor_xy:
                e.bind_to_anchor(e.anchor_id, anchor_xy[e.anchor_id])

        # G) Finally, move enzymes next to their reaction anchors
        self._place_enzymes_near_anchors()


    def to_etree(self):
        board_w, board_h = self._compute_board_size()
        root = self.xml_beginning(board_w, board_h)

        # A) Nodes
        for node in self.nodes.values():
            if node.x is None or node.y is None:
                raise ValueError(f"Nodo sin coordenadas: {node.label}")
            root.append(node.to_gpml())

        # B) Interactions
        for inter in self.interactions:
            root.append(inter.to_gpml())

        # C) xml finishing lines
        ET.SubElement(root, "InfoBox", {"CenterX": "0.0", "CenterY": "0.0"})
        ET.SubElement(root, "Biopax")
        return root


    def save(self, filename):
        tree = ET.ElementTree(self.to_etree()) 
        try:
            ET.indent(tree, space="  ") # Python ≥ 3.9
        except AttributeError:
            pass
        tree.write(filename, encoding="utf-8", xml_declaration=True)


############################## PARSING ########################################

class CSVPathwayParser:
    def __init__(self, csv_file, title="New Pathway"):
        self.csv_file = csv_file
        self.title = title
        self.pathway = Pathway(title)
        self.conversions = []                       # [(src_lbl, tgt_lbl)]
        self.pending_catalysis = defaultdict(list)  # (src_lbl, tgt_lbl) -> [enzyme_lbl] 
    
    def read(self):
        with open(self.csv_file, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=",")
            for row in reader:
                lbl = row["Node Label"]              
                if lbl not in self.pathway.nodes:    
                    self.pathway.add_node(Node(row["Node Type"], lbl, 
                                          row.get("Database",""), row.get("Database_ID","")))
                # conversion
                if row.get("Interaction Type") == "Conversion" and row.get("Interaction With"):
                    key = (lbl, row["Interaction With"])
                    self.conversions.append(key)
                    # Collect enzyme if it comes in the "Catalytic" column
                    enz = (row.get("Catalytic") or "").strip()
                    if enz:
                        self.pending_catalysis[key].append(enz)
                    # in case the target hasn't appeared as a separate row yet:
                    if row["Interaction With"] not in self.pathway.nodes:
                        self.pathway.add_node(Node("Metabolite", row["Interaction With"], "", ""))
        return self

    def build_interactions(self):
        # conversions (with anchor)
        for s_lbl, t_lbl in self.conversions:   # s = source, t = target
            s = self.pathway.nodes[s_lbl]
            t = self.pathway.nodes[t_lbl]
            conv = Interaction(s, t, "Conversion")
            self.pathway.add_interaction(conv)
            self.pathway._conv_key_to_inter[(s_lbl, t_lbl)] = conv # maps the source, target pair of each conversion to its object
        
        # catalysis (enzyme to anchor)
        for key, enz_list in self.pending_catalysis.items():
            conv = self.pathway._conv_key_to_inter.get(key)
            if not conv:
                continue
            if conv.anchor_id is None:
                conv.anchor_id = idgenerator.new("anchor")
            for enz_lbl in enz_list:
                if enz_lbl not in self.pathway.nodes:
                    self.pathway.add_node(Node("Enzyme", enz_lbl, "", ""))
                cat = Interaction(self.pathway.nodes[enz_lbl], None, "Catalysis")
                cat.bind_to_anchor(conv.anchor_id)
                self.pathway.add_interaction(cat)
            self.pathway._conv_key_to_catalysts[key] = enz_list
        return self

    def result(self):
        return self.pathway

############################## MAIN ###########################################

if __name__ == "__main__":
    csv_file = "ruta_facil.csv"
    pw = CSVPathwayParser(csv_file, "ruta_facil.csv").read().build_interactions().result()
    pw.assign_layout()
    pw.save("ruta_facil_8.gpml")

# cd C:\\Users\\deyan\\Desktop\\BIOINFORMÁTICA\\1TFM
# cd C:\\Users\\deyan\\GitHub\\Pathways-visualization-tool\\GPML
