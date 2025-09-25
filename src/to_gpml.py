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
ENZYME_OFFSET = 0.0     # distance from the anchor to the enzyme
ENZYME_STACK_GAP = 60.0 # separation between multiple enzymes on the same anchor

############################### ID GEN ########################################
## Description:
    # Generate unique id
## Attributes:
        # counters
## Methods:
        # _prefix: extract prefix according to node type
        # new: create a new id

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
## Description:
    # GPML DataNode (metabolite or enzyme). Contains the Node properties.
## Attributes:
        # node_type: Semantic type, e.g. "Metabolite" or "Enzyme".
        # label: Display text.
        # database: Xref DB name.
        # db_id: Xref ID.
        # graph_id: Unique GraphId in GPML.
        # x, y: Center coordinates (PathVisio uses centers).
        # width, height: Node box size in pixels.
## Methods:
        # coords: Set center coordinates
        # to_gpml: convert data to gpml

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
        self.width = 90.0 + len(label) * 4
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
## Description:
    # Represents conversions or catalysis from nodes
    # Computes anchors for catalytic reactions
## Attributes:
        # source: Source Node
        # target: Target Node: if target=conversion | None: if target=anchor
        # type: interaction type (conversion, catalysis...)
        # anchor_pos: (default_anchor_pos = 0.4) position for the anchor
        # anchor_id: set for conversion > used as a destination in catalysis
        # anchor_xy: calculate anchor coords (to be computed ONLY if coordinates exist)
## Methods:

class Interaction: 
    default_anchor_pos = 0.4
    def __init__(self, source, target, interaction_type, anchor_pos=None, anchor_id=None):
        self.source = source                
        self.target = target               
        self.type = (interaction_type or "").lower()
        self.anchor_pos = float(anchor_pos) if anchor_pos is not None else self.default_anchor_pos
        self.anchor_id = anchor_id          
        self._anchor_xy = None              # private: to be computed ONLY if coordinates exist

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
            axay = self.anchor_xy
            ax, ay = axay if axay else (self.source.x, self.source.y)


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
## Description:
    # Define the elements of "Pathway"
    # Computes canvas size
    # Adds nodes and interactions
    # Assigns layout and creates adequate coordinates according to pathway
    # Places enzymes
    # Serializes to GPML

## Attribute (see below)
## Methods: (see below)

class Pathway:  
    def __init__(self, title, organism=None):
        self.title = title
        self.organism = organism
        self.boardwidth = None
        self.boardheight = None
        self._laid_out = False 
        self.nodes = {}             # Node label
        self._nodes_by_id = {}      # Node graph_id
        self.interactions = []
        # conversion mappings -> pending anchor and catalysis
        self._conv_key_to_inter = {}     # (src_label, tgt_label) -> Interaction(conversion) 
        self._conv_key_to_catalysts = {} # (src_label, tgt_label) -> [enzyme_labels]

    def ensure_layout(self):
        """Run layout once if not already done."""
        if not self._laid_out:
            self.assign_layout()

    def add_node(self, node: Node):
        self.nodes[node.label] = node
        self._nodes_by_id[node.graph_id] = node
        # These indices avoid O(n) scans and keep parsing, layout, and GPML writing simple and reliable

    def add_interaction(self, inter: Interaction):
        self.interactions.append(inter)
        
    def _compute_board_size(self):
        xs = []
        ys = []
        # Using the node centre > compute node rectagle
        for n in self.nodes.values():                                
            if n.x is None: continue #skip nodes with no coords
            # collect all left/right edges into xs and all top/bottom edges into ys
            xs += [n.x - n.width/2.0, n.x + n.width/2.0]    # left/right edges
            ys += [n.y - n.height/2.0, n.y + n.height/2.0]  # top/bottom edges

        if not xs: # If nothing has coordinates yet, return a fallback board size
            return 500.0, 500.0  

        w = (max(xs) - min(xs)) + BOARD_MARGIN    
        h = (max(ys) - min(ys)) + 2*BOARD_MARGIN    
        self.boardwidth = max(w, 300.0)
        self.boardheight = max(h, 300.0)

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

    def _place_enzymes(self):
        """
        Place the enzymes at:
        - the range of their conversion if they only catalyze one,
        - the centroid (mean of X,Y) of all their ranges if they catalyze several.
        """       
        # 1) (src_label, tgt_label) -> (ax, ay) of the conversion
        conv_to_xy = {}   
        for (s_lbl, t_lbl), conv in self._conv_key_to_inter.items():
            axay = conv.anchor_xy
            if axay:
                conv_to_xy[(s_lbl, t_lbl)] = axay

        # 2) map: enzyme -> list (ax, ay) of the reactions it catalyzes
        enz_to_points = defaultdict(list)
        for key, enz_list in self._conv_key_to_catalysts.items():
            if key not in conv_to_xy:
                continue
            ax, ay = conv_to_xy[key]
            for enz_lbl in enz_list:
                if enz_lbl in self.nodes: #(defensive) ignore enzymes not created as a node 
                    enz_to_points[enz_lbl].append((ax, ay))
        
        # 3) placing on anchor or in centroid
        for enz_lbl, pts in enz_to_points.items():
            if not pts:
                continue
            if len(pts) == 1:
                ex, ey = pts[0]
            else:
                ex = sum(p[0] for p in pts) / len(pts)
                ey = sum(p[1] for p in pts) / len(pts)
            self.nodes[enz_lbl].coords(ex, ey)

    def assign_layout(self):
        """Layout BFS for metabolites/products; 
        then anchor & place enzymes; 
        finally compute board size."""
        
        if not self.nodes:
            return
            
        # A) Skeleton graph (metabolites only)
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
                                        # creates lists of nodes by depth from a root

        # C) Place metabolites row-by-row (grid)
        layers = list(nx.bfs_layers(G, [next(iter(G.nodes))])) if G.number_of_nodes() else []
        k_max = max((len(L) for L in layers), default=1)    # find widest row
        span_max = max(0, (k_max - 1) * COL_GAP)            # width occupied by the widest row(layer)

        placed = set()                                # save graph_ids already positioned (metabolites) 
        for ly, layer_nodes in enumerate(layers):
            k = len(layer_nodes)                # number of nodes in that row/layer
            span = max(0, (k-1) * COL_GAP)      # row width
            # 1) Calculate the initial offset to center this row within the width of the widest row
            start_x =  BOARD_MARGIN + (span_max - span) / 2.0
            y = BOARD_MARGIN + ly * LAYER_GAP
            # 2) Place the nodes in the row:
            for i, gid in enumerate(layer_nodes):
                node = self._nodes_by_id[gid]
                x = start_x + i * COL_GAP
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
                _ = e.anchor_xy # triggers compute if possible

        # F) Synchronize catalysis interactions with those anchor coordinates
        anchor_xy = {
            e.anchor_id: e.anchor_xy
            for e in self.interactions
            if e.type == "conversion" and e.anchor_id and e.anchor_xy
        }
        for e in self.interactions:
            if e.type == "catalysis"and e.anchor_id in anchor_xy:
                e.bind_to_anchor(e.anchor_id, anchor_xy[e.anchor_id])

        # G) Move enzymes next to their reaction anchors
        self._place_enzymes()
     
        # H) Compute the FINAL board size now that everyone has coords
        self._compute_board_size()

        self._laid_out = True

    def to_etree(self):
        self.ensure_layout()
        root = self.xml_beginning(self.boardwidth, self.boardheight)

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
## Description:
    # creates Pathway, empty structures
    # Reads csv 
    #    - creates Node objects (label → node, graph_id → node)
    #    - records conversions [(src_label, tgt_label)]
    #    - records pending catalysts { (src, tgt): [enz1, enz2, ...] }
    # Creates Interaction objects
    #    - adds conversion interactions; each gets (source, target)
    #    - ensures a shared anchor_id for each conversion
    #    - adds catalysis interactions enzyme → anchor_id (no coords yet)
    # Returns the Pathway instance
## Attributes:
        # csv_file name
        # pathway: pathway title
        # conversions: list of conversions tuples [(source_label, target_label)]
        # pending_catalysis: dict of catalytic reactions (src_label, tgt_label) -> [enzyme_label] 
## Methods:

class CSVPathwayParser:
    def __init__(self, csv_file, title="New Pathway", organism="Homo sapiens", delimiter=None):
        self.csv_file = csv_file
        self.pathway = Pathway(title, organism)
        self.conversions = []                       # [(src_lbl, tgt_lbl)]
        self.pending_catalysis = defaultdict(list)  # (src_lbl, tgt_lbl) -> [enzyme_lbl]
        self.delimiter = delimiter
    
    def read(self):
        with open(self.csv_file, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=self.delimiter)
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
        for s_lbl, t_lbl in self.conversions:   # s = source label, t = target label
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
    csv_file = "examples/data/ruta_media.csv"
    pw = CSVPathwayParser(csv_file, "Bile_Acids", "mus-musculus", ";").read().build_interactions().result()
    pw.assign_layout()
    pw.save("examples/gpml/ruta_media2.gpml")
