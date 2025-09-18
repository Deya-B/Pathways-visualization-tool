import csv
import xml.etree.ElementTree as ET  # create and manipulate XML structures
import xml.dom.minidom as minidom
import networkx as nx
import datetime
import uuid
from collections import defaultdict


# ----- layout / board -----
BOARD_MARGIN = 80.0     # margin around everything
LAYER_GAP = 140.0       # vertical separation between BFS layers
COL_GAP = 160.0         # approximate horizontal separation
ENZYME_OFFSET = 80.0    # perpendicular distance from the anchor to the enzyme
ENZYME_STACK_GAP = 60.0 # separation between multiple enzymes on the same anchor

############################### ID GEN ########################################
class IDGenerator: 
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
        return f"{p}{self.counters[p]:04d}"     

idgenerator = IDGenerator()

################################ NODE #########################################
# Metabolite(Node) and enzyme(Node)

class Node:
    """
    Contains the Node properties.
    
    Attributes:
        Node Type: Type of element (Metabolite, Enzyme...)
        Label: Name of the element represented in that node.
        Database: The database (DB) where is found (if applicable)
        Database ID (db_id): ID of that element in the given DB.
        Graph ID: Unique ID of the node in the pathway.
        Coordinates: x and y coordinates of the node.
        Width: the total length of the box with the node label.
        Height: the total height of the box with the node label.
    """ 

    def __init__(self, node_type, label, database, db_id):
        # Node properties:
        self.node_type = node_type      
        self.label = label
        self.database = database
        self.db_id = db_id
        self.graph_id = idgenerator.new(self.node_type)
        # self.graph_id = "id" + uuid.uuid5(
            # uuid.NAMESPACE_DNS, f"{node_type}|{label}|{database}|{db_id}").hex[:8]
        self.x, self.y = None, None
        self.width = 90.0 + len(label) * 2
        self.height = 25.0

    def coords(self, x, y):
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
    def __init__(self, source, target, interaction_type, anchor_pos=0.5, anchor_id=None):
        self.source = source                # Node
        self.target = target                # Node (for conversion) or None if target=anchor
        self.type = (interaction_type or "").lower()
        self.anchor_pos = float(anchor_pos)
        self.anchor_id = anchor_id          # set for conversion; used as a destination in catalysis
        self._anchor_xy = None              # (x,y) calculated after having source/target coords

    def compute_anchor_xy(self):
        """Coordenadas del ancla sobre la línea source-target."""
        if self.type != "conversion":
            return None
        if (self.source.x is None or self.source.y is None or
            self.target.x is None or self.target.y is None):
            return None
        x1 = self.source.x
        y1 = self.source.y + self.source.height / 2.0                          
        x2 = self.target.x
        y2 = self.target.y - self.target.height / 2.0
        ax = x1 + self.anchor_pos * (x2 - x1)
        ay = y1 + self.anchor_pos * (y2 - y1)
        self._anchor_xy = (ax, ay)
        return self._anchor_xy

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
            # target point
            target_point = {
                "X": str(self.target.x),
                "Y": str(self.target.y- self.target.height/2.0), # top edge of the destination node               
                "GraphRef": self.target.graph_id, "RelX": "0.0", "RelY": "-1.0",
                "ArrowHead": "mim-conversion"
            }
            ET.SubElement(graphics, "Point", target_point)
            
            # anchor with GraphId (to reference from Catalysis)
            if self.anchor_id is None:
                self.anchor_id = idgenerator.new("anchor")
            if self._anchor_xy is None:
                self.compute_anchor_xy()
            ax, ay = self._anchor_xy    # anchor point > linear                                           
            ET.SubElement(graphics, "Anchor", {
                "GraphId": self.anchor_id, "Position": str(self.anchor_pos), "Shape": "None"
            })

        elif self.type.lower() == "catalysis":
            # origin = enzyme (to nearest edge: use RelY=1 for simplicity)
            ET.SubElement(graphics, "Point", {    
                "X": str(self.source.x),
                "Y": str(self.source.y + self.source.height/2.0),          
                "GraphRef": self.source.graph_id, "RelX": "0.0", "RelY": "1.0"
            })
            # destination = ANCHOR of a conversion (GraphRef = anchor_id)
            # X,Y of the anchor are optional if there is a GraphRef, but we give them the same
            ax, ay = self._anchor_xy if self._anchor_xy else (self.source.x, self.source.y) 
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
        self.nodes = {}             # label -> Node
        self._nodes_by_id = {}      # graph_id -> Node
        self.interactions = []
        # conversion mappings -> pending anchor and catalysis
        self._conv_key_to_inter = {}     # (src_label, tgt_label) -> Interaction(conversion) 
        self._conv_key_to_catalysts = {} # (src_label, tgt_label) -> [enzyme_labels]


    def add_node(self, node: Node):
        self.nodes[node.label] = node
        self._nodes_by_id[node.graph_id] = node

    def add_interaction(self, inter: Interaction):
        self.interactions.append(inter)

    def _compute_board_size(self):
        xs = []
        ys = []
        for n in self.nodes.values():                                
            if n.x is None: continue
            xs += [n.x - n.width/2.0, n.x + n.width/2.0]
            ys += [n.y - n.height/2.0, n.y + n.height/2.0]
        if not xs:
            return 1000.0, 1000.0
        w = (max(xs) - min(xs)) + 2*BOARD_MARGIN
        h = (max(ys) - min(ys)) + 2*BOARD_MARGIN
        return max(w, 400.0), max(h, 400.0)

    def xml_beginning(self, board_w, board_h):
        root = ET.Element("Pathway", {
            "xmlns": "http://pathvisio.org/GPML/2013a",
            "Name": self.title,
            "Version": datetime.date.today().isoformat(),
            "Organism": self.organism
        })
        ET.SubElement(root, "Graphics", {"BoardWidth": f"{board_w:.1f}", "BoardHeight": f"{board_h:.1f}"})
        return root

    def _place_enzymes_near_anchors(self):
        """
        For each conversion with an anchor, place its enzymes perpendicular to the edge
        at ENZYME_OFFSET px of the anchor (stacked if there are multiple).
        """          
        for (s_lbl, t_lbl), conv in self._conv_key_to_inter.items():
            if conv._anchor_xy is None: # do not place enzymes if there's no anchor
                continue
            ax, ay = conv._anchor_xy
            x1 = conv.source.x; y1 = conv.source.y + conv.source.height/2.0
            x2 = conv.target.x; y2 = conv.target.y - conv.target.height/2.0
            dx, dy = (x2 - x1), (y2 - y1)
            L = (dx*dx + dy*dy) ** 0.5 or 1.0
            px, py = (-dy/L, dx/L)
            enzymes = self._conv_key_to_catalysts.get((s_lbl, t_lbl), [])
            for k, enz_lbl in enumerate(enzymes):
                n = self.nodes.get(enz_lbl)
                if not n: 
                    continue
                off = ENZYME_OFFSET + k*ENZYME_STACK_GAP
                n.coords(ax + px*off, ay + py*off)


    def assign_layout(self):
        """Layout BFS for metabolites/products; 
        The enzymes are then placed next to the anchor."""
        if not self.nodes:
            return
        # graph only with metabolites (not enzymes) so the skeleton is clean
        G = nx.Graph()
        for n in self.nodes.values():
            if n.node_type.lower() == "enzyme":
                continue
            G.add_node(n.graph_id)
        for e in self.interactions:
            if e.type == "conversion":
                G.add_edge(e.source.graph_id, e.target.graph_id)
        
        # BFS layers                            
        layers = list(nx.bfs_layers(G, [next(iter(G.nodes))])) if G.number_of_nodes() else []
                     # returns lists of nodes by depth from a root
        # place by layers
        placed = set()       # saves the graph_ids already positioned (metabolites)                                                    
        for ly, layer_nodes in enumerate(layers):
            for i, gid in enumerate(layer_nodes):
                node = self._nodes_by_id[gid]
                x = BOARD_MARGIN + i * COL_GAP
                y = BOARD_MARGIN + ly * LAYER_GAP # LAYER_GAP places each layer in a row separated by COL_GAP
                node.coords(x, y)
                placed.add(gid)

        #  remaining nodes (includes enzymes, isolates, etc.)
        rest = [n for gid, n in self._nodes_by_id.items() if gid not in placed] # remaining graph_ids (enzymes) are left aside
        for j, node in enumerate(rest):
            node.coords(BOARD_MARGIN + j * COL_GAP, 
                        BOARD_MARGIN + (len(layers) + 1) * LAYER_GAP)          
        
        # calculate conversion anchors (for catalysis) and place enzymes nearby
        for e in self.interactions:
            if e.type == "conversion":
                e.compute_anchor_xy()

        # synchronize catalysis with the anchor coordinates
        anchor_xy = {
            e.anchor_id: e._anchor_xy
            for e in self.interactions
            if e.type == "conversion" and e.anchor_id
        }
        for e in self.interactions:
            if e.type == "catalysis":
                e._anchor_xy = anchor_xy.get(e.anchor_id, e._anchor_xy)

        # reposition enzymes near the anchor
        self._place_enzymes_near_anchors()


    def to_etree(self):
        # recalculate board size once everything is positioned
        board_w, board_h = self._compute_board_size()
        root = self.xml_beginning(board_w, board_h)

        # nodos
        for node in self.nodes.values():
            if node.x is None or node.y is None:
                raise ValueError(f"Nodo sin coordenadas: {node.label}")
            root.append(node.to_gpml())

        # interacciones
        for inter in self.interactions:
            root.append(inter.to_gpml())

        # xml finishing lines
        ET.SubElement(root, "InfoBox", {"CenterX": "0.0", "CenterY": "0.0"})
        ET.SubElement(root, "Biopax")
        return root



    def save(self, filename):
        xml_str = ET.tostring(self.to_etree(), encoding="utf-8") # atención: el encoding no sale!
        pretty_xml = minidom.parseString(xml_str).toprettyxml(indent="  ")
        with open(filename, "w", encoding="utf-8") as f:
            f.write(pretty_xml)



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
            conv = Interaction(s, t, "Conversion", anchor_pos=0.4)
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
                cat = Interaction(self.pathway.nodes[enz_lbl], None, "Catalysis",
                                  anchor_pos=conv.anchor_pos, anchor_id=conv.anchor_id)
                # cat._anchor_xy = conv.compute_anchor_xy()
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
    pw.save("ruta_facil3.gpml")


# cd C:\\Users\\deyan\\Desktop\\BIOINFORMÁTICA\\1TFM


