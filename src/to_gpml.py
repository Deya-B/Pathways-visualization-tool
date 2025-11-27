# TODO: ADD Module docstring (what this file is about)
import csv
import datetime
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

# TODO: add DataNode "Pathway"  
# <DataNode TextLabel="Glucose metabolism" GraphId="f8db5" Type="Pathway">
#     <Graphics CenterX="178.5" CenterY="157.625" Width="151.0" Height="19.25" ZOrder="32768" FontWeight="Bold" FontSize="12" Valign="Middle" Color="14961e" />
#     <Xref Database="WikiPathways" ID="WP534" />
#   </DataNode>

# TODO: add elbows
    # In "mim-conversion" add ConnectorType = "Elbow"
    # Cambiar default_anchor_pos = 0.6

# TODO: add organs > make like groups

# TODO: add pathway info (id del pathway) en csv

# TODO: add checkup that all metabolite names/labels are unique


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

    def __init__(self, node_type, label, database, db_id):
        self.node_type = node_type  
        self.label = label
        self.database = database
        self.db_id = db_id
        self.graph_id = idgenerator.new(self.node_type)
        self.x, self.y = None, None
        self.width = 90.0 + len(self.label) * 4
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
            "FontSize": "12",
            "Valign": "Middle", 
            "Color": "0000ff" if self.node_type == "Metabolite" else "000000" 
        })
        ET.SubElement(datanode, "Xref", {
            "Database": self.database,
            "ID": self.db_id
        })
    

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
    def __init__(self, source, target, interaction_type):
        self.graph_id = idgenerator.new("interaction")
        self.source = source                
        self.target = target               
        self.type = (interaction_type or "").lower()

    def bind_to_anchor(self, anchor_id, xy=None):
        self.target = anchor_id
            

    def to_gpml(self):
    #     interaction = ET.Element("Interaction")
    #     graphics = ET.SubElement(interaction, "Graphics", {"LineThickness": "1.0"})

    #     if self.type.lower() == "conversion": 
    #         # source point       
    #         ET.SubElement(graphics, "Point", {    
    #             "X": str(self.source.x),
    #             "Y": str(self.source.y + self.source.height/2.0), # lower edge of the origin node              
    #             "GraphRef": self.source.graph_id, "RelX": "0.0", "RelY": "1.0"
    #         })
    #                                                     # NOTE: Attachment point: 
    #                                                         # left edge: (-1, 0) 
    #                                                         # right: (1, 0)
    #                                                         # top: (0, -1) 
    #                                                         # bottom: (0, 1)
    #         # target point
    #         ET.SubElement(graphics, "Point", {
    #             "X": str(self.target.x),
    #             "Y": str(self.target.y- self.target.height/2.0), # top edge of the destination node               
    #             "GraphRef": self.target.graph_id, "RelX": "0.0", "RelY": "-1.0",
    #             "ArrowHead": "mim-conversion"
    #         })
            
    #         # anchor with GraphId (to reference from Catalysis)
    #         if self.anchor_id is None:
    #             self.anchor_id = idgenerator.new("anchor")
    #         if self._anchor_xy is None:
    #             self.compute_anchor_xy()
    #         ax, ay = self._anchor_xy    # anchor point coords                                          
    #         ET.SubElement(graphics, "Anchor", {
    #             "GraphId": self.anchor_id, 
    #             "Position": str(self.anchor_pos), "Shape": "None"
    #         })

    #     elif self.type.lower() == "catalysis": 
    #         # destination=ANCHOR of a conversion
    #         axay = self.anchor_xy
    #         ax, ay = axay if axay else (self.source.x, self.source.y)


    #         def side_relxy(node, tx, ty):
    #             """Compute nearest side of enzyme to select the source for its arrow"""
    #             dx, dy = tx - node.x, ty - node.y
    #             if abs(dx) > abs(dy):
    #                 # left/right
    #                 return (-1.0, 0.0) if dx < 0 else (1.0, 0.0)  # left or right edge
    #             else:
    #                 # top/bottom
    #                 return (0.0, -1.0) if dy < 0 else (0.0, 1.0)  # top or bottom edge
        
    #         rx, ry = side_relxy(self.source, ax, ay) # enzyme side

    #         # origin at chosen side of enzyme
    #         ET.SubElement(graphics, "Point", {    
    #             "X": str(self.source.x), "Y": str(self.source.y),          
    #             "GraphRef": self.source.graph_id, 
    #             "RelX": str(rx), "RelY": str(ry)
    #         })
    #         # end at the anchor (GraphRef=anchor_id)
    #         ET.SubElement(graphics, "Point", {
    #             "X": str(ax), "Y": str(ay),
    #             "GraphRef": self.anchor_id, "RelX": "0.0", "RelY": "0.0",
    #             "ArrowHead": "mim-catalysis"
    #         })

    #     ET.SubElement(interaction, "Xref", {"Database": "", "ID":""})
    #     return interaction
        pass
    
class ConversionInteraction(Interaction):
    def __init__(self, source, target, anchor_pos=None, anchor_id=None):
        super().__init__(source, target, "mim-conversion")
        self.anchor_id = None
        self._anchor_xy = None # private: to be computed ONLY if coordinates exist
        self.anchor_pos = 0.4


###################### Layout/Board CONFIGURATION #############################

BOARD_MARGIN = 100.0    # margin around everything
LAYER_GAP = 100.0       # vertical separation between BFS layers
COL_GAP = 140.0         # approximate horizontal separation
ENZYME_OFFSET = 0.0     # distance from the anchor to the enzyme
ENZYME_STACK_GAP = 60.0 # separation between multiple enzymes on the same anchor

############################# LAYOUT CREATION #################################

class Layout:
    def __init__(self, nodes, interactions):
        self.nodes = nodes                # Dict of nodes, keyed by graph_id
        self.interactions = interactions  # List of Interactions
        self._laid_out = False


    def compute_anchor_xy(self, interactions):
        """Compute anchor (ax, ay) for a given Interaction. Based on the
        source and target nodes positions."""
        source = interactions.source
        target = interactions.target

        x1 = source.x; y1 = source.y + source.height / 2.0
        x2 = target.x; y2 = target.y - target.height / 2.0

        ax = x1 + interactions.anchor_pos * (x2 - x1)
        ay = y1 + interactions.anchor_pos * (y2 - y1)

        if interactions.anchor_id is None:
            interactions.anchor_id = idgenerator.new("anchor")
        return (ax, ay)


    def run(self):
        # A) Skeleton graph (metabolites only)
        G = nx.Graph()
        for node in self.nodes.values():
            if node.node_type == "GeneProduct": # Ignore enzymes for now
                continue
            G.add_node(node.graph_id)
        for inter in self.interactions:
            if isinstance(inter, ConversionInteraction):
                G.add_edge(inter.source.graph_id, inter.target.graph_id)

        # B) BFS layering (list of rows with the topology)           
        has_nodes = G.number_of_nodes() > 0
        if has_nodes:
            start_node = next(iter(G.nodes))
            bfs_layers_iterator = nx.bfs_layers(G, [start_node])
            layers = list(bfs_layers_iterator) 
            # creates lists of nodes by depth from a root (start node)
        else:
            layers = []

        # C) Setting X/Y coords
        # Step 1. Get maximal layer size and width
        k_max = max((len(L) for L in layers), default=1) # find widest row
        span_max = max(0, (k_max - 1) * COL_GAP) # width occupied by the widest row
        # Step 2. Loop through layers/rows
        placed = set() 
        for ly, layer_nodes in enumerate(layers): 
            k = len(layer_nodes) # number of nodes in this layer
            span = max(0, (k - 1) * COL_GAP) # row width for this layer
            start_x = BOARD_MARGIN + (span_max - span) / 2.0 # center this row  under the widest row
            y = BOARD_MARGIN + ly * LAYER_GAP # y-position for this row
        # Step 3. Assign positions to the nodes
            for i, graph_id in enumerate(layer_nodes):
                node = self.nodes[graph_id]
                x = start_x + i * COL_GAP # horizontal position in the row
                node.coords(x, y)
                placed.add(graph_id) # save graph_ids already positioned

        # D) Temporarily place enzymes on extra row
        rest = [n for graph_id, n in self.nodes.items() if graph_id not in placed]
        for j, node in enumerate(rest):
            node.coords(BOARD_MARGIN + j * COL_GAP, BOARD_MARGIN + (len(layers) + 1) * LAYER_GAP)

        # E) Compute anchors on conversions
        for inter in self.interactions:
            if isinstance(inter, ConversionInteraction):
                anchor_xy = self.compute_anchor_xy(inter) # calculate anchor coords
                inter._anchor_xy = self.compute_anchor_xy(inter) # update atribute
                
        # F) Synchronize catalysis interactions with those anchor coordinates
        # Build an anchor dictionary: anchor_id -> anchor_xy for conversions
        anchor_xy = {     
            inter.anchor_id: inter._anchor_xy
            for inter in self.interactions
            if isinstance(inter, ConversionInteraction) 
            and inter.anchor_id is not None
            and inter._anchor_xy is not None
        }

        # G) place enzymes near their anchors
        for inter in self.interactions:
            enz = inter.source          # Node (GeneProduct)
            a_id = getattr(inter, "anchor_id", None)
            if a_id in anchor_xy:
                ax, ay = anchor_xy[a_id]
                # simple placement: enzyme exactly at anchor
                enz.coords(ax, ay)  
                    
        # H) Compute the FINAL board size now that everyone has coords
        # self._compute_board_size()

        # self._laid_out = True


    def _place_enzymes(self): 
   
    #     """
    #     Place the enzymes at:
    #     - the range of their conversion if they only catalyze one,
    #     - the centroid (mean of X,Y) of all their ranges if they catalyze several.
    #     """       
    #     # 1) (src_label, tgt_label) -> (ax, ay) of the conversion
    #     conv_to_xy = {}   
    #     for (s_lbl, t_lbl), conv in self._conv_key_to_inter.items():
    #         axay = conv.anchor_xy
    #         if axay:
    #             conv_to_xy[(s_lbl, t_lbl)] = axay

    #     # 2) map: enzyme -> list (ax, ay) of the reactions it catalyzes
    #     enz_to_points = defaultdict(list)
    #     for key, enz_list in self._conv_key_to_catalysts.items():
    #         if key not in conv_to_xy:
    #             continue
    #         ax, ay = conv_to_xy[key]
    #         for enz_lbl in enz_list:
    #             if enz_lbl in self.nodes: #(defensive) ignore enzymes not created as a node 
    #                 enz_to_points[enz_lbl].append((ax, ay))
        
    #     # 3) placing on anchor or in centroid
    #     for enz_lbl, pts in enz_to_points.items():
    #         if not pts:
    #             continue
    #         if len(pts) == 1:
    #             ex, ey = pts[0]
    #         else:
    #             ex = sum(p[0] for p in pts) / len(pts)
    #             ey = sum(p[1] for p in pts) / len(pts)
    #         self.nodes[enz_lbl].coords(ex, ey)
        pass 


############################## BUILDER CLASSES ################################
################################## Parser #####################################
# Performs PIPELINE steps:
    # read TSV
    # match IDs
    # decide node types
    # filter rows
    # construct interactions
    # read databases

class Parser:
    def __init__(self, id_data_df, relations_df):
        self.id_data_df = id_data_df
        self.relations_df = relations_df
        self.db_index = {}
        self.nodes = {}
        self.interactions = []

        # Get column names
        self.id_col, self.db_col, self.name_col = (
            self.id_data_df.columns[0:3])
        self.source_col, self.target_col, self.catalyser_col = (
            self.relations_df.columns[0:3])


    def _get_create_node (self, node_id):
        """For source and target nodes"""
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
        is_pathway = is_wikipathways(node_id) or is_kegg(node_id)
        node_type = "Pathway" if is_pathway else "Metabolite"

        # Extract info from tabular data and transform into a Node object
        node = Node(node_type, info[self.name_col], info[self.db_col], node_id)
        self.nodes[node.graph_id] = node

        # Record that this metabolite is already mapped
        self.db_index[node_id] = node.graph_id        

        return node
 

    def build_conversions(self, row):
        # Obtain source and target ID
        source_id = getattr(row, self.source_col)
        target_id = getattr(row, self.target_col)
        # Build nodes
        source_node = self._get_create_node(source_id)
        target_node = self._get_create_node(target_id)

        # Conversion interaction (source -> target)
        if source_node and target_node:
            conversion = ConversionInteraction(source_node, target_node)
            self.interactions.append(conversion)
            return conversion
        
        return None
        

    def _create_catalyser_node(self, catal_id):
        """For catalyser nodes"""
        if pd.isnull(catal_id):
            return
        # Find metadata row
        match = self.id_data_df[self.id_data_df[self.id_col] == catal_id]
        if match.empty:
            return None
        
        info = match.iloc[0]
        node = Node("GeneProduct", info[self.name_col], info[self.db_col], catal_id)
        self.nodes[node.graph_id] = node

        return node
    

    def build_catalysis(self, row, conv_interaction):
        # Obtain catalyser ID
        catal_id = getattr(row, self.catalyser_col)
        if pd.isnull(catal_id):
            return None # no enzyme -> no interaction
        
        catal_node = self._create_catalyser_node(catal_id)
        if catal_node is None:
            return None
        
        anchor_id  = conv_interaction.anchor_id
        
        catal_inter = Interaction(
            source=catal_node, 
            target=anchor_id,
            interaction_type="mim-catalysis")

        self.interactions.append(catal_inter)
        return catal_inter


############################### XML build #####################################

# PASAR A AQUÍ todo LO RELATIVO AL XML...
    # Gather all your nodes/interactions from TSVParser
    # Build the XML tree using their to_gpml() methods
    # Optionally assign coordinates/layout here

    # header creation
    # mapping nodes → XML
    # mapping interactions → XML
    # handling attributes
    # board size
    # naming/organism metadata
    # consistency

    

################################ MAIN #########################################

def main(pathway_title, organism, ID_data_file, relations_file, delimiter="\t", ):
    # DataFrames from Data and Relations files
    id_data_df = pd.read_csv(ID_data_file, sep=delimiter, encoding="utf-8")
    relations_df = pd.read_csv(relations_file, sep=delimiter, encoding="utf-8")

    # Parse DF and get node and interaction objects
    builder = Parser(id_data_df, relations_df)
    
    # Build conversions
    conversions_list = []
    for row in relations_df.itertuples():
        conversion = builder.build_conversions(row)
        conversions_list.append(conversion)

    # Assign layout to compute x,y and anchors for catalysis
    layout = Layout(builder.nodes, builder.interactions)
    layout.run()  
    
    # Build catalytic reactions
    catalysis_list = []
    for row, conversion in zip(relations_df.itertuples(), conversions_list):
        catal = builder.build_catalysis(row, conversion)
        if catal is not None:     
            catalysis_list.append(catal)


    # Checkers
    for key, node in builder.nodes.items():
        print(f"{key}: {vars(node)}")  # vars() returns the attributes as a dict

    # for inter in conversions_list:
    #     source_label = inter.source.graph_id if inter.source else "None"
    #     target_label = inter.target.graph_id if inter.target else "None"
    #     print(f"{source_label} → {target_label}: {vars(inter)}")

    # for inter in catalysis_list:
    #     source_label = inter.source.graph_id if inter.source else "None"
    #     target_label = inter.target
    #     print(f"{source_label} → {target_label}: {vars(inter)}")


########################### HELPER FUNCIONS ###################################

def is_wikipathways(id_str):
    """Check if source/target ID match Pathway WikiPathways database, 
    which can have the form WP1234, WP12345_r2, for example"""
    return re.match(r'WP\d{1,5}(_r\d+)?$', str(id_str)) is not None

def is_kegg(id_str):
    """Check if source/target ID match Pathway WikiPathways database, 
    which has the form hsa00120, for example"""
    return re.match(r'^[a-zA-Z]{2,4}\d{5}$', str(id_str)) is not None


############################ ID GENERATOR #####################################
## Attributes:
        # counters
## Methods:
        # _prefix: extract prefix according to node type
        # new: create a new id
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
    ID_data_file = "c:/Users/dborrotoa/Desktop/TFM/src/examples/data/2-Secondary_BA_Synthesis_CA-Based_reactions_updated.tsv"
    relations_file = "c:/Users/dborrotoa/Desktop/TFM/src/examples/data/relationships.tsv"
    #home
    # ID_data_file = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/data/2-Secondary_BA_Synthesis_CA-Based_reactions_updated.tsv"
    # relations_file = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/data/relationships.tsv"
    
    pathway_title = "Secondary BA Synthesis (CA-Based reactions)"
    organism = "Homo sapiens, Mus-musculus"

    main(pathway_title, organism, ID_data_file, relations_file)



## NOTES:

# default delimiter="\t" 
    # can be changed as a 5th parameter. Example: 
        # main(pathway_title, organism, ID_data_file, relations_file, delimiter=";")

# Relations_file:
    # Must have the following columns in order:
    #   Source Node ID, Target Node ID, Catalyser Node ID
    # If the ID from the DataBase is not available (in cases where the  
    #   specific component is not found in a DataBase), then an ID such 
    #   as "metabo1", "enz1"... must be provided.

# Data_file:
    # Must have the following columns in order:
    #   DataBase ID, DataBase Name, Component Name
    # IMPORTANT: Every ID from the relations file must be also present in 
    #   the ID_data_file provided. Otherwise this wont appear in the pathway.