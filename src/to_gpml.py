# TODO: ADD Module docstring (what this file is about)
import csv
import datetime
import re

import pandas as pd
from collections import defaultdict
import xml.etree.ElementTree as ET  # create and manipulate XML structures
import networkx as nx


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


# ###################### Layout/Board CONFIGURATION #############################

# BOARD_MARGIN = 100.0    # margin around everything
# LAYER_GAP = 100.0       # vertical separation between BFS layers
# COL_GAP = 140.0         # approximate horizontal separation
# ENZYME_OFFSET = 0.0     # distance from the anchor to the enzyme
# ENZYME_STACK_GAP = 60.0 # separation between multiple enzymes on the same anchor

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
        # self.width = 90.0 + len(self.label) * 4
        # self.height = 25.0
    
    
    def obtain_coords(self, x: float, y: float) -> None: # making sure to get floats
        """Set positions as center coordinates in pixels.

        Args:
            x: Center X in pixels.
            y: Center Y in pixels.
        """
        self.x, self.y = float(x), float(y)
        pass

    def to_gpml(self):
        "Serialize to XML"
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
        self.graph_id = idgenerator.new("interaction")
        self.source = source                
        self.target = target               
        self.type = (interaction_type or "").lower()
        self.anchor_pos = float(anchor_pos) if anchor_pos is not None else self.default_anchor_pos
        self.anchor_id = anchor_id          
        self._anchor_xy = None              # private: to be computed ONLY if coordinates exist

        # self.ArrowHead = if enzyme is "mim-catalysis" else "mim-conversion"

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
    
    def to_gpml(self):
        pass


        
################################ MAIN #########################################

def main(pathway_title, organism, ID_data_file, relations_file, delimiter="\t", ):
    # DF from Relations file 
    relations_df = pd.read_csv(relations_file, sep=delimiter, encoding="utf-8")
    source_col, target_col, catalyser_col = relations_df.columns[0:3] # extract column names

    # DF from Data file
    id_data_df = pd.read_csv(ID_data_file, sep=delimiter, encoding="utf-8")
    id_col,db_col,name_col = id_data_df.columns[0:3] # extract column names
    
    nodes_dict_info = {} # Keep nodes info
    db_index = {} # Keep record of metabolite nodes already mapped (not to repeat)
    interactions_list = []

    # Iterate over the rows as relations in tuples of three
    #  print(relations_df.itertuples)
    for row in relations_df.itertuples():
        # Get IDs
        source_id = getattr(row, source_col)
        target_id = getattr(row, target_col)
        catal_id = getattr(row, catalyser_col)

        # Build SOURCE Nodes
        if pd.notnull(source_id) and source_id not in db_index:
            # Matching ID with that from id_data_df
            match_row = id_data_df[id_data_df[id_col] == source_id]
            if not match_row.empty:
                info = match_row.iloc[0]
                is_pathway = is_wikipathways(source_id) or is_kegg(source_id)
                node_type = "Pathway" if is_pathway else "Metabolite"
                # Extract info from tabular data and transform into a Node object
                node = Node(node_type, info[name_col], info[db_col], source_id)
                nodes_dict_info[node.graph_id] = node
                # record that this metabolite is already mapped
                db_index[source_id] = node.graph_id

        # Build TARGET Nodes
        if pd.notnull(target_id) and target_id not in db_index:
            match_row = id_data_df[id_data_df[id_col] == target_id]
            if not match_row.empty:
                info = match_row.iloc[0]
                is_pathway = is_wikipathways(target_id) or is_kegg(target_id)
                node_type = "Pathway" if is_pathway else "Metabolite"
                node = Node(
                    node_type, info[name_col], info[db_col], target_id
                    )
                nodes_dict_info[node.graph_id] = node
                # record that this metabolite is already mapped
                db_index[target_id] = node.graph_id

        # Interactions
        source_node = nodes_dict_info.get(db_index.get(source_id))
        target_node = nodes_dict_info.get(db_index.get(target_id))
        # print("source:", source_node.graph_id, "target:", target_node.graph_id)

        if source_node and target_node:
            conversion = Interaction(
                source=source_node.graph_id, target=target_node.graph_id, 
                interaction_type="mim-conversion"
                )
            interactions_list.append(conversion)
            

        # Build CATALYSER Nodes (store multiple nodes even if repeated)
        if pd.notnull(catal_id):
            catal_node = None
            match_row = id_data_df[id_data_df[id_col] == catal_id]
            if not match_row.empty:
                info = match_row.iloc[0]
                catal_node = Node(
                    "GeneProduct", info[name_col], info[db_col], catal_id
                    )
                nodes_dict_info[catal_node.graph_id] = catal_node
            catal_inter = Interaction(
                source=catal_node.graph_id, target=None, interaction_type="mim-catalysis"
                )  # target -> the anchor of a conversion must be assigned
            interactions_list.append(catal_inter)


    # for key, node in nodes_dict_info.items():
    #     print(f"{key}: {vars(node)}")  # vars() returns the attributes as a dict
        
    for inter in interactions_list:
        print(f"{inter.source} → {inter.target}: {vars(inter)}")



############################# LAYOUT CREATION #################################

## PASAR A AQUÍ todo LO RELATIVO A COORDENADAS Y LAYOUT...

class Layout:
    def __init__(self, nodes, interactions):
        self.nodes = nodes          # Dict of nodes, keyed by graph_id
        self.interactions = interactions  # List of Interaction objects
        self._laid_out = False

    # def assign_layout(self):
    #     # A) Skeleton graph (metabolites only)
    #     G = nx.Graph()
    #     for n in self.nodes.values():
    #         if n.node_type.lower() == "enzyme":
    #             continue
    #         G.add_node(n.graph_id)
    #     for e in self.interactions:
    #         if e.type == "conversion":
    #             G.add_edge(e.source.graph_id, e.target.graph_id)

    #     # B) BFS layering (topology -> rows)
    #     layers = list(nx.bfs_layers(G, [next(iter(G.nodes))])) if G.number_of_nodes() else []

    #     # C) Place metabolites row-by-row (grid)
    #     k_max = max((len(L) for L in layers), default=1)
    #     span_max = max(0, (k_max - 1) * COL_GAP)

    #     placed = set()
    #     for ly, layer_nodes in enumerate(layers):
    #         k = len(layer_nodes)
    #         span = max(0, (k - 1) * COL_GAP)
    #         start_x = BOARD_MARGIN + (span_max - span) / 2.0
    #         y = BOARD_MARGIN + ly * LAYER_GAP
    #         for i, gid in enumerate(layer_nodes):
    #             node = self.nodes[gid]
    #             x = start_x + i * COL_GAP
    #             node.coords(x, y)
    #             placed.add(gid)

    #     # D) Temporarily place rest (enzymes/others) on extra row
    #     rest = [n for gid, n in self.nodes.items() if gid not in placed]
    #     for j, node in enumerate(rest):
    #         node.coords(BOARD_MARGIN + j * COL_GAP, BOARD_MARGIN + (len(layers) + 1) * LAYER_GAP)

    #     self._laid_out = True


########################### HELPER FUNCIONS ###################################

def is_wikipathways(id_str):
    """Check if source/target ID match Pathway WikiPathways database, 
    which can have the form WP1234, WP12345_r2, for example"""
    return re.match(r'WP\d{1,5}(_r\d+)?$', str(id_str)) is not None

def is_kegg(id_str):
    """Check if source/target ID match Pathway WikiPathways database, 
    which has the form hsa00120, for example"""
    return re.match(r'^[a-zA-Z]{2,4}\d{5}$', str(id_str)) is not None



# def get_columns(df, keyword):
#     for col in df.columns:
#         if keyword in col.lower():
#             return col
#     raise KeyError(f"Column with keyword '{keyword}' not found")


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


############################## XML BUILD #####################################

# PASAR A AQUÍ todo LO RELATIVO AL XML...
    # Gather all your nodes/interactions from TSVParser
    # Build the XML tree using their to_gpml() methods
    # Optionally assign coordinates/layout here

def xml_build(nodes, interactions, title="New Pathway", organism=None):
    pass
    def __init__(self, title="New Pathway", organism=None):
        self.title = title
        self.organism = organism
        self._laid_out = False 
        self.nodes = {}             # Node label
        self.interactions = []


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

# relations_file 
    # must have an id for the NODE > if this specific component is not 
    # found in a DB, then an ID such as metabo1,enz1... must be provided.

# IMPORTANT: Every ID from the relations file must be also present in the 
#            ID_data_file provided, otherwise this wont appear in the pathway