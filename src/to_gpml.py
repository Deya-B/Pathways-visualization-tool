# Module docstring (what this file is about)
import csv
import datetime

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


################################ PARSER #######################################

class TSVParser:
    def __init__(self, ID_data_file, relations_file, delimiter="\t"):
        self.delimiter = delimiter
        self.ID_data_file = ID_data_file
        self.relations_file = relations_file
        self.id_data_df = self.read(self.ID_data_file)
        self.relations_df = self.read(self.relations_file)
        self.node_relations = []

    def read(self, file):
        """Read file and return DataFrame"""
        return pd.read_csv(file, sep=self.delimiter, encoding="utf-8")
    
    def find_column(self, df, keyword):
        for col in df.columns:
            if keyword in col.lower():
                return col
        raise KeyError(f"Column with keyword '{keyword}' not found")

    def get_node_relations(self):
        """Extracts relation tuples and saves in self.node_relations."""
        # find column names containing source, target and catalyser
        source_col = self.find_column(self.relations_df, "source")
        target_col = self.find_column(self.relations_df, "target")
        catalyser_col = self.find_column(self.relations_df, "catal")
        for _, row in self.relations_df.iterrows():
            source_db_id = row[source_col]
            target_db_id = row[target_col]
            catalyser_db_id = row[catalyser_col]
            # Append as tuple (source, target, catalyser)
            self.node_relations.append((source_db_id, target_db_id, catalyser_db_id))

    def get_node_ids(self):
        """Get IDs"""
        # get source node ids
        source_ids = [source_db_id for source_db_id,target_db_id,catalyser_db_id 
                      in self.node_relations]
        # removing the NaN
        source_ids_clean = [x for x in source_ids if pd.notnull(x)]

        # get catalyser ids
        catalyser_ids = [catalyser_db_id 
                         for source_db_id,target_db_id,catalyser_db_id
                         in self.node_relations]
        # removing the NaN
        catalyser_ids_clean = [x for x in catalyser_ids if pd.notnull(x)]

        return source_ids_clean, catalyser_ids_clean

    def get_node_info(self):
        pass
        # """Extract information from nodes"""
        # source_id, target_id = get_node_ids()
        # source_data = self.id_data_df[self.id_data_df["ID"].isin(source_ids_clean)]
        # print(source_data)
        # node_list = [Node.from_row(row) for _, row in source_data.iterrows()]

        # catalysis_data = self.id_data_df[self.id_data_df["ID"].isin(catalyser_ids_clean)]
        # print(catalysis_data)

    def get_interactions(self):
        pass

    def save(self):
        pass


################################ MAIN #########################################

def execute_main():
    # ID_data_file = "c:/Users/dborrotoa/Desktop/TFM/src/examples/data/2-Secondary_BA_Synthesis_CA-Based_reactions_updated.tsv"
    # relations_file = "c:/Users/dborrotoa/Desktop/TFM/src/examples/data/relationships.tsv"
    #home
    ID_data_file = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/data/2-Secondary_BA_Synthesis_CA-Based_reactions_updated.tsv"
    relations_file = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/data/relationships.tsv"
    
    pathway_title = "Secondary BA Synthesis (CA-Based reactions)"
    organism = "Homo sapiens, Mus-musculus"

    parser = TSVParser(ID_data_file, relations_file)
    parser.get_node_relations()
    parser.get_node_ids()


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
        pass


    # Extract each row from tabular data and transform into a Node
    @classmethod
    def from_row(cls, row):
        node_type, label, database, db_id = row[:4]  # Adjust columns as needed
        return cls(node_type, label, database, db_id)
    

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
        self.interaction_type = (interaction_type or "").lower()
        self.anchor_pos = float(anchor_pos) if anchor_pos is not None else self.default_anchor_pos
        self.anchor_id = anchor_id          
        self._anchor_xy = None              # private: to be computed ONLY if coordinates exist


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


############################# LAYOUT CREATION #################################

## PASAR A AQUÍ TODO LO RELATIVO A COORDENADAS Y LAYOUT...

class Layout:
    def __init__(self):
        pass


############################### XML BUILD #####################################

## PASAR A AQUÍ TODO LO RELATIVO AL XML...
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
    execute_main()

## NOTES
# default delimiter="\t" 
    # can be changed as a 3rd parameter. Example: 
        # TSVParser(ID_data_file, relations_file, delimiter=";")
# relations_file 
    # must have an id for the NODE > if this specific component is not 
    # found in a DB, then an ID such as metabo1,enz1... must be provided.