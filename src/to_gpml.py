# TODO: ADD Module docstring (what this file is about)

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

# TODO: add organs > make like groups
# TODO: incorporate a yaml config file
# TODO: incorporate hypothetical/multi-step reaction

################################ NODE #########################################
## Description:
    # GPML DataNode (metabolite, pathway or enzyme). 
    # Contains the Node properties.
## Attributes:
        # node_type: Semantic type. "Metabolite" "Pathway" or "GeneProduct"
        # label: Display text
        # database: Xref DB name
        # db_id: Xref ID
        # graph_id: Unique GraphId in GPML.
        # x, y: Center coordinates (PathVisio uses centers)
        # width, height: Node box size in pixels
## Methods:
        # coords: Set center coordinates
        # to_gpml: convert data to GPML

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
            "ShapeType": "None" if self.node_type == "Pathway" 
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
    def __init__(self, source, target, interaction_type):
        self.graph_id = idgenerator.new("interaction")
        self.source = source                
        self.target = target               
        self.type = (interaction_type or "").lower()
        self.anchor_xy = None


    def bind_to_anchor(self, anchor_id, xy=None):
        self.target = anchor_id
            

    def to_gpml(self):
        """GPML for catalysis interactions."""
        inter_el = ET.Element("Interaction", {
            "GraphId": self.graph_id
        })
        graphics = ET.SubElement(inter_el, "Graphics", {
            "LineThickness": "1.0"
        })

        if self.type == "mim-catalysis":
            enz = self.source       # Node = GeneProduct
            anchor_id = self.target
            ex, ey = enz.x, enz.y   # Enzyme coordinates
        
            # Origin point: enzyme (center)
            ET.SubElement(graphics, "Point", {
                "X": str(ex),
                "Y": str(ey),
                "GraphRef": enz.graph_id,
                "RelX": "0.0",
                "RelY": "0.0"
            })
            # Destination point: anchor (only GraphRef=anchor_id, without X/Y)
            ax, ay = self.anchor_xy if self.anchor_xy is not None else (ex, ey)
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
        self.anchor_pos = 0.5


    def to_gpml(self):
        """GPML for an anchored conversion."""
        inter_el = ET.Element("Interaction", {
            "GraphId": self.graph_id
        })
        graphics = ET.SubElement(inter_el, "Graphics", {
            "ConnectorType": "Elbow",
            "LineThickness": "1.0"
        })
        
        # Source point (Lower edge of origin node)
        # NOTE: Attachment point: 
                # left edge: (-1, 0) 
                # right: (1, 0)
                # top: (0, -1) 
                # bottom: (0, 1)
        sx = self.source.x
        sy = self.source.y + self.source.height / 2.0   
                                                         
        ET.SubElement(graphics, "Point", {    
                "X": str(sx),
                "Y": str(sy),              
                "GraphRef": self.source.graph_id, 
                "RelX": "0.0", 
                "RelY": "1.0"
        })                                              
        # Target point (Top edge of the destination node)
        tx = self.target.x
        ty = self.target.y - self.target.height / 2.0
        ET.SubElement(graphics, "Point", {
            "X": str(tx),
            "Y": str(ty), 
            "GraphRef": self.target.graph_id, 
            "RelX": "0.0", 
            "RelY": "-1.0",
            "ArrowHead": "mim-conversion"
        })
            
        # Anchor with GraphId: uses anchor_id and anchor_pos already 
        #   calculated in Layout                                         
        ET.SubElement(graphics, "Anchor", {
            "GraphId": self.anchor_id, 
            "Position": str(self.anchor_pos), 
            "Shape": "None"
        })

        ET.SubElement(inter_el, "Xref", {"Database": "", "ID":""})
        return inter_el


###################### Layout/Board CONFIGURATION #############################

BOARD_MARGIN = 200.0    # margin around everything
LAYER_GAP = 120.0       # vertical separation between BFS layers
COL_GAP = 250.0         # approximate horizontal separation

############################# LAYOUT CREATION #################################

class Layout:
    def __init__(self, nodes, interactions):
        self.nodes = nodes                # Dict of nodes, keyed by graph_id
        self.interactions = interactions  # List of Interactions
        self.boardwidth = None
        self.boardheight = None
        self._laid_out = False


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
        G = nx.Graph()

        for node in self.nodes.values():
            if node.node_type == "GeneProduct": # Ignore enzymes for now
                continue
            G.add_node(node.graph_id)

        for inter in self.interactions:
            if isinstance(inter, ConversionInteraction):
                G.add_edge(inter.source.graph_id, inter.target.graph_id)


        # See graph
        # labels = {n.graph_id: n.label for n in self.nodes.values()
        #   if n.node_type != "GeneProduct"}
        # pos = nx.spring_layout(G)  # o cualquier layout
        # nx.draw(G, pos, labels=labels, node_color="lightblue", node_size=800, font_size=8)
        # plt.show()


        # BFS layering, get a list of rows with the topology           
        has_nodes = G.number_of_nodes() > 0
        if has_nodes: # creates lists of nodes by depth from a root/start node
            start_node = next(iter(G.nodes))
            bfs_layers_iterator = nx.bfs_layers(G, [start_node])
            layers = list(bfs_layers_iterator) 
        else:
            layers = []

        # Setting x,y coords:
        # Step 1 - Get maximal layer size and width
        k_max = max((len(L) for L in layers), default=1) # find widest row
        span_max = max(0, (k_max - 1) * COL_GAP) # width occupied by the widest row
        
        # Step 2 - Loop through layers/rows
        placed = set() 
        for ly, layer_nodes in enumerate(layers): 
            k = len(layer_nodes) # number of nodes in this layer
            span = max(0, (k - 1) * COL_GAP) # row width for this layer
            start_x = BOARD_MARGIN + (span_max - span) / 2.0 # center this row  under the widest row
            y = 80 + ly * LAYER_GAP # y-position for this row
        
        # Step 3 - Assign positions to the nodes
            for i, graph_id in enumerate(layer_nodes):
                node = self.nodes[graph_id]
                x = start_x + i * COL_GAP # horizontal position in the row
                node.coords(x, y)
                placed.add(graph_id) # save graph_ids already positioned

        # Temporarily place enzymes on extra row
        rest = [n for graph_id, n in self.nodes.items() if graph_id not in placed]
        for j, node in enumerate(rest):
            node.coords(BOARD_MARGIN + j * COL_GAP, BOARD_MARGIN + (len(layers) + 1) * LAYER_GAP)

        self._laid_out = True
    

    def layout_anchors(self):
        """Compute conversion anchor coordinates. To obtain anchor ID and xy"""
        for inter in self.interactions:
            if isinstance(inter, ConversionInteraction):
                xy = self.compute_anchor_xy(inter) # calculate anchor coords
                inter._anchor_xy = xy              # update atribute
        

    def layout_catalysis(self):
        """Place enzyme nodes at their anchors once catalysis interactions exist."""
        # Build anchor lookup: anchor_id -> anchor_xy from conversion inter
        anchor_xy_dict = {     
            inter.anchor_id: inter._anchor_xy
            for inter in self.interactions
            if isinstance(inter, ConversionInteraction) 
                and inter.anchor_id is not None
                and inter._anchor_xy is not None}

        for inter in self.interactions:
            if inter.type == "mim-catalysis":
            # Catalysis: source = enzyme/GeneProduct, target = anchor_id
                enz = inter.source       
                anchor_id = inter.target
                if anchor_id in anchor_xy_dict and enz is not None:
                    ax, ay = anchor_xy_dict[anchor_id]
                    enz.coords(ax, ay) # place enzyme exactly at anchor
                    inter.anchor_xy = (ax, ay)


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
        
        # WikiPathways: WP + 4-5 dígitos + opcional _rN
        if re.fullmatch(r"WP\d{4,5}(?:_r\d+)?", str(node_id)):
            return True

        # KEGG Pathway típico: 3 letras (organismo) + 5 dígitos
        if re.fullmatch(r"[a-z]{3}\d{5}", str(node_id)):
            return True

        # Buscar el registro en el DataFrame
        match = self.id_data_df[self.id_data_df[self.id_col] == node_id]
        if match.empty:
            return False
        db_value = str(match.iloc[0][self.db_col])
        return "pathway" in db_value.lower() # Contiene "pathway"
    

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
        is_pathway_flag = self._is_pathway(node_id)
        node_type = "Pathway" if is_pathway_flag else "Metabolite"
        colour = self.get_node_colour(node_type)
        
        # Extract info from tabular data and transform into a Node object
        node = Node(node_type, info[self.name_col], info[self.db_col], node_id, colour)
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
        colour = self.get_node_colour("GeneProduct")
        node = Node("GeneProduct", info[self.name_col], info[self.db_col], catal_id, colour)
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

# PASAR A AQUÍ todo LO RELATIVO AL XML...
    # Gather all your nodes/interactions from TSVParser
    # Build the XML tree using their to_gpml() methods

    # header creation
    # mapping nodes → XML
    # mapping interactions → XML
    # handling attributes
    # board size
    # naming/organism metadata
    # consistency

class XMLBuilder:  
    def __init__(self, title, organism=None, nodes=None, interactions=None):
        self.title = title
        self.organism = organism
        self._laid_out = False 
        self.nodes = nodes or {}             # Node label
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
    # DataFrames from Data and Relations files
    id_data_df = pd.read_csv(ID_data_file, sep=delimiter, encoding="utf-8")
    relations_df = pd.read_csv(relations_file, sep=delimiter, encoding="utf-8")

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



    # Checkers
    # for key, node in builder.nodes.items():
    #     print(f"{key}: {vars(node)}")  # vars() returns the attributes as a dict

    # for inter in conversions_list:
    #     source_label = inter.source.graph_id if inter.source else "None"
    #     target_label = inter.target.graph_id if inter.target else "None"
    #     print(f"{source_label} → {target_label}: {vars(inter)}")

    # for inter in catalysis_list:
    #     source_label = inter.source.graph_id if inter.source else "None"
    #     target_label = inter.target
    #     print(f"{source_label} → {target_label}: {vars(inter)}")


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
    # ID_data_file = "c:/Users/dborrotoa/Desktop/TFM/src/examples/data/2-Secondary_BA_Synthesis_CA-Based_reactions_updated.tsv"
    # relations_file = "c:/Users/dborrotoa/Desktop/TFM/src/examples/data/relationships.tsv"
    # output_filename = "c:/Users/dborrotoa/Desktop/TFM/src/examples/gpml/ruta.gpml"
    #home
    ID_data_file = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/data/1-Alternative-Acidic_BA_updated.tsv"
    relations_file = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/data/relationships2.tsv"
    output_filename = "C:/Users/deyan/Desktop/BIOINFORMATICA/1TFM/src/examples/gpml/Alternative-NoDupsElbows.gpml"
    
    pathway_title = "Alternative/Acidic BA biosynthesis pathway"
    organism = "Homo sapiens, Mus-musculus"

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
    #   as "metabo1", "enz1"... must be provided.

# Data_file:
    # Must have the following columns in order:
    #   DataBase ID, DataBase Name, Component Name
    # IMPORTANT: Every ID from the relations file must be also present in 
    #   the ID_data_file provided. Otherwise this wont appear in the pathway.

# Pathway node must be defined somewhere

# ET.indent(tree, space="  "), part of xml.etree.ElementTree requires 
#   Python ≥ 3.9 to work