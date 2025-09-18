import csv
import xml.etree.ElementTree as ET  # create and manipulate XML structures
import xml.dom.minidom as minidom
import networkx as nx
import datetime
import uuid


# ----- layout / board -----
BOARD_MARGIN = 80.0     # margin around everything
LAYER_GAP = 140.0       # vertical separation between BFS layers
COL_GAP = 160.0         # approximate horizontal separation
ENZYME_OFFSET = 80.0    # perpendicular distance from the anchor to the enzyme
ENZYME_STACK_GAP = 60.0 # separation between multiple enzymes on the same anchor

############################### ID GEN ########################################
class IDGenerator:
    def __init__(self):
        self.counters = {}

    def generate_id(self, node_type):
        prefix = self.get_prefix(node_type)
        self.counters[node_type] += 1
        return f"{prefix}{self.counters[node_type]:04d}"

    def get_prefix(self, node_type):
        if node_type not in self.counters:
            self.counters(node_type) = 1
            prefix = node_type[0]
            return prefix.get(node_type)

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
        self.graph_id = idgenerator.generate_id(node_type)
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
        self.anchor_id = None  # only for conversions

    def coords(self):
        pass # Atención: para ser completado

    def to_gpml(self):
        interaction = ET.Element("Interaction")
        graphics = ET.SubElement(interaction, "Graphics", {"LineThickness": "1.0"})
        # source point         
        ET.SubElement(graphics, "Point", {    
            "X": str(self.source.x),
            "Y": str(self.source.y),
            "GraphRef": self.source.graph_id,
            "RelX": "0.0", "RelY": "1.0"
        })
        # target point
        target_point = {
            "X": str(self.target.x),
            "Y": str(self.target.y),
            "GraphRef": self.target.graph_id,
            "RelX": "0.0", "RelY": "-1.0",
            }
        if self.type.lower() == "conversion":
            target_point["ArrowHead"] = "mim-conversion"
            ET.SubElement(graphics, "Point", target_point)
            # add anchor in the middle
            self.anchor_id = uuid.uuid4().hex[:5]
            ET.SubElement(graphics, "Anchor", {
                "Position": "0.4",
                "Shape": "None"
            })
        elif self.type.lower() == "catalysis": # Atención: las flechas de catalisis no salen
            target_point["ArrowHead"] = "mim-catalysis"
            ET.SubElement(graphics, "Point", target_point)

        ET.SubElement(interaction, "Xref", {"Database": "", "ID":""})
        return interaction

############################## PATHWAY ########################################
# Define the elements of "Pathway"
# Adds nodes and interactions
# Creates and saves the pathway itself

class Pathway:  
    def __init__(self, title, organism="Homo sapiens"):
        self.title = title
        self.organism = organism
        self.nodes = {}
        self._nodes_by_id = {}
        self.interactions = []

    def add_node(self, node: Node):
        self.nodes[node.label] = node
        self._nodes_by_id[node.graph_id] = node

    def add_interaction(self, inter: Interaction):
        self.interactions.append(inter)

    def to_etree(self):
        """
        Define the elements of "Pathway" in GPML.
        
        Returns:
            ElementTree: GPML-based ElementTree
        """

        root = ET.Element("Pathway", {  
            "xmlns": "http://pathvisio.org/GPML/2013a",
            "Name": self.title,
            "Version": datetime.date.today().isoformat(),
            "Organism": self.organism
        })
        ET.SubElement(root, "Graphics", {"BoardWidth": "1000.0", "BoardHeight": "1000.0"}) 
                                      # ADD to change this according to total height/width
        # nodes
        for node in self.nodes.values():
            if node.x is None or node.y is None:
                raise ValueError(f"Nodo sin coordenadas: {node.label}")
            root.append(node.to_gpml())
        # interactions
        for inter in self.interactions:
            root.append(inter.to_gpml())
        # xml finishing lines
        ET.SubElement(root, "InfoBox", {"CenterX": "0.0", "CenterY": "0.0"})
        ET.SubElement(root, "Biopax")
        return root

    def assign_layout(self):
        """Layout BFS simple: capas en Y, orden por índice en X."""
        if not self.nodes:
            return
        G = nx.Graph()
        for n in self.nodes.values():
            G.add_node(n.graph_id)
        for e in self.interactions:
            G.add_edge(e.source.graph_id, e.target.graph_id)
        
        root_id = next(iter(self._nodes_by_id)) # first node
        layers = list(nx.bfs_layers(G,[root_id])) if G.number_of_nodes() else []

        # margin = 120.0
        # layer_gap = 140.0
        # board_w = 1000.0

        placed = set()
        for ly, layer_nodes in enumerate(layers):
            k = len(layer_nodes)
            for i, gid in enumerate(layer_nodes):
                node = self._nodes_by_id[gid]
                x = margin + (i + 1) * ((board_w - 2 * margin) / (k + 1))
                y = margin + ly * layer_gap
                node.coords(x, y)
                placed.add(gid)

        # nodos aislados (si los hay)
        rest = [n for gid, n in self._nodes_by_id.items() if gid not in placed]
        for j, node in enumerate(rest):
            node.coords(margin + (j + 1) * 160.0, margin + (len(layers) + 1) * layer_gap)

    def save(self, filename):
        xml_str = ET.tostring(self.to_etree(), encoding="utf-8") # atención: el encoding no sale!
        pretty_xml = minidom.parseString(xml_str).toprettyxml(indent="  ")
        with open(filename, "w", encoding="utf-8") as f:
            f.write(pretty_xml)



############################## PARSING ########################################

def parse_csv_to_pathway(csv_file, title="New Pathway"):
    pathway = Pathway(title)
    interactions_data = []

    with open(csv_file, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=",")
        for row in reader:
            node = Node(row["Node Type"], row["Node Label"], row.get("Database",""), row.get("Database_ID",""))
            if row["Node Label"] not in pathway.nodes:
                pathway.add_node(node)
            # acumular interacciones (se conectan por etiqueta)
            if row.get("Interaction Type") and row.get("Interaction With"):
                interactions_data.append((row["Node Label"], row["Interaction With"], row["Interaction Type"]))

    for source_label, target_label, inter_type in interactions_data:
        try:
            source = pathway.nodes[source_label]
            target = pathway.nodes[target_label]
        except KeyError as e:
            raise KeyError(f"Etiqueta de nodo no encontrada en CSV: {e}")
        pathway.add_interaction(Interaction(source, target, inter_type))

    return pathway


############################## MAIN ###########################################

if __name__ == "__main__":
    csv_file = "ruta_facil.csv"
    pathway = parse_csv_to_pathway(csv_file, "ruta_facil.csv")
    pathway.assign_layout()                                     # Node gets coordinates from this layout
    pathway.save("ruta_facil2.gpml")




# cd C:\\Users\\deyan\\Desktop\\BIOINFORMÁTICA\\1TFM


