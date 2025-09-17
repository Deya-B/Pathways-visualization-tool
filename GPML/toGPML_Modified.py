import csv
import xml.etree.ElementTree as ET  # create and manipulate XML structures
import xml.dom.minidom as minidom
import networkx as nx
import datetime
import uuid

################################ NODE #########################################
# Metabolite(Node) and enzyme(Node) > inherit from Node

class Node:
    """
    Contains all the Node properties.
    
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
        # all the node properties:
        self.node_type = node_type      
        self.label = label
        self.database = database
        self.db_id = db_id
        self.graph_id = f"id{hash(label) & 0xffff}"
        self.x, self.y = None, None
        self.width = 90.0 + len(label) * 2
        self.height = 25.0

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
    def __init__(self, source, target, interaction_type):
        self.source = source
        self.target = target
        self.type = interaction_type
        self.anchor_id = None  # only for conversions

    def to_gpml(self):
        interaction = ET.Element("Interaction")
                    # ET.Element(tag, attrib)
        graphics = ET.SubElement(interaction, "Graphics", {"LineThickness": "1.0"})
                 # ET.SubElement(parent, tag, attrib)
        
        # source point         
        ET.SubElement(graphics, "Point", {      # define arrow source
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

            # add anchor in middle
            self.anchor_id = uuid.uuid4().hex[:5]
            ET.SubElement(graphics, "Anchor", {
                "Position": "0.4",
                "Shape": "None"
            })

        elif self.type.lower() == "catalysis":
            target_point["ArrowHead"] = "mim-catalysis"
            ET.SubElement(graphics, "Point", target_point)

        ET.SubElement(interaction, "Xref", {"Database": "", "ID":""})
        return interaction

############################## PATHWAY ########################################
# Add everything to a pathway (contains nodes + interactions) 
# Export to GPML

class Pathway:  
    def __init__(self, title, organism="Homo sapiens"):
        self.title = title
        self.organism = organism
        self.nodes = {}
        self.interactions = []

    def add_node(self, node):
        self.nodes[node.label] = node

    def add_interaction(self, interaction):
        self.interactions.append(interaction)

    def assign_layout(self):
        G = nx.DiGraph()
        for node_id, node in self.nodes.items():
            G.add_node(node_id)
        for inter in self.interactions:
            G.add_edge(inter.source.label, inter.target.label)

        pos = nx.spring_layout(G, k=2, scale=400)
        for label, (x, y) in pos.items():
            node = self.nodes[label]
            node.x = float(x * 500 + 250)   # rescale & center
            node.y = float(y * 500 + 400)

    def to_gpml(self):
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
        
        ET.SubElement(root, "Graphics", {"BoardWidth": "1000.0", "BoardHeight": "1000.0"}) # ADD to change this according to total height/width                                                              
        for node in self.nodes.values():
            root.append(node.to_gpml())
        for inter in self.interactions:
            root.append(inter.to_gpml())
        ET.SubElement(root, "InfoBox", {"CenterX": "0.0", "CenterY":"0.0"})
        ET.SubElement(root, "Biopax")
        return root
        
    def save(self, filename):
        xml_str = ET.tostring(self.to_gpml(), encoding="utf-8")
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
            node = Node(row["Node Type"], row["Node Label"], row["Database"], row["Database_ID"])
            pathway.add_node(node)
            if row["Interaction Type"] and row["Interaction With"]:
                interactions_data.append((row["Node Label"], row["Interaction With"], row["Interaction Type"]))

    for source_label, target_label, inter_type in interactions_data:
        source = pathway.nodes[source_label]
        target = pathway.nodes[target_label]
        pathway.add_interaction(Interaction(source,target,inter_type))

    return pathway



if __name__ == "__main__":
    csv_file = "ruta_facil.csv"
    pathway = parse_csv_to_pathway(csv_file, "ruta_facil.csv")
    pathway.assign_layout()                                     # Node gets coordinates from this layout
    pathway.save("ruta_facil2.gpml")




# cd C:\\Users\\deyan\\Desktop\\BIOINFORMÁTICA\\1TFM


