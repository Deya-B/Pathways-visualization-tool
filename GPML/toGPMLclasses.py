import csv
import xml.etree.ElementTree as ET # create and manipulate XML structures
import xml.dom.minidom as minidom

################################ NODE #########################################
# Metabolite(Node) and enzyme(Node) > inherit from Node

class Node:             
    def __init__(self, node_type, label, database, db_id):
        self.node_type = node_type      # all the node properties
        self.label = label
        self.database = database
        self.db_id = db_id
        self.graph_id = f"id{hash(label) & 0xffff}"
        self.x = 100.0
        self.y = 100.0
        self.width = 90.0 + len(label) * 2
        self.height = 25.0

    def to_gpml(self):
        node_type = "GeneProduct" if self.node_type.lower() == "enzyme" else "Metabolite"
        datanode = ET.Element("DataNode", {
            "TextLabel": self.label,
            "GraphId": self.graph_id,
            "Type": node_type
        })
        graphics = ET.SubElement(datanode, "Graphics", {
            "CenterX": str(self.x),
            "CenterY": str(self.y),
            "Width": str(self.width),
            "Height": str(self.height),
            "FontSize": "12",
            "Valign": "Middle", 
            "Color": "0000ff" if node_type == "Metabolite" else "000000" 
        })
        xref = ET.SubElement(datanode, "Xref", {
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

    def to_gpml(self):
        interaction = ET.Element("Interaction")
                    # ET.Element(tag, attrib)
        graphics = ET.SubElement(interaction, "Graphics", {"LineThickness": "1.0"})
                 # ET.SubElement(parent, tag, attrib)
        ET.SubElement(graphics, "Point", {      # define arrow source
            "X": str(self.source.x),
            "Y": str(self.source.y),
            "GraphRef": self.source.graph_id,
            "RelX": "0.0", "RelY": "1.0"
        })
        ET.SubElement(graphics, "Point", {      # define arrow target
            "X": str(self.target.x),
            "Y": str(self.target.y),
            "GraphRef": self.target.graph_id,
            "RelX": "0.0", "RelY": "-1.0",
            "ArrowHead": "mim-conversion" if self.type == "Conversion" else "mim-catalysis"
        })
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

    def to_gpml(self):
        root = ET.Element("Pathway", {  # define the elements of the pathway
            "xmlns": "http://pathvisio.org/GPML/2013a",
            "Name": self.title,
            "Version": "Date",                               # ADD DATE CONTROL
            "Organism": self.organism
        })
        ET.SubElement(root, "Graphics", {"BoardWidth": "500.0", "BoardHeight": "800.0"})
                                        
                            # ADD to change this according to total height/width

        for node in self.nodes.values():
            root.append(node.to_gpml())
        for inter in self.interactions:
            root.append(inter.to_gpml())

        ET.SubElement(root, "InfoBox", {"CenterX": "0.0", "CenterY":"0.0"})
        ET.SubElement(root, "Biopax")
        return root
        
    def save(self, filename):
        tree = ET.ElementTree(self.to_gpml())
             # ET.ElementTree(element) > turns a root node into an entire tree
             #                         > we can then save the file with tree.write()
        tree.write(filename, encoding="UTF-8", xml_declaration=True)

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

############################### USAGE  ########################################

if __name__ == "__main__":
    csv_file = "ruta_facil.csv"
    pathway = parse_csv_to_pathway(csv_file, "ruta_facil.csv")
    pathway.save("ruta_facil.gpml")

xml_str = ET.tostring(pathway.to_gpml(), encoding ="utf-8")
pretty_xml = minidom.parseString(xml_str).toprettyxml(indent="  ")
with open ("ruta_facil.gpml", "w", encoding="utf-8") as f:
    f.write(pretty_xml)

# cd C:\\Users\\deyan\\Desktop\\BIOINFORMÁTICA\\1TFM