import datetime
import xml.etree.ElementTree as ET

class XMLBuilder:  
    def __init__(self, title, organism=None, nodes=None, interactions=None):
        self.title = title
        self.organism = organism
        self._laid_out = False 
        self.nodes = nodes or {}    # Node label
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