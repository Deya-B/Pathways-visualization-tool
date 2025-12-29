import datetime
import xml.etree.ElementTree as ET

class XMLBuilder:  
    """Build a GPML XML tree from nodes and interactions.

    This helper assembles a complete GPML ``Pathway`` element, including
    nodes, interactions, and minimal InfoBox/Biopax elements, ready to be
    written to disk.
    """

    def __init__(self, title, organism=None, nodes=None, interactions=None):
        """Initialize the XML builder with pathway metadata and graph data.

        Parameters
        ----------
        title : str
            Title of the pathway to store in the GPML ``Name`` attribute.
        organism : str, optional
            Organism name recorded in the GPML metadata, for example
            "Homo sapiens". If None, the attribute is left as None.
        nodes : dict, optional
            Mapping from graph ID to ``Node`` instances to be serialized
            as GPML DataNode elements.
        interactions : list, optional
            List of Interaction/ConversionInteraction instances to be
            serialized as GPML Interaction elements.
        """
        self.title = title
        self.organism = organism
        self._laid_out = False 
        self.nodes = nodes or {}    # Node label
        self._nodes_by_id = {}      # Node graph_id
        self.interactions = interactions or []


    def to_etree(self):
        """Create an ElementTree root for the GPML pathway.

        Returns
        -------
        xml.etree.ElementTree.Element
            Root ``Pathway`` element containing all nodes, interactions,
            and closing InfoBox/Biopax elements.
        """
        root = self.xml_header()

        # Nodes
        for node in self.nodes.values(): # iterate through the Node objects
            if node.x is None or node.y is None:
                raise ValueError(f"Nodo sin coordenadas: {node.label}")
            root.append(node.to_gpml())

        # Interactions
        for inter in self.interactions:
            root.append(inter.to_gpml())

        # InfoBox & Biopax xml finishing lines
        ET.SubElement(root, "InfoBox", {"CenterX": "0.0", "CenterY": "0.0"})
        ET.SubElement(root, "Biopax")
        return root
    

    def xml_header(self):
        """Create the GPML Pathway root element with header metadata.

        The header includes namespace, pathway title, current date as the
        version, organism, and a basic Graphics element with board size.

        Returns
        -------
        xml.etree.ElementTree.Element
            Newly created ``Pathway`` element ready to receive children.
        """

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