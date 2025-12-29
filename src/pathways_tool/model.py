import xml.etree.ElementTree as ET

from .idgen import idgenerator


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

        Parameters
        ----------
        x : float
            Center X in pixels.
        y : float
            Center Y in pixels.
        """
        self.x, self.y = float(x), float(y)
        pass


    def to_gpml(self):
        """Serialize the node to GPML XML element.

        Returns
        -------
        xml.etree.ElementTree.Element
            XML element representing the DataNode.
        """
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
            "ShapeType": "RoundedRectangle" if self.node_type == "Pathway" 
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


class Interaction:
    def __init__(self, source, target, interaction_type):
        self.graph_id = idgenerator.new("interaction")
        self.source = source                
        self.target = target               
        self.type = (interaction_type or "").lower()
        self.anchor_xy = None


    def bind_to_anchor(self, anchor_id, xy=None):
        """Bind the interaction to an anchor point.

        Parameters
        ----------
        anchor_id : str
            ID of the anchor.
        xy : tuple, optional
            Coordinates of the anchor.
        """
        self.target = anchor_id
        if xy is not None:
            self.anchor_xy = xy
            

    def to_gpml(self):
        """Serialize the interaction to GPML XML element.

        Returns
        -------
        xml.etree.ElementTree.Element
            XML element representing the Interaction.
        """
        inter_el = ET.Element("Interaction", {
            "GraphId": self.graph_id
        })
        graphics = ET.SubElement(inter_el, "Graphics", {
            "LineThickness": "1.0"
        })

        if self.type == "mim-catalysis":
            enz = self.source       # GeneProduct node
            anchor_id = self.target
            ex, ey = enz.x, enz.y   # Enzyme coordinates
        
            # Origin point: enzyme center
            ET.SubElement(graphics, "Point", {
                "X": str(ex),
                "Y": str(ey),
                "GraphRef": enz.graph_id,
                "RelX": "0.0",
                "RelY": "0.0"
            })
            # Destination point: anchor (use anchor_xy)
            if self.anchor_xy is None:
                # defensive: fall back to enzyme center
                ax, ay = ex, ey
            else:
                ax, ay = self.anchor_xy

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
        self.anchor_pos = 0.8
        self.attachment = "vertical"  # attachment mode: "vertical" or "side"

        # GPML point info (filled by Layout)
        self._src_point = None  # dict with X,Y,RelX,RelY
        self._tgt_point = None  # dict with X,Y,RelX,RelY


    def to_gpml(self):
        """Serialize the conversion interaction to GPML XML element.

        Returns
        -------
        xml.etree.ElementTree.Element
            XML element representing the Interaction.
        """
        inter_el = ET.Element("Interaction", {
            "GraphId": self.graph_id
        })
        graphics = ET.SubElement(inter_el, "Graphics", {
            "ConnectorType": "Elbow",
            "LineThickness": "1.0"
        })
        
        # Use points computed by Layout        
        ET.SubElement(graphics, "Point", self._src_point)
        ET.SubElement(graphics, "Point", self._tgt_point)

        # Anchor with GraphId: uses anchor_id and anchor_pos calculated in Layout                                        
        ET.SubElement(graphics, "Anchor", {
            "GraphId": self.anchor_id, 
            "Position": str(self.anchor_pos), 
            "Shape": "None"
        })

        ET.SubElement(inter_el, "Xref", {"Database": "", "ID":""})
        return inter_el
