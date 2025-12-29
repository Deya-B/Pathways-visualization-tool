import re
import pandas as pd
import logging

from .model import Node, Interaction, ConversionInteraction

# Global variables with accepted variants
PATHWAY_DB_NAMES = {"wikipathways", "reactome", "kegg pathway"}

class Parser:
    def __init__(self, id_data_df, relations_df):
        self.id_data_df = id_data_df
        self.relations_df = relations_df
        self.db_index = {}
        self.nodes = {}
        self.interactions = []
        self.conversion_by_pair = {}  # (source_id, target_id) -> ConvInter 
                                      # to create only one ConvInter per pair
        # Get column names
        self.id_col, self.db_col, self.name_col = (
            self.id_data_df.columns[0:3])
        self.source_col, self.target_col, self.catalyser_col = (
            self.relations_df.columns[0:3])
        self._clean_id_metadata()

    def _clean_field(self, value, *, empty_to_none=True):
        """Normalize string-like fields from TSV/CSV."""
        if isinstance(value, str):
            value = value.strip()
        if pd.isna(value):
            value = ""
        if empty_to_none and value == "":
            return None
        return value


    def _clean_id_metadata(self):
        """Normalize key columns in ID metadata (strip spaces, etc.)."""
        for col in (self.id_col, self.db_col, self.name_col):
            self.id_data_df[col] = self.id_data_df[col].apply(
                lambda v: self._clean_field(v, empty_to_none=False)
            )


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
        
        # WikiPathways: WP + 4-5 dígitos
        if re.fullmatch(r"WP\d{4,5}.*", str(node_id)):
            return True

        # KEGG Pathway típico: 3 letras (organismo) + 5 dígitos
        if re.fullmatch(r"[a-z]{3}\d{5}", str(node_id)):
            return True

        # Buscar el registro en el DataFrame
        match = self.id_data_df[self.id_data_df[self.id_col] == node_id]
        if match.empty:
            return False
        db_value = str(match.iloc[0][self.db_col]).lower()
        if any(name in db_value for name in PATHWAY_DB_NAMES):
            return True
        return False


    def _get_create_node (self, node_id):
        """For source and target nodes"""
        node_id = self._clean_field(node_id)  # normalize ID
        if pd.isnull(node_id):
            return None
        
        # If already created, return it
        if node_id in self.db_index:
            graph_id = self.db_index[node_id]
            return self.nodes[graph_id]
        
        # Find metadata row
        match = self.id_data_df[self.id_data_df[self.id_col] == node_id]
        if match.empty:
            logging.warning(f"ID '{node_id}' not found in ID metadata; skipping node")
            return None
        
        info = match.iloc[0]
        is_pathway_flag = self._is_pathway(node_id)
        node_type = "Pathway" if is_pathway_flag else "Metabolite"
        colour = self.get_node_colour(node_type)
        
        label = self._clean_field(info[self.name_col], empty_to_none=False)
        database = self._clean_field(info[self.db_col], empty_to_none=False)

        # Extract info from tabular data and transform into a Node object
        node = Node(node_type, label, database, node_id, colour)
        self.nodes[node.graph_id] = node

        # Record that this metabolite is already mapped
        self.db_index[node_id] = node.graph_id     
        return node


    def build_conversions(self, row):
        # Obtain source and target ID
        raw_source = getattr(row, self.source_col)
        raw_target = getattr(row, self.target_col)

        source_id = self._clean_field(raw_source)
        target_id = self._clean_field(raw_target)

        # Build nodes
        source_node = self._get_create_node(source_id)
        target_node = self._get_create_node(target_id)

        if not (source_node and target_node):
            return None

        key = (source_id, target_id)

        # Reuse existing ConversionInteraction if present
        if key in self.conversion_by_pair:
            return self.conversion_by_pair[key]
        
        conversion = ConversionInteraction(source_node, target_node)
        self.interactions.append(conversion)
        self.conversion_by_pair[key] = conversion
        return conversion
        

    def _create_catalyser_node(self, catal_id):
        """For catalyser nodes"""
        catal_id = self._clean_field(catal_id)
        if pd.isnull(catal_id):
            return None
        
        # Find metadata row
        match = self.id_data_df[self.id_data_df[self.id_col] == catal_id]
        if match.empty:
            logging.warning(f"ID '{catal_id}' not found in ID metadata; skipping node")
            return None
        
        info = match.iloc[0]
        colour = self.get_node_colour("GeneProduct")
            
        label = self._clean_field(info[self.name_col], empty_to_none=False)
        database = self._clean_field(info[self.db_col], empty_to_none=False)

        node = Node("GeneProduct", label, database, catal_id, colour)
        self.nodes[node.graph_id] = node
        return node
    

    def build_catalysis(self, row, conv_interaction):
        # Obtain catalyser ID
        raw_catal = getattr(row, self.catalyser_col)
        catal_id = self._clean_field(raw_catal)
        if pd.isnull(catal_id):
            return None # no enzyme -> no interaction
        
        catal_node = self._create_catalyser_node(catal_id)
        if catal_node is None:
            return None
        
        anchor_id  = conv_interaction.anchor_id # shared anchor
        
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