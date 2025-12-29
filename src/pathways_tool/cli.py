from .xml_builder import XMLBuilder
import xml.etree.ElementTree as ET  # create and manipulate XML structures
import logging

from .parser import Parser
from .layout import Layout
from .io_utils import read_csv
from .config import GPMLConfig, setup_logging


def run_from_config(config: GPMLConfig):
    setup_logging(config)
    main(
        pathway_title=config.pathway_title,
        organism=config.organism,
        ID_data_file=config.id_data_file,
        relations_file=config.relations_file,
        output_filename=config.output_filename,
        delimiter=config.delimiter,
    )


def main(pathway_title, organism, 
         ID_data_file, relations_file, output_filename, 
         delimiter="\t"):
    """Main function to build GPML from ID metadata and relations files."""
    id_data_df = read_csv(ID_data_file, sep=delimiter)
    relations_df = read_csv(relations_file, sep=delimiter)

    # Strip spaces in ID_metadata
    for col in id_data_df.select_dtypes(include="object").columns:
        id_data_df[col] = id_data_df[col].str.strip()

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

    try:
        tree.write(output_filename, encoding="utf-8", xml_declaration=True)
        logging.info(f"Successfully wrote GPML to...\n  {output_filename}")
    except Exception:
        logging.exception(f"Failed to write GPML file {output_filename}")
        raise