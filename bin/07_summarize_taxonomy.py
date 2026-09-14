#!/usr/bin/env python3
"""
Description: Summarizes ortholog hits across eukaryotic taxonomy using ETE3.
Features: 
- Filters for Eukaryota (TaxID 2759)
- Collapses major clades for better visualization
- Adds color-coded markers for different search sources (BLASTP, OrthoDB, etc.)
"""

import os
import argparse
import pandas as pd
from ete3 import NCBITaxa, TreeStyle, TextFace, CircleFace

def get_args():
    parser = argparse.ArgumentParser(description="Summarize taxonomic hits and render a tree.")
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing TaxID files.")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save the results.")
    return parser.parse_args()

def is_eukaryota(taxid, ncbi):
    try:
        lineage = ncbi.get_lineage(int(taxid))
        return 2759 in lineage
    except:
        return False


def main():
    args = get_args()
    ncbi = NCBITaxa()
    os.makedirs(args.output_dir, exist_ok=True)

    # merge all taxids from different sources
    files = {
        'BLASTP': 'BLASTP_taxids.txt',
        'OrthoDB': 'OrthoDB_taxids.txt',
        'InterPro': 'InterPro_taxids.txt',
        'TBLASTN': 'TBLASTN_taxids.txt'
    }

    # Load and combine data
    dfs = []
    for label, filename in files.items():
        path = os.path.join(args.input_dir, filename)
        if os.path.exists(path):
            tmp_df = pd.read_csv(path, sep='\t', header=None, names=['taxid'], dtype={'taxid': str})
            tmp_df['source'] = label
            dfs.append(tmp_df)
    
    if not dfs:
        print("[ERROR] No input TaxID files found in the specified directory.")
        return

    combined_df = pd.concat(dfs)
    combined_df['taxid'] = combined_df['taxid'].astype(str)

    # Filter for Eukaryota
    print("[INFO] Filtering for Eukaryotic lineages...")
    unique_ids = combined_df['taxid'].unique().tolist()
    filtered_ids = [tid for tid in unique_ids if is_eukaryota(tid, ncbi)]
    
    # Generate topology
    tree = ncbi.get_topology(filtered_ids, intermediate_nodes=True)
    name_map = ncbi.get_taxid_translator(filtered_ids)
    source_map = combined_df.groupby('taxid')['source'].apply(lambda x: set(x)).to_dict()

    # colors
    colors = {'BLASTP': 'red', 'OrthoDB': 'yellow', 'InterPro': 'blue', 'TBLASTN': 'green'}

    # Mark the sources
    for node in tree.traverse():
        if node.is_leaf():
            taxid = int(node.name)
            node.add_feature("original_taxid", str(taxid))
            node.name = name_map.get(int(taxid), taxid)
            
            sources = source_map.get(str(taxid), set())
            node.add_feature("original_taxid", str(taxid)) #keep original taxid for reference

            circle_faces =[] 
            if "BLASTP" in sources:
                circle_faces.append(CircleFace(10,"red"))
            if "OrthoDB" in sources: 
                circle_faces.append(CircleFace(10, "yellow"))
            if "InterPro" in sources:
                circle_faces.append(CircleFace(10, "blue"))
            if "TBLASTN" in sources:
                circle_faces.append(CircleFace(10, "green"))

            # Map sources to colors
            for idx, face in enumerate(circle_faces):
                node.add_face(face, column=idx + 1, position='aligned')

    # Collapse
    collapse_taxa = ["Bilateria", "Cnidaria", "Placozoa", "Ctenophora", "Porifera", "Choanoflagellata",
                     "Filasterea", "Ichthyosporea", "Aphelida", "Fungi", "Amoebozoa", "Streptophyta",
                     "Chlorophyta", "Prasinodermophyta", "Stramenopiles", "Alveolata", "Rhizaria",
                     "Haptista", "Metamonada", "Discoba", "Rhodophyta", "Cryptophyceae"]

    for taxon in collapse_taxa:
        taxid = ncbi.get_name_translator([taxon]).get(taxon)
        if taxid:
            node = tree.search_nodes(name=str(taxid[0]))
            if node:
                node = node[0]
                original_leaves = node.get_leaves()
                num_leaves = len(original_leaves)
            
            # Gather all sources from the collapsed leaves
                all_sources = set()
                for leaf in original_leaves:
                    taxid_str = getattr(leaf, 'original_taxid', leaf.name)
                    all_sources.update(source_map.get(taxid_str, set()))
            
            # Hide children to simulate collapse
                node.children = []
                node.name = f"{taxon} ({num_leaves} species)"

            # Add colored circles horizontally for collapsed nodes
                circle_faces = []
                if 'BLASTP' in all_sources:
                    circle_faces.append(CircleFace(5, 'red'))
                if 'OrthoDB' in all_sources:
                    circle_faces.append(CircleFace(5, 'yellow'))
                if 'InterPro' in all_sources:
                    circle_faces.append(CircleFace(5, 'blue'))
                if 'TBLASTN' in all_sources:
                    circle_faces.append(CircleFace(5, 'green'))

                for idx, face in enumerate(circle_faces):
                    node.add_face(face, column=idx + 1, position='aligned')

    # Render
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.scale = 50
    ts.title.add_face(TextFace("Phylogenetic Distribution Summary", fsize=20), column=0)
    
    output_png = os.path.join(args.output_dir, "taxonomy_summary.png")
    tree.render(output_png, w=2000, units='mm', tree_style=ts)
    tree.write(format=1, outfile=os.path.join(args.output_dir, "taxonomy_summary.nwk"))
    print("[SUCCESS] Results saved")

if __name__ == "__main__":
    main()
