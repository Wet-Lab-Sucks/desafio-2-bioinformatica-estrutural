#!/usr/bin/env python 

import numpy as np
import urllib.request
import io

def parse_atom_lines(pdb_lines, residue_number, atom_types=None):
    """
    Extract coordinates for specific atoms in a given residue.
    
    Parameters:
    pdb_lines (list): Lines from a PDB file
    residue_number (int): Residue number to extract
    atom_types (list): List of atom types to extract (e.g., ['N', 'CA', 'C', 'O'])
    
    Returns:
    dict: Dictionary mapping atom types to their coordinates
    """
    if atom_types is None:
        atom_types = ['N', 'CA', 'C', 'O']  # Main chain atoms
    
    coords = {}
    
    for line in pdb_lines:
        if line.startswith("ATOM"):
            # Parse PDB ATOM line
            try:
                atom_name = line[12:16].strip()
                #print(atom_name)
                res_name = line[17:20].strip()
                chain = line[21:22].strip()
                res_num = int(line[22:26].strip())
                
                if res_num == residue_number and atom_name in atom_types:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords[atom_name] = np.array([x, y, z])
            except (ValueError, IndexError):
                continue
    
    return coords

def calculate_rmsd(coords1, coords2):
    """
    Calculate RMSD between two sets of coordinates.
    
    Parameters:
    coords1, coords2 (dict): Dictionaries mapping atom types to coordinates
    
    Returns:
    float: RMSD value
    """
    if not coords1 or not coords2:
        raise ValueError("Empty coordinate sets")
    
    # Ensure we're comparing the same atoms
    common_atoms = set(coords1.keys()) & set(coords2.keys())
    if not common_atoms:
        raise ValueError("No common atoms to compare")
    
    # Calculate squared differences
    squared_diffs = []
    for atom in common_atoms:
        diff = coords1[atom] - coords2[atom]
        squared_diff = np.sum(diff * diff)
        squared_diffs.append(squared_diff)
    
    # Calculate RMSD
    rmsd = np.sqrt(np.mean(squared_diffs))
    return rmsd

def download_pdb(pdb_id):
    """Download PDB file content from the RCSB PDB server"""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        with urllib.request.urlopen(url) as response:
            return response.read().decode('utf-8').splitlines()
    except Exception as e:
        print(f"Error downloading PDB {pdb_id}: {e}")
        return []

def calculate_gly55_rmsd(experimental_pdb_lines, model_pdb_lines):
    """Calculate RMSD for Gly55 between experimental and model structures"""
    # Extract coordinates for Gly55 from both structures
    exp_coords = parse_atom_lines(experimental_pdb_lines, 55)
    print(f"Resíduo Gly55 do PDB ID 2V84.pdb: {exp_coords}")
    model_coords = parse_atom_lines(model_pdb_lines, 55)
    print(f"Resíduo Gly55 do PotD_Modelada_LBB.pdb: {model_coords}")

    # Check if we found all required atoms
    required_atoms = ['N', 'CA', 'C', 'O']
    for atom in required_atoms:
        if atom not in exp_coords:
            print(f"Warning: {atom} atom not found in experimental structure for Gly55")
        if atom not in model_coords:
            print(f"Warning: {atom} atom not found in model structure for Gly55")
    
    # Calculate RMSD
    rmsd = calculate_rmsd(exp_coords, model_coords)
    return rmsd

# Main execution
# Download experimental structure
experimental_pdb_lines = download_pdb("2V84")
#print(experimental_pdb_lines)

with open("dataset/PotD_Modelada_LBB.pdb", "r") as f:
    model_pdb_lines = f.readlines()

# compare both structures
rmsd = calculate_gly55_rmsd(experimental_pdb_lines, model_pdb_lines)
print(f"RMSD for Gly55: {rmsd:.2f} Å")

# output: 31.09
# seems ok.

# as we discussed, this will be delivered as a single output file, e.g. RMSD_residuo_Gly55.txt 
