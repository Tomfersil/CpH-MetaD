import re
import sys 

def parse_pdb(file_path):
    """Parse the PDB file and return a list of atoms."""
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:

            
            if line.startswith(('ATOM', 'HETATM')):
                atom_name = line[12:16].strip()
                if atom_name.startswith('H'):
                    continue
                else:
                    atom = {
                        'record': line[0:6].strip(),
                        'atom_num': int(line[6:11]),
                        'atom_name': line[12:16].strip(),
                        'res_name': line[17:20].strip(),
                        'chain_id': line[21],
                        'res_seq': int(line[22:26]),
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54]),
                        'element': line[76:78].strip(),
                        'line': line.strip()  # Save the original line for output
                    }
                    atoms.append(atom)

    return atoms

def update_atom_numbers(atoms):
    """Update atom numbers sequentially."""
    for idx, atom in enumerate(atoms):
        atom['atom_num'] = idx + 1
    return atoms

def remap_nucleotides(atoms, add_caps=True):
    """Remap the nucleotides into ROS, POX, UR0/AR0 blocks, with OH5 and OH3 caps, and renumber the res_seq."""
    remapped_atoms = []
    
    # Identifying the first O5' and last O3'
    first_O5_idx = None
    last_O3_idx  = None
    o3_p_pairs = {'O3':[],'P':[]}  



    for idx, atom in enumerate(atoms):
        if atom['atom_name'] == "O5'":
            if first_O5_idx is None and add_caps == True:
                first_O5_idx = idx
        if atom['atom_name'] == "O3'" and add_caps == True:
            last_O3_idx = idx
            for j in range(idx+1, len(atoms)):
                if atoms[j]['atom_name'] == "P":
                    o3_p_pairs['O3'].append(idx)
                    o3_p_pairs['P'].append(j)
                    break
    for idx, atom in enumerate(atoms):
        # Special cases for OH5 and OH3 caps
        if idx == first_O5_idx and add_caps is True:
            atom['res_name'] = 'OH5'
        elif idx == last_O3_idx and add_caps is True:
            atom['res_name'] = 'OH3'
        else:
            # Remap other residues based on the nucleotide name and atom name
            if atom['res_name'] == 'U':
                if atom['atom_name'] in ['C5\'', 'C4\'', 'C3\'', 'C2\'', 'C1\'', 'O4\'', 'O2\'']:
                    atom['res_name'] = 'ROS'
                elif atom['atom_name'] in ['P', 'O1P', 'O2P', 'O5\'', 'O3\'','OP1','OP2']:
                    atom['res_name'] = 'POX'
                else:
                    atom['res_name'] = 'UR0'
            
            elif atom['res_name'] == 'A':
                if atom['atom_name'] in ['C5\'', 'C4\'', 'C3\'', 'C2\'', 'C1\'', 'O4\'', 'O2\'']:
                    atom['res_name'] = 'ROS'
                elif atom['atom_name'] in ['P', 'O1P', 'O2P', 'O5\'', 'O3\'','OP1','OP2']:
                    atom['res_name'] = 'POX'
                else:
                    atom['res_name'] = 'AR0'
            
            elif atom['res_name'] == 'G':
                if atom['atom_name'] in ['C5\'', 'C4\'', 'C3\'', 'C2\'', 'C1\'', 'O4\'', 'O2\'']:
                    atom['res_name'] = 'ROS'
                elif atom['atom_name'] in ['P', 'O1P', 'O2P', 'O5\'', 'O3\'','OP1','OP2']:
                    atom['res_name'] = 'POX'
                else:
                    atom['res_name'] = 'GR0'
            
            elif atom['res_name'] == 'C':
                if atom['atom_name'] in ['C5\'', 'C4\'', 'C3\'', 'C2\'', 'C1\'', 'O4\'', 'O2\'']:
                    atom['res_name'] = 'ROS'
                elif atom['atom_name'] in ['P', 'O1P', 'O2P', 'O5\'', 'O3\'','OP1','OP2','OP3']:
                    atom['res_name'] = 'POX'
                else:
                    atom['res_name'] = 'CR0'
        remapped_atoms.append(atom)
        
            # Increment residue sequence number after processing each atom
            #current_residue_num += 1
    
    # Move the first O5' atom to the front of the list
    if first_O5_idx is not None:
        o5_atom = remapped_atoms.pop(first_O5_idx)
        remapped_atoms.insert(0, o5_atom)
    # Move every O3' before the following P atom, except the last O3'
    # Move the first O5' atom to the front of the list
    if last_O3_idx is not None:
        o3_atom = remapped_atoms.pop(last_O3_idx)
        remapped_atoms.append(o3_atom)
        # Move every O3' before the following P atom, except the last O3'
    
    # Update the atom numbering after remapping
    remapped_atoms = update_atom_numbers(remapped_atoms)

    rearranged_atoms = []
    for idx, atom in enumerate(remapped_atoms):
        rearranged_atoms.append(atom)  # Start by adding the atom normally
        # Find the next P atom and place the O3' before it
    
    # Move O3' (except the last one) to the next POX residue
    for idx, atom in enumerate(rearranged_atoms):
        if atom['atom_name'] == "O3'" and idx != last_O3_idx:
            if atom['atom_num']-1 in o3_p_pairs['O3']:
                O3_idx = o3_p_pairs['O3'].index(atom['atom_num']-1)
                P_idx  = o3_p_pairs['P'][O3_idx]

                o3_atom = rearranged_atoms.pop(idx)
                rearranged_atoms.insert(P_idx, o3_atom)
    
    # Update the atom numbering after moving atoms around
    rearranged_atoms = update_atom_numbers(rearranged_atoms)

    # Move every O5' (except the first one) to the previous POX residue
    for idx, atom in enumerate(rearranged_atoms):

        if atom['atom_name'] == "O5'" and idx != first_O5_idx:
            # Find the previous POX atom and move the O5' before it
            for i in range(idx - 1, -1, -1):
                if remapped_atoms[i]['res_name'] == 'POX':
                    # Move the O5' before the found POX atom
                    o5_atom = rearranged_atoms.pop(idx)
                    rearranged_atoms.insert(i+1, o5_atom)
                    # Update the atom numbering after moving atoms around
                    rearranged_atoms = update_atom_numbers(rearranged_atoms)
                    break


    # Move O2' following AR0 blocks to correct the ROS blocks
    for idx, atom in enumerate(rearranged_atoms):
        if atom['atom_name'] == "O2'" and idx > 0:
            # Check if the previous residue is AR0
            if atoms[idx - 1]['res_name'] == 'AR0' or atoms[idx - 1]['res_name'] == 'GR0':
                #print(idx)

                # Find the previous ROS residue
                for i in range(idx - 2, -1, -1):  # Searching for ROS block before AR0
                    if atoms[i]['res_name'] == 'ROS':
                        #print(i)
                        # Move O2' before the ROS atom
                        o2_atom = rearranged_atoms.pop(idx)
                        rearranged_atoms.insert(i, o2_atom)
                         # Update the atom numbering after moving atoms around
                        rearranged_atoms = update_atom_numbers(rearranged_atoms)
                        break



    # Move every O1P/OP1 and O2P/OP2 before the previous POX residue
    for idx, atom in enumerate(rearranged_atoms):
        if atom['atom_name'] in ["O1P", "OP1", "O2P", "OP2"]:
            # Find the previous POX atom
            for k in range(idx - 1, -1, -1):
                if rearranged_atoms[k]['res_name'] == 'POX':
                    #print(idx,k,atom['atom_name'])
                    # Move the current atom before the found POX atom
                    phosphate_atom = rearranged_atoms.pop(idx)
                    rearranged_atoms.insert(k + 1, phosphate_atom)
                    # Update the atom numbering after moving atoms around
                    rearranged_atoms = update_atom_numbers(rearranged_atoms)
                    break

    
    return rearranged_atoms

def renumber_residues(atoms):
    """Renumber the res_seq based on residue name changes."""
    current_residue_num = 0  # Start renumbering from 0
    previous_residue_name = None
    remapped_atoms = []
    
    for atom in atoms:
        # If the residue name changes, increment the residue number
        if atom['res_name'] != previous_residue_name:
            previous_residue_name = atom['res_name']
            current_residue_num += 1

        # Assign the current residue number
        atom['res_seq'] = current_residue_num
        remapped_atoms.append(atom)
    
    return remapped_atoms

def write_pdb(atoms, output_file):
    """Write atoms to a new PDB file."""
    with open(output_file, 'w') as file:
        for atom in atoms:
            file.write(f"{atom['record']:<6}{atom['atom_num']:>5} {atom['atom_name']:<4} {atom['res_name']:<4}{atom['chain_id']}{atom['res_seq']:>4}    {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}  1.00  0.00          {atom['element']:>2}\n")

# Input and output file paths
input_file = 'WW_bevilacqua_OH5.pdb'
output_file = 'WW_cph.pdb'

#input_file  = 'A1mer_nomenclature_pymol.pdb'
#output_file = 'converted_A3mer.pdb'

# Parse input PDB file
atoms = parse_pdb(input_file)

# Remap residues into ROS, POX, UR0, AR0 blocks with OH5 and OH3 caps
remapped_atoms = remap_nucleotides(atoms)

# Renumber the residues
remapped_atoms = renumber_residues(remapped_atoms)

# Write the output PDB
write_pdb(remapped_atoms, output_file)

print(f"Remapping completed. Output saved to {output_file}.")
