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
                atom = {
                    'record': line[0:6].strip(),
                    'atom_num': int(line[6:11]),
                    'atom_name': atom_name,
                    'res_name': line[17:20].strip(),
                    'res_renamed': "",
                    'chain_id': line[21],
                    'res_seq': int(line[22:26]),
                    'res_num': "",
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

def renumber_residues(atoms):
    """
    Renumber residues based solely on contiguous blocks of the same res_renamed.
    This ensures that all POX atoms (or any block) that are adjacent get the same final residue number.
    """
    if not atoms:
        return atoms

    current_residue_num = 1
    atoms[0]['res_new'] = current_residue_num
    for i in range(1, len(atoms)):
        # If the current atom has the same renamed residue as the previous, keep the same number.
        if atoms[i]['res_renamed'] == atoms[i-1]['res_renamed']:
            atoms[i]['res_new'] = current_residue_num
        else:
            current_residue_num += 1
            atoms[i]['res_new'] = current_residue_num
    return atoms

def remap_nucleotides(atoms, add_caps=True):
    """
    Remap the nucleotides into ROS, POX, UR0/AR0 blocks, with OH5 and OH3 caps.
    Then rearrange the atoms and update numbering.
    """
    # Work on a copy to avoid side effects.
    remapped_atoms = []
    for atom in atoms:
        # Convert OP1/OP2 to O1P/O2P
        if atom['atom_name'] == 'OP1':
            atom['atom_name'] = 'O1P'
        elif atom['atom_name'] == 'OP2':
            atom['atom_name'] = 'O2P'
        remapped_atoms.append(atom)

    # Remap residue names based on nucleotide and atom names.
    for atom in remapped_atoms:
        if atom['res_name'] == 'U':
            if atom['atom_name'] in ['C5\'', 'C4\'', 'C3\'', 'C2\'', 'C1\'', 'O4\'', 'O2\'']:
                atom['res_renamed'] = 'ROS'
            elif atom['atom_name'] in ['P', 'O1P', 'O2P', 'O5\'', 'O3\'', 'OP1', 'OP2']:
                atom['res_renamed'] = 'POX'
            else:
                atom['res_renamed'] = 'UR0'
        elif atom['res_name'] == 'A':
            if atom['atom_name'] in ['C5\'', 'C4\'', 'C3\'', 'C2\'', 'C1\'', 'O4\'', 'O2\'']:
                atom['res_renamed'] = 'ROS'
            elif atom['atom_name'] in ['P', 'O1P', 'O2P', 'O5\'', 'O3\'', 'OP1', 'OP2']:
                atom['res_renamed'] = 'POX'
            else:
                atom['res_renamed'] = 'AR0'
        elif atom['res_name'] == 'G':
            if atom['atom_name'] in ['C5\'', 'C4\'', 'C3\'', 'C2\'', 'C1\'', 'O4\'', 'O2\'']:
                atom['res_renamed'] = 'ROS'
            elif atom['atom_name'] in ['P', 'O1P', 'O2P', 'O5\'', 'O3\'', 'OP1', 'OP2']:
                atom['res_renamed'] = 'POX'
            else:
                atom['res_renamed'] = 'GR0'
        elif atom['res_name'] == 'C':
            if atom['atom_name'] in ['C5\'', 'C4\'', 'C3\'', 'C2\'', 'C1\'', 'O4\'', 'O2\'']:
                atom['res_renamed'] = 'ROS'
            elif atom['atom_name'] in ['P', 'O1P', 'O2P', 'O5\'', 'O3\'', 'OP1', 'OP2', 'OP3']:
                atom['res_renamed'] = 'POX'
            else:
                atom['res_renamed'] = 'CR0'

    # Identify and assign the cap atoms if add_caps is True.
    if add_caps:
        # First O5' becomes OH5.
        for atom in remapped_atoms:
            if atom['atom_name'] == "O5'":
                atom['res_renamed'] = 'OH5'
                break
        # Last O3' becomes OH3.
        for atom in reversed(remapped_atoms):
            if atom['atom_name'] == "O3'":
                atom['res_renamed'] = 'OH3'
                break

    # Renumber residues and update atom numbers.
    remapped_atoms = renumber_residues(remapped_atoms)
    remapped_atoms = update_atom_numbers(remapped_atoms)

    # Rearrangement steps:
    rearranged_atoms = remapped_atoms.copy()

    # Move O3' atoms (except the one capped as OH3) before the next POX residue.
    temp_atoms = rearranged_atoms.copy()
    for atom in temp_atoms:
        if atom['atom_name'] == "O3'" and atom['res_renamed'] != 'OH3':
            curr_idx = rearranged_atoms.index(atom)
            target_idx = None
            for later in rearranged_atoms[curr_idx+1:]:
                if later['res_renamed'] == 'POX':
                    target_idx = rearranged_atoms.index(later)
                    break
            if target_idx is not None:
                rearranged_atoms.remove(atom)
                rearranged_atoms.insert(target_idx, atom)
                rearranged_atoms = update_atom_numbers(rearranged_atoms)

    # Move every O5' (except the one capped as OH5) to after the previous POX residue.
    temp_atoms = rearranged_atoms.copy()
    for atom in temp_atoms:
        if atom['atom_name'] == "O5'" and atom['res_renamed'] != 'OH5':
            curr_idx = rearranged_atoms.index(atom)
            target_idx = None
            for i in range(curr_idx - 1, -1, -1):
                if rearranged_atoms[i]['res_renamed'] == 'POX':
                    target_idx = i + 1
                    break
            if target_idx is not None:
                rearranged_atoms.remove(atom)
                rearranged_atoms.insert(target_idx, atom)
                rearranged_atoms = update_atom_numbers(rearranged_atoms)

    # Move O2' atoms that follow AR0 or GR0 blocks to before the previous ROS residue.
    temp_atoms = rearranged_atoms.copy()
    for atom in temp_atoms:
        if atom['atom_name'] == "O2'" and rearranged_atoms.index(atom) > 0:
            curr_idx = rearranged_atoms.index(atom)
            prev_atom = rearranged_atoms[curr_idx - 1]
            if prev_atom['res_renamed'] in ('AR0', 'GR0'):
                target_idx = None
                for i in range(curr_idx - 2, -1, -1):
                    if rearranged_atoms[i]['res_renamed'] == 'ROS':
                        target_idx = i
                        break
                if target_idx is not None:
                    rearranged_atoms.remove(atom)
                    rearranged_atoms.insert(target_idx, atom)
                    rearranged_atoms = update_atom_numbers(rearranged_atoms)

    # Move every O1P/OP1 and O2P/OP2 atom before the previous POX residue.
    temp_atoms = rearranged_atoms.copy()
    for atom in temp_atoms:
        if atom['atom_name'] in ["O1P", "OP1", "O2P", "OP2"]:
            curr_idx = rearranged_atoms.index(atom)
            target_idx = None
            for i in range(curr_idx - 1, -1, -1):
                if rearranged_atoms[i]['res_renamed'] == 'POX':
                    target_idx = i + 1
                    break
            if target_idx is not None:
                rearranged_atoms.remove(atom)
                rearranged_atoms.insert(target_idx, atom)
                rearranged_atoms = update_atom_numbers(rearranged_atoms)

    return rearranged_atoms

def write_pdb(atoms, output_file):
    """Write atoms to a new PDB file."""
    with open(output_file, 'w') as file:
        for atom in atoms:
            file.write(f"{atom['record']:<6}{atom['atom_num']:>5} {atom['atom_name']:<4} {atom['res_renamed']:<4}{atom['chain_id']}{atom['res_new']:>4}    {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}  1.00  0.00          {atom['element']:>2}\n")

# Input and output file paths
input_file  = sys.argv[1]
output_file = sys.argv[2]

# Parse input PDB file
atoms = parse_pdb(input_file)

# Remap residues into ROS, POX, UR0, AR0 blocks with OH5 and OH3 caps and rearrange atoms.
remapped_atoms = remap_nucleotides(atoms)

# Renumber the residues (using our new renumbering that groups contiguous blocks)
remapped_atoms = renumber_residues(remapped_atoms)
remapped_atoms = update_atom_numbers(remapped_atoms)

# Write the output PDB file
write_pdb(remapped_atoms, output_file)

print(f"Remapping completed. Output saved to {output_file}.")
