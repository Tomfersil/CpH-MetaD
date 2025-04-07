import sys

def parse_pdb(file_path):
    """Parse the PDB file and return a list of atoms."""
    atoms = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(('ATOM', 'HETATM')):
                atom = {
                    'record': line[0:6].strip(),
                    'atom_num': int(line[6:11]),
                    'atom_name': line[12:16].strip(),
                    'res_renamed': line[17:20].strip(),
                    'res_name': "",
                    'chain_id': line[21],
                    'res_new': int(line[22:26]),
                    'x': float(line[30:38]),
                    'y': float(line[38:46]),
                    'z': float(line[46:54]),
                    'element': line[76:78].strip(),
                    'line': line.strip()  # Save the original line for output
                }
                atoms.append(atom)
    return atoms

def reverse_remap(atoms):
    """Reverse the residue renaming and atom rearrangements."""
    # Step 1: Rename GR0/GR1, AR0/AR1, CR0/CR1, UR0/UR1 back to G, A, C, U
    for atom in atoms:
        if atom['res_renamed'] in ['GR0', 'GR1']:
            atom['res_name'] = 'G'
        elif atom['res_renamed'] in ['AR0', 'AR1']:
            atom['res_name'] = 'A'
        elif atom['res_renamed'] in ['CR0', 'CR1']:
            atom['res_name'] = 'C'
        elif atom['res_renamed'] in ['UR0', 'UR1']:
            atom['res_name'] = 'U'

    # Step 2: Apply rules for OH5, OH3, ROS, and POX
    for i, atom in enumerate(atoms):
        if atom['res_renamed'] == 'OH5':
            # OH5 should belong to the same residue as the first A/C/G/U that appears next
            for next_atom in atoms[i + 1:]:
                if next_atom['res_name'] in ['A', 'C', 'G', 'U']:
                    atom['res_name'] = next_atom['res_name']
                    break
        elif atom['res_renamed'] == 'OH3':
            # OH3 should belong to the same residue as the last A/C/G/U that appears before
            for prev_atom in reversed(atoms[:i]):
                if prev_atom['res_name'] in ['A', 'C', 'G', 'U']:
                    atom['res_name'] = prev_atom['res_name']
                    break
        elif atom['res_renamed'] == 'ROS':
            # ROS should belong to the same residue as the first A/C/G/U that appears next
            for next_atom in atoms[i + 1:]:
                if next_atom['res_name'] in ['A', 'C', 'G', 'U']:
                    atom['res_name'] = next_atom['res_name']
                    break
        elif atom['res_renamed'] == 'POX':
            if atom['atom_name'] == "O3'":
                # O3' should belong to the previous A/C/G/U residue
                for prev_atom in reversed(atoms[:i]):
                    if prev_atom['res_name'] in ['A', 'C', 'G', 'U']:
                        atom['res_name'] = prev_atom['res_name']
                        break
            else:
                # Remaining POX atoms should belong to the next A/C/G/U residue
                for next_atom in atoms[i + 1:]:
                    if next_atom['res_name'] in ['A', 'C', 'G', 'U']:
                        atom['res_name'] = next_atom['res_name']
                        break

        # Reverse atom name changes
        if atom['atom_name'] == 'O1P':
            atom['atom_name'] = 'OP1'
        elif atom['atom_name'] == 'O2P':
            atom['atom_name'] = 'OP2'

    # Sort atoms by their original atom numbers to reverse rearrangements
    atoms.sort(key=lambda x: x['atom_num'])
    return atoms

def correct_residue_numbers(atoms):
    """Update residue numbers sequentially based on the residue name and the appearance of a new P atom."""
    current_residue_number = 0
    previous_residue_name = None
    previous_p_atom = False

    for atom in atoms:
        # Increment residue number if the residue name changes or a new O5' atom appears
        if atom['res_name'] != previous_residue_name or atom['atom_name'] == 'O5\'':
            current_residue_number += 1
            previous_residue_name = atom['res_name']
            previous_p_atom = atom['atom_name'] == 'P'

        atom['res_new'] = current_residue_number

    return atoms

def write_pdb(atoms, output_file):
    """Write atoms to a new PDB file."""
    with open(output_file, 'w') as file:
        for atom in atoms:
            file.write(f"{atom['record']:<6}{atom['atom_num']:>5} {atom['atom_name']:<4} {atom['res_name']:<4}{atom['chain_id']}{atom['res_new']:>4}    {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}  1.00  0.00          {atom['element']:>2}\n")

# Input and output file paths
input_file = sys.argv[1]
output_file = sys.argv[2]

# Parse the modified PDB file
atoms = parse_pdb(input_file)

# Reverse the remapping and rearrangements
restored_atoms = reverse_remap(atoms)

# Correct the residue numbers
corrected_atoms = correct_residue_numbers(restored_atoms)

# Write the restored PDB file
write_pdb(corrected_atoms, output_file)

print(f"Restoration completed. Output saved to {output_file}.")