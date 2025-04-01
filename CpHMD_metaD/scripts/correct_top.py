#!/bin/python3
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]
strings = ['Protein','; [ begin_tit_molecule ]',';[ end_tit_molecule ]']




with open(input_file, 'r') as file:
    lines = file.readlines()

    output_lines = []
    begins = False
    ends   = False

    for line in lines:

        if not begins and line.startswith('#include') and 'forcefield.itp' in line:
            begins = True
            output_lines.append(line)
            continue

        if begins:
            output_lines.append('\n' + strings[1] + '\n')
            begins = False
            
        if 'Other' in line:
            output_lines.append(line.replace('Other',strings[0]))
        elif 'RNA' in line:
            output_lines.append(line.replace('RNA',strings[0]))

        else:
            output_lines.append(line)


    fstring = output_lines.index('; Include Position restraint file\n')
    sstring = output_lines.index('#endif\n')
    tstring = output_lines.index('; Include water topology\n')

    if sstring > fstring and sstring < tstring:
        output_lines.insert(sstring+1,'\n' + strings[2] + '\n')

with open(output_file, 'w') as file:
    file.writelines(output_lines)