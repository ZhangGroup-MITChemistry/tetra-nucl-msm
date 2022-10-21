import os
import tools
import sys

nrl = int(sys.argv[1])
num_copy = int(sys.argv[2])

# continued
# plug in histones
# now let's collect all the histones into one file
output_file = open('histone-all.pdb', 'w')
input_list = [[] for i in range(num_copy)]

# Read each histone file
for i in range(1, num_copy + 1):
    histone_file = open('histone-%d.pdb' % (i), 'r')
    histone_file_lines = histone_file.readlines()
    for each_line in histone_file_lines:
        if each_line[0:4] == 'ATOM' and (each_line[17:20] not in ['HOH', 'CL ', 'MN ']) and (each_line[13:16] != 'OXT'):
            input_list[i - 1].append(each_line)
    histone_file.close()
    os.system('rm -f histone-%d.pdb' % (i))

# write the output
for i in range(0, num_copy):
    N = len(input_list[i])
    for j in range(N):
        output_file.write(input_list[i][j])
        if j < N - 1:
            if input_list[i][j][21] != input_list[i][j + 1][21]:
                output_file.write('TER\n')
    if i == num_copy - 1:
        output_file.write('END\n')
    else:
        output_file.write('TER\n')
output_file.close()
