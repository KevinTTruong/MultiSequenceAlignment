gap_char = '-'
filename = ""
instructions = "Instructions:\n" \
                "* = mutation\n- = gap\n" \
                "-The top line indicates gaps/mutation for the OVERALL alignment, while the gaps/mutation\n" \
                "between sequences for each PAIR.\n" \
                "-For the \"Potential mutated\" section, the \"# seq\" indicates how many sequences contain\n" \
                "each INDIVIDUAL function (listed just below) and the #s following indicate possible mutation\n" \
                "locations (the numbers header should help for reference).\n" \
                "-For example, \"[func] was found in [X] sequences and could either be the [nucleotide] at\n" \
                "[location] or [another nucleotide] at [another location]\".\n" \

# List of sequences
seq = []

# key = functionality; value = # of seqs containing func (frequency)
mut = {}

# key = # of seqs containing func (frequency); value = func (inverse of mut)
mut_by_freq = {}

# key = frequency; value = possible mutation locations for frequency
loc_by_freq = {}


# Take in the name of the file and attempt to open file.
# If fail, would ask for file name again.
# Filters all the functions in the header and adds it the 'mut', while counting repeats.
# Adds enter sequences to 'seq'
# Finally, set 'mut_by_freq' by inversing 'mut'.
def start_io():
    while True:  # Keep trying to open file until opened
        try:
            global filename
            filename = input("Enter file name (including .txt): ")
            with open(filename, "r") as file:
                line = file.readline()
                while line:  # while string not empty
                    if line[0] == '>':
                        # Filter functions to list all mutations in mut and count how many times they appear
                        key = line[1:-1].split(", ")
                        for k in range(len(key)):
                            try:
                                mut[key[k]] += 1
                            except KeyError as ke:
                                mut[key[k]] = 1
                        # #
                        line = file.readline()
                    sequence = ""
                    while line and line[0] != '>':
                        sequence += line[:-1]   # Remove last char (\n)
                        line = file.readline()

                    seq.append(sequence)
                file.close()
                # Creates a reverse version of mut, now that all the frequency is counted.
                keys = list(mut.values())
                values = list(mut.keys())
                for x in range(len(mut)):
                    try:
                        mut_by_freq[keys[x]].append(values[x])
                    except KeyError as ke:
                        mut_by_freq[keys[x]] = []
                        mut_by_freq[keys[x]].append(values[x])
                        loc_by_freq[keys[x]] = []
                break
        except IOError as e:
            print("Could not find file.")


# Needleman-Wunsch algorithm
# Match: +1
# Mismatch: -2
# Gap: -2
def align(index1, index2):
    matrix = [[0 for x in range(len(seq[index2])+1)] for y in range(len(seq[index1])+1)]
    # Initialize first row/column
    matrix[0][0] = 0
    for i in range(1, len(seq[index1])+1):
        matrix[i][0] = -2 * i
    for j in range(1, len(seq[index2])+1):
        matrix[0][j] = -2 * j
    # Fill matrix
    for x in range(1, max(len(seq[index1])+1, len(seq[index2])+1)):
        if len(seq[index1]) > x-1 and len(seq[index2]) > x-1:
            for i in range(x, len(seq[index1])+1):
                calc = -2
                if seq[index1][i-1] == seq[index2][x-1]:
                    calc = 1
                matrix[i][x] = max(matrix[i-1][x]-2, matrix[i][x-1]-2, matrix[i-1][x-1]+calc)
            for j in range(x, len(seq[index2])+1):
                calc = -2
                if seq[index1][x-1] == seq[index2][j-1]:
                    calc = 1
                matrix[x][j] = max(matrix[x-1][j]-2, matrix[x][j-1]-2, matrix[x-1][j-1]+calc)

    # Traverse and locate gaps
    i = 0
    j = 0
    i_gaps = []
    j_gaps = []
    for x in range(max(len(seq[index1]), len(seq[index2]))):
        if i+1 >= len(seq[index1]) and j+1 < len(seq[index2]):
            i_gaps.append(x+1)  # Recall: gap in matrix for 0
            j += 1
        elif j+1 >= len(seq[index2]) and i+1 < len(seq[index1]):
            j_gaps.append(x+1)
            i += 1
        else:
            k = max(matrix[i+1][j+1], matrix[i][j+1], matrix[i+1][j])
            if k == matrix[i+1][j+1]:
                i += 1
                j += 1
            # Priority route 1
            elif k == matrix[i+1][j]:
                j_gaps.append(x+1)
                i += 1
            elif k == matrix[i][j+1]:
                i_gaps.append(x+1)   # Recall: gap in matrix for 0
                j += 1
            # Priority route 2
            '''
            elif k == matrix[i][j+1]:
                i_gaps.append(x+1)   # Recall: gap in matrix for 0
                j += 1
            elif k == matrix[i+1][j]:
                j_gaps.append(x+1)
                i += 1
            '''

    # Insert gaps into sequences and updates previous ones if needed
    length = range(max(len(seq[index1]), len(seq[index2])))
    if len(i_gaps) > 0:
        i_gap_ptr = len(i_gaps) - 1
        ptr = index1
        original_len = len(seq[ptr])
        for i in reversed(length):
            if i_gaps[i_gap_ptr] > original_len:
                seq[ptr] += gap_char
                i_gap_ptr -= 1
                if i_gap_ptr == -1:
                    break
            elif i == i_gaps[i_gap_ptr]:
                seq[ptr] = seq[ptr][:i] + gap_char + seq[ptr][i:]
                # Exclusive to i: insert gaps to previous sequences
                if ptr > 0:
                    for p in range(ptr):
                        seq[p] = seq[p][:i] + gap_char + seq[p][i:]
                # -----
                i_gap_ptr -= 1
                if i_gap_ptr == -1:
                    break
    # Same algorithm as above, but for the 2nd/latest sequence
    if len(j_gaps) > 0:
        j_gap_ptr = len(j_gaps) - 1
        ptr = index2
        original_len = len(seq[ptr])
        for j in reversed(length):
            if j_gaps[j_gap_ptr] > original_len:
                seq[ptr] += gap_char
                j_gap_ptr -= 1
                if j_gap_ptr == -1:
                    break
            elif j == j_gaps[j_gap_ptr]:
                seq[ptr] = seq[ptr][:j] + gap_char + seq[ptr][j:]
                j_gap_ptr -= 1
                if j_gap_ptr == -1:
                    break
    # Fills any leading gaps to match sequence lengths
    fill_remaining_gaps(index1, index2)


# Leaves results on screen until user presses enter.
def end_io():
    file = filename[:-4]+"_results.txt"
    input("\nResults are also saved to "+file+". Press ENTER to exit.")


# Print out mutation locations heading for the overall alignment (the ---***-*-- thing)
def mutation_pointer():
    output = ""
    for x in range(len(seq[0])):
        clean = True
        for y in range(1, len(seq)):
            if seq[0][x] != seq[y][x]:
                output += "*"
                clean = False
                break
        if clean:
            output += '-'
    return output


# Print mutation locations heading for a selected pair of sequences (the ---***-*-- thing)
def pairwise_mutation_pointer(index1, index2):
    output = ""
    for x in range(len(seq[index1])):
        if seq[index1][x] != seq[index2][x]:
            output += "*"
        else:
            output += '-'
    return output


# Calls align() for all sequences in 'seq'
def multi_align():
    if len(seq) == 1:
        print("Only one sequence available")
    elif len(seq) == 0:
        print("No sequences available")
    else:
        for x in range(1, len(seq)):
            align(x-1, x)


# Print results.
def print_results():
    output = ""
    # Prints sequence and its mutation location heading
    for x in range(len(seq)-1):
        output += seq[x] + "\n"
        output += pairwise_mutation_pointer(x, x+1) + "\n"
    output += seq[len(seq)-1]
    output += "\n"
    # Prints Potential mutated locations
    output += "\nPotential mutated locations for each function (sorted by number of sequences with function):\n"
    mut_by_freq_keys = list(mut_by_freq.keys())
    mut_by_freq_values = list(mut_by_freq.values())
    for x in range(len(mut_by_freq)):
        output += str(mut_by_freq_keys[x]) + " seqs: "
        # Print locations
        values = loc_by_freq[mut_by_freq_keys[x]]
        if len(values) != 0:
            for y in range(len(values) - 1):
                output += str(values[y]) + ", "
            output += str(values[len(values) - 1])
        else:
            output += "None found"
        # Print functionality
        output += "\n-Possible functions: "
        values = mut_by_freq_values[x]
        for y in range(len(values)-1):
            output += str(values[y]) + ", "
        output += str(values[len(values)-1]) + "\n"
    print(output)
    return output


# Output the results onto a txt called [filename]_result.txt
def output_results():
    trim = filename[:-4]    # Removes ".txt"
    file = open(trim+"_results.txt", "w")
    file.write(print_header()+"\n")
    file.write(mutation_pointer()+"\n")
    print(mutation_pointer())
    file.write(print_results())
    file.write("\n"+instructions)
    print(instructions)
    file.close()


# Fills any leading gaps to match sequence lengths
def fill_remaining_gaps(index1, index2):
    while len(seq[index1]) > len(seq[index2]):
        seq[index2] += gap_char
    while len(seq[index1]) < len(seq[index2]):
        seq[index1] += gap_char


# Prints numerical header (starts at 1 NOT 0)
# Header repeatedly counts from 0 to 9
# Header above counts the 10s from reference
def print_header():
    x = 1
    top = ""
    output = ""
    while x <= len(seq[0]):
        output += str(x % 10)
        if x % 10 == 0:
            top += str(int(x/10))
        else:
            top += " "
        x += 1
    top = top[1:]
    top += '\n'
    top += output

    print(top)
    return top

# Finds possible mutation locations for each functionality.
# For each mutation, count amount of sequences with a nucleotide
# that match the number of sequences that contain the functionality.
# Append to 'loc_by_freq' for each matching frequency.
def match_mutations():
    mut_loc = mutation_pointer()
    for loc in range(len(mut_loc)):
        if mut_loc[loc] == "*":
            nuc = {
                'A': [],
                'T': [],
                'G': [],
                'C': [],
                '-': []
            }
            for s in range(len(seq)):
                nuc[seq[s][loc]].append(loc)
            for n in range(len(nuc)):
                for m in range(len(loc_by_freq)):
                    if list(loc_by_freq.keys())[m] == len(list(nuc.values())[n]):
                        result = str(loc+1) + " (" + list(nuc.keys())[n] + ")"
                        loc_by_freq[list(loc_by_freq.keys())[m]].append(result)


# Used for debugging
def print_matrix(matrix):
    for x in range(len(matrix)):
        print(matrix[x])


# "Main" starts here
start_io()
multi_align()

match_mutations()
output_results()

end_io()

# Test cases
# seq.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCACCGCACCCCATTTTGTCCGCTGTATTATCCCCAATGAGTTTAAGCAATCGG")
# seq.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCACCGCACCCCATTTTGTCCGCTGTATTATCCCCAATGAGTTTAAGCAATCGG")
# seq.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCACAGCACTGCACCCCATTTTGTCCGCTGTATTGTGCCCAATGAGTTTAAGCAGTCAG")
# seq.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCACCGCACCCCATTTTGTCCGCTGTATTATCCCCAATGAGTTTAAGCAATCGG")
# seq.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCACAGCACTGCACCCCATTTTGTCCGCTGTATTGTGCCCAATGAGTTTAAGCAGTCAG")
# seq.append("GAGCAGCTGAACAAGCTGATGACCACCCTCCATAGCCGCACCCCATTTTGTCCGCTGTATTATCCCCAATGAGTTTAAGCAATCGG")
# seq.append("CGCAC")
# seq.append("GGCTTC")
# print(mutation_pointer())
# print(print_results())
