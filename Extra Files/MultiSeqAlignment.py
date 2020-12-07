seq = []
gap_char = '-'
filename = ""
instructions = "* = mutation\n" \
               "The top line indicates gaps for the OVERALL alignment, " \
               "while the gaps between sequences indicates mutations for each PAIR.\n\n "


# IO
def start_io():
    while True:  # Keep trying to open file until opened
        try:
            global filename
            filename = input("Enter file name (including .txt): ")
            with open(filename, "r") as file:
                line = file.readline()
                while line:  # while string not empty
                    if line[0] == '>':
                        # Do stuff potentially
                        line = file.readline()
                    sequence = ""
                    while line and line[0] != '>':
                        sequence += line[:-1]   # Remove last char (\n)
                        line = file.readline()

                    seq.append(sequence)
                file.close()
                break
        except IOError as e:
            print("Could not find file.")


# Needleman-Wunsch algorithm
def align(index1, index2):
    matrix = [[0 for x in range(len(seq[index2])+1)] for y in range(len(seq[index1])+1)]
    # Initialize first row/column
    matrix[0][0] = 0
    for i in range(1, len(seq[index1])+1):
        matrix[i][0] = -2 * i
    for j in range(1, len(seq[index2])+1):
        matrix[0][j] = -2 * j
    # Fill matrix (WOO it works!)
    for x in range(1, max(len(seq[index1])+1, len(seq[index2])+1)):
        if len(seq[index1]) > x-1 and len(seq[index2]) > x-1:
            for i in range(x, len(seq[index1])+1):
                calc = -1
                if seq[index1][i-1] == seq[index2][x-1]:
                    calc = 1
                matrix[i][x] = max(matrix[i-1][x]-2, matrix[i][x-1]-2, matrix[i-1][x-1]+calc)
            for j in range(x, len(seq[index2])+1):
                calc = -1
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

    # Insert gaps into sequences
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
                # Exclusive to i
                if ptr > 0:
                    for p in range(ptr):
                        seq[p] = seq[p][:i] + gap_char + seq[p][i:]
                # -----
                i_gap_ptr -= 1
                if i_gap_ptr == -1:
                    break
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
    # Change bandaid to add remaining gaps to the shorter sequence (count from opposing gaps size)
    fill_remaining_gaps(index1, index2)


def end_io():
    file = filename[:-4]+"_results.txt"
    input("\n\nResults are also saved to "+file+". Press ENTER to exit.")


# Expose mutation locations
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


# Expose mutation locations for a selected pair of sequences
def pairwise_mutation_pointer(index1, index2):
    output = ""
    for x in range(len(seq[index1])):
        if seq[index1][x] != seq[index2][x]:
            output += "*"
        else:
            output += '-'
    return output


# Multi-align
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
    for x in range(len(seq)-1):
        output += seq[x] + "\n"
        output += pairwise_mutation_pointer(x, x+1) + "\n"
    output += seq[len(seq)-1]
    return output


# Output the results onto a txt called [filename]_result.txt
def output_results():
    trim = filename[:-4]    # Removes ".txt"
    file = open(trim+"_results.txt", "w")
    file.write(instructions)
    file.write(mutation_pointer()+"\n")
    file.write(print_results())
    file.close()


def fill_remaining_gaps(index1, index2):
    while len(seq[index1]) > len(seq[index2]):
        seq[index2] += gap_char
    while len(seq[index1]) < len(seq[index2]):
        seq[index1] += gap_char


# Used for debugging
def print_matrix(matrix):
    for x in range(len(matrix)):
        print(matrix[x])


# "Main" starts here
start_io()
multi_align()
output_results()
print()
print(mutation_pointer())
print(print_results())
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
