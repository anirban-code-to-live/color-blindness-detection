
# coding: utf-8

# In[1]:

import numpy as np
import pickle
from bitarray import bitarray


# In[2]:

def generate_first_column(bwt_data):
    first_col_dict = {"A": 0, "G": 0, "C": 0, "T": 0, "$": 0, "N": 0, "U": 0}
    counts = [0, 0, 0, 0, 0, 0, 0]
    for line in bwt_data:
        for token in line:
            if token == "A":
                counts[0] += 1
            elif token == "C":
                counts[1] += 1
            elif token == "T":
                counts[2] += 1
            elif token == "G":
                counts[3] += 1
            elif token == "$":
                counts[4] += 1
            elif token == "N":
                counts[5] += 1
            else:
                counts[6] += 1

    first_col_dict["A"] = counts[0]
    first_col_dict["C"] = counts[1]
    first_col_dict["T"] = counts[2]
    first_col_dict["G"] = counts[3]
    first_col_dict["$"] = counts[4]
    first_col_dict["N"] = counts[5]
    first_col_dict["U"] = counts[6]
#     print(counts)
#     print(first_col_dict)
    return first_col_dict


# In[3]:

def generate_binary_and_count_arrays(bwt_full, save=True):
    bit_arr_length = len(bwt_full)
    print(bit_arr_length)
    G_bit_arr = bitarray(bit_arr_length)
    A_bit_arr = bitarray(bit_arr_length)
    C_bit_arr = bitarray(bit_arr_length)
    T_bit_arr = bitarray(bit_arr_length)
    
    G_bit_arr.setall(False)
    A_bit_arr.setall(False)
    C_bit_arr.setall(False)
    T_bit_arr.setall(False)
    
    delta = 100
    count_arr_length = bit_arr_length/delta + 1
    print(count_arr_length)
    
    G_count_arr = []
    A_count_arr = []
    C_count_arr = []
    T_count_arr = []
            
    g_count = 0
    c_count = 0
    a_count = 0
    t_count = 0
    for i in range(len(bwt_full)):
        token = bwt_full[i]
        if token == "A":
            A_bit_arr[i] = True
            a_count += 1
        elif token == "C":
            C_bit_arr[i] = True
            c_count += 1
        elif token == "T":
            T_bit_arr[i] = True
            t_count += 1
        elif token == "G":
            G_bit_arr[i] = True
            g_count += 1
        else:
            print('No match')

        if i % delta == 0:
            G_count_arr.append(g_count)
            A_count_arr.append(a_count)
            T_count_arr.append(t_count)
            C_count_arr.append(c_count)
            
    if save == True:       
        g_bit_arr_filename = '../processed_data/g_bit_arr.pkl'
        c_bit_arr_filename = '../processed_data/c_bit_arr.pkl'
        a_bit_arr_filename = '../processed_data/a_bit_arr.pkl'
        t_bit_arr_filename = '../processed_data/t_bit_arr.pkl'
        with open(g_bit_arr_filename, 'wb') as g:
            pickle.dump(G_bit_arr, g)
        with open(c_bit_arr_filename, 'wb') as c:
            pickle.dump(C_bit_arr, c)
        with open(a_bit_arr_filename, 'wb') as a:
            pickle.dump(A_bit_arr, a)
        with open(t_bit_arr_filename, 'wb') as t:
            pickle.dump(T_bit_arr, t)

        g_count_arr_filename = '../processed_data/g_count_arr.pkl'
        c_count_arr_filename = '../processed_data/c_count_arr.pkl'
        a_count_arr_filename = '../processed_data/a_count_arr.pkl'
        t_count_arr_filename = '../processed_data/t_count_arr.pkl'
        with open(g_count_arr_filename, 'wb') as g1:
            pickle.dump(G_count_arr, g1)
        with open(c_count_arr_filename, 'wb') as c1:
            pickle.dump(C_count_arr, c1)
        with open(a_count_arr_filename, 'wb') as a1:
            pickle.dump(A_count_arr, a1)
        with open(t_count_arr_filename, 'wb') as t1:
            pickle.dump(T_count_arr, t1)
    
    return G_bit_arr, A_bit_arr, C_bit_arr, T_bit_arr, G_count_arr, A_count_arr, C_count_arr, T_count_arr


# In[4]:

def rank_query(token, index):
    if index == 0:
        return 0
    rank = 0
    if(token == 'G'):
        while (index - 1) % delta != 0:
            index -= 1
            if G_bit_arr[index]:
                rank += 1
        rank += G_count_arr[(index-1) / delta]
        return rank
    
    elif(token == 'A'):
        while (index - 1) % delta != 0:
            index -= 1
            if A_bit_arr[index]:
                rank += 1
        rank += A_count_arr[(index-1) / delta]
        return rank
    
    elif(token == 'C'):
        while (index - 1) % delta != 0:
            index -= 1
            if C_bit_arr[index]:
                rank += 1
        rank += C_count_arr[(index-1) / delta]
        return rank
                
    elif(token == 'T'):
        while (index - 1) % delta != 0:
            index -= 1
            if T_bit_arr[index]:
                rank += 1
        rank += T_count_arr[(index-1) / delta]
        return rank
                
    else:
        print(token)
        print("Token mismatch!")


# In[5]:

def get_band_index(token, rank):
    if token == 'A':
        return band_start_index_a + rank
    elif token == 'C':
        return band_start_index_c + rank
    elif token == 'G':
        return band_start_index_g + rank
    elif token == 'T':
        return band_start_index_t + rank
    else:
        print("No match for token!!")


# In[6]:

red_exon_1_start_index = 149249757
red_exon_1_end_index = 149249868
red_exon_2_start_index = 149256127
red_exon_2_end_index = 149256423
red_exon_3_start_index = 149258412
red_exon_3_end_index = 149258580
red_exon_4_start_index = 149260048
red_exon_4_end_index = 149260213
red_exon_5_start_index = 149261768
red_exon_5_end_index = 149262007
red_exon_6_start_index = 149264290
red_exon_6_end_index = 149264400

green_exon_1_start_index = 149288166
green_exon_1_end_index = 149288277
green_exon_2_start_index = 149293258
green_exon_2_end_index = 149293554
green_exon_3_start_index = 149295542
green_exon_3_end_index = 149295710
green_exon_4_start_index = 149297178
green_exon_4_end_index = 149297343
green_exon_5_start_index = 149298898
green_exon_5_end_index = 149299137
green_exon_6_start_index = 149301420
green_exon_6_end_index = 149301530

def is_matching_red(index, length):
    if index >= red_exon_1_start_index and index <= red_exon_1_end_index:
        return True, 0
    elif index >= red_exon_2_start_index and index <= red_exon_2_end_index:
        return True, 1
    elif index >= red_exon_3_start_index and index <= red_exon_3_end_index:
        return True, 2
    elif index >= red_exon_4_start_index and index <= red_exon_4_end_index:
        return True, 3
    elif index >= red_exon_5_start_index and index <= red_exon_5_end_index:
        return True, 4
    elif index >= red_exon_6_start_index and index <= red_exon_6_end_index:
        return True, 5
    else:
        return False, -1
    
def is_matching_green(index, length):
    if index >= green_exon_1_start_index and index <= green_exon_1_end_index:
        return True, 0
    elif index >= green_exon_2_start_index and index <= green_exon_2_end_index:
        return True, 1
    elif index >= green_exon_3_start_index and index <= green_exon_3_end_index:
        return True, 2
    elif index >= green_exon_4_start_index and index <= green_exon_4_end_index:
        return True, 3
    elif index >= green_exon_5_start_index and index <= green_exon_5_end_index:
        return True, 4
    elif index >= green_exon_6_start_index and index <= green_exon_6_end_index:
        return True, 5
    else:
        return False, -1


# In[7]:

bwt_data = np.loadtxt('../data/chrX_last_col.txt', dtype=str)
# print(bwt_data)
# print(len(bwt_data))


# In[ ]:

reference_read = np.loadtxt("../data/chrX_map.txt", dtype=str)
print("Length of the reference sequence :: " + str(len(reference_read)))


# In[9]:

ref_sequence = np.loadtxt("../data/chrX.fa", dtype=str)
# print(len(ref_sequence))
merged_sequence = ''.join(ref_sequence[1:])
# print(len(merged_sequence))


# In[10]:

bwt_full = ''.join(bwt_data)
# print(len(bwt_full))


# In[11]:

# test_read = 'GAGGACAGCACCCAGTCCAGCATCTTCACCTACACCAACAGCAACTCCACCAGAGGTGAGCCAGCAGGCCCGTGGAGGCTGGGTGGCTGCACTGGGGGCCA'
# print(len(test_read))
# delta = 100


# In[12]:

first_col_dict = generate_first_column(bwt_data)
a_count = first_col_dict['A']
g_count = first_col_dict['G']
c_count = first_col_dict['C']
t_count = first_col_dict['T']

band_start_index_a = 0
band_end_index_a = a_count - 1
band_start_index_c = a_count
band_end_index_c = a_count + c_count - 1
band_start_index_g = a_count + c_count
band_end_index_g = a_count + c_count + g_count - 1
band_start_index_t = a_count + c_count + g_count
band_end_index_t = a_count + c_count + g_count + t_count - 1

# print(band_start_index_a, band_end_index_a)
# print(band_start_index_c, band_end_index_c)
# print(band_start_index_g, band_end_index_g)
# print(band_start_index_t, band_end_index_t)

# first_token = test_read[-1]
# print(first_token)

# if first_token == 'A':
#     band_start_index = 0
#     band_end_index = a_count - 1
# elif first_token == 'C':
#     band_start_index = a_count
#     band_end_index = a_count + c_count - 1
# elif first_token == 'G':
#     band_start_index = a_count + c_count
#     band_end_index = a_count + c_count + g_count - 1
# elif first_token == 'T':
#     band_start_index = a_count + c_count + g_count
#     band_end_index = a_count + c_count + g_count + t_count - 1

# print(band_start_index)
# print(band_end_index)

G_bit_arr, A_bit_arr, C_bit_arr, T_bit_arr, G_count_arr, A_count_arr, C_count_arr, T_count_arr = generate_binary_and_count_arrays(bwt_full, False)


# In[ ]:

# for i in range(2, len(test_read)+1):
#     token = test_read[-1*i]
#     start_token_rank = rank_query(token, band_start_index)
#     end_token_rank = rank_query(token, band_end_index)
#     if start_token_rank == end_token_rank:
#         break;
#     band_start_index = get_band_index(token, start_token_rank)
#     band_end_index = get_band_index(token, end_token_rank)
# #     print(start_token_rank)
# #     print(end_token_rank)

# print(band_start_index)
# print(band_end_index)


# In[ ]:

# print(reference_read[band_end_index-2])
# start_seq_index = reference_read[band_end_index-2]


# In[ ]:

# read_from_seq = merged_sequence[int(start_seq_index):int(start_seq_index)+101]
# print(merged_sequence[int(start_seq_index):int(start_seq_index)+101])


# In[ ]:

# assert test_read == read_from_seq


# In[ ]:

## DOnt go beyond this


# In[13]:

red_exon_match_count = [0, 0, 0, 0, 0, 0]
green_exon_match_count = [0, 0, 0, 0, 0, 0]


# In[14]:

all_reads = np.loadtxt('../data/reads', dtype=str)
# print(len(all_reads))
# print(len(all_reads[0]))
max_length_read = max(len(read) for read in all_reads)
print("Maximum length read" + str(max_length_read))
min_length_read = min(len(read) for read in all_reads)
print("Minimum length read" + str(min_length_read))
count_read_containing_N = len([read for read in all_reads if 'N' in read])
print("Count of reads containing N" + str(count_read_containing_N))

all_reads_removing_N = []
for read in all_reads:
    if 'N' in read:
        new_read = read.replace('N', 'A')
        all_reads_removing_N.append(new_read)
    else:
        all_reads_removing_N.append(read)

new_count_read_containing_N = len([read for read in all_reads_removing_N if 'N' in read])
# print(new_count_read_containing_N)

# with open('../processed_data/processed_reads.pkl', 'wb') as f:
#     pickle.dump(all_reads_removing_N, f)
    
all_reads_reverse_compliment = []
for read in all_reads_removing_N:
    reverse_read = read[::-1]
    temp_M = reverse_read.replace('G', 'M')
    temp_G = temp_M.replace('C', 'G')
    temp_C = temp_G.replace('M', 'C')
    temp_N = temp_C.replace('A', 'N')
    temp_T = temp_N.replace('T', 'A')
    temp_A = temp_T.replace('N', 'T')
    all_reads_reverse_compliment.append(temp_A)
    
# with open('../processed_data/processed_reverse_compliment_reads.pkl', 'wb') as f:
#     pickle.dump(all_reads_reverse_compliment, f)


# In[15]:

def get_match_indices(read):
    first_token = read[-1]
    if first_token == 'A':
        band_start_index = 0
        band_end_index = a_count - 1
    elif first_token == 'C':
        band_start_index = a_count
        band_end_index = a_count + c_count - 1
    elif first_token == 'G':
        band_start_index = a_count + c_count
        band_end_index = a_count + c_count + g_count - 1
    elif first_token == 'T':
        band_start_index = a_count + c_count + g_count
        band_end_index = a_count + c_count + g_count + t_count - 1
    
    for i in range(2, len(read)+1):
        token = read[-1*i]
        start_token_rank = rank_query(token, band_start_index)
        end_token_rank = rank_query(token, band_end_index)
        if start_token_rank == end_token_rank:
            return None, None;
        band_start_index = get_band_index(token, start_token_rank)
        band_end_index = get_band_index(token, end_token_rank)
    return band_start_index, band_end_index


# In[16]:

def count_exon_matches(start_index, end_index, read_length):
    ref_seq_index_matches = []
    for i in range(end_index - start_index):
        ref_seq_index_matches.append(int(reference_read[start_index+i]))
    red_exon_numbers = []
    green_exon_numbers = []
    for ref_seq_index in ref_seq_index_matches:
        is_red_match, red_exon_number = is_matching_red(ref_seq_index, read_length)
        is_green_match, green_exon_number = is_matching_green(ref_seq_index, read_length)
        if is_red_match:
            red_exon_numbers.append(red_exon_number)
        if is_green_match:
            green_exon_numbers.append(green_exon_number)
            
    if len(red_exon_numbers) > 0 and len(green_exon_numbers) > 0:
        for exon_number in red_exon_numbers:
            red_exon_match_count[exon_number] += 0.5
        for exon_number in green_exon_numbers:
            green_exon_match_count[exon_number] += 0.5
    elif len(red_exon_numbers) > 0:
        for exon_number in red_exon_numbers:
            red_exon_match_count[exon_number] += 1
    elif len(green_exon_numbers) > 0:
        for exon_number in green_exon_numbers:
            green_exon_match_count[exon_number] += 0.5
    else:
        pass


# In[18]:

delta = 100
for read in all_reads_removing_N:
    band_first_match_index, band_last_match_index = get_match_indices(read)
    if band_first_match_index is not None and band_last_match_index is not None:
#         print(band_first_match_index)
#         print(band_last_match_index)
        count_exon_matches(band_first_match_index, band_last_match_index, len(read))


# In[19]:

print("Exact match counts for red gene")
print(red_exon_match_count)
print("Exact match count for green gene")
print(green_exon_match_count)


# In[ ]:

# for read in all_reads_reverse_compliment:
#     band_first_match_index, band_last_match_index = get_match_indices(read)
#     if band_first_match_index is not None and band_last_match_index is not None:
#         count_exon_matches(band_first_match_index, band_last_match_index, len(read))


# In[ ]:

# print(red_exon_match_count)
# print(green_exon_match_count)


# In[20]:

def get_matching_ref_indices(start_index, end_index):
    ref_seq_index_matches = []
    for i in range(end_index - start_index):
        ref_seq_index_matches.append(int(reference_read[start_index+i]))
    return ref_seq_index_matches


# In[21]:

len_ref_seq = len(merged_sequence)
def get_mismatch_count_ref_seq(index, read):
    mismatch_count = 0
    if len_ref_seq <= index + len(read):
#         print('Length mismatch')
        return -1
    for i in range(len(read)):
        if(merged_sequence[index+i] != read[i]):
            mismatch_count += 1
    return mismatch_count


# In[22]:

## Check for mismatches
for read in all_reads_removing_N:
    read_length = len(read)
    first_part = read[0 : read_length/3]
    second_part = read[read_length/3 : 2*read_length/3]
    third_part = read[2*read_length/3 : read_length]
    first_part_band_match_start_index, first_part_band_match_end_index = get_match_indices(first_part)
    second_part_band_match_start_index, second_part_band_match_end_index = get_match_indices(second_part)
    third_part_band_match_start_index, third_part_band_match_end_index = get_match_indices(third_part)
#     print(first_part_band_match_start_index, first_part_band_match_end_index)
#     print(second_part_band_match_start_index, second_part_band_match_end_index)
#     print(third_part_band_match_start_index, third_part_band_match_end_index)
    
    first_ref_matches = []
    second_ref_matches = []
    third_ref_matches = []
    
    if first_part_band_match_start_index is not None and first_part_band_match_end_index is not None:
        first_ref_matches = get_matching_ref_indices(first_part_band_match_start_index, first_part_band_match_end_index)
    
    if second_part_band_match_start_index is not None and second_part_band_match_end_index is not None:
        second_ref_matches = get_matching_ref_indices(second_part_band_match_start_index, second_part_band_match_end_index)
    
    if third_part_band_match_start_index is not None and third_part_band_match_end_index is not None:
        third_ref_matches = get_matching_ref_indices(third_part_band_match_start_index, third_part_band_match_end_index)
#     print(first_ref_matches)
#     print(second_ref_matches)
#     print(third_ref_matches)
    
    len_first_part = len(first_part)
    len_second_part = len(second_part)
    len_third_part = len(third_part)
#     print(len_first_part, len_second_part, len_third_part)
    second_ref_matches = [val - len_first_part for val in second_ref_matches]
#     print(second_ref_matches)
    third_ref_matches = [val -len_first_part - len_second_part for val in third_ref_matches]
#     print(third_ref_matches)
    
    common_matches = set(first_ref_matches) & set(second_ref_matches) & set(third_ref_matches)
#     print(common_matches)
    first_remaining_matches = set(first_ref_matches).difference(common_matches)
    second_remaining_matches = set(second_ref_matches).difference(common_matches)
    third_remaining_matches = set(third_ref_matches).difference(common_matches)
#     print(first_remaining_matches, second_remaining_matches, third_remaining_matches)
    all_remaining_matches = list(first_remaining_matches.union(second_remaining_matches.union(third_remaining_matches)))
#     print(all_remaining_matches)
    
    if(len(all_remaining_matches) > 0):
        red_exon_numbers = []
        green_exon_numbers = []
        for i in range(len(all_remaining_matches)):
            mismatch_count = get_mismatch_count_ref_seq(all_remaining_matches[i], read)
#             print(mismatch_count)
            if mismatch_count > 0 and mismatch_count < 3:
                is_red_match, red_exon_number = is_matching_red(all_remaining_matches[i], read_length)
                is_green_match, green_exon_number = is_matching_green(all_remaining_matches[i], read_length)
                if is_red_match:
                    red_exon_numbers.append(red_exon_number)
                if is_green_match:
                    green_exon_numbers.append(green_exon_number)
        
        if len(red_exon_numbers) > 0 and len(green_exon_numbers) > 0:
            for exon_number in red_exon_numbers:
                red_exon_match_count[exon_number] += 0.5
            for exon_number in green_exon_numbers:
                green_exon_match_count[exon_number] += 0.5
        elif len(red_exon_numbers) > 0:
            for exon_number in red_exon_numbers:
                red_exon_match_count[exon_number] += 1
        elif len(green_exon_numbers) > 0:
            for exon_number in green_exon_numbers:
                green_exon_match_count[exon_number] += 0.5
        else:
            pass


# In[23]:

print("Red exon count considering up to two mismatches :: ")
print(red_exon_match_count)
print("Grren exon count considering up to two mismatches :: ")
print(green_exon_match_count)
# print(read)
# print i
# print(len(merged_sequence))
# print(all_remaining_matches[i])
# print(len(read))


# In[ ]:

# ## Check for mismatches
# for read in all_reads_reverse_compliment:
#     read_length = len(read)
#     first_part = read[0 : read_length/3]
#     second_part = read[read_length/3 : 2*read_length/3]
#     third_part = read[2*read_length/3 : read_length]
#     first_part_band_match_start_index, first_part_band_match_end_index = get_match_indices(first_part)
#     second_part_band_match_start_index, second_part_band_match_end_index = get_match_indices(second_part)
#     third_part_band_match_start_index, third_part_band_match_end_index = get_match_indices(third_part)
# #     print(first_part_band_match_start_index, first_part_band_match_end_index)
# #     print(second_part_band_match_start_index, second_part_band_match_end_index)
# #     print(third_part_band_match_start_index, third_part_band_match_end_index)
    
#     first_ref_matches = []
#     second_ref_matches = []
#     third_ref_matches = []
    
#     if first_part_band_match_start_index is not None and first_part_band_match_end_index is not None:
#         first_ref_matches = get_matching_ref_indices(first_part_band_match_start_index, first_part_band_match_end_index)
    
#     if second_part_band_match_start_index is not None and second_part_band_match_end_index is not None:
#         second_ref_matches = get_matching_ref_indices(second_part_band_match_start_index, second_part_band_match_end_index)
    
#     if third_part_band_match_start_index is not None and third_part_band_match_end_index is not None:
#         third_ref_matches = get_matching_ref_indices(third_part_band_match_start_index, third_part_band_match_end_index)
# #     print(first_ref_matches)
# #     print(second_ref_matches)
# #     print(third_ref_matches)
    
#     len_first_part = len(first_part)
#     len_second_part = len(second_part)
#     len_third_part = len(third_part)
# #     print(len_first_part, len_second_part, len_third_part)
#     second_ref_matches = [val - len_first_part for val in second_ref_matches]
# #     print(second_ref_matches)
#     third_ref_matches = [val -len_first_part - len_second_part for val in third_ref_matches]
# #     print(third_ref_matches)
    
#     common_matches = set(first_ref_matches) & set(second_ref_matches) & set(third_ref_matches)
# #     print(common_matches)
#     first_remaining_matches = set(first_ref_matches).difference(common_matches)
#     second_remaining_matches = set(second_ref_matches).difference(common_matches)
#     third_remaining_matches = set(third_ref_matches).difference(common_matches)
# #     print(first_remaining_matches, second_remaining_matches, third_remaining_matches)
#     all_remaining_matches = list(first_remaining_matches.union(second_remaining_matches.union(third_remaining_matches)))
# #     print(all_remaining_matches)
    
#     if(len(all_remaining_matches) > 0):
#         for i in range(len(all_remaining_matches)):
#             mismatch_count = get_mismatch_count_ref_seq(all_remaining_matches[i], read)
# #             print(mismatch_count)
#             if mismatch_count > 0 and mismatch_count < 3:
#                 is_red_match, red_exon_number = is_matching_red(all_remaining_matches[i], read_length)
#                 is_green_match, green_exon_number = is_matching_green(all_remaining_matches[i], read_length)
#                 if is_red_match and is_green_match:
#                     red_exon_match_count[red_exon_number] += 0.5
#                     green_exon_match_count[green_exon_number] += 0.5
#                 elif is_red_match:
#                     red_exon_match_count[red_exon_number] += 1
#                 elif is_green_match:
#                     green_exon_match_count[green_exon_number] += 1
#                 else:
#                     assert True == True


# In[ ]:

# print(red_exon_match_count)
# print(green_exon_match_count)

