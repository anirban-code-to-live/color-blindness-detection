import numpy as np
import pickle


if __name__ == '__main__':
    print("Welcome to the world of color-blindness!!")
    # bwt_data = np.loadtxt('../data/chrX_last_col.txt', dtype=str)
    # print(bwt_data)
    # print(len(bwt_data))
    # # print(list(set(bwt_data)))
    # first_col_dict = {"A": 0, "G": 0, "C": 0, "T": 0, "$": 0, "N": 0, "U": 0}
    # counts = [0, 0, 0, 0, 0, 0, 0]
    # for line in bwt_data:
    #     for token in line:
    #         if token == "A":
    #             counts[0] += 1
    #         elif token == "C":
    #             counts[1] += 1
    #         elif token == "T":
    #             counts[2] += 1
    #         elif token == "G":
    #             counts[3] += 1
    #         elif token == "$":
    #             counts[4] += 1
    #         elif token == "N":
    #             counts[5] += 1
    #         else:
    #             counts[6] += 1
    #
    # first_col_dict["A"] = counts[0]
    # first_col_dict["C"] = counts[1]
    # first_col_dict["T"] = counts[2]
    # first_col_dict["G"] = counts[3]
    # first_col_dict["$"] = counts[4]
    # first_col_dict["N"] = counts[5]
    # first_col_dict["U"] = counts[6]
    # print(counts)
    # print(first_col_dict)
    #
    # filename = '../processed_data/first_col.pkl'
    # with open(filename, 'wb') as f:
    #     pickle.dump(first_col_dict, f)

    reference_read = np.loadtxt("../data/chrX_map.txt", dtype=str)
    print(len(reference_read))
    test_read = "GAGGACAGCACCCAGTCCAGCATCTTCACCTACACCAACAGCAACTCCACCAGAGGTGAGCCAGCAGGCCCGTGGAGGCTGGGTGGCTGCACTGGGGGCCA"
