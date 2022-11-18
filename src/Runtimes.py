################################################################
# libraries:

# import time
# import random
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns

################################################################
# functions:

# from Readmapper_Cython_Convert import SuffixArray
# from Readmapper_Cython_Convert import approx_search
# from SEQsimulator import simulate_string
# from SEQsimulator import get_exact_read
# from naive import naive_algorithm

# # Usage:
# string  = 'aaaaaaaaaaaaaaa'
# pattern = 'ba'
# SA = SuffixArray(string)
# print(approx_search(string, pattern, SA, 2))

################################################################
# tests:
 
# # 1mio test SA:
# string = simulate_string(1000000)
# start_time = time.time()
# SuffixArray(string)
# end_time = time.time()
# print(end_time-start_time)


# # Test suffixtree-algorithm vs naive-algorithm for same result:
# for i in range(50000):
#     print('Iteration nr: ', i+1)
#     ref = simulate_string(random.randint(40,100))
#     read = get_exact_read(ref, random.randint(0,30))
#     SA = SuffixArray(ref)
#     offsets = approx_search(ref, read, SA, 0)
#     if sorted(offsets) != sorted(naive_algorithm(ref,read)):
#         print('Algorithm mistake!')
#         print(ref)
#         print(read)
#         print(sorted(offsets))
#         print(sorted(naive_algorithm(ref,read)))
#         break
#     if i+1 == 50000:
#         print('DONE')


# # Runtimes for the array construction (varying ref lengths):
# ref_lengths = [25000,50000,75000,100000,125000,150000,175000]
# runtimes = []
# for idx in range(7):
#     print('Iteration nr: ', idx+1) 
#     replicate = []
#     for j in range(10):
#         ref = simulate_string(ref_lengths[idx])
#         start_time = time.time()
#         SuffixArray(ref)
#         end_time = time.time()
#         replicate.append(end_time-start_time)
#     runtimes.append(np.mean(replicate))
# # plot running times:
# fig, ax = plt.subplots()
# sns.lineplot(x=ref_lengths, y=runtimes, ax=ax)
# plt.xlabel('ref length')
# plt.ylabel('runtime (s)')
# plt.tight_layout()
# plt.show()


# # Runtimes for mapping (varying read lengths):
# ref_lengths = [2500,5000,7500,10000,12500]
# read_lengths_10 = [100]*5
# read_lengths_20 = [200]*5
# read_lengths_30 = [300]*5
# read_lengths_40 = [400]*5
# read_lengths_50 = [500]*5
# runtimes_10 = []
# runtimes_20 = []
# runtimes_30 = []
# runtimes_40 = []
# runtimes_50 = []
# for idx in range(5):
#     print('Iteration nr: ', idx+1)
#     runtimes_10_replicate = []
#     runtimes_20_replicate = []
#     runtimes_30_replicate = []
#     runtimes_40_replicate = []
#     runtimes_50_replicate = []
    
#     for i in range(100):
#         ref = simulate_string(ref_lengths[idx])
#         SA = SuffixArray(ref)
#         mm_allowed = 5
        
#         read_10 = get_exact_read(ref, read_lengths_10[idx])
#         read_20 = get_exact_read(ref, read_lengths_20[idx])
#         read_30 = get_exact_read(ref, read_lengths_30[idx])
#         read_40 = get_exact_read(ref, read_lengths_40[idx])
#         read_50 = get_exact_read(ref, read_lengths_50[idx])
        
#         # Dont know why it is nessesary to read run this chunk of code, but 
#         # if i dont the first the first runtime is affected. Mayde something 
#         # to do with loading the modules??
#         offsets = approx_search(ref, read_10, SA, mm_allowed)
        
#         start_time = time.time()
#         offsets = approx_search(ref, read_10, SA, mm_allowed)
#         end_time = time.time()
#         runtimes_10_replicate.append(end_time-start_time)
        
#         start_time = time.time()
#         offsets = approx_search(ref, read_20, SA, mm_allowed)
#         end_time = time.time()
#         runtimes_20_replicate.append(end_time-start_time)
        
#         start_time = time.time()
#         offsets = approx_search(ref, read_30, SA, mm_allowed)
#         end_time = time.time()
#         runtimes_30_replicate.append(end_time-start_time)
        
#         start_time = time.time()
#         offsets = approx_search(ref, read_40, SA, mm_allowed)
#         end_time = time.time()
#         runtimes_40_replicate.append(end_time-start_time)
        
#         start_time = time.time()
#         offsets = approx_search(ref, read_50, SA, mm_allowed)
#         end_time = time.time()
#         runtimes_50_replicate.append(end_time-start_time)
        
#     runtimes_10.append(np.mean(runtimes_10_replicate))
#     runtimes_20.append(np.mean(runtimes_20_replicate))
#     runtimes_30.append(np.mean(runtimes_30_replicate))
#     runtimes_40.append(np.mean(runtimes_40_replicate))
#     runtimes_50.append(np.mean(runtimes_50_replicate))

# # plot running times:
# fig, ax = plt.subplots()
# sns.lineplot(x=ref_lengths, y=runtimes_50, ax=ax, label='read length = 500')
# sns.lineplot(x=ref_lengths, y=runtimes_40, ax=ax, label='read length = 400')
# sns.lineplot(x=ref_lengths, y=runtimes_30, ax=ax, label='read length = 300')
# sns.lineplot(x=ref_lengths, y=runtimes_20, ax=ax, label='read length = 200')
# sns.lineplot(x=ref_lengths, y=runtimes_10, ax=ax, label='read length = 100')
# plt.xlabel('ref length')
# plt.ylabel('runtime (s)')
# plt.tight_layout()
# plt.show()

################################################################
