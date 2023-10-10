
import numpy as np
from random import randrange

problem_size_lst = [6, 10, 20, 40, 60]

for i in range(len(problem_size_lst)):
    
    f = open(f'{problem_size_lst[i]}.txt', 'w')
    

    distance = int(problem_size_lst[i] / 2 + 1)

    f.write(f'{problem_size_lst[i]}\n')

    for j in range(problem_size_lst[i]):

        if j+distance > problem_size_lst[i] - 1:
            f.write(f'{j} {j + distance - problem_size_lst[i]} {distance}\n')

        else:
            f.write(f'{j} {j + distance} {distance}\n')






    
    # constraint_count = (problem_size_lst[i] - 1) // 2
    # f.write(f'{constraint_count}\n')
    # f.write(f'{problem_size_lst[i]}\n')

    # standard_permutation = list(range(problem_size_lst[i]))
    # random_permutation = np.random.permutation(problem_size_lst[i])

    # print("random_permutation: ", random_permutation)
    # print("standard_permutation: ", list(range(problem_size_lst[i])))
    # print()

    # f.write(f'{random_permutation}\n')

    # for j in range(constraint_count):
    #     # f.write(f'2 * j: {2 * j}, 2 * j + 1: {2 * j+1}\n')
    #     idx1 = random_permutation[2 * j]
    #     idx2 = random_permutation[2 * j + 1]
    #     # print(j)
    #     # print(f'idx1: ', idx1)
    #     # print(f'idx2: ', idx2)
        
    #     if idx1 < idx2:
    #         # print("distance: ", idx2 - idx1)
    #         f.write(f'{random_permutation[2 * j]} {random_permutation[2 * j + 1]} {idx2 - idx1}\n')
    #     else:
    #         # print("distance: ", idx1 - idx2)
    #         f.write(f'{random_permutation[2 * j + 1]} {random_permutation[2 * j]} {idx1 - idx2}\n')

    # print()

    f.close()