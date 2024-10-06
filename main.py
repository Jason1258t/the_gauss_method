import math
from functools import reduce
from enum import Enum

A1 = [
    [2, 2, 1, 9],  # x + 2y + z = 9
    [2, 3, 3, 21],  # 2x + 3y + 3z = 21
    [2, 4, 2, 18]  # 3x + 4y + 2z = 18
]

A2 = [
    [1, 2, 1, 3],  # x + 2y + z = 3
    [2, 4, 2, 6],  # 2x + 4y + 2z = 6
    [1, 1, 1, 2]  # x + y + z = 2
]

A3 = [
    [1, 2, 5],  # x + 2y = 5
    [2, 4, 10],  # 2x + 4y = 10
    [1, 2, 7]  # x + 2y = 7
]


def gaussian_elimination(m):
    stepwise_m = get_stepwise(m)

    # print_m(stepwise_m)
    matrix_type = get_matrix_type(m)
    if matrix_type == MatrixTypes.determine:
        print(calculate_determined_matrix(stepwise_m))


def calculate_determined_matrix(m):
    variables_length = len(m)
    variables = dict(zip(['x' + str(x) for x in range(1, variables_length + 1)], [0] * variables_length))

    for row in m[::-1]:
        ind, main_ratio = get_non_nullable_index(row)
        if ind == len(row) - 2:
            variables[f'x{ind + 1}'] = row[-1] / main_ratio
        else:
            s = row[-1]
            for j in range(ind + 1, len(row) - 1):
                s -= row[j] * variables[f'x{j + 1}']

            variables[f'x{ind + 1}'] = s / main_ratio

    return variables



class MatrixTypes(Enum):
    determine = 1
    not_determine = 2
    incompatible = 3


def get_matrix_type(m):

    variables_count = len(m[0]) - 1
    extended_rank = get_extended_matrix_rank(m)
    base_rank = get_base_matrix_rank(m)

    if base_rank == extended_rank == variables_count:
        return MatrixTypes.determine
    elif base_rank == extended_rank:
        return MatrixTypes.not_determine
    else:
        return MatrixTypes.incompatible

def get_extended_matrix_rank(stepwise_matrix):
    c = 0
    for i in stepwise_matrix:
        f = False
        for j in i:
            if j != f:
                f = True
                break

        if f: c += 1

    return c

def get_base_matrix_rank(stepwise_matrix):
    c = 0
    for i in stepwise_matrix:
        f = False
        for j in i[:-1]:
            if j != f:
                f = True
                break

        if f: c += 1

    return c

def get_stepwise(m):
    m = sort_matrix(m)
    column = 0

    while not is_stepwise(m):
        for i in range(len(m) - 1, column, -1):
            row = m[i]
            if row[column] == 0: continue
            variants = [x[column] for x in m[column:i]]
            best_row_index = find_best_match_element(row[column], variants)

            if best_row_index is not None:
                another_row = m[column + best_row_index]
                ratio = round(row[column] / another_row[column])
                for j in range(len(row)):
                    m[i][j] -= ratio * another_row[j]
            else:
                m = math_rows_using_lcm(m, i, column, column)

        m = sort_matrix(m)
        column += 1

    for i in range(len(m)):
        gcd = find_gcd(m[i])
        for j in range(len(m[i])):
            m[i][j] //= gcd

    return m

def is_stepwise(m):
    prev_index = get_non_nullable_index(m[0])[0]
    for i in range(1, len(m)):
        current_ind = get_non_nullable_index(m[i])[0]
        if current_ind > prev_index:
            prev_index = current_ind
        else:
            return False

    return True

def math_rows_using_lcm(m, main_index, another_index, column):
    another_row = m[another_index]
    row = m[main_index]

    nok = lcm(row[column], another_row[column])
    main_ratio = nok // row[column]
    another_ratio = nok // another_row[column]
    for j in range(len(row)):
        m[main_index][j] *= main_ratio
        m[another_index][j] *= another_ratio

    ratio = round(m[main_index][column] / m[another_index][column])
    for j in range(len(row)):
        m[main_index][j] -= ratio * another_row[j]

    return m

def sort_matrix(m):
    for i in range(len(m)):
        start = find_optimal_start(m, i)
        m[i], m[start] = m[start], m[i]

    return m

def find_optimal_start(m, starts_at=0):
    optimal_row_start, start_element = get_non_nullable_index(m[starts_at])
    optimal_row_index = starts_at
    for i in range(starts_at, len(m)):
        row = m[i]
        ind, new_start_element = get_non_nullable_index(row)
        if ind < optimal_row_start or ind == optimal_row_index and (abs(new_start_element) < abs(start_element)):
            optimal_row_index = i
            optimal_row_start = ind
            start_element = new_start_element

    return optimal_row_index

def get_non_nullable_index(row):
    for i in range(len(row)):
        if row[i] != 0:
            return i, row[i]

    return len(row) + 1

def find_best_match_element(item, variants):
    for i in range(len(variants)):
        v = variants[i]
        if item % v == 0: return i

def lcm(a, b):
    return abs(a * b) // math.gcd(a, b)

def find_gcd(numbers):
    return reduce(math.gcd, numbers)

def print_m(m):
    for i in m:
        print(*i)


gaussian_elimination(A1)
