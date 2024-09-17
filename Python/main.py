import numpy as np
import time
import copy

def gaussian_elimination(mat):
    n = len(mat)
    sign = 1
    row = 0

    for col in range(n):
        # Find Pivot
        pivot_row = None
        for r in range(row, n):
            if mat[r][col] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            sign = 0
            continue

        # Swap rows
        if pivot_row != row:
            mat[row], mat[pivot_row] = mat[pivot_row], mat[row]
            sign *= -1

        # Elimination
        for r in range(row + 1, n):
            if mat[r][col] != 0:
                factor = mat[r][col] / mat[row][col]
                for c in range(col, n):
                    mat[r][c] -= factor * mat[row][c]

        row += 1

    return mat, sign

def Determinant(mat):
    matrix,det = gaussian_elimination(mat)
    if det == 0:
        return 0
    for i in range(len(matrix)):
        det *= matrix[i][i]

    return det

def rank(mat):
    matrix, _ = gaussian_elimination(mat)
    rank = 0
    for row in matrix:
        if any(abs(x) > 1e-10 for x in row):
            rank += 1
    return rank

def Test_Ele():
    test = [[1,2,3,4]
           ,[0,4,3,4]
           ,[0,8,6,9]
           ,[0,0,0,2]]
    result ,s = gaussian_elimination(test)
    for i in result:
        print(i)
    print(f"rank        = > result = {rank(test)} , expected = {np.linalg.matrix_rank(test)} ")
    print(f"determinant = > result = {Determinant(test)} , expected= {np.linalg.det(test)}")




def matrix_inverse(matrix):
    n = len(matrix)
    # Check if the matrix is square
    if n != len(matrix[0]):
        raise ValueError("Input matrix must be square")

    identity_matrix = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    augmented_matrix = [row[:] + identity_row for row, identity_row in zip(matrix, identity_matrix)]

    # Gaussian elimination with partial pivoting
    for k in range(n):
        # Find pivot row
        pivot_row = None
        for i in range(k, n):
            if augmented_matrix[i][k] != 0:
                pivot_row = i
                break
        if pivot_row is None:
            return None

        # Swap rows
        augmented_matrix[k], augmented_matrix[pivot_row] = augmented_matrix[pivot_row], augmented_matrix[k]

        # Normalize pivot row
        pivot = augmented_matrix[k][k]
        for j in range(k, 2 * n):
            augmented_matrix[k][j] /= pivot

        # Eliminate other rows
        for i in range(n):
            if i != k:
                factor = augmented_matrix[i][k]
                for j in range(k, 2 * n):
                    augmented_matrix[i][j] -= factor * augmented_matrix[k][j]

    # Extract inverse matrix
    inverse_matrix = [[row[n + i] for i in range(n)] for row in augmented_matrix]

    return inverse_matrix

"""                   Tests                                 """

def assert_matrix_equal(matrix1, matrix2, tol=1e-9):
    assert len(matrix1) == len(matrix2) and len(matrix1[0]) == len(matrix2[0]), "Matrices have different dimensions"
    for i in range(len(matrix1)):
        for j in range(len(matrix1[0])):
            assert abs(matrix1[i][j] - matrix2[i][j]) < tol, f"Matrices are not equal at element ({i}, {j})"
def test_inverse():
    # Test case 1: Inverse of identity matrix
    matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    expected = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    result = matrix_inverse(matrix)
    assert_matrix_equal(result, expected)
    # Test case 2: Inverse of singular matrix
    matrix = [[1, 2], [2, 4]]
    assert matrix_inverse(matrix) is None, "Expected None for singular matrix"
    # Test case 3: Inverse of non-singular matrix
    matrix = [[4, 7], [2, 6]]
    expected = [[0.6, -0.7], [-0.2, 0.4]]
    result = matrix_inverse(matrix)
    assert_matrix_equal(result, expected)

    # Test case 4: Inverse of another non-singular matrix
    matrix = [[1, 2, 3], [0, 1, 4], [5, 6, 0]]
    expected = [[-24, 18, 5], [20, -15, -4], [-5, 4, 1]]
    result = matrix_inverse(matrix)
    assert_matrix_equal(result, expected)



"""      Final Test                      """
def read_matrix_from_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    new_lines = [s.replace(',', '') for s in lines]

    matrix = []
    for line in new_lines:
        row = [float(num) for num in line.split()]
        matrix.append(row)

    return matrix

def Final_Test():
    st = time.time()
    filename1 = 'Inverse.txt'
    matrix = read_matrix_from_file(filename1)
    result = matrix_inverse(matrix)
    expected = np.linalg.inv(matrix)
    assert_matrix_equal(expected,result)
    en = time.time()
    print(f"Calculate Inverse time : {en-st} sec")
    for i in range(5):
        print(result[0][i],end = " ")
    print()
    for i in range(5):
        print(result[1][i],end = " ")
    print()
    for i in range(5):
        print(result[-2][i],end = " ")
    print()
    for i in range(5):
        print(result[-1][i],end = " ")
    print()
    print()

    st = time.time()
    filename2 = 'Determinant.txt'
    matrix2 = read_matrix_from_file(filename2)
    det_matrix2 = Determinant(matrix2)
    assert np.isclose(np.linalg.det(matrix2), det_matrix2), f"{det_matrix2} must be : Determinant is not as expected: {np.linalg.det(matrix2)}\n"
    en = time.time()
    print(f"Determinant => result =  {det_matrix2}, expected = {np.linalg.det(matrix2)}")
    print(f"Calculate Determinant time : {en-st} sec")
    print()

    st = time.time()
    filename3 = 'Rank.txt'
    matrix3 = read_matrix_from_file(filename3)
    mcopy = copy.deepcopy(matrix3)
    result_rank = rank(mcopy)
    expected_rank = np.linalg.matrix_rank(matrix3)
    assert result_rank == expected_rank , f"get {result_rank} but must get {expected_rank}"
    en = time.time()
    print(f"Rank = {result_rank} , expected : {expected_rank}")
    print(f"Calculate Rank time : {en-st}")


if __name__ == '__main__':
    #test_inverse()
    #Test_Ele()
    Final_Test()
