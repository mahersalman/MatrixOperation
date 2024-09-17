# Matrix Operations Calculator (Rank, Determinant, and Inverse)

This project implements algorithms in both Haskell and Python to calculate the rank, determinant, and inverse of an NxN matrix using Gaussian Elimination and Gauss-Jordan Elimination.

## Algorithms and Explanation

### 1. **Gaussian Elimination (for Determinant and Rank Calculation)**
Gaussian Elimination transforms a matrix into an upper triangular form through row operations (row swapping, scaling, and adding multiples of rows). 
- **Determinant**: After reducing to upper triangular form, the determinant is the product of the diagonal elements. Adjust the sign for each row swap.
- **Rank**: The rank is determined by counting the non-zero rows after transformation.

#### Complexity:
- **Time Complexity**: O(N³) 
- **Space Complexity**: O(N²)

### 2. **Gauss-Jordan Elimination (for Inverse Calculation)**
Gauss-Jordan extends Gaussian Elimination to convert a matrix into the identity matrix while transforming the identity matrix into the inverse.
- **Inverse**: The augmented matrix [A|I] is transformed into [I|A⁻¹] through row operations.

#### Complexity:
- **Time Complexity**: O(N³)
- **Space Complexity**: O(N²)

## Implementations

### **Python**:
The Python version uses procedural programming to implement matrix operations with nested loops and conditionals, ensuring ease of use and flexibility.

### **Haskell**:
The Haskell version leverages functional programming paradigms such as lazy evaluation, offering an efficient and concise approach to matrix operations.

Both implementations ensure optimized matrix processing for large NxN matrices.

## Results:
- For matrix inverse, both Python and Haskell implementations print the first 5 elements from the first and last two rows.