# Matrix

This repo provides a library for matrix operations. It is intended for use in Arduino projects, but can be used in other projects as well.

## Array

Class for array manipulations (coefficient wise):

- `+`, `-`, `*`, `/`: operations per element in the array with another array or with a scalar
- `+=`, `-=`, `*=`: updates element by the operation (argument is array or scalar) 			
- `/=`: divide elements by scalar
- `==`: test for equality of all elements
- `rowCount`, `columnCount`: number of rows/columns in the array
- `getColumn`, `getRow`: get one row or column
- `isSizeEqual`: do the arrays have the same number of rows and columns
- `isSquare`: is the row count equal to the column count
- `me`, `()`: get an element indicated bu row and column
- `pow2`: multiply elements by themselves
- `swapRows`, `swapColumns`: swap two rows or columns
- `setColumn`: set all elements of a column to a value, or to a vector.

## Matrix

Class for basic matrix operations

- Multiply: matrix multiplication
- `getAdjoint`: cofactor matrix
- `getAdjugate`: transpose of adjoint
- `getCofactor`: product of the minor of the element and -1^(positional value of element)
- `getDeterminant`: the [determinant](https://en.wikipedia.org/wiki/Determinant) of the matrix 
- `getTrace`: sum of the elements on the main diagonal
- `getMinor`: determinant of matrix not including the indicated row/column
- `getIdentity`: identity matrix
- `inverted`: the matrix returning the indentity matrix when multiplied by the original matrix
- `isInvertible`: whether or not a matrix is invertible
- `normalized`: each element divided by the square root of the sum of the squared elements
- `squared`: matrix multiplied by itself
- `toArray`: convert matrix to array
- `transposed`: swap rows and columns

## SolverMatrix

- `getEigenvales`: determine the [eigenvalues](https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors). Works up to 3x3 matrices
- `getEigenvectors`: determine the [eigenvectors](https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors). Works up to 3x3 matrices
- `getEigenvectorFor`: determine the eigenvector belonging to an eigenvalue. Expects an earlier calculated eigenvalue.
- `getFreeVariables`: determine the [free variables](https://en.wikipedia.org/wiki/Free_variables_and_bound_variables). Expects a matrix in reduced row echelon form.
- `getNullSpace`: determine the [Null space](https://en.wikipedia.org/wiki/Kernel_(linear_algebra)). Expects a matrix in reduced row echelon form.
- `toReducedRowEchelonForm`: determine the [Reduced Row Echelon Form](https://en.wikipedia.org/wiki/Row_echelon_form#rref). Doing this makes finding eigenvalues much simpler.

## Structure

The repo was setup such that it can be built both with CMake and with MSBuild. There is a Visual Studio solution that uses MSBuild.
In Visual Studio, the release build also creates a ZIP file that can be imported into the Arduino IDE as a library.
This is why we also have `keywords.txt` and `library.properties`.

Google Test is used as the testing framework. Asserts are used to ensure that preconditions are met. They have as consequence that code 
coverage gets a bit lower (now ~97% overall), but they are useful for debugging. 
The asserts are only used in the debug build. In the release build, they are replaced by empty macros.