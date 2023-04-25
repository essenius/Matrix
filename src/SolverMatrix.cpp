#include "SolverMatrix.h"

#include <stdexcept>
#include <cmath>
#include <iostream>

double SolverMatrix::_eigenEpsilon = 1e-6;

SolverMatrix::SolverMatrix(Dimension rows, Dimension columns): Matrix(rows, columns) {}

SolverMatrix::SolverMatrix(std::initializer_list<std::initializer_list<double>> list) : Matrix(list) {}

SolverMatrix::SolverMatrix(const Matrix &other) : Matrix(other) {}

double SolverMatrix::getEigenEpsilon() {
    return _eigenEpsilon;
}

Matrix SolverMatrix::getEigenvalues() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
    if (_rows > 3) {
        throw std::invalid_argument("Only works for 3d matrices");
    }
    if (_rows == 1) {
        return Matrix({ {me(0, 0)} });
    }
    if (_rows == 2) {
        const double determinant = getDeterminant();
        const double trace = getTrace();
        const double discriminant = trace * trace - 4 * determinant;
        if (discriminant < 0) {
            // no real eigenvalues
            return Matrix(0, 0);
        }
        const auto rootDiscriminant = sqrt(discriminant);
        return Matrix{
            {(trace + rootDiscriminant) / 2},
            {(trace - rootDiscriminant) / 2}
        };
    }
    // 3x3 matrix
    // Characteristic polynomial, see https://mathworld.wolfram.com/CharacteristicPolynomial.html, but swapping signs
    const double trace = getTrace();
    const double traceSquared = squared().getTrace();
    const double determinant = getDeterminant();
    const double a2 = -trace;
    const double a1 = (trace * trace - traceSquared) / 2.0;
    //const double a0 = (-trace * trace * trace - 2 * cubic().getTrace() + 3 * trace * traceSquared) / 6.0; 
    const double a0 = -determinant;

    // see https://mathworld.wolfram.com/CubicFormula.html
    double q = (3 * a1 - a2 * a2) / 9.0;
    const double r = (9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54.0;

    // x^3 + 3qx -2r = 0

    const double discriminant = q * q * q + r * r ;

    if (discriminant < -_epsilon) {
        // Three distinct real roots (or >1 coinciding roots if discriminant = 0)
        q = -q;
        const double theta = acos(r / sqrt(q * q * q));
        const double qFactor = 2 * sqrt(q);
        return Matrix({
            { qFactor * cos(theta / 3.0) - a2 / 3.0 },
            { qFactor * cos((theta + 2 * M_PI) / 3.0) - a2 / 3.0 },
            { qFactor * cos((theta + 4 * M_PI) / 3.0) - a2 / 3.0 }
        });
    }
    if (discriminant > _epsilon) {
        // One real root and two complex conjugate roots. Return just the real one.
        const double s = cbrt(r + sqrt(discriminant));
        const double t = cbrt(r - sqrt(discriminant));
        return Matrix({ {-a2 / 3.0 + s + t} });
    }
    // discriminant is 0
    const double alpha = cbrt(r);
    const double root1 = -a2 / 3 + 2 * alpha;
    const double root2 = -alpha - a2 / 3;
    return fabs(root1 - root2) < _epsilon ? Matrix({ { root1 } }) : Matrix({ {root1}, {root2} });
}


/// @brief get the null space of the matrix. This should return at least 1 vector for the eigenvalue matrix (A - lambda * I).
/// @note The matrix is already expected to be in row echelon form, so we can just look at the free variables.
/// @return the vectors in the null space
Matrix SolverMatrix::getNullSpace() const {
    const auto freeVariables = getFreeVariables();
    if (freeVariables.empty()) {
        // no free variables, so no null space
        return Matrix(0, 0);
    }
    Matrix result(_rows, freeVariables.size());
    auto resultColumn = 0;
    for (auto freeVariable : freeVariables) {
        result.setColumn(resultColumn, getColumn(freeVariable) * -1);
        result(freeVariable, resultColumn) = 1;
        for(auto otherFreeVariable: freeVariables) {
            if (otherFreeVariable != freeVariable) {
                result(otherFreeVariable, resultColumn) = 0;
            }
        }
        resultColumn++;
    }
    return result;
}

void SolverMatrix::setEigenEpsilon(double epsilon) {
    _eigenEpsilon = epsilon;
}

Matrix SolverMatrix::getEigenvectorFor(const double lambda) const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }

    auto beta =  SolverMatrix(*this - identity(_rows) * lambda);
    auto permutation = beta.toReducedRowEchelonFormWithPivot();
    return permutation * beta.getNullSpace();
}

Matrix SolverMatrix::getEigenvectors() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
    const auto eigenvalues = getEigenvalues();
    Matrix result(_rows, _rows);
    Dimension currentRow = 0;
    for (Dimension eigenValueIndex = 0; eigenValueIndex < eigenvalues.rows(); eigenValueIndex++) {
        auto eigenvectors = getEigenvectorFor(eigenvalues(eigenValueIndex, 0));
        for (Dimension vectorIndex = 0; vectorIndex < eigenvectors.columns(); vectorIndex++) {
            auto eigenvector = Matrix(eigenvectors.getColumn(vectorIndex)).normalize();
            result.setColumn(currentRow, eigenvector);
            currentRow++;
        }
    }
    result.setColumnCount(currentRow);
    return result;
}

/// @brief Finds the free variables in the matrix by searching for non-pivot columns.
/// This must already be a matrix in row echelon form
/// @return a vector of the free variable columns
std::vector<Dimension> SolverMatrix::getFreeVariables() const {
    std::vector<Dimension> result;
    Dimension row = 0;
    Dimension column = 0;
    Dimension resultsFound = 0;
    while (resultsFound < _rows) {
        if (abs(me(row, column)) <= _eigenEpsilon) {
            result.push_back(column);
            resultsFound++;
        } else {
            row++;
            resultsFound++;
        }
        column++;
        if (column >= _columns) {
            column = 0;
            row++;
        }
    }
    return result;
}

void SolverMatrix::eliminatePivotValueInRow(Dimension pivot, Dimension row) {
    if (row == pivot) return;
    if (me(pivot, pivot) < _eigenEpsilon) return;
    const double valueToEliminate = me(row, pivot);               
    if (abs(valueToEliminate) < _eigenEpsilon) return;
    const double compensationFactor = -valueToEliminate / me(pivot, pivot);
    for (Dimension column = 0; column < _columns; column++) {
        (*this)(row, column) += compensationFactor * me(pivot, column);
    }
}

void SolverMatrix::multiplyRow(Dimension row, double factor) {
    for (Dimension column = 0; column < _columns; column++) {
        (*this)(row, column) *= factor;
    }
}

Matrix SolverMatrix::toReducedRowEchelonFormWithPivot() {
    auto permutation = Matrix::identity(_columns);
    const auto maxPivot = std::min(_rows, _columns);
    
    for (Dimension pivot = 0; pivot < maxPivot; pivot++) {
        Dimension maxRow = pivot;
        Dimension maxColumn = pivot;
        double maxValue = abs(me(pivot, pivot));

        for(Dimension searchRow = pivot; searchRow < _rows; searchRow++) {
            for(Dimension searchColumn = pivot; searchColumn < _columns; searchColumn++) {
                if(abs(me(searchRow, searchColumn)) > maxValue) {
                    maxRow = searchRow;
                    maxColumn = searchColumn;
                    maxValue = abs(me(searchRow, searchColumn));
                }
            }
        }

        // swap rows and/or columns to bring pivot to (row, row)
        // if we do a column swap, we also need to swap the permutation matrix
        swapRows(pivot, maxRow);
        swapColumns(pivot, maxColumn);
        permutation.swapColumns(pivot, maxColumn);

        // make pivot element equal to 1

        auto pivotValue = me(pivot, pivot);
        if (abs(pivotValue) > _eigenEpsilon) {
            multiplyRow(pivot, 1.0 / pivotValue);
        }

        // eliminate all other elements below the pivot

        for (Dimension subRow = pivot + 1; subRow < _rows; subRow++) {
            eliminatePivotValueInRow(pivot, subRow);                
        }    
    }

    // back-subsitution    

    // using int instead of Dimension as Dimension is never negative

    for (int pivot = maxPivot - 1; pivot >= 0; pivot--) {
        if (abs(me(pivot, pivot)) > _eigenEpsilon) {
            multiplyRow(pivot, 1.0 / me(pivot, pivot));
        }

        // eliminate all entries above pivot
        for (int subRow = pivot - 1; subRow >= 0; subRow--) {
            eliminatePivotValueInRow(pivot, subRow);
        }
    }
    return permutation;
}