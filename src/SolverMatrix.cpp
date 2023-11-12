// Copyright 2023 Rik Essenius
// 
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file
// except in compliance with the License. You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software distributed under the License
// is distributed on an "AS IS" BASIS WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and limitations under the License.

#include "SolverMatrix.h"
#include <cassert>
#include <cmath>
#include <iostream>
#ifdef _WIN32
#include <corecrt_math_defines.h>
#endif
SolverMatrix::SolverMatrix(const std::initializer_list<std::initializer_list<double>> list) : Matrix(list) {}

SolverMatrix::SolverMatrix(const Matrix &other) : Matrix(other) {}

Matrix SolverMatrix::getEigenvalues() const {
    assert(isSquare() && rowCount() < 4);
    if (rowCount() == 1) {
        return Matrix({ {me(0, 0)} });
    }
    if (rowCount() == 2) {
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
    // another way for calculating a0 is (-trace * trace * trace - 2 * cubic().getTrace() + 3 * trace * traceSquared) / 6.0
    const double a0 = -determinant;

    // see https://mathworld.wolfram.com/CubicFormula.html
    double q = (3 * a1 - a2 * a2) / 9.0;
    const double r = (9 * a2 * a1 - 27 * a0 - 2 * a2 * a2 * a2) / 54.0;

    // x^3 + 3qx - 2r = 0

    const double discriminant = q * q * q + r * r ;

    if (discriminant < -Epsilon) {
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
    if (discriminant > Epsilon) {
        // One real root and two complex conjugate roots. Return just the real one.
        const double s = cbrt(r + sqrt(discriminant));
        const double t = cbrt(r - sqrt(discriminant));
        return Matrix({ {-a2 / 3.0 + s + t} });
    }
    // discriminant is 0
    const double alpha = cbrt(r);
    const double root1 = -a2 / 3 + 2 * alpha;
    const double root2 = -alpha - a2 / 3;
    return fabs(root1 - root2) < Epsilon ? Matrix({ { root1 } }) : Matrix({ {root1}, {root2} });
}


/// @brief get the null space of the matrix. This should return at least 1 vector for the eigenvalue matrix (A - lambda * I).
/// @note The matrix is already expected to be in row echelon form, so we can just look at the free variables.
/// @return the vectors in the null space
Matrix SolverMatrix::getNullSpace() const {
    const auto freeVariables = getFreeVariables();
    if (freeVariables.empty()) {
        // no free variables, so no null space. 
        // We need the rows to have dimension 3 to be able to multiply with the permutation matrix
        return Matrix(3, 0);
    }
    Matrix result(rowCount(), static_cast<Dimension>(freeVariables.size()));
    auto resultColumn = 0;
    for (const auto freeVariable : freeVariables) {
        result.setColumn(resultColumn, getColumn(freeVariable) * -1);
        result(freeVariable, resultColumn) = 1;
        for(const auto otherFreeVariable: freeVariables) {
            if (otherFreeVariable != freeVariable) {
                result(otherFreeVariable, resultColumn) = 0;
            }
        }
        resultColumn++;
    }
    return result;
}


Matrix SolverMatrix::getEigenvectorFor(const double lambda) const {
    assert(isSquare());

    auto beta =  SolverMatrix(*this - getIdentity(rowCount()) * lambda);
    const auto permutation = beta.toReducedRowEchelonFormWithPivot();
    return permutation * beta.getNullSpace();
}

Matrix SolverMatrix::getEigenvectors() const {
    const auto eigenvalues = getEigenvalues();
    Matrix result(rowCount(), rowCount());
    Dimension currentRow = 0;
    for (Dimension eigenValueIndex = 0; eigenValueIndex < eigenvalues.rowCount(); eigenValueIndex++) {
        auto eigenvectors = getEigenvectorFor(eigenvalues(eigenValueIndex, 0));
        for (Dimension vectorIndex = 0; vectorIndex < eigenvectors.columnCount(); vectorIndex++) {
            auto eigenvector = Matrix(eigenvectors.getColumn(vectorIndex)).normalized();
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
    while (resultsFound < rowCount()) {
        if (abs(me(row, column)) <= EigenEpsilon) {
            result.push_back(column);
            resultsFound++;
        } else {
            row++;
            resultsFound++;
        }
        column++;
        if (column >= columnCount()) {
            column = 0;
            row++;
        }
    }
    return result;
}

void SolverMatrix::eliminatePivotValueInRow(const Dimension pivot, const Dimension row) {
    assert(row != pivot);
    if (me(pivot, pivot) < EigenEpsilon) return;
    const double valueToEliminate = me(row, pivot);               
    if (abs(valueToEliminate) < EigenEpsilon) return;
    const double compensationFactor = -valueToEliminate / me(pivot, pivot);
    for (Dimension column = 0; column < columnCount(); column++) {
        (*this)(row, column) += compensationFactor * me(pivot, column);
    }
}

void SolverMatrix::multiplyRow(const Dimension row, const double factor) {
    for (Dimension column = 0; column < columnCount(); column++) {
        (*this)(row, column) *= factor;
    }
}

void SolverMatrix::findMaxPivot(const Dimension& pivot, Dimension& maxRow, Dimension& maxColumn) const {
    maxRow = pivot;
    maxColumn = pivot;
    double maxValue = 0;

    for(Dimension searchRow = pivot; searchRow < rowCount(); searchRow++) {
        for(Dimension searchColumn = pivot; searchColumn < columnCount(); searchColumn++) {
            if(abs(me(searchRow, searchColumn)) > maxValue) {
                maxRow = searchRow;
                maxColumn = searchColumn;
                maxValue = abs(me(searchRow, searchColumn));
            }
        }
    }
}

Matrix SolverMatrix::toReducedRowEchelonFormWithPivot() {
    auto permutation = getIdentity(columnCount());
    const auto maxPivot = std::min(rowCount(), columnCount());
    
    for (Dimension pivot = 0; pivot < maxPivot; pivot++) {
        Dimension maxRow = pivot;
        Dimension maxColumn = pivot;
        findMaxPivot(pivot, maxRow, maxColumn);

        // swap rows and/or columns to bring pivot to (row, row)
        // if we do a column swap, we also need to swap the permutation matrix
        swapRows(pivot, maxRow);
        swapColumns(pivot, maxColumn);
        permutation.swapColumns(pivot, maxColumn);

        // make pivot element equal to 1

        const auto pivotValue = me(pivot, pivot);
        if (abs(pivotValue) > EigenEpsilon) {
            multiplyRow(pivot, 1.0 / pivotValue);
        }

        // eliminate all other elements below the pivot

        for (Dimension subRow = pivot + 1; subRow < rowCount(); subRow++) {
            eliminatePivotValueInRow(pivot, subRow);                
        }    
    }

    // back-subsitution    

    // using int instead of Dimension as Dimension is never negative

    for (int pivot = maxPivot - 1; pivot >= 0; pivot--) {
        if (abs(me(pivot, pivot)) > EigenEpsilon) {
            multiplyRow(pivot, 1.0 / me(pivot, pivot));
        }

        // eliminate all entries above pivot
        for (int subRow = pivot - 1; subRow >= 0; subRow--) {
            eliminatePivotValueInRow(pivot, subRow);
        }
    }
    return permutation;
}