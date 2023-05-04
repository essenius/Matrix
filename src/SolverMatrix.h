#ifndef SOLVERMATRIX_H
#define SOLVERMATRIX_H

#include <vector>
#include "Matrix.h"

typedef unsigned int Dimension;

/// Class for more complex matrix manipulations
class SolverMatrix: public Matrix {
public:
    explicit SolverMatrix(std::initializer_list<std::initializer_list<double>> list);
    SolverMatrix(const Matrix& other);
    Matrix getEigenvalues() const;
    Matrix getEigenvectorFor(double lambda) const;
    Matrix getEigenvectors() const;
    std::vector<Dimension> getFreeVariables() const;
    Matrix getNullSpace() const;

    // converts itself to RREF and returns the permutation matrix. 
    Matrix toReducedRowEchelonFormWithPivot();

    static constexpr double EIGEN_EPSILON = 1e-6;
protected:
    void eliminatePivotValueInRow(Dimension pivot, Dimension row);
    void findMaxPivot(const Dimension &pivot, Dimension& maxRow, Dimension& maxColumn) const;
    void multiplyRow(Dimension row, double factor);
};

#endif // MATRIX_H
