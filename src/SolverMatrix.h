#ifndef SOLVERMATRIX_H
#define SOLVERMATRIX_H

#include <vector>
#include "Matrix.h"

typedef unsigned int Dimension;

/// Class for more complex matrix manipulations
class SolverMatrix: public Matrix {
public:
    SolverMatrix(Dimension rows, Dimension columns);
    SolverMatrix(std::initializer_list<std::initializer_list<double>> list);
    SolverMatrix(const Matrix& other);
    Matrix getEigenvalues() const;
    Matrix getEigenvectorFor(double lambda) const;
    Matrix getEigenvectors() const;
    std::vector<Dimension> getFreeVariables() const;
    Matrix getNullSpace() const;
    // void toRowEchelonForm();

    // converts itself to RREF and returns the permutation matrix. 
    Matrix toReducedRowEchelonFormWithPivot();

protected:
    void eliminatePivotValueInRow(Dimension pivot, Dimension row, double epsilon);
    void multiplyRow(Dimension row, double factor);
};

#endif // MATRIX_H
