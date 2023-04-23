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
    void toRowEchelonForm();

protected:
    void swapRows(Dimension row1, Dimension row2);
};

#endif // MATRIX_H
