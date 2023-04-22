#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include "BasicMatrix.h"

typedef unsigned int Dimension;

/// Class for more complex matrix manipulations
class Matrix: public BasicMatrix {
public:
    Matrix(Dimension rows, Dimension columns);
    Matrix(std::initializer_list<std::initializer_list<double>> list);
    Matrix(const BasicMatrix& other);
    BasicMatrix getEigenvalues() const;
    BasicMatrix getEigenvectorFor(double lambda) const;
    BasicMatrix getEigenvectors() const;
    std::vector<Dimension> getFreeVariables() const;
    BasicMatrix getNullSpace() const;
    void toRowEchelonForm();

protected:
    void swapRows(Dimension row1, Dimension row2);
};

#endif // MATRIX_H
