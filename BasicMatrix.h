#ifndef BASICMATRIX_H
#define BASICMATRIX_H

#include <vector>
#include "Array.h"

/// Class for basic matrix manipulations
class BasicMatrix : public Array {
public:
    /// Constructors
    BasicMatrix(Dimension rows, Dimension columns);
    BasicMatrix(std::initializer_list<std::initializer_list<double>> list);
    BasicMatrix(const BasicMatrix& other);
    BasicMatrix(const Array other);

    void operator*=(const BasicMatrix& other);

    BasicMatrix adjoint() const;
    BasicMatrix adjugate() const;
    double cofactor(Dimension row, Dimension column) const;
    BasicMatrix cubic() const;
    double getDeterminant() const;
    static BasicMatrix diagonal(std::initializer_list<double> list);
    double getTrace() const;
    static BasicMatrix identity(Dimension size);
    BasicMatrix inverse() const;
    bool isInvertible() const;
    BasicMatrix normalize() const;
    BasicMatrix squared() const;
    BasicMatrix transpose() const;

    friend BasicMatrix operator*(BasicMatrix left, const BasicMatrix& right);
};

#endif // MATRIX_H
