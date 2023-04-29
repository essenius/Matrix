#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include "Array.h"

/// Class for basic matrix manipulations
class Matrix : public Array {
public:
    /// Constructors
    Matrix(Dimension rows, Dimension columns);
    explicit Matrix(std::initializer_list<std::initializer_list<double>> list);
    explicit Matrix(const Array& other);

    // *= works differently in matrices
    void operator*=(const Matrix& other);
    // this one doesn't, but it's still needed
    void operator*=(double other);

    Matrix adjoint() const;
    Matrix adjugate() const;
    double getCofactor(Dimension row, Dimension column) const;
    double getDeterminant() const;
    double getTrace() const;
    Matrix getMinor(Dimension row, Dimension column) const;
    static Matrix identity(Dimension size);
    Matrix inverse() const;
    bool isInvertible() const;
    Matrix normalize() const;
    Matrix squared() const;
    Array toArray() const;
    Matrix transpose() const;

    // need to redefine these as well, didn't find a better way
    friend Matrix operator+(Matrix left, const Matrix& right);
    friend Matrix operator-(Matrix left, const Matrix& right); 
    friend Matrix operator*(Matrix left, const Matrix& right);
    friend Matrix operator*(Matrix left, double right);
    friend Matrix operator*(double left, Matrix right);
};

#endif
