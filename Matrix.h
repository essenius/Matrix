#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

typedef unsigned int Dimension;
/// Class for matrix manipulations
class Matrix {
public:
    /// Constructors
    Matrix(Dimension rows, Dimension columns);
    Matrix(std::initializer_list<std::initializer_list<double>> list);
    Matrix(const Matrix& other);


    double& operator()(Dimension row, Dimension column);
    const double& operator()(Dimension row, Dimension column) const;

    void operator+=(const Matrix& other);
    void operator+=(double other);

    void operator-=(const Matrix& other);
    void operator-=(double other);

    void operator*=(const Matrix& other);
    void operator*=(double other);

    void operator/=(double other);

    std::vector<Dimension> getFreeVariables() const;
    bool operator==(const Matrix &other) const;

    Matrix adjoint() const;
    Matrix adjugate() const;
    double cofactor(Dimension row, Dimension column) const;
    Dimension columns() const;
    Matrix cubic() const;
    double getDeterminant() const;
    Matrix getNullSpace() const;
    bool contains(double value, double epsilon = 1e-10) const;
    static Matrix diagonal(std::initializer_list<double> list);
    Matrix eigenvalues() const;
    Matrix eigenvectorFor(double lambda, double epsilon = 1e-10) const;
    bool equalSize(const Matrix& other) const;
    Matrix getColumn(Dimension column) const;
    Dimension getRank() const;
    Matrix getRow(Dimension row) const;
    double getTrace() const;
    static Matrix identity(Dimension size);
    Matrix inverse() const;
    bool isInvertible() const;
    bool isSquare() const;
    double me(Dimension row, Dimension column) const;
    //const double me(Dimension row, Dimension column) const;
    Matrix normalize() const;
    Dimension rows() const;
    Matrix squared() const;
    void setAt(const Matrix& input, Dimension row, Dimension column);
    void setAll(double value);
    void setColumn(Dimension column, double value);
    void setColumn(Dimension column, const Matrix& input);
    void setRow(Dimension row, const Matrix& input);
    void toRowEchelonForm();
    Matrix transpose() const;

    friend Matrix operator+(Matrix left, const Matrix& right);
    friend Matrix operator-(Matrix left, const Matrix& right);
    friend Matrix operator*(Matrix left, const Matrix& right);
    friend Matrix operator*(Matrix left, double right);
    friend Matrix operator/(Matrix left, double right);


private:
    std::vector<double> _data;
    static const double _epsilon;
    Dimension _rows;
    Dimension _columns;
    Dimension _arraySize;
    /*Dimension _count;
    void setValues() {}
    template<typename T, typename... Ts>
    void setValues(T value, Ts... values) {
        if (_count < _rows * _columns) {
            _data[_count] = value;
            _count++;
            setValues(values...);
        }
    } */

    void swapRows(Dimension row1, Dimension row2);
};

#endif // MATRIX_H
