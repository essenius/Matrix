#include "Matrix.h"

#include <stdexcept>
#include <cmath>
#include <iostream>

const double Matrix::_epsilon = 1e-12; //std::numeric_limits<double>::epsilon();

void Matrix::swapRows(Dimension row1, Dimension row2) {
    if (row1 == row2) return;
    for (Dimension column = 0; column < _columns; ++column) {
        double temp = me(row1, column);
        (*this)(row1, column) = me(row2, column);
        (*this)(row2, column) = temp;
    }
}

Matrix::Matrix(const Dimension rows, const Dimension columns) : _data(rows* columns) {
    _rows = rows;
    _columns = columns;
    _arraySize = rows * columns;
}

Matrix::Matrix(const Matrix& other) {
    _rows = other._rows;
    _columns = other._columns;
    _arraySize = _rows * _columns;
    _data = other._data;
}

Matrix::Matrix(const std::initializer_list<std::initializer_list<double>> list) {
    _rows = static_cast<Dimension>(list.size());
    _columns = static_cast<Dimension>(list.begin()->size());
    _arraySize = _rows * _columns;
    _data = std::vector<double>(_rows * _columns);
    Dimension row = 0;
    for (auto rowList : list) {
        Dimension column = 0;
        for (const auto value : rowList) {
            _data[row * _columns + column] = value;
            column++;
        }
        row++;
    }
}

double& Matrix::operator()(Dimension row, Dimension column) {
    return _data[row * _columns + column];
}

const double& Matrix::operator()(Dimension row, Dimension column) const {
    return _data[row * _columns + column];
}

void Matrix::operator+=(const Matrix& other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
        _data[i] += other._data[i];
    }
}

void Matrix::operator+=(double other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
        _data[i] += other;
    }
}

void Matrix::operator-=(const Matrix& other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
        _data[i] -= other._data[i];
    }
}

void Matrix::operator/=(const double other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
        _data[i] /= other;
    }
}

void Matrix::operator*=(const Matrix& other) {
    if (_columns != other._rows) {
        throw std::invalid_argument("Matrix dimensions do not match");
    }
    Matrix result(_rows, other._columns);
    for (Dimension row = 0; row < _rows; ++row) {
        for (Dimension otherColumn = 0; otherColumn < other._columns; ++otherColumn) {
            double sum = 0;
            for (Dimension column = 0; column < _columns; ++column) {
                sum += me(row, column) * other(column, otherColumn);
            }
            result(row, otherColumn) = sum;
        }
    }
    *this = result;
}

void Matrix::operator-=(const double other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
        _data[i] -= other;
    }
}

void Matrix::operator*=(const double other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
        _data[i] *= other;
    }
}

bool Matrix::operator==(const Matrix& other) const {
    if (!equalSize(other)) {
        return false;
    }
    for (Dimension i = 0; i < _arraySize; ++i) {
        if (abs(_data[i]- other._data[i]) > _epsilon) {
            return false;
        }
    }
    return true;
}

Matrix Matrix::adjoint() const {
    Matrix result(_rows, _columns);
    for (Dimension row = 0; row < _rows; ++row) {
        for (Dimension column = 0; column < _columns; ++column) {
            result(row, column) = cofactor(row, column);
        }
    }
    return result;
}

Matrix Matrix::adjugate() const {
    return adjoint().transpose();
}

double Matrix::cofactor(Dimension row, Dimension column) const {
    Matrix minor(_rows - 1, _columns - 1);
    for (Dimension subRow = 0; subRow < _rows; ++subRow) {
        for (Dimension subColumn = 0; subColumn < _columns; ++subColumn) {
            if (subRow < row && subColumn < column) {
                minor(subRow, subColumn) = me(subRow, subColumn);
            }
            if (subRow < row && subColumn > column) {
                minor(subRow, subColumn - 1) = me(subRow, subColumn);
            }
            if (subRow > row && subColumn < column) {
                minor(subRow - 1, subColumn) = me(subRow, subColumn);
            }
            if (subRow > row && subColumn > column) {
                minor(subRow - 1, subColumn - 1) = me(subRow, subColumn);
            }
        }
    }
    return minor.getDeterminant() * pow(-1, row + column);
}

Dimension Matrix::columns() const {
    return _columns;
}

Matrix Matrix::cubic() const {
    auto result = squared();
    result *= *this;
    return result;
}

double Matrix::getDeterminant() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
    if (_rows == 1) {
        return _data[0];
    }
    if (_rows == 2) {
        return _data[0] * _data[3] - _data[1] * _data[2];
    }
    double result = 0;
    for (Dimension row = 0; row < _rows; ++row) {
        Matrix minor(_rows - 1, _columns - 1);
        for (Dimension subRow = 1; subRow < _rows; ++subRow) {
            for (Dimension column = 0; column < _columns; ++column) {
                if (column < row) {
                    minor(subRow - 1, column) = me(subRow, column);  //_data[subRow * _columns + column];
                }
                else if (column > row) {
                    minor(subRow - 1, column - 1) = me(subRow, column);
                }
            }
        }
        result += _data[row] * minor.getDeterminant() * (row % 2 == 0 ? 1 : -1);
    }
    return result;
}

bool Matrix::contains(const double value, const double epsilon) const {
    for (const auto entry: _data) {
        if (fabs(entry - value) <= epsilon) return true;
    }
    return false;
}

Matrix Matrix::diagonal(const std::initializer_list<double> list) {
    const auto size = static_cast<Dimension>(list.size());
    Matrix result(size, size);
    Dimension i = 0;
    for (const auto value : list) {
        result(i, i) = value;
        i++;
    }
    return result;
}

Matrix Matrix::eigenvalues() const {
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
    // const double a3 = 1;
    const double a2 = -trace;
    const double a1 = (trace * trace - traceSquared) / 2.0;
    const double a0 = (-trace * trace * trace - 2 * cubic().getTrace() + 3 * trace * traceSquared) / 6.0; // = -determinant

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

Matrix Matrix::eigenvectorFor(const double lambda, const double epsilon) const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }

    Matrix aMinusLabdaI = *this - identity(_rows) * lambda;
    /* if (abs(aMinusLabdaI.getDeterminant()) < _epsilon) {
        std::cout << "Non invertible matrix";
    } */
    Matrix beta =  *this - identity(_rows) * lambda;

    // perform Gaussian elimination

    beta.toRowEchelonForm();

    /*for (Dimension pivotRow = 0; pivotRow < _rows; pivotRow++) {
        Dimension maxRow = pivotRow;
        for (Dimension subRow = pivotRow + 1; subRow < _rows; subRow++) {
            if (abs(beta(subRow, pivotRow)) > abs(beta(maxRow, pivotRow))) {
                maxRow = subRow;
            }
        }

        if (abs(beta(maxRow, pivotRow)) < epsilon) {
            continue;
        }

        beta.swapRows(pivotRow, maxRow);
            
        // normalize row
        const double pivot = beta(pivotRow, pivotRow);
        // <= is on purpose - we have an extra column in beta
        for (Dimension column = pivotRow; column <= _columns; column++) {
            beta(pivotRow, column) /= pivot;
        }

        // Set columns to 0 for the pivot by subtracting pivot row (with factor)
        for (Dimension row = 0; row < _rows; row++) {
            if (row == pivotRow) {
                continue;
            }
            const double factor = beta(row, pivotRow);
            for (Dimension column = pivotRow; column <= _columns; column++) {
                beta(row, column) -= factor * beta(pivotRow, column);
            }
        }
    } */

    auto result = beta.getColumn(_columns-1) * -1;
    result(_rows - 1, 0) = 1;
    return result; //.normalize()
}

bool Matrix::equalSize(const Matrix& other) const {
    return _rows == other._rows && _columns == other._columns;
}

Matrix Matrix::getColumn(const Dimension column) const {
    Matrix result(_rows, 1);
    for (Dimension row = 0; row < _rows; row++) {
        result(row, 0) = me(row, column);
    }
    return result;
}

Dimension Matrix::getRank() const {
    Matrix result = *this;
    result.toRowEchelonForm();
    Dimension rank = 0;
    for (Dimension row = 0; row < _rows; row++) {
        for (Dimension column = 0; column < _columns; column++) {
            if (abs(result(row, column)) > _epsilon) {
                rank++;
                break;
            }
        }
    }
    return rank;
}

Matrix Matrix::getRow(const Dimension row) const {
    Matrix result(1, _columns);
    for (Dimension column = 0; column < _columns; column++) {
        result(0, column) = me(row, column);
    }
    return result;
}

double Matrix::getTrace() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
    double result = 0;
    for (Dimension row = 0; row < _rows; ++row) {
        result += me(row, row);
    }
    return result;
}

Matrix Matrix::identity(Dimension size) {
    Matrix result(size, size);
    for (Dimension i = 0; i < size; ++i) {
        result(i, i) = 1;
    }
    return result;
}

Matrix Matrix::inverse() const {
    if (!isInvertible()) {
        throw std::invalid_argument("Matrix is not invertible");
    }
    return adjugate() / getDeterminant();
}

bool Matrix::isInvertible() const {
    return isSquare() && abs(getDeterminant()) > _epsilon;
}

bool Matrix::isSquare() const {
    return _rows == _columns;
}

double Matrix::me(const Dimension row, const Dimension column) const {
    return _data[row * _columns + column];
}

Matrix Matrix::normalize() const {
    double norm = 0;
    for (const auto entry : _data) {
        norm += entry * entry;
    }
    norm = sqrt(norm);
    if (abs(norm) > _epsilon) return *this / norm;
    return *this;
}

Dimension Matrix::rows() const {
    return _rows;
}

void Matrix::setAt(const Matrix& input, const Dimension startRow, const Dimension startColumn) {
    for (Dimension row = 0; row < input._rows; row++) {
        for (Dimension column = 0; column < input._columns; column++) {
            (*this)(startRow + row, startColumn + column) = input(row, column);
        }
    }
}

void Matrix::setAll(const double value) {
    for (double& entry : _data) {
        entry = value;
    }
}

void Matrix::setColumn(const Dimension column, const double value) {

    for (Dimension row = 0; row < _rows; row++) {
        (*this)(row, column) = value;
    }
}

void Matrix::setColumn(Dimension column, const Matrix& input) {
    for (Dimension row = 0; row < _rows; row++) {
        (*this)(row, column) = input(row, 0);
    }
}

void Matrix::setRow(const Dimension row, const Matrix& input) {
    for (Dimension column = 0; column < _columns; column++) {
        (*this)(row, column) = input(0, column);
    }
}

Matrix Matrix::squared() const {
    return *this * *this;
}

void Matrix::toRowEchelonForm() {
    Dimension lead = 0;
    for (Dimension row = 0; row < _rows; row++) {
        if (lead > _columns) {
            return;
        }
        Dimension pivotRow = row;
        while (abs(me(pivotRow, lead)) < _epsilon) {
            pivotRow++;
            if (pivotRow == _rows) {
                pivotRow = row;
                lead++;
                if (lead == _rows) {
                    return;
                }
            }
        }

        swapRows(pivotRow, row);

        const double leadingValue = me(row, lead);
        if (abs(leadingValue) < _epsilon) {
            // we have a zero row, so we need to swap it with a non-zero row
            for (Dimension nextRow = row + 1; nextRow < _rows; row++) {
                if (abs(me(nextRow, lead)) > _epsilon) {
                    swapRows(row, nextRow);
                    break;
                }
                // none found, move to next column
                if (nextRow == _rows - 1) {
                    lead++;
                    row--;
                    break;
                }
            }
         } else {
            // normalize the row using leading value
            for (Dimension column = 0; column < _columns; column++) {
                (*this)(row, column) /= leadingValue;
            }
         }
        // subtract the row from all other rows to get zeros in the column
        for (Dimension subRow = 0; subRow < _rows; subRow++) {
            if (subRow != row) {
                const double leadingValue2 = me(subRow, lead);
                for (Dimension column = 0; column < _columns; column++) {
                    (*this)(subRow, column) -= leadingValue2 * me(row, column);
                }
            }
        }
        lead++;
    }
}

Matrix Matrix::transpose() const {
    Matrix result(_columns, _rows);
    for (Dimension row = 0; row < _rows; ++row) {
        for (Dimension column = 0; column < _columns; ++column) {
            result(column, row) = me(row, column);
        }
    }
    return result;
}


Matrix operator/(Matrix left, const double right) {
    left /= right;
    return left;
}

Matrix operator*(Matrix left, const Matrix& right) {
    left *= right;
    return left;
}

Matrix operator*(Matrix left, const double right) {
    left *= right;
    return left;
}

Matrix operator+(Matrix left, const Matrix& right) {
    left += right;
    return left;
}

Matrix operator-(Matrix left, const Matrix& right) {
    left -= right;
    return left;
}