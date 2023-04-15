#include "Matrix.h"

Matrix::Matrix(Dimension rows, Dimension columns) : _data(rows * columns) {
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

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> list) {
    _rows = list.size();
    _columns = list.begin()->size();
    _arraySize = _rows * _columns;
    _data = std::vector<double>(_rows * _columns);
    Dimension row = 0;
    for (auto rowList : list) {
        Dimension column = 0;
        for (auto value : rowList) {
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

void Matrix::operator+=(double other){
    for (Dimension i = 0; i < _arraySize; ++i) {
            _data[i] += other;
    }
}

void Matrix::operator-=(const Matrix &other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
            _data[i] -= other._data[i];
    }
}

void Matrix::operator/=(double other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
            _data[i] /= other;
    }
}

void Matrix::operator*=(const Matrix &other) {
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

void Matrix::operator-=(double other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
            _data[i] -= other;
    }
}

void Matrix::operator*=(double other) {
    for (Dimension i = 0; i < _arraySize; ++i) {
            _data[i] *= other;
    }
}

bool Matrix::operator==(const Matrix &other) const {
    if (!equalSize(other)) {
        return false;
    }
    for (Dimension i = 0; i < _arraySize; ++i) {
        if (_data[i] != other._data[i]) {
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
    return minor.determinant() * pow(-1, row + column);
}

Dimension Matrix::columns() const {
    return _columns;
}

double Matrix::determinant() const {
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
                } else if (column > row) {
                    minor(subRow - 1, column - 1) = me(subRow, column);
                }
            }
        }
        result += _data[row] * minor.determinant() * (row % 2 == 0 ? 1 : -1);
    }
    return result;
}

Matrix Matrix::diagonal(std::initializer_list<double> list) {
    Matrix result(list.size(), list.size());
    Dimension i = 0;
    for (auto value : list) {
        result(i, i) = value;
        i++;
    }
    return result;
}

Matrix Matrix::eigenvectorFor(double lambda) const {
    const double EPSILON = 1e-10;
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
    if (_rows > 3) {
        throw std::invalid_argument("Only works for 3d matrices");
    }
    Matrix beta = *this - identity(_rows) * lambda;

    // perform Gaussian elimination
    for (int pivotRow = 0; pivotRow < 3; pivotRow++) {
        double pivot = beta(pivotRow, pivotRow);
        if (fabs(pivot) < EPSILON) {
            // matrix is singular, can't find eigenvectors
            return Matrix(0, 0);
        }
        for (int row = pivotRow + 1; row < 3; row++) {
            double factor = beta(row, pivotRow) / pivot;
            for (int column = pivotRow; column < 3; column++) {
                beta(row, column) -= factor * beta(pivotRow, column);
            }
        }
    }

    // back-subsitution
    Matrix result(3, 1);
    for (int row = 2; row >= 0; row--) {
        double sum = 0;
        for (int column = row + 1; column < 3; column++) {
            sum += beta(row, column) * result(column, 0);
        }
        result(row, 0) = (beta(row, 3) - sum) / beta(row, row);
    }
    return result;
}


bool Matrix::equalSize(const Matrix& other) const {
    return _rows == other._rows && _columns == other._columns;
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
    return adjugate() / determinant();
}

bool Matrix::isInvertible() const {
    return isSquare() && determinant() != 0;
}

Matrix Matrix::eigenvalues() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
    if (_rows > 3) {
        throw std::invalid_argument("Only works for 3d matrices");
    }
    if (_rows == 1) {
        return Matrix({{me(0, 0)}});
    } 
    if (_rows == 2) {
        double a = _data[0];
        double b = _data[1];
        double c = _data[2];
        double d = _data[3];
        double discriminant = b * b - 4 * a * d;
        if (discriminant < 0) {
            // no real eigenvalues
            return Matrix(0, 0);
        }
        return Matrix{
            {(-b + sqrt(discriminant)) / (2 * a)},
            {(-b - sqrt(discriminant)) / (2 * a)}};
    }
    // 3x3 matrix
    double b = trace();
    double c = (1/2.0) * (pow(b, 2.0) - squared().trace());
    double d = determinant();
    double q = (pow(b, 2.0) - 3*c) / 9.0;
    double r = (2*pow(b, 3.0) - 9*b*c + 27*d) / 54.0;
    double s = pow(q, 3.0) - pow(r, 2.0);

    if (s > 0) {
        // Three distinct real roots
        double theta = acos(r / sqrt(pow(q, 3.0)));
        return Matrix({
            { -2 * sqrt(q) * cos(theta / 3.0) - b / 3.0 },
            { -2 * sqrt(q) * cos((theta + 2 * M_PI) / 3.0) - b / 3.0 },
            { -2 * sqrt(q) * cos((theta- 2 * M_PI) / 3.0) - b / 3.0 }
        });
    } 
    if (s < 0) {
        // One real root and two complex conjugate roots
        double rho = sqrt(pow(q, 3.0));
        double theta = acos(r / rho);
        return Matrix({
            { -2 * pow(rho, 1/3.0) * cos(theta/3.0) - b / 3.0 },
            { pow(rho, 1/3.0) * (cos(theta/3.0) + sqrt(3)*sin(theta/3.0)) - b / 3.0 },
            { pow(rho, 1/3.0) * (cos(theta/3.0) - sqrt(3)*sin(theta/3.0)) - b / 3.0 }
        });
    }
    // One real root and two identical real roots
    auto identicalRoot = pow(r, 1/3.0) - b / 3.0;
    return Matrix({
        { -2 * pow(r, 1/3.0) - b / 3.0 },
        { identicalRoot }, 
        { identicalRoot }
    });
}

bool Matrix::isSquare() const {
    return _rows == _columns;
}

double Matrix::me(Dimension row, Dimension column) const {
    return _data[row * _columns + column];
}

//const double Matrix::me(Dimension row, Dimension column) const {
//    return _data[row * _columns + column];
//}

Dimension Matrix::rows() const {
    return _rows;
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

Matrix Matrix::squared() const{
    return *this * *this;
}

double Matrix::trace() const {
    if (!isSquare()) {
        throw std::invalid_argument("Matrix is not square");
    }
    double result = 0;
    for (Dimension row = 0; row < _rows; ++row) {
        result += me(row, row);
    }
    return result;
}

Matrix operator/(Matrix left, double right) {
    left /= right;
    return left;
}

Matrix operator*(Matrix left, const Matrix &right) {
    left *= right;
    return left;    
}

Matrix operator*(Matrix left, double right) {
    left *= right;
    return left;
}

Matrix operator+(Matrix left, const Matrix &right) {
    left += right;
    return left;
}

Matrix operator-(Matrix left, const Matrix& right) {
    left -= right;
    return left;
}