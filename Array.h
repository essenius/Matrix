#ifndef ARRAY_H
#define ARRAY_H

#include <vector>

typedef unsigned int Dimension;

/// Class for array manipulations (coefficient wise)
class Array {
public:
    Array(Dimension rows, Dimension columns);
    Array(const std::initializer_list<std::initializer_list<double>> list);
    Array(const Array& other);

    double& operator()(Dimension row, Dimension column);
    const double& operator()(Dimension row, Dimension column) const;

    void operator+=(const Array& other);
    void operator+=(double other);

    void operator-=(const Array& other);
    void operator-=(double other);

    void operator*=(const Array& other);
    void operator*=(double other);

    void operator/=(double other);

    bool operator==(const Array &other) const;

    void assertMatchingSize(const Array& other) const;
    Dimension columns() const;
    bool equalSize(const Array& other) const;
    Array getColumn(const Dimension column) const;
    Array getRow(const Dimension row) const;
    bool isSquare() const;
    double me(Dimension row, Dimension column) const;

    Array pow2() const;
    Dimension rows() const;

    void setColumn(Dimension column, const Array & input);
    void setColumn(const Dimension column, const double value);
    void setColumnCount(const Dimension columns);

    friend Array operator+(Array left, const Array& right);
    friend Array operator-(Array left, const Array& right);
    friend Array operator*(Array left, const Array& right);
    friend Array operator*(Array left, double right);
    friend Array operator/(Array left, double right);

protected:
    std::vector<double> _data;
    static const double _epsilon;
    Dimension _rows;
    Dimension _columns;
    Dimension _arraySize;
};
#endif