//
// Created by ilia on 17.09.23.
//

#include "../headers/Matrix.hpp"

using namespace std;

Matrix::Matrix() {
    rows_ = 0;
    columns_ = 0;
    data_ = vector<vector<float>>();
}


Matrix::Matrix(int rows, int columns) {
    rows_ = rows;
    columns_ = columns;

    data_ = vector<vector<float>>(rows, vector<float>(columns, 0.0f));
}


Matrix::Matrix(const Matrix &other) {
    rows_ = other.rows_;
    columns_ = other.columns_;

    data_ = vector<vector<float>>(rows_, vector<float>(columns_, 0.0f));

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++) {
            data_[i][j] = other.data_[i][j];
        }
    }
}

Matrix::~Matrix() {
    rows_ = 0;
    columns_ = 0;
    data_ = vector<vector<float>>(0);
}

[[nodiscard]] int Matrix::rows() const {
    return rows_;
}

[[nodiscard]] int Matrix::columns() const {
    return columns_;
}

float &Matrix::operator()(int row, int column) {
    if (row < 0 || row >= rows_ || column < 0 || column >= columns_) {
        throw out_of_range("Index out of range. Matrix &operator()");
    }
    return data_[row][column];
}

float Matrix::operator()(int row, int column) const {
    if (row < 0 || row >= rows_ || column < 0 || column >= columns_) {
        throw out_of_range("Index out of range. Matrix &operator()");
    }
    return data_[row][column];
}

Matrix &Matrix::operator=(const Matrix &other) {
    if (this == &other) {
        return *this;
    }
    Matrix temp(other);
    swap(*this, temp);

    return *this;
}


void Matrix::swap(Matrix &first, Matrix &second) {
    std::swap(first.rows_, second.rows_);
    std::swap(first.columns_, second.columns_);
    std::swap(first.data_, second.data_);
}


Vector Matrix::operator*(Vector vector) const {
    if (columns_ != vector.size()) {
        throw invalid_argument("Matrix and vector dimensions are incompatible for multiplication.");
    }

    Vector result(rows_);

    for (int i = 0; i < rows_; ++i) {
        float sum = 0.0;
        for (int j = 0; j < columns_; ++j) {
            sum += data_[i][j] * vector[j];
        }
        result[i] = sum;
    }

    return result;
}


Vector Matrix::getRow(int index) const {
    if (index < 0 || index >= rows_) {
        throw out_of_range("Row index is out of range.");
    }

    Vector rowVector(columns_);

    for (int i = 0; i < columns_; ++i) {
        rowVector[i] = data_[index][i];
    }

    return rowVector;
}

void Matrix::setRow(int index, Vector &vector) {
    if (index < 0 || index >= rows_ || vector.size() != columns_) {
        throw invalid_argument("Invalid row index or vector dimensions.");
    }

    for (int i = 0; i < columns_; ++i) {
        data_[index][i] = vector[i];
    }
}

Vector Matrix::getCol(int index) const {
    if (index < 0 || index >= columns_) {
        throw out_of_range("Column index is out of range.");
    }

    Vector colVector(rows_);

    for (int i = 0; i < rows_; ++i) {
        colVector[i] = data_[i][index];
    }

    return colVector;
}


void Matrix::setCol(int index, Vector &vector) {
    if (index < 0 || index >= columns_ || vector.size() != rows_) {
        throw invalid_argument("Invalid column index or vector dimensions.");
    }

    for (int i = 0; i < rows_; ++i) {
        data_[i][index] = vector[i];
    }
}


Matrix Matrix::transpose() const {
    Matrix transposed(columns_, rows_);

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < columns_; ++j) {
            transposed(j, i) = data_[i][j];
        }
    }

    return transposed;
}
