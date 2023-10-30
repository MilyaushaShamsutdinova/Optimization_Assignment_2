//
// Created by ilia on 17.09.23.
//

#include "../headers/Matrix.hpp"

using namespace std;

Matrix::Matrix() {
    rows_ = 0;
    columns_ = 0;
    data_ = vector<vector<double>>();
}


Matrix::Matrix(int rows, int columns) {
    rows_ = rows;
    columns_ = columns;

    data_ = vector<vector<double>>(rows, vector<double>(columns, 0.00000000000f));
}


Matrix::Matrix(const Matrix &other) {
    rows_ = other.rows_;
    columns_ = other.columns_;

    data_ = vector<vector<double>>(rows_, vector<double>(columns_, 0.00000000000f));

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++) {
            data_[i][j] = other.data_[i][j];
        }
    }
}

Matrix::~Matrix() {
    rows_ = 0;
    columns_ = 0;
    data_ = vector<vector<double>>(0);
}

[[nodiscard]] int Matrix::rows() const {
    return rows_;
}

[[nodiscard]] int Matrix::columns() const {
    return columns_;
}

double &Matrix::operator()(int row, int column) {
    if (row < 0 || row >= rows_ || column < 0 || column >= columns_) {
        throw out_of_range("Index out of range. Matrix &operator()");
    }
    return data_[row][column];
}

double Matrix::operator()(int row, int column) const {
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


Vector Matrix::operator*(Vector vector) {
    if (columns_ != vector.size()) {
        throw invalid_argument("Matrix and vector dimensions are incompatible for multiplication.");
    }

    Vector result(rows_);

    for (int i = 0; i < rows_; ++i) {
        double sum = 0.0;
        for (int j = 0; j < columns_; ++j) {
            sum += data_[i][j] * vector[j];
        }
        result[i] = sum;
    }

    return result;
}


Matrix& Matrix::operator*=(const Matrix& item) {
    if (columns_ == item.rows_) {
        Matrix tmp_matrix_(rows_, item.columns_);
        for (int i = 0; i < rows_; i++) {
            for (int k = 0; k < item.rows_; k++) {
                for (int j = 0; j < item.columns_; j++) {
                    tmp_matrix_(i, j) += (data_[i][k] * item(k, j));
                }
            }
        }
        swap(tmp_matrix_, *this);
    } else {
        cout << "Operation * of wrong dimensions for matrix multiplication" << endl;
    }
    return *this;
}

Matrix Matrix::operator-(const Matrix& item) const {
    Matrix tmp_matrix_(*this);
    tmp_matrix_ -= item;
    return tmp_matrix_;
}

Matrix& Matrix::operator-=(const Matrix& item) {
    if (rows_ == item.rows_ && columns_ == item.columns_) {
        Matrix tmp_matrix_(rows_, item.columns_);
        for (int i = 0; i < rows_; i++){
            for (int j = 0; j < rows_; ++j) {
                tmp_matrix_(i, j) = data_[i][j] - item(i, j);
            }
        }
        swap(tmp_matrix_, *this);
    } else {
        cout << "Operation - of wrong dimensions for matrix subtraction" << endl;
    }
    return *this;
}

Matrix Matrix::operator*(const Matrix& item) const {
    Matrix tmp_matrix_(*this);
    tmp_matrix_ *= item;
    return tmp_matrix_;
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


Matrix Matrix::transpose() {
    Matrix transposed(columns_, rows_);

    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < columns_; ++j) {
            transposed(j, i) = data_[i][j];
        }
    }

    return transposed;
}

Matrix Matrix::inverse() {
    Matrix inverse_matrix(rows_, columns_);
    Matrix tmp_matrix(rows_, columns_ * 2);

    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < columns_; j++) {
            tmp_matrix(i, j) = data_[i][j];
        }

    }
    for (int i = 0; i < rows_; ++i) {
        for (int j = columns_; j < columns_*2; ++j) {
            tmp_matrix(i, j) = 0;
        }
        tmp_matrix(i, columns_ + i) = 1;
    }

    for (int i = 0; i < rows_; i++) {
        double tmp = tmp_matrix(i, i);
        for (int j = 0; j < columns_*2; j++) {
            tmp_matrix(i, j) /= tmp;
        }
        for (int j = 0; j < rows_; j++) {
            if (i != j) {
                tmp = tmp_matrix(j, i);
                for (int k = 0; k < columns_*2; k++) {
                    tmp_matrix(j, k) -= tmp_matrix(i, k) * tmp;
                }
            }
        }
    }

    for (int i = columns_; i < columns_*2; ++i) {
        for (int j = 0; j < rows_; ++j) {
            inverse_matrix(j, i - columns_) = tmp_matrix(j,i);
        }
    }

    return inverse_matrix;
}
