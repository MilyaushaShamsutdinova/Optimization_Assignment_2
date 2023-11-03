//
// Created by ilia on 17.09.23.
//

#include "../headers/Vector.hpp"

#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;

Vector::Vector() {
    vector_ = vector<double>();
}

Vector::Vector(int vector_size, double initial_value) {
    vector_ = vector<double>(vector_size, initial_value);
    for (int i = 0; i < vector_size; i++) {
        vector_[i] = initial_value;
    }
}

Vector::Vector(const Vector &other) {
    if (this == &other) {
        return;
    }

    vector_ = vector<double>(other.size(), 0.0f);

    for (int i = 0; i < vector_.size(); i++) {
        vector_[i] = other.vector_[i];
    }
}

Vector::~Vector() {
    vector_ = vector<double>(0);
}

int Vector::size() const {
    return (int) vector_.size();
}

double &Vector::operator[](int index) {
    if (index < 0 || index >= vector_.size()) {
        throw out_of_range("Index out of range. Vector &operator[]");
    }
    return vector_[index];
}

double Vector::operator[](int index) const {
    if (index < 0 || index >= vector_.size()) {
        throw out_of_range("Index out of range. Vector operator[]");
    }
    return vector_[index];
}

Vector Vector::operator-(Vector &other) const {
    Vector tmp(*this);
    tmp -= other;
    return tmp;
}

Vector &Vector::operator-=(Vector &other) {
    if (other.size() == vector_.size()) {
        Vector tmp(vector_.size());
        for (int i = 0; i < other.size(); i++) {
            tmp[i] = vector_[i] - other[i];
        }
        swap(tmp, *this);
    } else {
        int size = min((long long)vector_.size(), (long long)other.size());
        Vector tmp(size);
        for (int i = 0; i < size; i++) {
            tmp[i] = vector_[i] - other[i];
        }
        swap(tmp, *this);
    }
    return *this;
}

Vector Vector::operator*(const double scalar) const {
    Vector tmp(*this);
    tmp *= scalar;
    return tmp;
}

Vector &Vector::operator*=(const double scalar) {
    for (double &i: vector_) {
        i *= scalar;
    }
    return *this;
}

Vector Vector::operator/(const double scalar) const {
    Vector tmp(*this);
    tmp /= scalar;
    return tmp;
}

Vector &Vector::operator/=(const double scalar) {
    if (scalar == 0) {
        cerr << "Division by zero is not allowed." << endl;
        return *this;
    }

    for (double &i: vector_) {
        i /= scalar;
    }
    return *this;
}

void Vector::resize() {
    vector_.resize(vector_.size() * 2);
}

Vector &Vector::operator=(Vector other) {
    if (this == &other) {
        return *this;
    }
    Vector temp(other);
    swap(*this, temp);

    return *this;
}

void Vector::swap(Vector &first, Vector &second) {
    std::swap(first.vector_, second.vector_);
}

Vector Vector::operator+(Vector& other) const {
    Vector tmp(*this);
    tmp += other;
    return tmp;
}

Vector& Vector::operator+=(Vector& other) {
    if (other.size() == vector_.size()) {
        Vector tmp(*this);
        for (int i = 0; i < other.size(); ++i) {
            tmp.vector_[i] += other.vector_[i];
        }
        swap(tmp, *this);
    } else {
        cout << "Vectors must be of the same size." << endl;
        //cout << other.size() << " " << vector_.size() << endl;
    }
    return *this;
}
