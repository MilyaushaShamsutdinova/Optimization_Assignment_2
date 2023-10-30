//
// Created by ilia on 17.09.23.
//

#include "../headers/Vector.hpp"

Vector::Vector() {
    vector_ = vector<float>();
}

Vector::Vector(int vector_size, float initial_value) {
    vector_ = vector<float>(vector_size, initial_value);
    for (int i = 0; i < vector_size; i++) {
        vector_[i] = initial_value;
    }
}

Vector::Vector(const Vector &other) {
    if (this == &other) {
        return;
    }

    vector_ = vector<float>(other.size(), 0.0f);

    for (int i = 0; i < vector_.size(); i++) {
        vector_[i] = other.vector_[i];
    }
}

Vector::~Vector() {
    vector_ = vector<float>(0);
}

int Vector::size() const {
    return (int) vector_.size();
}

float &Vector::operator[](int index) {
    if (index < 0 || index >= vector_.size()) {
        throw out_of_range("Index out of range. Vector &operator[]");
    }
    return vector_[index];
}

float Vector::operator[](int index) const {
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
        for (int i = 0; i < other.size(); i++) {
            vector_[i] -= other[i];
        }
    } else {
        cout << "Vectors must be of the same size." << endl;
    }
    return *this;
}

Vector Vector::operator*(const float scalar) const {
    Vector tmp(*this);
    tmp *= scalar;
    return tmp;
}

Vector &Vector::operator*=(const float scalar) {
    for (float &i: vector_) {
        i *= scalar;
    }
    return *this;
}

Vector Vector::operator/(const float scalar) const {
    Vector tmp(*this);
    tmp /= scalar;
    return tmp;
}

Vector &Vector::operator/=(const float scalar) {
    if (scalar == 0) {
        cerr << "Division by zero is not allowed." << endl;
        return *this;
    }

    for (float &i: vector_) {
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
        Vector tmp(other);
        for (int i = 0; i < other.size(); ++i) {
            tmp[i] += other[i];
        }
        swap(*this, tmp);
    } else {
        cout << "Vectors must be of the same size." << endl;
    }
    return *this;
}
