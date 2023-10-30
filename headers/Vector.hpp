//
// Created by ilia on 17.09.23.
//

#ifndef OPTIMIZATION_ASSIGNMENT_2_VECTOR_HPP
#define OPTIMIZATION_ASSIGNMENT_2_VECTOR_HPP


#include <iostream>
#include <vector>


using namespace std;


class Vector {
public:
    Vector();

    explicit Vector(int vector_size, double initial_value = 0.0000000000f);

    Vector(const Vector &other);

    ~Vector();

    [[nodiscard]] int size() const;

    double &operator[](int index);

    double operator[](int index) const;

    Vector operator-(Vector &other) const;

    Vector &operator-=(Vector &other);

    Vector operator*(double scalar) const;

    Vector &operator*=(double scalar);

    Vector operator/(double scalar) const;

    Vector &operator/=(double scalar);

    Vector &operator=(Vector other);

    Vector operator+(Vector &other) const;

    Vector &operator+=(Vector &other);

    //overloaded ostream operator
    friend std::ostream &operator<<(std::ostream &os, const Vector &item) {
        for (int i = 0; i < item.size(); ++i) {
            os << item[i] << " ";
        }
        cout << "\n";
        return os;
    }

    //overloaded istream operator
    friend std::istream &operator>>(std::istream &input_stream, Vector &item) {
        for (int i = 0; i < item.size(); ++i) {
            input_stream >> item[i];
        }
        if (!input_stream) {
            throw std::runtime_error("Invalid input (matrix input error)\n");
        }
        return input_stream;
    }

private:
    vector<double> vector_;

    // resize function since the standard arrays used for storing data
    void resize();


    static void swap(Vector &first, Vector &second);
};


#endif //OPTIMIZATION_ASSIGNMENT_2_VECTOR_HPP