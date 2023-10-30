
#ifndef OPTIMIZATION_ASSIGNMENT_2_INTERIORPOINT_HPP
#define OPTIMIZATION_ASSIGNMENT_2_INTERIORPOINT_HPP

#include <iostream>

#include "Matrix.hpp"
#include "Vector.hpp"

using namespace std;

class InteriorPoint {
public:
    static void start_interior_point(const Matrix& A, const Vector& b, const Vector& c,
                                     const Vector& trial_solution, float alpha, float accuracy);

private:
    static bool check_data(const Matrix &A, const Vector &b, const Vector &c, float alpha, float accuracy);

    static void initialize_algorithm_data(const Matrix &A, const Vector &B, const Vector &C, Matrix &main_matrix,
                                          Vector &func_coefficients);
};

#endif //OPTIMIZATION_ASSIGNMENT_2_INTERIORPOINT_HPP