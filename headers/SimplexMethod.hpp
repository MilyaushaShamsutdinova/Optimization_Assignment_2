//
// Created by ilia on 17.09.23.
//

#ifndef OPTIMIZATION_COURSE_ASSIGNMENT_1_SIMPLEXMETHOD_HPP
#define OPTIMIZATION_COURSE_ASSIGNMENT_1_SIMPLEXMETHOD_HPP


#include "Vector.hpp"
#include "Matrix.hpp"

using namespace std;

class SimplexMethod {
public:

    static void start_simplex(const Matrix &A, const Vector &B, const Vector &C, double accuracy);

private:

    static void initialize_algorithm_data(const Matrix &A, const Vector &B, const Vector &C, Matrix &main_matrix,
                                          Vector &func_coefficients, Vector &net_eval);

    static int define_pivot_col(Vector &net_eval);

    static Vector calculate_ratio(Matrix &main_matrix, int &pivot_col, double accuracy);

    static int define_pivot_row(Vector &ratio);

    static double define_pivot_element(Matrix &main_matrix, int &pivot_col, int &pivot_row);

    static Vector define_basis(Vector &basis, int &pivot_row, int &pivot_col, Vector &function_coefficients);

    static Matrix
    update_main_matrix(Matrix &main_matrix, int &pivot_row, int &pivot_col, double pivot_element, double accuracy);

    static Vector define_basis_element(Vector &basis, Vector &basis_el, int &pivot_row, int &pivot_col);

    static Vector calculate_profit(Matrix &main_matrix, Vector &basis, double accuracy);

    static Vector calculate_net_evaluation(Vector &function_coefficients, Vector &profit, double accuracy);

    static bool check_net_evaluation(Vector &net_eval);

    // methods for rounding
    static double rounding(double epsilon, double variable);

    static Vector rounding(double epsilon, Vector &variable);

    static Matrix rounding(double epsilon, Matrix &variable);

    static bool check_data(const Matrix &A, const Vector &B, const Vector &C, double epsilon);

    static Vector form_x_vector(Vector &x, Vector &basis_el, Matrix &main_matrix);
};

#endif //OPTIMIZATION_COURSE_ASSIGNMENT_1_SIMPLEXMETHOD_HPP