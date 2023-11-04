#include "../headers/InteriorPoint.hpp"

#include <cmath>

void InteriorPoint::start_interior_point(const Matrix& A, const Vector& b, const Vector& c,
                                         const Vector& trial_solution, double alpha, double accuracy) {
    bool lpp_is_solvable = check_data(A, b, c, alpha, accuracy);

    if (!lpp_is_solvable){
        return;
    }

    Matrix main_matrix;
    Vector func_coefficients;
    Matrix D;

    initialize_algorithm_data(A, b, c, main_matrix, func_coefficients);

    D = Matrix(main_matrix.columns(), main_matrix.columns());

    Vector x = trial_solution;

    Matrix A_tilda;
    Vector c_tilda, c_p, x_tilda, x_unit(func_coefficients.size(), 1.000000f);
    Matrix P;
    Matrix I(D.rows(), D.columns());

    for (int i = 0; i < min(I.rows(),I.columns()); ++i) {
        I(i, i) = 1.000000f;
    }

    Vector x_last;
    int iter = 0;
    while (true) {
        ++iter;
        for (int i = 0; i < D.rows(); ++i) {
            for (int j = 0; j < D.columns(); ++j) {
                D(i, j) = 0.000000f;
            }
            D(i, i) = x[i];
        }

        A_tilda = main_matrix * D;
        c_tilda = D * func_coefficients;

        Matrix A_tilda_transpose = A_tilda.transpose();
        Matrix A_square = A_tilda * A_tilda_transpose;
        Matrix A_inverse = A_square.inverse();

        Matrix A_tmp = A_tilda_transpose * A_inverse;
        A_tmp *= A_tilda;

        P = I - A_tmp;

        c_p = P * c_tilda;

        double v, m = 1.000000f;
        for (int i = 0; i < c_p.size(); ++i) {
            if (c_p[i] < 0.000000f) {
                m = min(m, c_p[i]);
                v = abs(m);
            }
        }

        if (m > 0.000000f) {
            break;
        }

        Vector tmp = c_p * (alpha / v);
        x_tilda = x_unit + tmp;

        x = D * x_tilda;

        if (iter == 1){
            x_last = Vector(x);
            continue;
        }

        double dist = 0;
        for (int i = 0; i < x.size(); ++i) {
            dist += (x_last[i] - x[i]) * (x_last[i] - x[i]);
        }
        if (abs(sqrt(dist) - accuracy) < accuracy) {
            break;
        }

        x_last = Vector(x);
        double profit = 0.000000f;

        for (int i = 0; i < A.columns(); ++i) {
            profit += c[i] * x[i];
        }
    }

    cout << iter << " iterations for alpha=" << alpha << "\n";
    cout << "Answer:\n";
    for (int i = 0; i < c.size(); ++i) {
        cout << "x_" << i+1 << " = " << x[i] << "\n";
    }
    double profit = 0.000000f;

    for (int i = 0; i < A.columns(); ++i) {
        profit += c[i] * x[i];
    }

    cout << "Profit = " << profit << endl;

}

void InteriorPoint::initialize_algorithm_data(const Matrix &A, const Vector &B, const Vector &C, Matrix &main_matrix,
                                              Vector &func_coefficients) {
    main_matrix = Matrix(A.rows(), A.columns() + B.size());

    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.columns(); ++j) {
            main_matrix(i, j) = A(i, j);
        }
        main_matrix(i, A.columns() + i) = 1.000000f;
    }

    func_coefficients = Vector(main_matrix.columns(), 0.000000);
    for (int i = 0; i < C.size(); ++i) {
        func_coefficients[i] = C[i];
    }

}


bool InteriorPoint::check_data(const Matrix &A, const Vector &b,const Vector &c, double alpha, double accuracy){

    if (alpha > 1.0 || alpha < 0.0){
        cout << "Incorrect input for alpha!\n";
        return false;
    }

    for (int i = 0; i < b.size(); ++i) {
        if (isless(b[i], accuracy - accuracy)) {
            cout << "The problem does not have solution!\n";
            return false;
        }
    }

    //Check for sparse constraint matrix
    int nonZeroCount = 0;
    double threshold = 0.6;
    int totalElements = A.rows() * A.columns();

    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.columns(); ++j) {
            if (A(i, j) != 0.0) {
                nonZeroCount++;
            }
        }
    }
    double sparsityRatio = static_cast<double>(nonZeroCount) / totalElements;

    if (sparsityRatio < threshold){
        cout << "The Interior-Point Method is not applicable!\n";
        return false;
    }

    return true;
}