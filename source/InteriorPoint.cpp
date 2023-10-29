#include "../headers/InteriorPoint.hpp"

#include <cmath>

void InteriorPoint::start_interior_point(const Matrix& A, const Vector& b, const Vector& c,
                                         const Vector& trial_solution, float alpha, float accuracy) {
    bool is_feasible = check_data(A, b, c, alpha, accuracy);

    if (!is_feasible){
        cout << "There is no feasible solution!\n";
        return;
    }

    Matrix main_matrix;
    Vector func_coefficients;
    Matrix D;

    initialize_algorithm_data(A, b, c, main_matrix, func_coefficients, D, trial_solution);

    Vector x(func_coefficients.size(), accuracy);

    Matrix A_tilda;
    Vector c_tilda, c_p, x_tilda, x_iter, x_unit(func_coefficients.size(), 1.000000f);
    Matrix P;
    Matrix I(D.rows(), D.columns());

    for (int i = 0; i < I.rows(); ++i) {
        I(i, i) = 1.000000f;
    }

    int iter = 0;
    while (true) {
        ++iter;
        if (iter != 1) {
            for (int i = 0; i < D.rows(); ++i) {
                for (int j = 0; j < D.columns(); ++j) {
                    D(i, j) = 0.000000f;
                }
                D(i, i) = x[i];
            }
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

        float v, m;
        for (int i = 0; i < c_p.size(); ++i) {
            m = min(c_p[0], c_p[i]);
            v = abs(m);
        }

        if (m > 0.000000f) {
            cout << "In iteration " << iter << " no negative values V\n";
            break;
        }

        c_p *= alpha / v;
        x_tilda = x_unit + c_p;

        x = D * x_tilda;

        cout << "Iteration: " << iter << " x= " << x << endl;

        bool flag = false;
        for (int i = 0; i < x.size(); ++i) {
            if (x[i] < accuracy - accuracy) {
                flag = true;
                break;
            }
        }
        if (!flag) {
            break;
        }
    }

    cout << "The answer is:\n";
    for (int i = 0; i < x.size(); ++i) {
        cout << "x" << i << " = " << x[i] << " ";
    }
    float profit = 0.000000f;

    for (int i = 0; i < A.columns(); ++i) {
        profit += c[i] * x[i];
    }

    cout << "\n\nProfit of these x values: " << profit << endl;

}

void InteriorPoint::initialize_algorithm_data(const Matrix &A, const Vector &B, const Vector &C, Matrix &main_matrix,
                                              Vector &func_coefficients, Matrix &D, const Vector &trivial_solution) {


    main_matrix = Matrix(A.rows(), A.columns() + B.size());
    D = Matrix(main_matrix.columns(), main_matrix.columns());

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

    for (int i = 0; i < A.rows(); ++i) {
        D(i, i) = 1.000000f;
    }

    for (int i  = A.rows(); i < D.rows(); ++i) {
        D(i, i) = trivial_solution[i - A.rows()];
    }
}


bool InteriorPoint::check_data(const Matrix &A, const Vector &b,const Vector &c, float alpha, float accuracy){

    for (int i = 0; i < b.size(); ++i) {
        if (isless(b[i], accuracy - accuracy)) {
            cout << "There is no feasible solution!";
            return false;
        }
    }
    return true;
}