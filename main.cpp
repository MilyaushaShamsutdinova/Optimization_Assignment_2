#include <iostream>

#include "headers/Matrix.hpp"
#include "headers/Vector.hpp"
#include "headers/InteriorPoint.hpp"
#include "headers/SimplexMethod.hpp"

using namespace std;

void perform() {
    int n, m;
    cin >> n >> m;

    Vector b(m), c(n), trial_solution(m+n);
    Matrix A(m,n);
    double alpha, accuracy;

    cout << "Enter coefficients of the objective function:\n";
    cin >> c;
    cout <<"Enter the coefficients matrix of the constraints:\n";
    cin >> A;
    cout << "Enter the right hand side of the constraints:\n";
    cin >> b;
    cout << "Enter the trial solution:\n";
    //cin >> trial_solution;
//    cout << "Enter the value of alpha:\n";
//    cin>> alpha;
    cout << "Enter the value of accuracy needed in calculations:\n";
    cin >> accuracy;

//    cout << "\n\n------------------------------------------------";
//    cout << "\n----- Interior Point Method with alpha=0.5 -----\n";
//    cout << "------------------------------------------------\n\n";
//    alpha = 0.5;
//    InteriorPoint::start_interior_point(A, b, c, trial_solution, alpha, accuracy);
//
//    cout << "\n\n------------------------------------------------";
//    cout << "\n----- Interior Point Method with alpha=0.9 -----\n";
//    cout << "------------------------------------------------\n\n";
//    alpha = 0.9;
//    InteriorPoint::start_interior_point(A, b, c, trial_solution, alpha, accuracy);

    cout << "\n\n------------------------------------------------";
    cout << "\n---------------- Simplex Method ----------------\n";
    cout << "------------------------------------------------\n\n";
    SimplexMethod::start_simplex(A, b, c, accuracy);

}

int main() {
    perform();
    return EXIT_SUCCESS;
}