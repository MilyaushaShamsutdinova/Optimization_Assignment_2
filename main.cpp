#include <iostream>

#include "headers/Matrix.hpp"
#include "headers/Vector.hpp"
//#include "headers/SimplexMethod.hpp"

using namespace std;

void perform() {
    int n, m;
    cin >> n >> m;

    Vector b(m), c(n);
    Matrix A(m,n);
    float alpha, accuracy;
    cin >> c >> A >> b >> alpha >> accuracy;

    //SimplexMethod::start_simplex(A, b, C, accuracy);
    // run the interior point algorithm
}

int main() {
    perform();
    return EXIT_SUCCESS;
}