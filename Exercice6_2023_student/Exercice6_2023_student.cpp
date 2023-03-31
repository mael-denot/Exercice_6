#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ConfigFile.tpp"


using namespace std;

// Résolution d'un système d'équations linéaires par élimination de
// Gauss-Jordan:
template<class T>
vector<T>
solve(const vector<T>& diag,
      const vector<T>& lower,
      const vector<T>& upper,
      const vector<T>& rhs)
{
    vector<T> solution(diag.size());
    vector<T> new_diag(diag);
    vector<T> new_rhs(rhs);

    for (int i = 1; i < diag.size(); ++i) {
        if (new_diag[i - 1] == 0.0) {
            throw std::runtime_error("solve: division by zero");
        }
        double pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    if (new_diag[diag.size() - 1] == 0.0) {
        throw std::runtime_error("solve: division by zero");
    }
    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    for (int i(diag.size() - 2); i >= 0; --i)
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];

    //test if solution is nan
    for (int i = 0; i < solution.size(); ++i) {
        if (solution[i] != solution[i]) {
            throw std::runtime_error("solve: solution is nan");
        }
    }

    return solution;
}

// done: Computation of epsilon relative (vacuum + medium)
double
epsilon(double r, 
    double ra,
	double rb, 
	double R, 
	double epsilon_a, 
	double epsilon_b, 
	double epsilon_R)
{
    if ((ra <= r) && (r <= rb)) {
        return 1.0;
    } else if ((rb <= r) && (r <= R)) {
        return epsilon_b + (epsilon_R - epsilon_b) * (r - rb) / (R - rb);
    } else {
        std::cout << "r = " << r << std::endl;
        throw std::runtime_error("epsilon: r is out of bounds");
    }
}

// done: Computation of rho_eps=rho_free(r)/eps0

double
rho_epsilon(double r,
	double ra,
	double rb,
    double R,
	double A)
{
    if ((ra <= r) && (r <= rb)) {
        return 4*A*(r - ra)*(rb - r) / ((rb - ra)*(rb - ra));
    } else if ((rb <= r) && (r <= R)) {
        return 0.0;
    } else {
        std::cout << "r = " << r << std::endl;
        throw std::runtime_error("rho_epsilon: r is out of bounds");
    }
}


double f_A(  double r, 
            double ra, 
            double rb, 
            double R, 
            double A, 
            double epsilon_a, 
            double epsilon_b, 
            double epsilon_R, 
            double h_i,
            double h_j
            ) {
    if (h_i == 0 || h_j == 0) {
        throw std::runtime_error("f_A: h_i or h_j is zero");
    }
    return epsilon(r, ra, rb, R, epsilon_a, epsilon_b, epsilon_R) * r / (h_i * h_j);
}

double f_b( double r, 
            double ra, 
            double rb, 
            double R, 
            double A, 
            double r_i,
            double h_i) {
    if (h_i == 0) {
        throw std::runtime_error("f_b: h_i is zero");
    }
    return rho_epsilon(r, ra, rb, R, A) * r*(r-r_i)/h_i;
}



int
main(int argc, char* argv[])
{
    // USAGE: Exercise6 [configuration-file] [<settings-to-overwrite> ...]

    // Read the default input
    string inputPath = "configuration.in.example";
    // Optionally override configuration file.
    if (argc > 1)
        inputPath = argv[1];

    ConfigFile configFile(inputPath);
    // Override settings
    for (int i = 2; i < argc; i++)
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Read geometrical inputs
    const double R  = configFile.get<double>("R");
    const double ra = configFile.get<double>("ra");
    const double rb = configFile.get<double>("rb");
    
    // Free charge source
    const double A = configFile.get<double>("A");
    
    // Dielectric relative permittivity
    const double epsilon_a = configFile.get<double>("epsilon_a");
    const double epsilon_b = configFile.get<double>("epsilon_b");
    const double epsilon_R = configFile.get<double>("epsilon_R");
    
    // Boundary conditions
    const double Va = configFile.get<double>("Va");
    const double VR = configFile.get<double>("VR");
    
    // Discretization
    const int N1 = configFile.get<int>("N1");
    const int N2 = configFile.get<int>("N2");
    
    // Fichiers de sortie:
    string fichier = configFile.get<string>("output");
    string fichier_phi = fichier+"_phi.out";
    string fichier_E   = fichier+"_E.out";
    string fichier_D   = fichier+"_D.out";

    // done: Create the finite elements
    const int pointCount = N1 + N2 + 1;

    // done: Initialize position of elements 
    vector<double> r(pointCount);
    for (int i = 0; i < N1; ++i){
        r[i] = ra + (rb - ra) * i / (N1);
    }
    for (size_t i = N1; i < pointCount; i++)
    {
        r[i] = rb + (R - rb) * (i - N1) / (N2);
    }

    std::cout << "r = ";
    for (int i = 0; i < pointCount; ++i)
        std::cout << r[i] << " ";

    std::cout << std::endl;
    
    
    // done: Calculate distance between elements
    vector<double> h(pointCount-1);
    vector<double> midPoint(pointCount-1);

    for (int i = 0; i < pointCount-1; ++i){
        h[i] = r[i+1] - r[i];
        midPoint[i] = (r[i] + h[i]/2.0);
    }
    std::cout << "h = ";
    for (int i = 0; i < h.size(); ++i)
        std::cout << h[i] << " ";
    
    std::cout << std::endl;

    // done: Construct the matrices
    vector<double> diagonal(pointCount, 1.0);  // Diagonal
    vector<double> lower(pointCount - 1, 0.0); // Lower diagonal
    vector<double> upper(pointCount - 1, 0.0); // Upper diagonal
    vector<double> rhs(pointCount, 0.0);       // Right-hand-side
    
    double p(1.0); // p = 1.0 for the trapezoidal rule

    for (int i = 0; i < pointCount-2; ++i) { 
        // use the trapezoidal rule to approximate the integral
        upper[i] = h[i] * (p * (f_A(r[i], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i+1]) + 
                                f_A(r[i+1], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i+1])) + 
                     (1.0-p) *  f_A(midPoint[i], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i+1]));

        lower[i] = upper[i];

        diagonal[i+1] = -( h[i+1] * (p * (f_A(r[i], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i]) + 
                                   f_A(r[i+1], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i])) + 
                        (1.0-p) *  f_A(midPoint[i], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i])) +
                    h[i+1] * (p * (f_A(r[i+1], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i]) + 
                                   f_A(r[i+2], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i])) + 
                        (1.0-p) *  f_A(midPoint[i], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i])));
  

        rhs[i] = h[i]*(p * f_b(r[i+1], ra, rb, R, A, r[i], h[i])) + 
                            (1.0-p)*f_b(midPoint[i], ra, rb, R, A, r[i], h[i]);
                            
    }

    // done: Set boundary conditions
    diagonal[0] = 1.0;
    upper[0] = 0.0;
    upper[pointCount-2] = 39.0;
    lower[pointCount-2] = 0.0;
    diagonal[pointCount-1] = 1.0;
    rhs[0] = Va;
    rhs[pointCount-1] = VR;

//    diagonal[1] = h[1] * (p * (f_A(r[i], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i+1]) + 
//                                    f_A(r[i+1], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i])) + 
//                         (1.0-p) *  f_A(midPoint[i], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i])) +
//                       h[i] * (p * (f_A(r[i], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i]) + 
//                                    f_A(r[i+1], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i])) + 
//                         (1.0-p) *  f_A(midPoint[i], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[i], h[i]));
    //upper[pointCount-2] = h[pointCount-2] * (p * (f_A(r[pointCount-2], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[pointCount-2], h[pointCount-1]) + 
                    //             f_A(r[pointCount-1], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[pointCount-2], h[pointCount]-1)) + 
                    //  (1.0-p) *  f_A(midPoint[pointCount-2], ra, rb, R, A, epsilon_a, epsilon_b, epsilon_R, h[pointCount-2], h[pointCount-1]));;

    // cout the elements of the matrix
    std::cout << "diagonal = ";
    for (int i = 0; i < pointCount; ++i) {
        std::cout << diagonal[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "lower = ";
    for (int i = 0; i < pointCount-1; ++i) {
        std::cout << lower[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "upper = ";
    for (int i = 0; i < pointCount-1; ++i) {
        std::cout << upper[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "rhs = ";
    for (int i = 0; i < pointCount; ++i) {
        std::cout << rhs[i] << " ";
    }
    std::cout << std::endl;

    
    // Solve the system of equations
    vector<double> phi = solve(diagonal, lower, upper, rhs);

    std::cout << "phi = ";
    for (int i = 0; i < pointCount; ++i) {
        std::cout << phi[i] << " ";
    }
    std::cout << std::endl;

    // done but unsure: Calculate electric field E and displacement vector D
    vector<double> E(pointCount - 1, 0);
    vector<double> D(pointCount - 1, 0);
    double dphidr=0.0;
    for (int i = 0; i < E.size(); ++i) {
        dphidr = (phi[i+1] - phi[i])/h[i];
        E[i] = -dphidr;
        D[i] = -dphidr*epsilon(r[i], ra, rb, R, epsilon_a, epsilon_b, epsilon_R);
    }

    // Export data
    {
        // Electric potential phi
        ofstream ofs(fichier_phi);
        ofs.precision(15);

        if (r.size() != phi.size())
            throw std::runtime_error("error when writing potential: r and "
                                     "phi does not have size");

        for (int i = 0; i < phi.size(); ++i) {
            ofs << r[i] << " " << phi[i] << endl;
        }
    }

    {
        // Electric field E
        ofstream ofs(fichier_E);
        ofs.precision(15);

        if (r.size() != (E.size() + 1))
            throw std::runtime_error("error when writing electric field: size of "
                                     "E should be 1 less than r");

        for (int i = 0; i < E.size(); ++i) {
            ofs << midPoint[i] << " " << E[i] << endl;
        }
    }
    {
        // Displacement field D
        ofstream ofs(fichier_D);
        ofs.precision(15);

        if (E.size() != D.size())
            throw std::runtime_error("error when writing displacement field: size of "
                                     "D should be equal than E");

        for (int i = 0; i < D.size(); ++i) {
            ofs << midPoint[i] << " " << D[i] << endl;
        }
    }

    return 0;
}

