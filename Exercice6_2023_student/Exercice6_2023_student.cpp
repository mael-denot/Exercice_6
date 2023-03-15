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
        double pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution[diag.size() - 1] =
      new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    for (int i(diag.size() - 2); i >= 0; --i)
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];

    return solution;
}

// TODO: Computation of epsilon relative (vacuum + medium)
double
epsilon(double r, 
	double R, 
	double rb, 
	double epsilon_a, 
	double epsilon_b, 
	double epsilon_R)
{
     return 1.0;
}

// TODO: Computation of rho_eps=rho_free(r)/eps0

double
rho_epsilon(double r,
	double ra,
	double rb,
	double A)
{
    return 0.0;
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

    // TODO: Create the finite elements
    const int pointCount = 10;

    // TODO: Initialize position of elements
    vector<double> r(pointCount);
    for (int i = 0; i < r.size(); ++i)
        r[i] = i;
    
    // TODO: Calculate distance between elements
    vector<double> h(pointCount-1);
     vector<double> midPoint(pointCount-1);
    for (int i = 0; i < h.size(); ++i){
        h[i] = 1;
        midPoint[i] = i;
    }
    // TODO: Construct the matrices
    vector<double> diagonal(pointCount, 1.0);  // Diagonal
    vector<double> lower(pointCount - 1, 0.0); // Lower diagonal
    vector<double> upper(pointCount - 1, 0.0); // Upper diagonal
    vector<double> rhs(pointCount, 0.0);       // Right-hand-side
    
    for (int i = 0; i < pointCount; ++i) { 
        // TODO: Matrix entries
        // upper[i] = ...
        // lower[i] = ...
        // diagonal[i] = ...
        // rhs[i] = ...
    }
    
    // TODO: Set boundary conditions
    // ...
    
    // Solve the system of equations
    vector<double> phi = solve(diagonal, lower, upper, rhs);

    // TODO: Calculate electric field E and displacement vector D
    vector<double> E(pointCount - 1, 0);
    vector<double> D(pointCount - 1, 0);
    double dphidr=0.0;
    for (int i = 0; i < E.size(); ++i) {
        //E[i] = ...;
        //D[i] = ...; //NB: Normalized D (factor of eps_0)
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

