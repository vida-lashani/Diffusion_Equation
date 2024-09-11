#include <fstream> 
#include <iostream> 
#include <cmath>
#include <iomanip>
using namespace std;

int main() 
{ 
    double r_values[] = {0.25, 0.5, 1, 2.5, 10, 100}; // r values
    double theta_values[] = {0, 0.5, 1}; // theta values
    int num_r = sizeof(r_values) / sizeof(r_values[0]);
    int num_theta = sizeof(theta_values) / sizeof(theta_values[0]);

    // Arrays for storing the solution and errors
    double y[32]; // location vector 
    double x[1000][32]; // u values vector
    double a[32][32]; // coefficient matrix
    double b[1000][32]; // known vector
    double Linf[10000]; // error history

    // Iterate through all combinations of r and theta
    for (int r_index = 0; r_index < num_r; r_index++) {
        for (int theta_index = 0; theta_index < num_theta; theta_index++) {
            
            double r = r_values[r_index];
            double theta = theta_values[theta_index];
            int n = 0; // iteration number
            int k = 0; // auxiliary iteration number

            // Initial value setup
            for (int i = 0; i <= 31; i++) {
                if (i <= 15) 
                    x[0][i] = 0; 
                else 
                    x[0][i] = 1;
            }

            // Create file names with r and theta values
            string base_filename = "solution_r_" + to_string(r) + "_theta_" + to_string(theta) + ".txt";
            ofstream solution_file(base_filename);
            solution_file << "x t u" << endl;

            // Saving initial condition
            for (int i = 0; i <= 31; i++) {
                y[i] = i / 31.0;
                solution_file << y[i] << " " << r * n / 32.0 << " " << x[k][i] << endl;
            }

            // Initialize coefficient matrix a
            for (int i = 0; i <= 31; i++) 
                for (int j = 0; j <= 31; j++) 
                    a[i][j] = 0;
            
            a[0][0] = 1;
            a[31][31] = 1;
            for (int i = 1; i <= 30; i++) {
                a[i-1][i] = -1 * r * theta;
                a[i][i] = 1 + 2 * r * theta;
                a[i][i+1] = -1 * r * theta;
            }

            // Modify matrix to account for Thomas method
            for (int i = 1; i <= 30; i++) {
                double m = a[i][i-1] / a[i-1][i-1];
                a[i][i] = a[i][i] - m * a[i][i+1];
            }

            // Set initial known vector b
            b[0][0] = 0;
            b[0][31] = 1;
            for (int i = 1; i <= 30; i++) {
                b[0][i] = x[0][i] + r * (1 - theta) * (x[0][i-1] - 2 * x[0][i] + x[0][i+1]);
            }

            double max_error = 1;
            ofstream error_file("error_r_" + to_string(r) + "_theta_" + to_string(theta) + ".txt");
            error_file << "time step  Linf" << endl;

            // Main loop for solving the system
            while (max_error > 0.00001) {
                max_error = 0;
                // Update known vector
                for (int i = 1; i <= 30; i++) {
                    b[k][0] = 0;
                    b[k][31] = 1;
                    double m = -1.0 * r * theta / (1 + 2 * r * theta);
                    b[k][i] = b[k][i] - m * b[k][i-1];
                }

                // Thomas method to solve system
                x[k+1][0] = 0;
                x[k+1][31] = 1;
                for (int i = 30; i >= 1; i--) {
                    x[k+1][i] = (b[k][i] - a[i][i+1] * x[k+1][i+1]) / a[i][i];
                }

                // Calculate error and update known vector
                for (int i = 1; i <= 30; i++) {
                    b[k+1][i] = x[k+1][i] + r * (1 - theta) * (x[k+1][i-1] - 2 * x[k+1][i] + x[k+1][i+1]);
                    double current_error = fabs(x[k+1][i] - x[k][i]);
                    if (current_error > max_error) {
                        max_error = current_error;
                    }
                }
                Linf[n] = max_error;
                error_file << n << " " << Linf[n] << endl;

                // Save solution at each time step
                for (int i = 0; i <= 31; i++) {
                    y[i] = i / 31.0;
                    solution_file << y[i] << " " << r * n / 32.0 << " " << x[k+1][i] << endl;
                }

                n++;
                k++;
                if (k == 999) {
                    for (int i = 0; i <= 31; i++) {
                        x[1][i] = x[1000][i];
                        b[0][i] = b[999][i];
                    }
                    k = 0;
                }
            }

            // Close files after finishing for this combination of r and theta
            solution_file.close();
            error_file.close();
            cout << "Finished for r = " << r << " and theta = " << theta << endl;
        }
    }

    return 0;
}
