#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include "include.h"
#include "newmat.h"
#include "newmatio.h"
#include "newmatap.h"


using namespace std;

double Frobenious_Norm(Matrix Data)
{
    Data = Data*Data.t();
    double F = 0;
    for (int l=1; l<=Data.Nrows(); l++){
        for (int m=1; m<=Data.Ncols(); m++) {
            F += pow(Data(l,m),2);
        }
    }
    return pow(F,0.5);
}

Matrix Matrix_Sketch(Matrix Data, double epsilon)
{
    int cols_of_sketch = ceil(2.0/epsilon);
    cout << "l= " << cols_of_sketch << endl;
    IdentityMatrix I(cols_of_sketch);
    Matrix A, B;
    DiagonalMatrix C;
    double d;
    
    
    if (cols_of_sketch < Data.Nrows())
    {
        Matrix R(Data.Nrows(), cols_of_sketch), Dn(cols_of_sketch,cols_of_sketch);
        R = 0.0;
        Dn= 0.0;
        Matrix CH(Data.Nrows(),1);
        CH = 0.0;
        
        for (int l=1; l<=Data.Ncols(); l++){
            for (int m=1; m<=cols_of_sketch; m++ )
            {
                if(R.Column(m)==CH)
                {
                    R.Column(m) = Data.Column(l);
                    break;
                }
                else
                {
                    
                    
                    SVD(R, C, A, B);
                    d = pow(C(cols_of_sketch),2);
                    
                    for (int l=1; l<=cols_of_sketch; l++)
                    {
                        Dn(l,l) = pow(C(l)*C(l) - d,0.5);
                        
                    }
                    R = A*Dn;
                    
                }
                
            }
        }
        return R;
    }
    
    else
    {
        Matrix R(Data.Nrows(), Data.Nrows()), Dn(Data.Nrows(),Data.Nrows());
        R = 0.0;
        Dn= 0.0;
        Matrix CH(Data.Nrows(),1);
        CH = 0.0;
        
        SVD(Data*Data.t(), C, A, B);
        
        for (int m=1; m<=Data.Nrows(); m++)
        {
            Dn(m,m) = pow(C(m),0.5);
            
        }
        R = A*Dn*B.t();
        
        return R;
    }
    
    
}

int main (int argc, char* argv[])
{
    int dimension, no_of_data_points;
    double epsilon;
    
    sscanf (argv[1], "%d", &dimension);
    sscanf (argv[2], "%d", &no_of_data_points);
    sscanf (argv[3], "%lf", &epsilon);
    ifstream input_file(argv[4]);
    ofstream output_file(argv[5]);
    
    Matrix Data(dimension, no_of_data_points);
    
    cout << "Edo Liberty's Matrix Sketching Algorithm" << endl;
    cout << "----------------------------------------" << endl;
    cout << "Original Data-Matrix has " << dimension << "-rows & " << no_of_data_points << "-cols" << endl;
    cout << "Epsilon = " << epsilon << " (i.e. max. of " << 100*epsilon << "% reduction of  Frobenius-Norm of the Sketch Matrix)"<< endl;
    cout << "Input File = " << argv[4] << endl;
    
    // Read the Data
    for (int i = 1; i <= dimension; i++)
        for (int j = 1; j <= no_of_data_points; j++)
        {
            double x;
            input_file >> x;
            Data(i,j) = x;
        }
    
    // Compute the Frobenius-Norm of the original Data-Matrix
    double Data_Forbenius_Norm = Frobenious_Norm(Data);
    cout << "Frobenius Norm of the (" << Data.Nrows() << " x " << Data.Ncols() << ") Data Matrix = ";
    cout << Data_Forbenius_Norm << endl;
    
    Matrix Sketch(dimension, min(dimension, (int) ceil(2/epsilon)));
    Sketch = Matrix_Sketch(Data, epsilon);
    double Sketch_Forbenius_Norm = Frobenious_Norm(Sketch);
    cout << "Frobenius Norm of the (" << Sketch.Nrows() << " x " << Sketch.Ncols() << ") Sketch Matrix = ";
    cout << Sketch_Forbenius_Norm << endl;
    cout << "Change in Frobenius-Norm between Sketch & Original  = ";
    cout << setprecision(3) << 100*(Sketch_Forbenius_Norm - Data_Forbenius_Norm)/Data_Forbenius_Norm << "%" << endl;
    
    output_file << Sketch;
    cout << "File `" << argv[5] << "' contains a (" << Sketch.Nrows() << " x " << Sketch.Ncols();
    cout << ") Matrix-Sketch" << endl;
    
    
}
