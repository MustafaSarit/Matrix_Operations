///\author Mustafa SARITEMUR
///\date December 2017

/*! \mainpage Introduction
 *
 * \section intro_sec Overview
 *
 * The program is intended for Matrix Operations.
 * The input to the program have to be .txt file format that stores an augmented matrix.
 * Then, the matrix is processed by the program, which gives determinant, solution with gauss elimination or prints the matrix.
 */

/**
 *  @file main.cpp
 *  @brief The main file. Contains the whole program.
 */

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <cmath>
#include <vector>


using namespace std;                        ///\namespace std

int g_row;                                  /**< global storage of the given matrix's length. */
vector<double> g_line;                      /**< global storage of the given matrix's width. */
vector< vector<double> > g_input;           /**< global storage of the given matrix itself. */

/** \brief The function controls if the given input files name format is correct.
 *
 * \param txts The local storage of the given input file name.
 * \return True if input file name format correct, false if input file name format wrong.
 *
 */
bool controler(string txts){
    string tx = ".txt";
    int a = 3;
    for(unsigned int i = txts.size() - 1; a >= 0; i--, a--){
        if(txts[i] == tx[a])
            continue;
        else
            return false;
    }
    return true;
}


/** \brief The function calculate the determinant of the given matrix with recursion.
 *
 * \param row The local storage of the given matrix's length.
 * \param input The local storage of the given matrix itself.
 * \param line The local storage of the given matrix's width.
 * \return Double data type of Determinant.
 *
 */
double determinant(int row, vector< vector<double> > input, vector<double> line){
    int c, subi, i, j, subj;
    vector< vector<double> > submat(row, line);
    double d = 0;           /**< Local storage of the determinant. */
    if (row == 2)
    {
        return( (input[0][0] * input[1][1]) - (input[1][0] * input[0][1]));
    }
    else
    {
        for(c = 0; c < row; c++)
        {
            subi = 0;
            for(i = 1; i < row; i++)
            {
                subj = 0;
                for(j = 0; j < row; j++)
                {
                    if (j == c)
                    {
                        continue;
                    }
                    submat[subi][subj] = input[i][j];
                    subj++;
                }
                subi++;
            }
        d = d + (pow(-1 ,c) * input[0][c] * determinant(row - 1 ,submat, line));
        }
    }
    return d;
}


/** \brief The function do the Gaussian elimination to find the solution.
 *
 * \param input The local storage of the given matrix itself.
 * \return The solution with vector<double> type.
 *
 */
vector<double> gauss(vector< vector<double> > input) {
    int n = input.size();

    for (int i = 0; i < n; i++) {
        // Search for maximum in this column
        double maxEl = abs(input[i][i]);
        int maxRow = i;
        for (int k = i+1; k < n; k++) {
            if (abs(input[k][i]) > maxEl) {
                maxEl = abs(input[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k = i; k < n+1; k++) {
            double tmp = input[maxRow][k];
            input[maxRow][k] = input[i][k];
            input[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k = i+1; k < n; k++) {
            double c = -input[k][i] / input[i][i];
            for (int j = i; j < n+1; j++) {
                if (i == j) {
                    input[k][j] = 0;
                } else {
                    input[k][j] += c * input[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);            /**< Local storage of the solution */
    for (int i = n-1; i >= 0; i--) {
        x[i] = input[i][n] / input[i][i];
        for (int k = i-1; k >= 0; k--) {
            input[k][n] -= input[k][i] * x[i];
        }
    }
    return x;
}


/** \brief The function Hails the user and get the input file's name.
 *
 * \return input file's name with string data type.
 *
 */
string welcome(){
    string txts;            /**< Local storage of the input file's name */
    cout << "Welcome to MATRIX OPERATIONS \t\t\t\t | Mustafa SARITEMUR | \n\n"
    << "This program can do the Gaussian elimination to find solution of given square matrix,\n"
    << "find out the value of determinant or  simply print the given matrix\n\n" << endl;
    cout << "Enter input file's name: \t!!(format of input have to be '.txt')" << endl;
    cin  >> txts;
    cout << "\n";
    if(controler(txts) == false){
        cout << "Format of given input file is invalid!!" << endl;
        exit(0);
    }
    return txts;
}


/** \brief The function print the matrix, determinant and (if rank allows) solution to console also create an output file with .html format and copy all to it. If encounter any error print the unique error message and exit program.
 *
 * \param row The local storage of the given matrix's length.
 * \param input The local storage of the given matrix itself.
 * \param line The local storage of the given matrix's width.
 * \param determ The local storage of the determinant.
 *
 */
void print(vector< vector<double> > input, int row, vector<double> line, double determ) {
    ofstream output;
    output.open("output.html");
    if(output.bad() || output.fail()){
        cout << "The output file is not working properly!!" << endl;
        output.close();
        exit(0);
    }

    output << "<!DOCTYPE html><html lang=\"en\"><head> <meta charset=\"utf-8\"> <title>Semester Project</title> <style>body{background: #fff; padding: 0px; margin: 0px;}table.given td:hover{background-color: aqua; -moz-transition: all 0.6s ease-in; -webkit-transition: all 0.6s ease-in; -o-transition: all 0.6s ease-in; transition: all 0.6s ease-in;}table.given{width: 300px; font-family: arial, sans-serif; border-collapse: collapse; margin: auto; margin-top: 25px;}table.given td, th{text-align: center; padding: 8px;}table.given td:last-child{border-left: 1px solid #000; border-right: 1px solid #000;}table.given td:first-child{border-left: 1px solid #000}table.result{width: 160px; font-family: arial, sans-serif; border-collapse: collapse; margin: auto;}table.result td, th{text-align: center; padding: 8px; border: 1px solid #d1d1d1;}table.result tr:nth-child(even){background-color: #e3e3e3;}div#text{margin: 25px auto; width: 300px; text-align: center;}</style></head><body>";
    output << "<table class=\"given\">";
    int n = input.size();
    for (int i = 0; i < n; i++) {
        output << "<tr>";
        for (int j = 0; j < n+1; j++) {
            output << "<td>" << input[i][j] <<"</td>";
            cout << input[i][j] << "\t";
            if (j == n-1) {
                cout << "| ";
            }
        }
        output << "</tr>";
        cout << "\n";
    }
    output << "</table>";

    cout << endl;

    // Determinant of the matrix
    output << "<div id=\"text\"> Determinant of the given matrix: " << determ << "</div>";
    cout << "Determinant of the given matrix: " << determ << endl;

        if(determ == 0){
            output << " <div id=\"text\"> Rank of this matrix does not allow for Gaussian elimination!! </div>";
            cout << "Rank of this matrix does not allow for Gaussian elimination!!" << endl;
        } else {
            // Doing the Gaussian elimination
            vector<double> resultG(row);
            resultG = gauss(input);

            // Printing the Result of Gaussian elimination
            cout << "Results:\n";

            output << " <table class=\"result\"> <tr> <th>Results</th> </tr></table> <table class=\"result\">";
            string parameter = "Variable";
            for (int i = 0; i < row; i++) {
                output << " <tr> <td>" << parameter << i+1 << "</td><td>" << resultG[i] << "</td></tr>";
                cout << parameter << i+1 << " = " <<resultG[i] << endl;
            }
        }
    output << "</table></body></html>";
    output.close();
}


/** \brief The function reads the input file and store matrix to g_input. If encounter any error print the unique error message and exit program.
 *
 * \param txts Local storage of the input file's name.
 *
 */
void readInput(string txts){

    ifstream txt;

    txt.open(txts);
    if(txt.bad() || txt.fail()){
        cout << "The file you have entered either is not working or does not exist!" << endl;
        txt.close();
        exit(0);
    }
    int row;
    txt >> row;
    if(txt.fail()){
        cout << "The contents of the given file are not correct!!" << endl;
        txt.close();
        exit(0);
    }

    vector<double> line(row + 1, 0);
    vector< vector<double> > input(row, line);

    // Read input data
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < row; j++) {
            txt >> input[i][j];
            if(txt.fail()){
                cout << "The contents of the given file are not correct!!" << endl;
                txt.close();
                exit(0);
            }
        }
    }

    for (int i = 0; i < row; i++) {
        txt >> input[i][row];
        if(txt.fail()){
            cout << "The contents of the given file are not correct!!" << endl;
            txt.close();
            exit(0);
        }
    }

    g_row = row;
    g_line = line;
    g_input = input;

    txt.close();

}


/** \brief The main function. Calls other functions.
 *
 * \return 1 if program worked correctly, 0 if program encounter an error.
 *
 */
int main() {
    string txts = welcome();
    readInput(txts);
    double determ = determinant(g_row, g_input, g_line);
    print(g_input, g_row, g_line, determ);
    return 1;
}
