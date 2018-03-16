#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include "jama_svd.h"
#include "tnt.h"
#include "half.hpp"

using namespace TNT;
using namespace JAMA;
using half_float::half;
using namespace half_float::detail;
using namespace std;

void helpInfo();
bool checkFilePath(string);
Array2D< double > fileToMatrix(string, int&, int&, int&); // Read the pgm file into a 2D Matrix
void pgma2b(const string&);
void pgmb2a(const string&);
void compressPgm(string, int);
void compressPgma(string, int);
void SVD2Pgm(string);
void SVDB2Pgm(string);
Array2D< double > product(Array2D< double >, Array2D< double >);
Array2D< double > transpose2D(Array2D<double>);
double norm2(Array2D< double >, Array2D< double >);
long fileSize(string);

int main(int argc, char** argv)
{
    // Check whether there are two arguments, if not exit and return to the command line
    if(argc != 3 && argc != 4){
        cout << "The number of inputs is wrong. Please correct your input." << endl;
        helpInfo();
        return -1;
    }
    // Check whether can successfully open the input file, if not, return to the command line
    string filePath = argv[2];
    if(!checkFilePath(filePath)){
        cout << "Cannot correctly open the file. Please check your file path." << endl;
        return -1;
    }

    char option = argv[1][0];
    switch(option){
        case '1': pgma2b(filePath);
                break;
        case '2': pgmb2a(filePath);
                break;
        case '3': {
                istringstream ss(argv[3]);
                int k;
                ss >> k;
                compressPgma(filePath, k);
                break;}
        case '4': SVD2Pgm(filePath);
                break;
        default: cout << "Something wrong!" << endl;
                break;
    }

//    compressPgma("./bin/debug/baboon.pgm", 32);
//    SVG2Pgm("./bin/debug/baboon.pgm.svd");
    return 0;
}

void helpInfo(){
    cout << "Usage 1: SVD.exe 1 image.pgm" << endl;
    cout << "\t Convert the ASCII pgm file to a binary file." << endl;
    cout << "Usage 2: SVD.exe 2 image_b.pgm" << endl;
    cout << "\t Convert the binary file generated in Usage 1 back to an ASCII file" << endl;
    cout << "Usage 3: SVD.exe 3 image.pgm" << endl;
    cout << "\t Compress the image.pgm via SVD and save it into a binary file called image_b.pgm.SVD" << endl;
    cout << "Usage 4: SVD.exe 4 image_b.pgm.SVD" << endl;
    cout << "\t Convert the binary file from Usage 3 back to an ASCII file that readable by the pgm viewer" << endl;
}

bool checkFilePath(string filePath){
    ifstream infile(filePath);
    if(!infile.good()){
        infile.close();
        return false;
    }
    infile.close();
    return true;
}

Array2D< double > fileToMatrix(string filePath, int &row, int &col, int &grayScale){
    ifstream infile(filePath);
    string line;
    getline(infile, line);
    getline(infile, line);
    infile >> row >> col >> grayScale;
    cout << row << " " << col << " " << grayScale << endl;
    Array2D< double > A(row, col, 0.0);
    for (int i=0; i < row; i++)
        for (int j=0; j < col; j++)
            infile >> A[i][j];
//    // Test the 2d Array
//    for (int i=0; i < row; i++){
//        for (int j=0; j < col; j++)
//            cout << A[i][j] << "\t";
//        cout << "\n";
//    }
    return A;
}

void pgma2b(const string& filename) {
      string new_name(filename.begin(), filename.end() - 4);
      new_name += "_b.pgm";
      ifstream pgm_ascii;
      ofstream pgm_bin;
      pgm_ascii.open(filename.c_str());
      pgm_bin.open(new_name.c_str(), ios_base::out | ios_base::binary);

      unsigned int height, width, buffer;
      unsigned char greyscale;
      string line;
      getline(pgm_ascii, line);
      getline(pgm_ascii, line);
      pgm_ascii >>height >> width >> buffer;
      greyscale = buffer;
      vector<unsigned char> image(height * width, 0);
      vector<unsigned char>::iterator pixel = image.begin();
      for (pixel = image.begin(); pixel != image.end(); ++pixel) {
            pgm_ascii >> buffer;
            *pixel = (unsigned char)buffer;
            cout << (int)(*pixel) << ' ';
      }
      pgm_bin.write((reinterpret_cast<char*>(&height)), 1).
              write((reinterpret_cast<char*>(&height))+ 1, 1).
              write((reinterpret_cast<char*>(&width)), 1).
              write((reinterpret_cast<char*>(&width)) + 1, 1).
              write(reinterpret_cast<char*>(&greyscale), 1);
      for (pixel = image.begin(); pixel != image.end(); ++pixel) {
            pgm_bin.write(reinterpret_cast<char*>(&(*pixel)), 1);
      }
      pgm_ascii.close();
      pgm_bin.close();
}

void pgmb2a(const string& filename) {
      string new_name(filename.begin(), filename.end() - 6);
      new_name += "2.pgm";
      ofstream pgm_ascii;
      ifstream pgm_bin;
      pgm_ascii.open(new_name.c_str());
      pgm_bin.open(filename.c_str(), ios_base::in | ios_base::binary);

      unsigned int height = 0, width = 0;
      unsigned char greyscale;
      pgm_bin.read(reinterpret_cast<char*>(&height), 1).
              read(reinterpret_cast<char*>(&height) + 1, 1).
              read(reinterpret_cast<char*>(&width), 1).
              read(reinterpret_cast<char*>(&width) + 1, 1).
              read(reinterpret_cast<char*>(&greyscale), 1);
      vector<unsigned char> image(height * width);
      vector<unsigned char>::iterator pixel;
      for (pixel = image.begin(); pixel != image.end(); ++pixel) {
            cout << *pixel << ' ';
            pgm_bin.read(reinterpret_cast<char*>(&(*pixel)), 1);
      }
      pgm_ascii << height << ' ' << width << '\n'
                << (unsigned int)greyscale << '\n';
      for (pixel = image.begin(); pixel != image.end(); ++pixel) {
            pgm_ascii << (unsigned int)(*pixel) << '\t';
      }
      pgm_ascii.close();
      pgm_bin.close();
}

void compressPgm(string filePath, int k){
    int row, col, grayScale;
    Array2D< double > pgmMatrix = fileToMatrix(filePath, row, col, grayScale);
    Array2D< double > U(pgmMatrix), S(pgmMatrix), V(pgmMatrix);
    SVD< double > SU(U), SS(S), SV(V);
    SU.getU(U);
    SS.getS(S);
    SV.getV(V);
    int rankS = SU.rank();
    if(rankS < k){
        cout << "k is greater than the rank, " << rankS << ". Please give a smaller k.";
        return;
    }
    if(k == 0)
        k = rankS;
//    for(int i=0; i < U.dim1(); i++){
//        for(int j=0; j < U.dim2(); j++)
//            cout << U[i][j] << '\t';
//        cout << '\n';
//    }
//    cout << '\n';
//    for(int i=0; i < S.dim1(); i++){
//        for(int j=0; j < S.dim2(); j++)
//            cout << S[i][j] << '\t';
//        cout << '\n';
//    }
//    cout << '\n';
//    for(int i=0; i < V.dim1(); i++){
//        for(int j=0; j < V.dim2(); j++)
//            cout << V[i][j] << '\t';
//        cout << '\n';
//    }
    // Write to the binary file
    string new_name(filePath.begin(), filePath.end() - 4);
    new_name += ".pgm.SVD";
    ofstream pgm_asc;
    pgm_asc.open(new_name.c_str(), ios_base::in);
    pgm_asc << row << col << grayScale;
    cout << row << " " << col << " " << grayScale << endl;
    pgm_asc.close();
//    ofstream pgm_bin;
//    pgm_bin.open(new_name.c_str(), ios_base::out | ios_base::binary);
//    unsigned int height, width, greyscale;
//    pgm_bin.write((reinterpret_cast<char*>(&height)), 1).
//            write((reinterpret_cast<char*>(&height)) + 1, 1).
//            write((reinterpret_cast<char*>(&width)), 1).
//            write((reinterpret_cast<char*>(&width)) + 1, 1).
//            write(reinterpret_cast<char*>(&greyscale), 1).
//            write(reinterpret_cast<char*>(&k), 1);
//    for(int i=0; i < U.dim1(); i++){
//        for(int j=0; j < k; j++)
//              cout << U[i][j];
//
//    }
    for(int i=0; i < k; i++){
//        half hNum = half_cast<half>(S[k][k]);
//        cout << hNum << '\t';
        cout << S[i][i] << '\t';
//        pgm_bin.write((reinterpret_cast<char*>(&hNum)), 1).
//                write((reinterpret_cast<char*>(&hNum)) + 1, 1);
    }
//    for(int i=0; i < k; i++){
//        for(int j=0; j < V.dim2(); j++)
//            half hNum = half_cast<half>(V[i][j]);
////            pgm_bin.write((reinterpret_cast<char*>(&hNum)), 1).
////                    write((reinterpret_cast<char*>(&hNum)) + 1, 1);
//    }
}

void compressPgma(string filePath, int k){
    int row, col, grayScale;
    Array2D< double > pgmMatrix = fileToMatrix(filePath, row, col, grayScale);
    Array2D< double > U(pgmMatrix), S(pgmMatrix), V(pgmMatrix);
    SVD< double > SU(U), SS(S), SV(V);
    SU.getU(U);
    SS.getS(S);
    SV.getV(V);
    Array2D< double > VT = transpose2D(V);
    int rankS = SU.rank();
    if(rankS < k){
        cout << "k is greater than the rank, " << rankS << ". Please give a smaller k.";
        return;
    }
    if(k == 0)
        k = rankS;
    string new_name(filePath.begin(), filePath.end() - 4);
    new_name += ".pgm.SVD";
    ofstream svd_asc;
    svd_asc.open(new_name.c_str());
    svd_asc << row << " " << col << " " << k << "\n" << grayScale << "\n";
    for(int i = 0; i < U.dim1(); i++){
        for(int j=0; j < k; j++)
            svd_asc << U[i][j] << " ";
        svd_asc << "\n";
    }
    for(int i = 0; i < k; i++){
        svd_asc << S[i][i] << " ";
    }
    svd_asc << "\n";
    for(int i = 0; i < k; i++){
        for(int j=0; j < VT.dim2(); j++)
            svd_asc << VT[i][j] << " ";
        svd_asc << "\n";
    }
    svd_asc.close();
}

void SVD2Pgm(string filePath){
    string new_name(filePath.begin(), filePath.end() - 8);
    new_name += "_asc.pgm";
    ifstream svd_asc;
    ofstream pgm_asc;
    svd_asc.open(filePath);
    pgm_asc.open(new_name);
    int row, col, k, grayScale;
    svd_asc >> row >> col >> k >> grayScale;
    cout << row << " " << col << " " << k << " " << grayScale << endl;
    pgm_asc << "P2\n";
    pgm_asc << row << " " << col << "\n";
    pgm_asc << grayScale << "\n";
    Array2D< double > U(row, k, 0.0);
    Array2D< double > S(k, k, 0.0);
    Array2D< double > VT(k, col, 0.0);
    for(int i = 0; i < row; i++)
        for(int j = 0; j < k; j++)
            svd_asc >> U[i][j];
    for(int i = 0; i < k; i++)
        svd_asc >> S[i][i];
    for(int i = 0; i < k; i++)
        for(int j = 0; j < col; j++)
            svd_asc >> VT[i][j];
    // Display U,S,V
//    for(int i = 0; i < row; i++){
//        for(int j = 0; j < k; j++)
//            cout << U[i][j] << "\t";
//        cout << "\n";
//    }
//    cout << endl;
//    for(int i = 0; i < k; i++){
//         for(int j = 0; j < k; j++)
//            cout << S[i][j] << "\t";
//         cout << "\n";
//    }
//    cout << endl;
//    for(int i = 0; i < k; i++){
//        for(int j = 0; j < col; j++)
//            cout << V[i][j] << "\t";
//        cout << "\n";
//    }
    // Display USV
    Array2D< double > US(product(U, S));
    Array2D< double > USV(row, col, 0.0);
    try{
        USV = product(US, VT);
    } catch (string e){
        cout << e << endl;
    }
    for(int i = 0; i < USV.dim1(); i++){
        for(int j = 0; j < USV.dim2(); j++)
            pgm_asc << USV[i][j] << " ";
        pgm_asc << "\n";
    }
//    for(int i = 0; i < USV.dim1(); i++){
//        for(int j = 0; j < USV.dim2(); j++)
//            cout << USV[i][j] << "\t";
//        cout << "\n";
//    }
    svd_asc.close();
    pgm_asc.close();
}

Array2D< double > product(Array2D< double > A, Array2D< double > B){
    int rowA = A.dim1();
    int colA = A.dim2();
    int rowB = B.dim1();
    int colB = B.dim2();
//    cout << "rowA" << rowA << endl;
//    cout << "colA" << colA << endl;
//    cout << "rowB" << rowB << endl;
//    cout << "colB" << colB << endl;
    if(colA != rowB){
        throw string("Cannot multiply these two matrices.");
    }
    Array2D< double > P(rowA, colB, 0.0);
    for(int i = 0; i < rowA; i++){
        for(int j = 0; j < colB; j++){
            for(int k = 0; k < colA; k++){
                P[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return P;
}

Array2D< double > transpose2D(Array2D< double > A){
    Array2D< double > T(A.dim2(), A.dim1(), 0.0);
    for(int i = 0; i < T.dim1(); i++)
        for(int j = 0; j < T.dim2(); j++)
            T[i][j] = A[j][i];
    return T;
}

double norm2(Array2D< double > A, Array2D< double > B){
    double n2 = 0.0;
    if (A.dim1() != B.dim1() || A.dim2() != B.dim2()){
        throw string("The two matrices do not have the same dimension.");
    }
    for(int i = 0; i < A.dim1(); i++)
        for(int j = 0; j < A.dim2(); j++)
            n2 += pow(abs(A[i][j] - B[i][j]), 2);
    return sqrt(n2);
}

long fileSize(string filePath){
    ifstream in(filePath, ios::binary|ios::ate);
    long fs = in.tellg();
    in.close();
    return fs;
}


