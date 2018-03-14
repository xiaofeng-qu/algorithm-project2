#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include "jama_svd.h"
#include "tnt.h"

using namespace TNT;
using namespace JAMA;
using namespace std;

void helpInfo();
bool checkFilePath(string);
Array2D<double> fileToMatrix(string); // Read the pgm file into a 2D Matrix
void pgma2b(const string&);
void pgmb2a(const string&);
void convertInitialPgmToASCII(string);
void compressInitialPgm(string);
void convertCompressedPgmToASCII(string);

int main(int argc, char** argv)
{
    // Check whether there are two arguments, if not exit and return to the command line
    if(argc != 3){
        cout << "There should be two arguments. Please correct your input." << endl;
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
 //   cout << option << endl;
    switch(option){
        case '1': pgma2b(filePath);
                break;
        case '2': pgmb2a(filePath);
                break;
        case '3': break;
        case '4': break;
        default: cout << "Something wrong!" << endl;
                break;
    }


//    ifstream infile(filePath);
//    string line;
//    char firstChar;
//    getline(infile, line);
//    cout << line[0] << endl;
    // Read the initial pgm file into a 2D matrix
//    ifstream infile(filePath);
//    string line;
//    while(true){
//        getline(inputFile, line);
//        istringstream ss(line);
//        char firstChar
//        ss >> firstChar;
//
//    }
//
//    int M = 5, N = 5;
//    Array2D< double > A(M,N, 0.0);    /* create MxN array; all zeros */
//    for (int i=0; i < M; i++)
//        for (int j=0; j < N; j++)
//            A[i][j] = rand() % 100;              /* initalize array values */
//    Array2D< double > B = A.copy();       /* create a new copy */
//    Array2D< double > C(B);               /* create a new view of B */
//                                         /* Both arrays (B & C) share data */
//    for (int i=0; i < M; i++){
//        for (int j=0; j < N; j++){
//            cout << A[i][j] << "\t";
//        }
//        cout << endl;
//    }
//    cout<<endl;
    Array2D< double > A = fileToMatrix(filePath);
    cout << endl;
    SVD< double > D(A);
    D.getV(A);
    for (int i=0; i < 4; i++){
        for (int j=0; j < 5; j++){
            cout << A[i][j] << "\t";
        }
        cout << endl;
    }
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

Array2D<double> fileToMatrix(string filePath){
    ifstream infile(filePath);
    string line;
    int row, col, grayScale;
    getline(infile, line);
    getline(infile, line);
    infile >> row >> col >> grayScale;
    cout << row << " " << col << " " << grayScale << endl;
    Array2D< double > A(row, col, 0.0);
    for (int i=0; i < row; i++)
        for (int j=0; j < col; j++)
            infile >> A[i][j];
    // Test the 2d Array
    for (int i=0; i < row; i++){
        for (int j=0; j < col; j++)
            cout << A[i][j] << "\t";
        cout << "\n";
    }
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
      char magic_number[2];
      pgm_ascii >> magic_number[0] >> magic_number[1] >> height >> width >> buffer;
      greyscale = buffer;
      vector<unsigned char> image(height * width, 0);
      vector<unsigned char>::iterator pixel = image.begin();
      for (pixel = image.begin(); pixel != image.end(); ++pixel) {
            pgm_ascii >> buffer;
            *pixel = (unsigned char)buffer;
            cout << (int)(*pixel) << ' ';
      }
      magic_number[1] = '5';
      pgm_bin.write(magic_number, 2);
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
      char magic_number[2];
      pgm_bin.read(magic_number, 2).
              read(reinterpret_cast<char*>(&height), 1).
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
      magic_number[1] = '2';
      pgm_ascii << magic_number[0] << magic_number[1] << '\n' << height << ' ' << width << '\n'
                << (unsigned int)greyscale << '\n';
      for (pixel = image.begin(); pixel != image.end(); ++pixel) {
            pgm_ascii << ' ' << (unsigned int)(*pixel);
      }
      pgm_ascii.close();
      pgm_bin.close();
}
