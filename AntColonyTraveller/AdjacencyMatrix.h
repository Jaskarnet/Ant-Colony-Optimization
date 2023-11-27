#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "AdjacencyMatrix.h"

using namespace std;

template <typename T>
class adjacencyMatrix {
public:
    adjacencyMatrix() {
        lowestWeight = 9999999;
    }

    adjacencyMatrix(int vertexNumber) {
        startingVertex = 0;
        lowestWeight = 9999999;
        routeLength = new int[vertexNumber];
        w = new int[vertexNumber];
        v = vertexNumber;
        e = (v * (v - 1)) / 2;
        A = new T* [vertexNumber];

        for (int i = 0; i < vertexNumber; i++) {
            A[i] = new T[vertexNumber];
        }

        for (int i = 0; i < vertexNumber; i++) {
            for (int j = 0; j < vertexNumber; j++) {
                A[i][j] = -1;
            }
        }
    }

    ~adjacencyMatrix() {

    }

	void showMatrix() {
		cout << endl << "    ";
		for (int i = 0; i < v; i++) cout << setw(4) << i;
		cout << endl << endl;
		for (int i = 0; i < v; i++)
		{
			cout << setw(4) << left << i;
			for (int j = 0; j < v; j++) {
				if (A[i][j] == -1) cout << setw(4) << right << "-";
				else cout << setw(4) << right << fixed << setprecision(0) << A[i][j];
			}
			cout << endl;
		}
		cout << endl;
	}

	void deleteMatrix() {
		for (int i = 0; i < v; i++) delete[] A[i];
		delete[] A;
	}

	void insert(int row, int column, int value) {
		A[row][column] = value;
	}

	void insertBidirected(int row, int column, int value) {
		A[row][column] = value;
		A[column][row] = value;
	}

	int getVertexNumber() {
		return v;
	}

	int getEdgeNumber() {
		return e;
	}

	T** getA() {
		return A;
	}

	int getWeight(int i) {
		return w[i];
	}

	int getRoadLength(int i) {
		return routeLength[i];
	}

	int getStartingVertex() {
		return startingVertex;
	}

	int getLowestWeight() {
		return lowestWeight;
	}

	bool getTest() {
		return test;
	}

	void setTestTrue() {
		test = true;
	}


	bool file_read_line(ifstream& file, int tab[], int size)
	{
		string s;
		getline(file, s);

		if (file.fail() || s.empty())
			return(false);

		istringstream in_ss(s);

		for (int i = 0; i < size; i++)
		{
			in_ss >> tab[i];
			if (in_ss.fail())
				return(false);
		}
		return(true);
	}


	adjacencyMatrix* file_read_graphDirected(std::string file_name) {
		ifstream file;
		int tab[500];
		file.open(file_name.c_str());

		if (file.is_open())
		{
			if (file_read_line(file, tab, 1))
			{
				int vertexNumber = tab[0];
				adjacencyMatrix* adjMatrix = new adjacencyMatrix(vertexNumber);
				adjMatrix->e = (vertexNumber * (vertexNumber - 1)) / 2;



				for (int i = 0; i < vertexNumber; i++)
					if (file_read_line(file, tab, vertexNumber)) {
						for (int j = 0; j < vertexNumber; j++) {
							// wczytywanie opisu kraw�dzi: tab[0] - wierzcho�ek pocz�tkowy, tab[1] - wierzcho�ek ko�cowy, tab[2] - waga
							adjMatrix->insert(i, j, tab[j]);
							if (i != j && adjMatrix->A[i][j] < adjMatrix->lowestWeight) adjMatrix->lowestWeight = adjMatrix->A[i][j];
						}
					}
					else
					{
						cout << "File error - READ DATA" << endl;
						break;
					}
				return adjMatrix;
			}
			else
				cout << "File error - READ INFO" << endl;
			file.close();
		}
		else
			cout << "File error - OPEN" << endl;
	}

private:
    T** A;
    int v, e;
    int* w, * routeLength;
    int startingVertex;
    int lowestWeight;
    bool test;
};