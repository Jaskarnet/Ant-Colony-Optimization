#include "adjacencyMatrix.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>


using namespace std;

void adjacencyMatrix::showMatrix() {
	cout << endl << "    ";
	for (int i = 0; i < v; i++) cout << setw(4) << i;
	cout << endl << endl;
	for (int i = 0; i < v; i++)
	{
		cout << setw(4) << left << i;
		for (int j = 0; j < v; j++) {
			if (A[i][j] == -1) cout << setw(4) << right << "-";
			else cout << setw(4) << right << A[i][j];
		}
		cout << endl;
	}
	cout << endl;
}

void adjacencyMatrix::deleteMatrix() {
	for (int i = 0; i < v; i++) delete[] A[i];
	delete[] A;
}

void adjacencyMatrix::insert(int row, int column, T value) {
	A[row][column] = value;
}

void adjacencyMatrix::insertBidirected(int row, int column, T value) {
	A[row][column] = value;
	A[column][row] = value;
}

int adjacencyMatrix::getVertexNumber() {
	return v;
}

int adjacencyMatrix::getEdgeNumber() {
	return e;
}

int** adjacencyMatrix::getA() {
	return A;
}

int adjacencyMatrix::getWeight(int i) {
	return w[i];
}

int adjacencyMatrix::getRoadLength(int i) {
	return routeLength[i];
}

int adjacencyMatrix::getStartingVertex() {
	return startingVertex;
}

int adjacencyMatrix::getLowestWeight() {
	return lowestWeight;
}

bool adjacencyMatrix::getTest() {
	return test;
}

void adjacencyMatrix::setTestTrue() {
	test = true;
}


bool adjacencyMatrix::file_read_line(ifstream& file, int tab[], int size)
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


adjacencyMatrix* adjacencyMatrix::file_read_graphDirected(std::string file_name) {
	ifstream file;
	int tab[500];
	file.open(file_name.c_str());

	if (file.is_open())
	{
		if (adjacencyMatrix::file_read_line(file, tab, 1))
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