#include <iostream>
#include <fstream>
#include <list>
#include <mpi.h>
#include <math.h>

struct imageHeader
{
	std::string p; //tipul imaginii: P5 sau P6
	std::string comment; //comentariu
	int width; //numarul de coloane din matricea de pixeli
	int height; //numarul de linii din matricea de pixeli
	int maxVal; //valoarea maxima a unui pixel
};

unsigned char** allocMatrix(int rows, int cols);
void readImageHeader(std::ifstream& inputStream, imageHeader& header);
unsigned char** readMatrixFromFile(std::ifstream& inputStream, imageHeader& header, int multiplier);
void writeImageFile(std::string outputFileName, imageHeader header, unsigned char** imageMatrix, int multiplier);
void freeImageMatrix(unsigned char** imageMatrix, imageHeader header);
void createFilter(float matrix[3][3], std::initializer_list<int> values, int factor);
void rotateMatrix(float matrix[3][3], int rows, int columns);
void setValueIfValid(int x, int y, int n, int m, imageHeader header, int multiplier, float elemMatrix[3][3], unsigned char** imageMatrix);
unsigned char** applyFilter(unsigned char** imageMatrix, imageHeader header, float filterMatrix[3][3], int multiplier, int startPoint, int endPoint);
unsigned char** readImageFromFile(char* fileName, imageHeader& header, int& multiplier);
void initializeAllFilters(float smoothingFilter[3][3], float gaussianFilter[3][3], float sharpenFilter[3][3], float meanFilter[3][3], float embossFilter[3][3]);
void movePointsToFinalMatrix(unsigned char** tempMatrix, unsigned char** finalMatrix, int start, int stop, int width);
unsigned char** applyFiltersToImage(int rank, int nProcesses, std::list<std::string>& filters, imageHeader header, int multiplier, unsigned char** imageMatrix);

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		std::cout << "Este necesar cel putin un filtru" << std::endl;
		return 0;
	}

	int rank; // indexul procesului curent
	int nProcesses; //numarul total de procese
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcesses);
	
	int multiplier; //1 pt pozele alb-negru si 3 pt color(1 pixel are 3 canale de culoare)
	imageHeader header;
	unsigned char** imageMatrix = nullptr;

	if (rank == 0)
	{
		//procesul 0 face citirea din fisier
		imageMatrix = readImageFromFile(argv[1], header, multiplier);
	}
	
	std::list<std::string> filters; //lista filtrelor
	for (int i = 3; i < argc; i++)
	{
		//pune fiecare filtru in lista
		filters.push_back(std::string{ argv[i] });
	}

	//procesul 0 trimite multiplicatorul, numarul de linii, de coloane si dimensiunea maxima a unei pixel celorlalte procese
	MPI_Bcast(&multiplier, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&header.width, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&header.height, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&header.maxVal, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//aplica filtrele pe imagine
	imageMatrix = applyFiltersToImage(rank, nProcesses, filters, header, multiplier, imageMatrix);

	if (rank == 0)
	{
		//procesul 0 o sa scrie rezultatul in fisier
		std::string outputFileName{ argv[2] };
		writeImageFile(outputFileName, header, imageMatrix, multiplier);
	}

	//elibereaza memoria
	freeImageMatrix(imageMatrix, header);

	MPI_Finalize();
}

unsigned char** allocMatrix(int rows, int cols)
{
	try
	{
		//aloca matricea ca un vector(matrice liniarizata) pentru
		//a avea o zona continua in memorie - se poate trimite cu un singur send
		unsigned char* liniarMatrix = new unsigned char[rows * cols];

		//vrem sa putem indexa totusi matricea ca m[i][j], deci alocam un vector de pointeri
		//si fiecare m[i] puncteaza unde ar incepe coloana lui in matricea liniara
		unsigned char** matrix2d = new unsigned char* [rows];
		for (int r = 0; r < rows; r++)
		{
			matrix2d[r] = liniarMatrix + r * cols;
		}

		return matrix2d;
	}
	catch (std::bad_alloc&)
	{
		std::cout << "Failed to allocate memory" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 0);
		MPI_Finalize();
	}
}

void readImageHeader(std::ifstream& inputStream, imageHeader& header)
{
	std::string newLine;
	std::getline(inputStream, header.p);
	std::getline(inputStream, header.comment);
	inputStream >> header.width;
	inputStream >> header.height;
	inputStream >> header.maxVal;
	std::getline(inputStream, newLine);
}

unsigned char** readImageFromFile(char* fileName, imageHeader& header, int& multiplier)
{
	std::string inputFileName{ fileName };
	std::ifstream inputStream{ inputFileName, std::ios::in | std::ios::binary };

	if (inputStream.fail())
	{
		//fisierul de intrare nu exista, inchide toate procesele
		std::cout << "Imaginea de intare nu exista" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 0);
		MPI_Finalize();
		return 0;
	}

	readImageHeader(inputStream, header);

	if (header.p.compare("P5") == 0)
	{
		multiplier = 1;
	}
	else if (header.p.compare("P6") == 0)
	{
		multiplier = 3;
	}

	unsigned char** imageMatrix = readMatrixFromFile(inputStream, header, multiplier);
	inputStream.close();

	return imageMatrix;
}

unsigned char** readMatrixFromFile(std::ifstream& inputStream, imageHeader& header, int multiplier)
{
	//citeste matricea de pixeli din fisier
	int width = header.width * multiplier;
	unsigned char** imageMatrix = allocMatrix(header.height, width);

	for (int i = 0; i < header.height; i++)
	{
		inputStream.read(reinterpret_cast<char*>(imageMatrix[i]), sizeof(unsigned char) * width);
	}

	return imageMatrix;
}

void initializeAllFilters(float smoothingFilter[3][3], float gaussianFilter[3][3], float sharpenFilter[3][3], float meanFilter[3][3], float embossFilter[3][3])
{
	//scrie valorile pentru fiecare filtru in matricele corespunzatoare
	createFilter(smoothingFilter, { 1, 1, 1, 1, 1, 1, 1, 1, 1 }, 9);
	rotateMatrix(smoothingFilter, 3, 3);

	createFilter(gaussianFilter, { 1, 2, 1, 2, 4, 2, 1, 2, 1 }, 16);
	rotateMatrix(gaussianFilter, 3, 3);

	createFilter(sharpenFilter, { 0, -2, 0, -2, 11, -2, 0, -2, 0 }, 3);
	rotateMatrix(sharpenFilter, 3, 3);

	createFilter(meanFilter, { -1, -1, -1, -1, 9, -1, -1, -1, -1 }, 1);
	rotateMatrix(meanFilter, 3, 3);

	createFilter(embossFilter, { 0, 1, 0, 0, 0, 0, 0, -1, 0 }, 1);
	rotateMatrix(embossFilter, 3, 3);
}

void movePointsToFinalMatrix(unsigned char** tempMatrix, unsigned char** finalMatrix, int start, int stop, int width)
{
	//copiaza punctele dintr-o matrice in matricea finala
	int line = 0;
	int column = 0;

	for (int x = start; x < stop; x++)
	{
		//calculam unde trebuie pus punctul in matricea finala
		int j = x % width;
		int i = x / width;

		finalMatrix[i][j] = tempMatrix[line][column];

		//luam punctele in ordine din matricea temporara
		column = column + 1;
		if (column == width)
		{
			line = line + 1;
			column = 0;
		}
	}
}

unsigned char** applyFiltersToImage(int rank, int nProcesses, std::list<std::string>& filters, imageHeader header, int multiplier, unsigned char** imageMatrix)
{
	//fiecare proces isi initializeaza filtrele
	float smoothingFilter[3][3], gaussianFilter[3][3], sharpenFilter[3][3], meanFilter[3][3], embossFilter[3][3];
	initializeAllFilters(smoothingFilter, gaussianFilter, sharpenFilter, meanFilter, embossFilter);

	//calculeaza numarul efectiv de coloane. pt pozele color este 3 * nr de coloane din fisier
	int width = header.width * multiplier;
	//calculeaza numarul de pixeli(elemente) din matrice
	int numElements = width * header.height;

	if (rank != 0)
	{
		//procesul 0 aloca matricea la citire, dar si restul proceselor trebuie sa o aloce
		imageMatrix = allocMatrix(header.height, width);
	}

	//calculeaza cate elemente ii revin procesului curent
	int chunkSize = numElements / nProcesses;
	int start = chunkSize * rank;
	int stop = chunkSize * (rank + 1);

	//in cazul in care numarul de pixeli nu se imparte direct la numarul de procese, ultimul proces va face calculele pt punctele ramase
	//nu are rost sa impartim egal restul de puncte pentru ca este foarte mic(maxim 7 pt 8 procese)
	if (rank == nProcesses - 1)
	{
		if (stop < numElements)
		{
			stop = numElements;
		}
	}

	for (auto const& filter : filters) //pentru fiecare filtru din lista
	{
		//procesul 0 trimite intreaga matrice de pixeli(matricea din fisier) catre toate celelalte procese
		//trimitem intreaga matrice pentru a avea acces usor la vecini, daca am trimite doar
		//portiunea corespunzatoare din matrice atunci nu am avea toti vecinii
		MPI_Bcast(&(imageMatrix[0][0]), width * header.height, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

		//fiecare filtru aplica produsul de convolutie doar pe punctele repartizate lui
		unsigned char** partialImageMatrix;

		if (filter.compare("smooth") == 0)
		{
			partialImageMatrix = applyFilter(imageMatrix, header, smoothingFilter, multiplier, start, stop);
		}
		else if (filter.compare("blur") == 0)
		{
			partialImageMatrix = applyFilter(imageMatrix, header, gaussianFilter, multiplier, start, stop);
		}
		else if (filter.compare("sharpen") == 0)
		{
			partialImageMatrix = applyFilter(imageMatrix, header, sharpenFilter, multiplier, start, stop);
		}
		else if (filter.compare("mean") == 0)
		{
			partialImageMatrix = applyFilter(imageMatrix, header, meanFilter, multiplier, start, stop);
		}
		else if (filter.compare("emboss") == 0)
		{
			partialImageMatrix = applyFilter(imageMatrix, header, embossFilter, multiplier, start, stop);
		}

		if (rank != 0) //trimite rezultatul catre 0
		{
			//deoarece fiecare proces calculeaza si trimite catre 0 doar portiunea lui de matrice trebuie sa determinam
			//dimensiunea matricei trimise de procesul rank catre procesul 0

			//numarul de puncte pe care a lucrat procesul rank
			int numberOfPoints = stop - start + 1;
			//calculam un nou numar de linii, tinand numarul de coloane constant, pentru simplitate
			int partialMatrixHeight = ceil((float)numberOfPoints / (float)width);

			//procesele 1 .. n vor trimite catre 0 bucata lor de matrice de dimensiune width * partialMatrixHeight
			MPI_Send(&(partialImageMatrix[0][0]), width * partialMatrixHeight, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
			//eliberam memoria ocupata de matricea partiala a fiecarui proces
			freeImageMatrix(partialImageMatrix, header);
		}

		if (rank == 0) //0 primeste toate rezultate
		{
			//procesul 0 primeste de la procesele 1...n bucatile lor de matrice si le asambleaza intr-o matrice completa
			unsigned char** assembledMatrix = allocMatrix(header.height, width);

			//prima data mutam datele din matricea partiala a procesul 0 in matricea completa calculata de 0
			movePointsToFinalMatrix(partialImageMatrix, assembledMatrix, start, stop, width);

			//pentru fiecare proces k de la 1 la n de la care 0 trebuie sa primeasca bucati de matrice
			for (int k = 1; k < nProcesses; k++)
			{
				//recalculam intervalul de puncte asociat fiecare proces, altfel ar trebui sa il trimitem catre 0
				//probabil este mai rapid sa facem iar operatiile decat sa facem 2 MPI_Send/MPI_Recv
				int chunkSize = numElements / nProcesses;
				int start = chunkSize * k;
				int stop = chunkSize * (k + 1);

				if (k == nProcesses - 1)
				{
					if (stop < numElements)
					{
						stop = numElements;
					}
				}

				int numberOfPoints = stop - start + 1;
				int partialMatrixHeight = ceil((float)numberOfPoints / (float)width);

				//aloca spatiu in procesul 0 pentru matricea temporara primita de la procesul k
				unsigned char** tempMatrix = allocMatrix(partialMatrixHeight, width);
				MPI_Recv(&(tempMatrix[0][0]), width * partialMatrixHeight, MPI_UNSIGNED_CHAR, k, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				//muta punctele calculate de procesul k in matricea finala
				movePointsToFinalMatrix(tempMatrix, assembledMatrix, start, stop, width);

				//elibereaza memoria matricei temporare, o folosim doar pt MPI_Recv
				freeImageMatrix(tempMatrix, header);
			}

			//in cazul in care rulam pe mai multe filtre, outputul de la pasul k este inputul pt pasul k + 1
			imageMatrix = assembledMatrix;
		}
	}

	return imageMatrix;
}

void writeImageFile(std::string outputFileName, imageHeader header, unsigned char** imageMatrix, int multiplier)
{
	//scrie headerul
	std::ofstream outputStream{ outputFileName, std::ios::out | std::ios::binary };
	outputStream << header.p << std::endl;
	outputStream << header.width << " " << header.height << std::endl;
	outputStream << header.maxVal << std::endl;

	//scrie matricea de pixeli
	for (int i = 0; i < header.height; i++)
	{
		outputStream.write(reinterpret_cast<const char*>(imageMatrix[i]), sizeof(unsigned char)* header.width * multiplier);
	}

	outputStream.close();
}


void freeImageMatrix(unsigned char** imageMatrix, imageHeader header)
{
	delete[] imageMatrix[0]; //elibereza matricea liniara
	delete[] imageMatrix;    //elibereaza vectorul de pointeri in matricea liniara
}

void createFilter(float matrix[3][3], std::initializer_list<int> values, int factor)
{
	//copiaza datele din initialializer_list in matrice si imparte la factor
	int k = 0;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matrix[i][j] = *(values.begin() + k);
			k = k + 1;
			matrix[i][j] = matrix[i][j] * (1.0 / factor);
		}
	}
}

void rotateMatrix(float matrix[3][3], int rows, int columns)
{
	//schimba liniile din matrice astfel incat sa realizam rotatia cu 180 de grade
	for (int i = 0; i < rows / 2; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			std::swap(matrix[i][j], matrix[rows - i - 1][columns - j - 1]);
		}
	}

	//cazul in care numarul de linii este impar
	if (rows % 2 == 1)
	{
		for (int j = 0; j < rows / 2; j++)
		{
			std::swap(matrix[rows / 2][j], matrix[rows / 2][rows - j - 1]);
		}
	}
}

void setValueIfValid(int x, int y, int n, int m, imageHeader header, int multiplier, float elemMatrix[3][3], unsigned char** imageMatrix)
{
	int height = header.height;
	int width = header.width * multiplier;

	if (x < 0 || x >= height || y < 0 || y >= width)
	{
		elemMatrix[n][m] = 0;
	}
	else
	{
		elemMatrix[n][m] = (float)imageMatrix[x][y];
	}
}

unsigned char** applyFilter(unsigned char** imageMatrix, imageHeader header, float filterMatrix[3][3], int multiplier, int startPoint, int endPoint)
{
	int numberOfPoints = endPoint - startPoint + 1;
	int width = header.width * multiplier;
	int height = ceil((float)numberOfPoints / (float)width);

	//aloca o matrice care retine exact punctele pe care procesul curent trebuie sa lucreze
	unsigned char** newImageMatrix = allocMatrix(height, width);

	//matricea in care extragem vecinii unui pixel pt a aplica filtrul
	float elemMatrix[3][3];

	int line = 0;
	int column = 0;

	for (int x = startPoint; x < endPoint; x++)
	{
		int j = x % width;
		int i = x / width;

		if (multiplier == 1)
		{
			//extrage vecinii pt o poza alb negru
			for (int n = 0; n < 3; n++)
			{
				for (int m = 0; m < 3; m++)
				{
					if (i - 1 + n < 0 || i - 1 + n >= header.height || j - 1 + m < 0 || j - 1 + m >= width)
					{
						//nu bordam efectiv matricea cu 0, ar creste dimensiunea
						//daca vrem sa accesam un vecin care nu exista, punem 0 in matricea rezultat
						elemMatrix[n][m] = 0;
					}
					else
					{
						elemMatrix[n][m] = (float)imageMatrix[i - 1 + n][j - 1 + m];
					}
				}
			}
		}
		else
		{
			//extrage vecinii pentru o poza color(tine cont ca vecinul pt canalul red este 3 puncte in stanga/dreapta)
			for (int n = 0; n < 3; n++)
			{
				for (int m = 0; m < 3; m++)
				{
					if (i - 1 + n < 0 || i - 1 + n >= header.height || j - 3 + m * 3 < 0 || j - 3 + m * 3 >= width)
					{
						//nu bordam efectiv matricea cu 0, ar creste dimensiunea
						//daca vrem sa accesam un vecin care nu exista, punem 0 in matricea rezultat
						elemMatrix[n][m] = 0;
					}
					else
					{
						elemMatrix[n][m] = (float)imageMatrix[i - 1 + n][j - 3 + m * 3];
					}
				}
			}
		}

		//calculeaza noua valoare
		float product = filterMatrix[0][0] * elemMatrix[0][0] + filterMatrix[0][1] * elemMatrix[0][1] + filterMatrix[0][2] * elemMatrix[0][2] +
			filterMatrix[1][0] * elemMatrix[1][0] + filterMatrix[1][1] * elemMatrix[1][1] + filterMatrix[1][2] * elemMatrix[1][2] +
			filterMatrix[2][0] * elemMatrix[2][0] + filterMatrix[2][1] * elemMatrix[2][1] + filterMatrix[2][2] * elemMatrix[2][2];

		//clamp
		if (product > header.maxVal)
		{
			product = header.maxVal;
		}

		if (product < 0)
		{
			product = 0;
		}

		newImageMatrix[line][column++] = (unsigned char)product;

		if (column == width)
		{
			line = line + 1;
			column = 0;
		}
	}

	return newImageMatrix;
}
