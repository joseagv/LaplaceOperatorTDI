/***** Librerias  *****/
#include <stdio.h>
#include <C_General.hpp>
#include <C_Trace.hpp>
#include <C_File.hpp>
#include <C_Arguments.hpp>
#include <C_Matrix.hpp>
#include <C_Image.hpp>
#include <iostream>
#include <cmath>
#include <conio.h>
#include <string>
#include <limits>
#include <wchar.h>
#include <locale.h>

/***** Variables *****/
//Variables de posicion y manejo de la matriz
int firstRow, lastRow, firstCol, lastCol, nRows, nCols, ruido;
//Variables imagen
C_Image imagen, imagenGauss, imagenLaplace, imagenPasoPorCero, imagenFinal;
//Variable del umbral
int MaxUmbral;

/***** Metodos *****/
//Init carga los valores de la imagen en las variables
void init(int mainMaxUmbral) {
	//Calculo de las primeras y ultimas filas y columnas de la imagen
	firstRow = imagen.FirstRow();
	lastRow = imagen.LastRow();
	firstCol = imagen.FirstCol();
	lastCol = imagen.LastCol();
	//Total de filas y columnas de la imagen
	nRows = imagen.RowN();
	nCols = imagen.ColN();
	MaxUmbral = mainMaxUmbral;
}

//Filtro Gauss 3x3
// |1|2|1|
// |2|4|2|
// |1|2|1|
C_Matrix filtroGauss3x3() {
	C_Matrix filtro = C_Matrix(0, 2, 0, 2, 0);
	filtro(0, 0) = 1; filtro(0, 1) = 2; filtro(0, 2) = 1;
	filtro(1, 0) = 2; filtro(1, 1) = 4; filtro(1, 2) = 2;
	filtro(2, 0) = 1; filtro(2, 1) = 2; filtro(1, 2) = 1;
	return filtro;
}

//Filtro Gauss 5x5
// |1| 4| 7| 4|1|
// |4|20|33|20|4|
// |7|33|55|33|7|
// |4|20|33|20|4|
// |1| 4| 7| 4|1|
C_Matrix filtroGauss5x5() {
	C_Matrix filtro = C_Matrix(0, 4, 0, 4, 0);
	filtro(0, 0) = 1; filtro(0, 1) = 4; filtro(0, 2) = 7; filtro(0, 3) = 4; filtro(0, 4) = 1;
	filtro(1, 0) = 4; filtro(1, 1) = 20; filtro(1, 2) = 33; filtro(1, 3) = 20; filtro(1, 4) = 4;
	filtro(2, 0) = 7; filtro(2, 1) = 33; filtro(2, 2) = 55; filtro(2, 3) = 33; filtro(2, 4) = 7;
	filtro(1, 0) = 4; filtro(1, 1) = 20; filtro(1, 2) = 33; filtro(1, 3) = 20; filtro(1, 4) = 4;
	filtro(0, 0) = 1; filtro(0, 1) = 4; filtro(0, 2) = 7; filtro(0, 3) = 4; filtro(0, 4) = 1;
	return filtro;
}

//Implementación del filtro Laplaciano3x3 con máscara4
// |-1|-1|-1|
// |-1| 8|-1|
// |-1|-1|-1|
C_Matrix filtroLaplace() {
	C_Matrix filtro = C_Matrix(0, 2, 0, 2, 0);
	filtro(0, 0) = -1; filtro(0, 1) = -1; filtro(0, 2) = -1;
	filtro(1, 0) = -1; filtro(1, 1) = 8; filtro(1, 2) = -1;
	filtro(2, 0) = -1; filtro(2, 1) = -1; filtro(2, 2) = -1;
	return filtro;
}
//Reflect, se utiliza para reposicionar los indices de la matriz de convolucion para tratar los bordes
int reflect(int M, int x)
{
	if (x < 0) //Si la posicion se sale del rango minimo recolocamos
	{
		return -x - 1;
	}
	if (x >= M)//Si la posicion se sale del rango maximo recolocamos
	{
		return 2 * M - x - 1;
	}
	return x;
}

//Aplicar Filtro
C_Matrix AplicarFiltro(C_Image imagen, C_Matrix filtro) {
	//Variable de salida
	C_Matrix salida = imagen;
	//Ponemos sus valores a cero que es el color negro
	salida.SetValue(0);
	//Sumatoria
	int sum;
	//Guardamos el tamano del filtro.
	int fSizeX = filtro.ColN();
	int fSizeY = filtro.RowN();
	//Guardamos el centro de los ejes del filtro
	int fCenterX = fSizeX / 2;
	int fCenterY = fSizeY / 2;
	//Guardamos donde se tiene que empezar aplicar el filtro
	int rowAux, colAux;

	//Recorremos la imagen
	for (int i = firstRow; i < nRows; i++) {
		for (int j = firstCol; j < nCols; j++) {
			sum = 0;
			//Recorremos el filtro
			for (int m = 0; m < fSizeY; m++) {
				for (int n = 0; n < fSizeX; n++) {
					//Posicion donde empieza la convolucion
					rowAux = i + m - fCenterY;
					colAux = j + n - fCenterX;
					//Comprobamos que los pixeles estan dentro de los limites
					if (rowAux >= firstRow && rowAux < nRows && colAux >= firstCol && colAux < nCols) {
						int a = imagen(rowAux, colAux);
						int b = filtro(m, n);
						sum += a * b;
					}//Como la posicion donde se realizaria la convolucion esta fuera del limite buscamos el pixel mas cercano que si lo este
					else {
						colAux = reflect(nCols, j - n);
						rowAux = reflect(nRows, i - m);
						//Comprobamos que la nueva posicion esta dentro de los limites
						if (rowAux >= firstRow && rowAux < nRows && colAux >= firstCol && colAux < nCols) {
							int a = imagen(rowAux, colAux);
							int b = filtro(m, n);
							sum += a * b;
						}
						else {//En caso de que se de la situacion de que sigue sin entrar en las condiciones ponemos la suma del pixel a 0 
							//y salimos de la iteracion actual
							sum = 0;
							m = fSizeY;
							n = fSizeX;
						}
					}
				}
			}//Sacamos el valor absoluto de la suma
			salida(i, j) = sum;
		}
	}
	return salida;
}

//Paso por cero
// |p1|p2|p3|
// |p4|p0|p5|
// |p6|p7|p8|
void PasoPorCero(C_Image imagenLaplace, int maxTreshold) {
	//Valores de los pixeles a tratar
	double p0, p1, p2, p3, p4, p5, p6, p7, p8;
	//Copiamos la imagen de laplace y le damos valor de 255
	imagenPasoPorCero = imagenLaplace;
	imagenPasoPorCero.SetValue(255);
	//Recorremos la imagen
	for (int i = firstRow + 1; i < lastRow; i++) {
		for (int j = firstCol + 1; j < lastCol; j++) {
			//Guardamos los valores
			p0 = imagenLaplace(i, j);
			p1 = imagenLaplace(i - 1, j - 1);
			p2 = imagenLaplace(i - 1, j);
			p3 = imagenLaplace(i - 1, j + 1);
			p4 = imagenLaplace(i, j - 1);
			p5 = imagenLaplace(i, j + 1);
			p6 = imagenLaplace(i + 1, j - 1);
			p7 = imagenLaplace(i + 1, j);
			p8 = imagenLaplace(i + 1, j + 1);
			// Comprobamos los pixeles candidatos
			if (p4 - p5 > maxTreshold) {
				imagenPasoPorCero(i, j) = 0;
			}
			else if (p2 - p7 > maxTreshold) {
				imagenPasoPorCero(i, j) = 0;
			}
			else if (p1 - p8 > maxTreshold) {
				imagenPasoPorCero(i, j) = 0;
			}
			else if (p3 - p6 > maxTreshold) {
				imagenPasoPorCero(i, j) = 0;
			}
			else {
				imagenPasoPorCero(i, j) = 255;
			}
		}
	}
}

//Resaltar Bordes
void ResaltarBordes(C_Image imagen, C_Image imagenPasoPorcero) {
	//Creamos una copia de la imagen original
	imagenFinal = imagen;
	//Modificamos la paleta de la imagen haciendo que el valor 0 corresponda a un color rojo
	imagenFinal.palette(0, C_RED) = 255;
	imagenFinal.palette(0, C_GREEN) = 0;
	imagenFinal.palette(0, C_BLUE) = 0;
	//Recorremos la imagen de paso por cero y los pixeles que sean borde los reflejamos en la imagen final de color rojo
	for (int i = firstRow; i < lastRow; i++) {
		for (int j = firstCol; j < lastCol; j++) {
			if (imagenPasoPorcero(i, j)  == 0) {
				imagenFinal(i, j) = 0;
			}
		}
	}
}

/***** Main *****/
//Main
int main(int argc, char **argv)
{
	int filGauss, umbral;
	char nombreImagen[60];
	cout << "Introduzca el nombre de la imagen.\n";
	cin.getline(nombreImagen, 60);
	imagen.ReadBMP(nombreImagen);
	if (imagen.Fail()) {
		cout << "No se ha encontrado la imagen.\n";
		system("pause");
		exit;
	}
	else {
		cout << "La imagen ha cargado.\n";
		system("pause");
		do {
			cout << "Indique la mascara Gauss a utilizar:\n";
			cout << "1 - Mascara 3x3\n2 - Mascara 5x5\n";
			cin >> filGauss;
			cout << "Indique el Umbral Maximo(0, 255): \n";
			cin >> umbral;
			if (filGauss < 1 && filGauss > 2) {
				cout << "Indique una mascara valida\n";
			}
			if (umbral < 0 || umbral > 255) {
				cout << "Indique un umbral valido\n";
			}
		} while (umbral < 0 || umbral > 255 || filGauss < 1 && filGauss > 2);

	}
	init(umbral);
	cout << "\nGenerando Imagen Gauss\n";
	if (filGauss == 1) {
		imagenGauss = C_Image(AplicarFiltro(imagen, filtroGauss3x3()));
		imagenGauss.WriteBMP("ImagenFiltroGaussiano.bmp");
	}
	if (filGauss == 2) {
		imagenGauss = C_Image(AplicarFiltro(imagen, filtroGauss5x5()));
		imagenGauss.WriteBMP("ImagenFiltroGaussiano.bmp");
	}
	imagenGauss = C_Image(AplicarFiltro(imagen, filtroGauss3x3()));
	imagenGauss.WriteBMP("ImagenFiltroGaussiano.bmp");
	cout << "\nGenerando Imagen Laplace\n";
	imagenLaplace = C_Image(AplicarFiltro(imagenGauss, filtroLaplace()));
	imagenLaplace.WriteBMP("ImagenFiltroLaplace.bmp");
	cout << "\nGenerando Imagen Paso por Cero\n";
	PasoPorCero(imagenLaplace, umbral);
	imagenPasoPorCero.WriteBMP("ImagenPasoPorCero.bmp");
	cout << "\nGenerando Imagen Final\n";
	ResaltarBordes(imagen, imagenPasoPorCero);
	imagenFinal.WriteBMP("ImagenFinal.bmp");
	cout << "\nEl proceso se ha realizado con exito\n";
	system("pause");
	exit;
}