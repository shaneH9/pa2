#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void rowMult(double **mtrx, int row, int colLength, double mult, char operation);
void augment(double **mtrx, int rowL, int colL);
void transposeMatrix(double **matrix, int rows, int cols, double ***transposedMatrix);
void invert(double **m, int row, double ***inversionResult);
void multiplyMatrices(double **matrixA, int rowsA, int colsA, double **matrixB, int rowsB, int colsB, double ***resultMatrix);
double **initializeMatrix(int rows, int cols);
void freeMtrx(double **mtrx, int rows);
void printMatrix(int rows, int cols, double **matrix);
void divRows(int rows, int currRow, double **m, double ***n);

int main(int argc, char *argv[])
{
    char name[1024];
    FILE *tFile = fopen(argv[1], "r");
    fscanf(tFile, "%s", name);
    double **w;
    int wRows;
    int wCols;

    if (tFile == NULL)
    {
        printf("null file");
        return EXIT_FAILURE;
    }

    if (strcmp(name, "train") == 0)
    {
        double c1;
        double r1;
        fscanf(tFile, "%lf", &c1);
        fscanf(tFile, "%lf", &r1);
        c1 = c1 + 1;

        double **mtrxY = malloc(r1 * sizeof(double *));
        for (int i = 0; i < r1; i++)
        {
            mtrxY[i] = malloc(sizeof(double));
        }

        double **mtrX = initializeMatrix(r1, c1);

        for (int x = 0; x < r1; x++)
        {
            mtrX[x][0] = 1.0;
            for (int y = 1; y < c1; y++)
            {
                fscanf(tFile, "%lf", &mtrX[x][y]);
            }
            fscanf(tFile, "%lf", &mtrxY[x][0]);
        }

        double **transposedMX; // mtrX^T
        transposeMatrix(mtrX, r1, c1, &transposedMX);

        double **transposeTimesX; // x^T * X
        multiplyMatrices(transposedMX, c1, r1, mtrX, r1, c1, &transposeTimesX);

        freeMtrx(mtrX, r1);

        double **invertedMatrix; //(x^T*X)^-1
        invert(transposeTimesX, c1, &invertedMatrix);

        freeMtrx(transposeTimesX, c1);

        double **invertTimesTranspose; //(x^T*X)^-1 * X^T
        multiplyMatrices(invertedMatrix, c1, c1, transposedMX, c1, r1, &invertTimesTranspose);

        freeMtrx(invertedMatrix, c1);

        // compute weight
        multiplyMatrices(invertTimesTranspose, c1, r1, mtrxY, r1, 1, &w);
        wRows = c1;
        wCols = 1;

        freeMtrx(invertTimesTranspose, c1);
        freeMtrx(transposedMX, c1);
        freeMtrx(mtrxY, r1);
    }
    else
    {
        printf("%s\n", "file read error silly billy");
        fclose(tFile);
        return EXIT_FAILURE;
    }
    fclose(tFile);
    FILE *dFile = fopen(argv[2], "r");
    fscanf(dFile, "%s", name);
    if (strcmp(name, "data") == 0)
    {
        double c1;
        double r1;

        fscanf(dFile, "%lf", &c1);
        fscanf(dFile, "%lf", &r1);
        c1 = c1 + 1;

        double **mtrX = initializeMatrix(r1, c1);

        for (int x = 0; x < r1; x++)
        {
            mtrX[x][0] = 1;
            for (int y = 1; y < c1; y++)
            {
                fscanf(dFile, "%lf", &mtrX[x][y]);
            }
        }

        double **pwices;
        multiplyMatrices(mtrX, r1, c1, w, wRows, wCols, &pwices);
        freeMtrx(mtrX, r1);
        for (int x = 0; x < r1; x++)
        {
            printf("%.0f\n", pwices[x][0]);
        }

        freeMtrx(pwices, r1);
        freeMtrx(w, wRows);
    }
    else
    {
        printf("%s\n", "file read error silly billy");
        fclose(dFile);
        return EXIT_FAILURE;
    }
    fclose(dFile);
    return EXIT_SUCCESS;
}

void transposeMatrix(double **matrix, int rows, int cols, double ***transposedMatrix)
{
    // Allocate memory for the transposed matrix
    *transposedMatrix = initializeMatrix(cols, rows);

    // Fill the transposed matrix
    for (int x = 0; x < rows; x++)
    {
        for (int y = 0; y < cols; y++)
        {
            (*transposedMatrix)[y][x] = matrix[x][y];
        }
    }
}

void multiplyMatrices(double **matrix1, int row1, int col1, double **matrix2, int row2, int col2, double ***product)
{

    *product = initializeMatrix(row1, col2);

    for (int x = 0; x < row1; x++)
    {
        for (int y = 0; y < col2; y++)
        {
            (*product)[x][y] = 0;
            for (int k = 0; k < col1; k++)
            {
                (*product)[x][y] += matrix1[x][k] * matrix2[k][y];
            }
        }
    }
}

void invert(double **m, int row, double ***n)
{

    *n = initializeMatrix(row, row);
    augment(*n, row, row); // transforms into 1s and 0s assuming square matrix

    for (int x = 0; x < row; x++)
    {
        divRows(row,x,m,n);
        // subtracting the row from the other rows
        for (int k = 0; k < row; k++)
        {
            if (k != x)
            {
                double pivot2 = m[k][x];
                for (int y = 0; y < row; y++)
                {
                    m[k][y] -= pivot2 * m[x][y];
                    (*n)[k][y] -= pivot2 * (*n)[x][y];
                }
            }
        }
    }
}

void divRows(int rows, int currRow, double **m, double ***n)
{
    double pivot1 = m[currRow][currRow];
        for (int z = 0; z < rows; z++)
        {
            m[currRow][z] /= pivot1;
            (*n)[currRow][z] /= pivot1;
        }
}

void augment(double **mtrx, int rowL, int colL)
{
    for (int x = 0; x < rowL; x++)
    {
        for (int y = 0; y < colL; y++)
        {
            if (x == y)
            {
                mtrx[x][y] = 1.0;
            }
            else
            {
                mtrx[x][y] = 0.0;
            }
        }
    }
}

double **initializeMatrix(int rows, int cols)
{
    double **matrix = (double **)malloc(rows * sizeof(double *));
    for (int x = 0; x < rows; x++)
    {
        matrix[x] = (double *)malloc(cols * sizeof(double));
    }

    for (int x = 0; x < rows; x++)
    {
        for (int y = 0; y < cols; y++)
        {
            matrix[x][y] = 0.0;
        }
    }
    return matrix;
}

void freeMtrx(double **mtrx, int rows)
{
    for (int x = 0; x < rows; x++)
    {
        free(mtrx[x]);
    }
    free(mtrx);
}

void printMatrix(int rows, int cols, double **matrix)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf("%lf ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n\n\n");
}
