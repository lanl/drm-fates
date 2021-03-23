#include <string>
using namespace std;


int** allocate_2d_int(int dim1, int dim2)
{
    //dimension 1
    int **arr = new int*[dim1];
    //dimension 2
    for (int i=0; i<dim1; i++) 
        arr[i] = new int[dim2];
    return arr;
}

double** allocate_2d_double(int dim1, int dim2)
{
    //dimension 1
    double **arr = new double*[dim1];
    //dimension 2
    for (int i=0; i<dim1; i++) 
        arr[i] = new double[dim2];
    return arr;
}

float** allocate_2d_float(int dim1, int dim2)
{
    //dimension 1
    float **arr = new float*[dim1];
    //dimension 2
    for (int i=0; i<dim1; i++) 
        arr[i] = new float[dim2];
    return arr;
}

int*** allocate_3d_int(int dim1, int dim2, int dim3)
{
    //dimension 1
    int ***arr = new int**[dim1];
    //dimension 2
    for (int i=0; i<dim1; i++) 
        arr[i] = new int*[dim2];
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            arr[i][j] = new int[dim3];
    return arr;
}

double*** allocate_3d_double(int dim1, int dim2, int dim3)
{
    //dimension 1
    double ***arr = new double**[dim1];
    //dimension 2
    for (int i=0; i<dim1; i++) 
        arr[i] = new double*[dim2];
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            arr[i][j] = new double[dim3];
    return arr;
}

float*** allocate_3d_float(int dim1, int dim2, int dim3)
{
    //dimension 1
    float ***arr = new float**[dim1];
    //dimension 2
    for (int i=0; i<dim1; i++) 
        arr[i] = new float*[dim2];
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            arr[i][j] = new float[dim3];
    return arr;
}

int**** allocate_4d_int(int dim1, int dim2, int dim3, int dim4)
{
    //dimension 1
    int ****arr = new int***[dim1];
    //dimension 2
    for (int i=0; i<dim1; i++) 
        arr[i] = new int**[dim2];
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            arr[i][j] = new int*[dim3];
    //dimension 4
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
                arr[i][j][k] = new int[dim4];
    return arr;
}

float**** allocate_4d_float(int dim1, int dim2, int dim3, int dim4)
{
    //dimension 1
    float ****arr = new float***[dim1];
    //dimension 2
    for (int i=0; i<dim1; i++) 
        arr[i] = new float**[dim2];
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            arr[i][j] = new float*[dim3];
    //dimension 4
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
                arr[i][j][k] = new float[dim4];
    return arr;
}

double**** allocate_4d_double(int dim1, int dim2, int dim3, int dim4)
{
    //dimension 1
    double ****arr = new double***[dim1];
    //dimension 2
    for (int i=0; i<dim1; i++) 
        arr[i] = new double**[dim2];
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            arr[i][j] = new double*[dim3];
    //dimension 4
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
                arr[i][j][k] = new double[dim4];
    return arr;
}

void deallocate_2d_int(int **(&arr), int dim1, int dim2)
{
    //dimension 2
    for (int i=0; i<dim1; i++) 
    {
        delete[] arr[i];
        arr[i] = NULL;
    }
    //dimension 1
    delete[] arr;
    arr = NULL;

    return;
}

void deallocate_2d_double(double **(&arr), int dim1, int dim2)
{
    //dimension 2
    for (int i=0; i<dim1; i++) 
    {
        delete[] arr[i];
        arr[i] = NULL;
    }
    //dimension 1
    delete[] arr;
    arr = NULL;

    return;
}

void deallocate_2d_float(float **(&arr), int dim1, int dim2)
{
    //dimension 2
    for (int i=0; i<dim1; i++) 
    {
        delete[] arr[i];
        arr[i] = NULL;
    }
    //dimension 1
    delete[] arr;
    arr = NULL;

    return;
}

void deallocate_3d_int(int ***(&arr), int dim1, int dim2, int dim3)
{
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
        {
            delete[] arr[i][j];
            arr[i][j] = NULL;
        }
    //dimension 2
    for (int i=0; i<dim1; i++) 
    {
        delete[] arr[i];
        arr[i] = NULL;
    }
    //dimension 1
    delete[] arr;
    arr = NULL;

    return;
}

void deallocate_3d_double(double ***(&arr), int dim1, int dim2, int dim3)
{
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
        {
            delete[] arr[i][j];
            arr[i][j] = NULL;
        }
    //dimension 2
    for (int i=0; i<dim1; i++) 
    {
        delete[] arr[i];
        arr[i] = NULL;
    }
    //dimension 1
    delete[] arr;
    arr = NULL;

    return;
}

void deallocate_3d_float(float ***(&arr), int dim1, int dim2, int dim3)
{
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
        {
            delete[] arr[i][j];
            arr[i][j] = NULL;
        }
    //dimension 2
    for (int i=0; i<dim1; i++) 
    {
        delete[] arr[i];
        arr[i] = NULL;
    }
    //dimension 1
    delete[] arr;
    arr = NULL;

    return;
}

void deallocate_4d_int(int ****(&arr), int dim1, int dim2, int dim3, int dim4)
{
    //dimension 4
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
            {
                delete arr[i][j][k];
                arr[i][j][k] = NULL;
            }
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
        {
            delete arr[i][j];
            arr[i][j] = NULL;
        }
    //dimension 2
    for (int i=0; i<dim1; i++) 
    {
        delete arr[i];
        arr[i] = NULL;
    }
    //dimension 1
    delete arr;
    arr = NULL;

    return;
}

void deallocate_4d_float(float ****(&arr), int dim1, int dim2, int dim3, int dim4)
{
    //dimension 4
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
            {
                delete arr[i][j][k];
                arr[i][j][k] = NULL;
            }
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
        {
            delete arr[i][j];
            arr[i][j] = NULL;
        }
    //dimension 2
    for (int i=0; i<dim1; i++) 
    {
        delete arr[i];
        arr[i] = NULL;
    }
    //dimension 1
    delete arr;
    arr = NULL;

    return;
}

void deallocate_4d_double(double ****(&arr), int dim1, int dim2, int dim3, int dim4)
{
    //dimension 4
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
            for (int k=0; k<dim3; k++)
            {
                delete arr[i][j][k];
                arr[i][j][k] = NULL;
            }
    //dimension 3
    for (int i=0; i<dim1; i++)
        for (int j=0; j<dim2; j++)
        {
            delete arr[i][j];
            arr[i][j] = NULL;
        }
    //dimension 2
    for (int i=0; i<dim1; i++) 
    {
        delete arr[i];
        arr[i] = NULL;
    }
    //dimension 1
    delete arr;
    arr = NULL;

    return;
}





