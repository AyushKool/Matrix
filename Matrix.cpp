#include <iostream>
#include <math.h>

using namespace std;

class Matrix
{
private:
    int **matrix, n, m;

public:
    //     Matrix() : Matrix(1, 1) {}
    Matrix(int n = 2) : Matrix(n, n) {}
    Matrix(int n, int m)
    {
        this->n = n;
        this->m = m;
        matrix = (int **)calloc(sizeof(int *), n);
        for (int i = 0; i < n; i++)
            matrix[i] = (int *)calloc(sizeof(int), m);
    }
    Matrix(int **arr, int n, int m)
    {
        this->n = n;
        this->m = m;
        matrix = (int **)calloc(sizeof(int *), n);
        for (int i = 0; i < n; i++)
            matrix[i] = (int *)calloc(sizeof(int), m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                matrix[i][j] = arr[i][j];
    }
    int row()
    {
        return n;
    }
    int col()
    {
        return m;
    }
    void init()
    {
        cout << "Enter elements: ";
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                cin >> matrix[i][j];
    }
    void display()
    {
        cout << endl;
        bool isUpper = true, isLower = true;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                printf("%5d", matrix[i][j]); //put -%5d for right align and no minus for left align
                if (i < j && matrix[i][j] != 0)
                    isLower = false;
                if (i > j && matrix[i][j] != 0)
                    isUpper = false;
            }

            cout << endl;
        }
        if (isUpper)
            cout << "This is Upper Triangular\n";
        if (isLower)
            cout << "This is Lower Triangular\n";
        if (!isUpper && !isLower)
            cout << "This is neither Upper nor Lower\n";
        cout << endl;
    }

    int getElem(int i, int j)
    {
        return matrix[i][j];
    }

    //Returns 1 if equal, 0 if not
    int operator==(Matrix other)
    {
        if (row() != other.row() || col() != other.col())
            return 0;
        for (int i = 0; i < row(); i++)
            for (int j = 0; j < col(); j++)
                if (matrix[i][j] != other.getElem(i, j))
                    return 0;
        return 1;
    }

    int operator!=(Matrix other)
    {
        return !(*this == other);
    }

    Matrix operator*(Matrix other)
    {
        if (m != other.row())
            throw "Error. Cannot multiply matrices.";
        else
        {
            int **ans;
            ans = (int **)calloc(sizeof(int *), n);
            for (int i = 0; i < n; i++)
                ans[i] = (int *)calloc(sizeof(int), other.col());

            for (int i = 0; i < n; i++)
                for (int j = 0; j < other.col(); j++)
                    for (int k = 0; k < m; k++)
                        ans[i][j] += matrix[i][k] * other.getElem(k, j);
            return Matrix(ans, n, other.col());
        }
    }

    // int operator*(Matrix other)
    // {
    //     return other * *this
    // }

    Matrix minor(Matrix mat, int r, int c)
    {
        int minor_row = mat.row() - 1, minor_col = mat.col() - 1;
        int **minor = (int **)calloc(sizeof(int *), minor_row);
        for (int i = 0; i < minor_row; i++)
            minor[i] = (int *)calloc(sizeof(int), minor_col);

        int r2, c2;
        for (int i = 0, r2 = 0; i < mat.row(); i++)
        {
            if (i == r)
                continue;
            for (int j = 0, c2 = 0; j < mat.col(); j++)
            {
                if (j == c)
                    continue;
                minor[r2][c2] = mat.getElem(i, j);
                c2++;
            }
            r2++;
        }
        return Matrix(minor, minor_row, minor_col);
    }

    void swapRow(int a, int b)
    {
        int *tmp = matrix[a];
        matrix[a] = matrix[b];
        matrix[b] = tmp;
    }

    void toUpperTri()
    {
        for (int j = 0, piv = 0; j < m; ++j, piv = j)
        {
            while (piv < n && matrix[piv][j] == 0)
                piv++;
            if (piv >= n)
                continue;
            else if (piv != j)
                swapRow(piv, j);
            piv = matrix[j][j];
            for (int i = j + 1; i < n; ++i)
            {
                float factor = (float)matrix[i][j] / piv;
                //printf("Row %d - Row %d\n", i, j);
                for (int k = j; k < m; ++k)
                    matrix[i][k] -= factor * matrix[j][k];
            }
        }
    }

    void transpose()
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = i + 1; j < m; ++j)
            {
                matrix[i][j] ^= matrix[j][i];
                matrix[j][i] ^= matrix[i][j];
                matrix[i][j] ^= matrix[j][i];
            }
        }
    }

    void toLowerTri()
    {
        transpose();
        toUpperTri();
        transpose();
    }

    //find it using upper triangular matrix
    float detCalc(Matrix mat)
    {
        //static int calls = 1;
        //     cout << "Call " << calls << endl;
        //     calls++;
        int rows = mat.row();
        if (rows == 1)
            return mat.getElem(0, 0);
        else if (rows == 2)
        {
            float calc = mat.getElem(0, 0) * mat.getElem(1, 1) - mat.getElem(0, 1) * mat.getElem(1, 0);
            cout << "Det: " << calc << endl;
            return calc;
        }
        float sum = 0;
        int sign = 1;
        Matrix min;
        for (int j = 0; j < mat.col(); j++, sign *= -1)
        {
            min = minor(mat, 0, j);
            cout << "Minor of: (0, " << j << ")" << endl;
            min.display();
            sum += sign * mat.getElem(0, j) * detCalc(min);
        }
        return sum;
    }

    float det()
    {
        if (n != m || n == 0)
            throw "Cannot determine Determinant.";
        else
            return detCalc(*this);
    }

    Matrix operator+(Matrix other)
    {
        if (n != other.row() || m != other.col())
            throw "Error. Cannot add matrices.";
        else
        {
            int **ans;
            ans = (int **)calloc(sizeof(int *), n);
            for (int i = 0; i < n; i++)
                ans[i] = (int *)calloc(sizeof(int), m);

            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                    ans[i][j] = matrix[i][j] + other.getElem(i, j);
            return Matrix(ans, n, other.col());
        }
    }

    Matrix operator-()
    {
        return *this * -1;
    }

    Matrix operator-(Matrix other)
    {
        return *this + (-other);
    }
    // Matrix uppertri()
    // {
    //     int ** ans = (int **)calloc(sizeof(int *), n);
    //     for (int i = 0; i < n; i++)
    //         ans[i] = (int *)calloc(sizeof(int), m);

    //     int piv, factor = 1;
    //     for(piv = 0; ans[0][piv] == 0; piv++);
    //     for(int i=1; i<n; i++)
    //     {
    //         factor = ans[i][piv]/ans[0][piv];
    //     }
    // }

    Matrix operator*(int scalar)
    {
        int **ans;
        ans = (int **)calloc(sizeof(int *), n);
        for (int i = 0; i < n; i++)
            ans[i] = (int *)calloc(sizeof(int), m);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                ans[i][j] = matrix[i][j];
                ans[i][j] *= scalar;
            }
        }
        return Matrix(ans, n, m);
    }

    Matrix friend operator*(int scalar, Matrix other);
};

Matrix operator*(int scalar, Matrix other)
{
    return other.operator*(scalar);
}

int main()
{
    //Matrix mat()//declares function with no arguments and Matrix return type
    // Matrix mat1(3), mat2(2);
    // mat1.init();
    // (2 * mat1).display();
    // mat2.init();

    // (2 * mat1 + mat2).display();
    // (mat1 - mat2).display();

    // if(mat1 == mat2)
    //     cout << "Equal";
    // else
    //     cout << "Not Equal";
    //cout << endl << "Det: " <<  mat1.det();
    // Matrix ans;
    // try
    // {
    //     ans = mat1 * mat2;
    // }
    // catch (const char *msg)
    // {
    //     cerr << msg << endl;
    // }
    // ans.display();
    // cout << endl << "Det:" << mat1.det() << endl;
    // cout << endl << "Det:" << (2*mat1).det() << endl;

    Matrix mat(3);
    mat.init();
    mat.display();

    mat.toUpperTri();
    //mat.toLowerTri();
    mat.display();
}