#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"
#include <stdio.h>
#include <stdlib.h>


/* Convert a Python list into the Matrix struct */
Matrix *listMatrixConverter(PyObject *data_list) {
    if (!PyList_Check(data_list)) {
        PyErr_SetString(PyExc_ValueError, "Input data must be a 2D list");
        return NULL;
    }

    Py_ssize_t rows = PyList_Size(data_list);

    /* Check the first row to determine the number of columns */
    PyObject *first_row = PyList_GetItem(data_list, 0);
    Py_ssize_t cols = PyList_Size(first_row);

    /* Create the Matrix struct */
    Matrix *mat = (Matrix*) malloc(sizeof(Matrix));
    *mat = createMatrix((int)rows, (int)cols);
    if (mat->mat == NULL) {
        fprintf(stdout, "An Error Has Occurred\n");
        free(mat);
        return NULL;
    }
    /* Fill the matrix from the Python list data */
    for (Py_ssize_t i = 0; i < rows; ++i) {
        PyObject *row = PyList_GetItem(data_list, i);
        if (!PyList_Check(row) || PyList_Size(row) != cols) {
            PyErr_SetString(PyExc_ValueError, "All rows must be lists of equal length");
            free(mat); /* Clean up in case of error */
            return NULL;
        }
        for (Py_ssize_t j = 0; j < cols; ++j) {
            PyObject *item = PyList_GetItem(row, j);
            if (!PyFloat_Check(item) && !PyLong_Check(item)) {
                PyErr_SetString(PyExc_TypeError, "All elements in the list must be numeric (float or int)");
                free(mat);
                return NULL;
            }
            mat->mat[i][j] = (float) PyFloat_AsDouble(item); /* Cast to float for storage */
        }
    }

    return mat;
}

/* Convert the Matrix struct to a Python list (2D list) */
static PyObject *matrixListConverter(Matrix result) {

    /* Get matrix dimensions */
    int rows = result.rows;
    int cols = result.cols;

    /* Create the top-level list (for rows) */
    PyObject *result_list = PyList_New(rows);
    if (result_list == NULL) {
        return NULL; /* Return NULL if memory allocation fails */
    }

    /* Fill the Python list with values from the C matrix */
    for (int i = 0; i < rows; ++i) {
        /* Create a list for each row */
        PyObject *row_list = PyList_New(cols);
        if (row_list == NULL) {
            Py_DECREF(result_list); /* Clean up in case of error */
            return NULL;
        }

        for (int j = 0; j < cols; ++j) {
            /* Convert each element to a Python float and set it in the row list */
            PyObject *value = PyFloat_FromDouble((double) result.mat[i][j]);
            if (value == NULL) {
                Py_DECREF(row_list);
                Py_DECREF(result_list);
                return NULL;
            }
            PyList_SetItem(row_list, j, value);  /* No need to `DECREF` value as PyList_SetItem steals the reference */
        }

        /* Add the row to the main list */
        PyList_SetItem(result_list, i, row_list);  /* PyList_SetItem steals the reference to row_list */
    }

    return result_list;
}

/* A mediator between C and Python, using a Python list for similarity matrix */
static PyObject *Csymnmf(PyObject *self, PyObject *args){
    PyObject *data_list_obj1;
    PyObject *data_list_obj2;

    /* Parse the Python argument to ensure we received a Python object */
    if (!PyArg_ParseTuple(args, "OO", &data_list_obj1,&data_list_obj2)) {
        return NULL;
    }

    /* Ensure the input is a 2D Python list */
    if (!PyList_Check(data_list_obj1)||(!PyList_Check(data_list_obj2))) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D list as input");
        return NULL;
    }

    /* Convert the Python list to the Matrix struct */
    Matrix *init_H = listMatrixConverter(data_list_obj1);
    Matrix *W = listMatrixConverter(data_list_obj2);
    if (init_H == NULL || W == NULL) {
        return NULL;
    }

    /* Call the sym function, which operates on the Matrix struct */
    Matrix result = converge(*init_H,*W);
    free(init_H); /* Free memory allocated for converted_mat */
    free(W); /* Free memory allocated for converted_mat */
    /* Convert the result Matrix back to a Python 2D list */

    PyObject *result_py = matrixListConverter(result);
    if (result_py == NULL) {
        return NULL;
    }

    return result_py;
}

/* A mediator between C and Python, using a Python list for similarity matrix */
static PyObject *Csym(PyObject *self, PyObject *args){
    PyObject *data_list_obj;

    /* Parse the Python argument to ensure we received a Python object */
    if (!PyArg_ParseTuple(args, "O", &data_list_obj)) {
        return NULL;
    }

    /* Ensure the input is a 2D Python list */
    if (!PyList_Check(data_list_obj)) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D list as input");
        return NULL;
    }

    /* Convert the Python list to the Matrix struct */
    Matrix *converted_mat = listMatrixConverter(data_list_obj);
    if (converted_mat == NULL) {
        return NULL;
    }

    /* Call the sym function, which operates on the Matrix struct */
    Matrix result = sym(*converted_mat);
    free(converted_mat); /* Free memory allocated for converted_mat */

    /* Convert the result Matrix back to a Python 2D list */
    PyObject *result_py = matrixListConverter(result);
    if (result_py == NULL) {
        return NULL;
    }

    return result_py;
}

/* A mediator between C and python diagonal degree Matrix */
static PyObject *Cddg(PyObject *self, PyObject *args){
    PyObject *data_list_obj;

    /* Parse the Python argument to ensure we received a Python object */
    if (!PyArg_ParseTuple(args, "O", &data_list_obj)) {
        return NULL;
    }

    /* Ensure the input is a 2D Python list */
    if (!PyList_Check(data_list_obj)) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D list as input");
        return NULL;
    }

    /* Convert the Python list to the Matrix struct */
    Matrix *converted_mat = listMatrixConverter(data_list_obj);
    if (converted_mat == NULL) {
        return NULL;
    }

    /* Call the sym function, which operates on the Matrix struct */
    Matrix result = ddg(*converted_mat);
    free(converted_mat); /* Free memory allocated for converted_mat */

    /* Convert the result Matrix back to a Python 2D list */
    PyObject *result_py = matrixListConverter(result);
    if (result_py == NULL) {
        return NULL;
    }

    return result_py;
}

/* A mediator between C and python normalized similarity Matrix */
static PyObject *Cnorm(PyObject *self, PyObject *args){
    PyObject *data_list_obj;

    /* Parse the Python argument to ensure we received a Python object */
    if (!PyArg_ParseTuple(args, "O", &data_list_obj)) {
        return NULL;
    }

    /* Ensure the input is a 2D Python list */
    if (!PyList_Check(data_list_obj)) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D list as input");
        return NULL;
    }

    /* Convert the Python list to the Matrix struct */
    Matrix *converted_mat = listMatrixConverter(data_list_obj);
    if (converted_mat == NULL) {
        return NULL;
    }

    /* Call the sym function, which operates on the Matrix struct */
    Matrix result = norm(*converted_mat);
    free(converted_mat); /* Free memory allocated for converted_mat */

    /* Convert the result Matrix back to a Python 2D list */
    PyObject *result_py = matrix_list_converter(result);
    if (result_py == NULL) {
        return NULL;
    }

    return result_py;
}


static PyMethodDef SymnmfLib_FunctionsTable[] ={
    {
        "Csym",
        Csym,
        METH_VARARGS,
        "returns the similarity matrix using python API" 
    },{
        "Cddg",
        Cddg,
        METH_VARARGS,
        "returns the diagonal degree matrix using python API" 
    },{
        "Cnorm",
        Cnorm,
        METH_VARARGS,
        "returns the  normalized similarity  matrix using python API" 
    },{
        "Csymnmf",
        Csymnmf,
        METH_VARARGS,
        "returns the symnmf matrix using python API" 
    },{
        NULL, NULL, 0, NULL
        }
};

static struct PyModuleDef SymnmfLib_Module = {
    PyModuleDef_HEAD_INIT,
    "mysymnmf",
    "mysymnmf Python wrapper for custom C extension library.",
    -1,
    SymnmfLib_FunctionsTable
};

PyMODINIT_FUNC PyInit_mysymnmf(void) {
    return PyModule_Create(&SymnmfLib_Module);
}

