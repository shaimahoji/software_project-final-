#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"
#include <numpy/arrayobject.h>
#include <stdio.h>
#include <stdlib.h>


// Convert a NumPy array into the Matrix struct
Matrix *array_matrix_converter(PyArrayObject *data_array) {
    if (data_array == NULL) {
        return NULL;
    }

    int nd = PyArray_NDIM(data_array);//number of dim
    npy_intp *dim = PyArray_SHAPE(data_array); // the shape of matrix 

    if (nd != 2) {
        PyErr_SetString(PyExc_ValueError, "Input array must be 2D");
        return NULL;
    }

    npy_intp rows = dim[0];
    npy_intp cols = dim[1];


    double *data = (double*) PyArray_DATA(data_array);

    Matrix *mat = (Matrix*) malloc(sizeof(Matrix));
    *mat = create_matrix((int)rows, (int)cols);
    if (mat->mat == NULL) {
        fprintf(stdout, "An Error Has Occurred\n");
        free(mat);
        return NULL;
    }

    // Fill the matrix
    for (npy_intp i = 0; i < rows; ++i) {
        for (npy_intp j = 0; j < cols; ++j) {
            mat->mat[i][j] = (float)data[i * cols + j]; // Cast to float for storage
        }
    }
        return mat;
    }

// convert the Matrix struct we created in C to a numpy array (python) 
static PyObject *matrix_array_converter(Matrix result){
    // take the matrix dimensions so we can create numpy array that fits
    npy_intp rows = result.rows;
    npy_intp cols = result.cols;
    npy_intp dims[2] = {rows, cols};

    // creat numpy array so we can return the results and print it in python
    PyObject *result_array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    if (result_array == NULL) {
        return NULL;}

    double *result_data = (double*) PyArray_DATA((PyArrayObject*) result_array);

    // fill the result array with the results in the "C" matrix
    for (npy_intp i = 0; i < rows; ++i) {
        for (npy_intp j = 0; j < cols; ++j){
            result_data[i * cols + j] = result.mat[i][j];
        }
    }
    return result_array;
}

//a mediator between C and python similarity matrix
static PyObject *Csym(PyObject *self, PyObject *args){
    PyObject *data_array_obj;
    // check if we got the right variables
    if (!PyArg_ParseTuple(args, "O" ,&data_array_obj)) {
        return NULL;
    }
   
    //ensure that this is a numpy array
    PyArrayObject *data_array = (PyArrayObject*) PyArray_FROM_OTF(data_array_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
   
    //clear ifits not the right input
    if (data_array == NULL ) {
        Py_XDECREF(data_array);
        return NULL;
    }
    
    // convert the numpy array to a list so we can pass it as an arguments for the sym function
    Matrix *converted_mat= array_matrix_converter(data_array);
    if (converted_mat == NULL) {
        Py_DECREF(data_array);
        return NULL;
    }

    Matrix result = sym(converted_mat[0]);
    PyObject *result_py = matrix_array_converter(result);

    if (result_py == NULL){
        Py_XDECREF(data_array);
        return NULL;}    
    
    // clear the lists we got as inputs
    Py_DECREF(data_array);

    return result_py;
}

//a mediator between C and python diagonal degree Matrix
static PyObject *Cddg(PyObject *self, PyObject *args){
    PyObject *data_array_obj;

    // check if we got the right variables
    if (!PyArg_ParseTuple(args, "O" ,&data_array_obj)) {
        return NULL;
    }
    //ensure that this is a numpy array
    PyArrayObject *data_array = (PyArrayObject*) PyArray_FROM_OTF(data_array_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
   
    //clear ifits not the right input
    if (data_array == NULL ) {
        Py_XDECREF(data_array);
        return NULL;
    }
    
    // convert the numpy array to a list so we can pass it as an arguments for the sym function
    Matrix *converted_mat= array_matrix_converter(data_array);
    if (converted_mat->mat == NULL) {
        Py_DECREF(data_array);
        return NULL;
    }

    Matrix result = ddg(converted_mat[0]);
    PyObject *result_py = matrix_array_converter(result);
    
    if (result_py == NULL){
        Py_XDECREF(data_array);
        return NULL;}

    // clear the lists we got as inputs
    Py_DECREF(data_array);

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
    import_array();  // Initialize NumPy C API
    return PyModule_Create(&SymnmfLib_Module);
}

