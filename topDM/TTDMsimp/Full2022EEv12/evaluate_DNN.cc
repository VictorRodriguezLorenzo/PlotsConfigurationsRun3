#ifndef EVALUATE_DNN_DARK_HIGGS
#define EVALUATE_DNN_DARK_HIGGS

#include <vector>
#include <iostream>
#include <TMath.h>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "ROOT/RVec.hxx"

#include <Python.h>


using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps;

float evaluate_dnn(
    float dphill, 
    float PuppiMET_pt,  
    float mT2, 
    float pdark, 
    float chel, 
    float dphi_ttbar,
    float dphi_met_llb)
{
    Py_Initialize();

    double result = -1;

    // Import the module
    PyObject* sysPath = PySys_GetObject("path");
    PyList_Append(sysPath, PyUnicode_DecodeFSDefault("/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp/Full2022EEv12"));
    PyObject* pModule = PyImport_ImportModule("EvaluateDNN");

    if (pModule == NULL) {
        printf("ERROR importing module \n");
        exit(-1);
    } 

    if (pModule != NULL) {
        // Retrieve the function
        PyObject* pFunction = PyObject_GetAttrString(pModule, "load_neural_network");
        if (pFunction == NULL) {
            printf("ERROR getting function");
            exit(-1);
        }
	if (pFunction != NULL) {
            // Prepare arguments
	    std::vector<float> input;

	    input.push_back(dphill);
            input.push_back(PuppiMET_pt);
            input.push_back(mT2);
            input.push_back(pdark);
            input.push_back(chel);
            input.push_back(dphi_ttbar);
            input.push_back(dphi_met_llb);

            // Input
            PyObject* pList = PyList_New(7);
            for (int i = 0; i < 7; ++i) {
                PyList_SetItem(pList, i, PyFloat_FromDouble((double)input[i]));
            }
            PyObject* pArgs = PyTuple_Pack(1, pList);
            if (pArgs != NULL) {
		    // Call the function
		    PyObject* pValue = PyObject_CallObject(pFunction, pArgs);
		    if (pValue != NULL) {
			    if (PyList_Check(pValue)) {
				    Py_ssize_t listSize = PyList_Size(pValue);
 
				    for (Py_ssize_t i = 0; i < listSize; i++) {
					    PyObject* listItem = PyList_GetItem(pValue, i);
					    result = PyFloat_AsDouble(listItem);
				    }
			    } else {
				    PyErr_Print();
			    } 
		    } else {
			    PyErr_Print();
		    }

		    Py_DECREF(pArgs);
	    } else {
		    PyErr_Print();
	    }

	    Py_DECREF(pFunction);
	} else {
		PyErr_Print();
	}

	Py_DECREF(pModule);
    } else {
	    PyErr_Print();
    }

//    cout << "Returning result DNN: " << result << endl;
    return (float)result;
    Py_Finalize();
}

#endif
