#ifndef EVALUATE_DNN_TWDM
#define EVALUATE_DNN_TWDM

#include <vector>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <unordered_map>
#include <TMath.h>
#include <math.h>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "ROOT/RVec.hxx"

#include <Python.h>


using namespace std;
using namespace ROOT;
using namespace ROOT::VecOps;


RVecF evaluate_dnn_tWDM(
    const RVecF& input,
    const RVecF& mPhi,
    const std::string& model_type
)
{
    // Each alias calls this evaluator with either the scalar ("s") or
    // pseudoscalar ("ps") model.  Cache one Python callable per model type:
    // a single static pointer would make whichever alias runs first select the
    // model for every subsequent call in the same mkShapesRDF process.
    static std::unordered_map<std::string, PyObject*> functions;

    if (!Py_IsInitialized()) {
        Py_Initialize();
    }

    PyGILState_STATE gstate = PyGILState_Ensure();

    PyObject* pFunctionAll = nullptr;
    const auto cachedFunction = functions.find(model_type);

    if (cachedFunction == functions.end()) {
        PyObject* sysPath = PySys_GetObject("path");
        PyList_Append(sysPath, PyUnicode_DecodeFSDefault("/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2023v12"));

        std::string moduleName = "EvaluateDNN_tWDM_" + model_type;
        PyObject* pModule = PyImport_ImportModule(moduleName.c_str());

        if (pModule == NULL) {
            PyErr_Print();
            throw std::runtime_error("ERROR importing EvaluateDNN module");
        }

        pFunctionAll = PyObject_GetAttrString(pModule, "load_neural_network_all_mPhi");
        Py_DECREF(pModule);

        if (pFunctionAll == NULL || !PyCallable_Check(pFunctionAll)) {
            PyErr_Print();
            Py_XDECREF(pFunctionAll);
            throw std::runtime_error("ERROR getting callable load_neural_network_all_mPhi");
        }

        functions.emplace(model_type, pFunctionAll);
    } else {
        pFunctionAll = cachedFunction->second;
    }

    RVecF result(mPhi.size(), -99.9);

    PyObject* pList = PyList_New(input.size());

    for (size_t i = 0; i < input.size(); ++i) {
        PyList_SetItem(pList, i, PyFloat_FromDouble((double)input[i]));
    }

    PyObject* pMassList = PyList_New(mPhi.size());

    for (size_t i = 0; i < mPhi.size(); ++i) {
        PyList_SetItem(pMassList, i, PyFloat_FromDouble((double)mPhi[i]));
    }

    PyObject* pArgs = PyTuple_Pack(2, pList, pMassList);

    Py_DECREF(pList);
    Py_DECREF(pMassList);

    if (pArgs != NULL) {
        PyObject* pValue = PyObject_CallObject(pFunctionAll, pArgs);

        if (pValue != NULL) {
            if (PyList_Check(pValue)) {
                Py_ssize_t listSize = PyList_Size(pValue);
                Py_ssize_t nValues = std::min<Py_ssize_t>(listSize, static_cast<Py_ssize_t>(result.size()));

                for (Py_ssize_t i = 0; i < nValues; i++) {
                    PyObject* listItem = PyList_GetItem(pValue, i);
                    result[i] = (float)PyFloat_AsDouble(listItem);
                }
            } else {
                PyErr_Print();
            }

            Py_DECREF(pValue);
        } else {
            PyErr_Print();
        }

        Py_DECREF(pArgs);
    } else {
        PyErr_Print();
    }

    PyGILState_Release(gstate);

    return result;
}

#endif
