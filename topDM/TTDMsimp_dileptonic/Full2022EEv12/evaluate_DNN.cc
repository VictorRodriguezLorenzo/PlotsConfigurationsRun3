#ifndef EVALUATE_DNN_TTDM
#define EVALUATE_DNN_TTDM

#include <vector>
#include <iostream>
#include <stdexcept>
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
//    float lep_pt1,
//    float lep_pt2,
//    float lep_eta1,
//    float lep_eta2,
//    float mll,
//    float ptll,
//    float drll,
//    float detall,
    float dphill,
//    float yll,
    float PuppiMET_pt,
//    float PuppiMET_phi,
//    float dphilmet,
//    float dphilmet1,
//    float dphilmet2,
//    float dphillmet,
//    float mtw1,
//    float mtw2,
//    float mth,
//    float mTi,
//    float mR,
    float mT2,
//    float mTe,
//    float recoil,
//    float upara,
//    float uperp,
//    float pTWW,
//    float mcoll,
//    float mcollWW,
//    float choiMass,
//    float nbjet_jet_ratio,
//    int njet,
//    float ht,
//    float vht_pt,
//    float dphijet1met,
//    float dphijet2met,
//    float dphijjmet,
    float chel,
    float pdark,
    float dphi_ttbar,
    float dphi_met_llb,
    int mPhi
                )
{
    static PyObject* pFunction = nullptr;
    if (!Py_IsInitialized()) {
        Py_Initialize();
    }
    PyGILState_STATE gstate = PyGILState_Ensure();
    if (pFunction == nullptr) {
        PyObject* sysPath = PySys_GetObject("path");
        PyList_Append(sysPath, PyUnicode_DecodeFSDefault("/afs/cern.ch/user/v/victorr/private/PlotsConfigurationsRun3/topDM/TTDMsimp_dileptonic/Full2022EEv12"));

        PyObject* pModule = PyImport_ImportModule("EvaluateDNN");
        if (pModule == NULL) {
            PyErr_Print();
            throw std::runtime_error("ERROR importing EvaluateDNN module");
        }

        pFunction = PyObject_GetAttrString(pModule, "load_neural_network");
        Py_DECREF(pModule);
        if (pFunction == NULL) {
            PyErr_Print();
            throw std::runtime_error("ERROR getting load_neural_network");
        }
    }

    double result = -1;
    // Prepare arguments
    std::vector<float> input;

//            input.push_back(lep_pt1);
//            input.push_back(lep_pt2);
//            input.push_back(lep_eta1);
//            input.push_back(lep_eta2);
// 
//            input.push_back(mll);
//            input.push_back(ptll);
//            input.push_back(drll);
//            input.push_back(detall);
            input.push_back(dphill);
//            input.push_back(yll);
 
            input.push_back(PuppiMET_pt);
//            input.push_back(PuppiMET_phi);
//            input.push_back(dphilmet);
//            input.push_back(dphilmet1);
//            input.push_back(dphilmet2);
//            input.push_back(dphillmet);
// 
//            input.push_back(mtw1);
//            input.push_back(mtw2);
//            input.push_back(mth);
//            input.push_back(mTi);
//            input.push_back(mR);
            input.push_back(mT2);
//            input.push_back(mTe);
 
//            input.push_back(recoil);
//            input.push_back(upara);
//            input.push_back(uperp);
//            input.push_back(pTWW);
// 
//            input.push_back(mcoll);
//            input.push_back(mcollWW);
//            input.push_back(choiMass);
// 
//            input.push_back(nbjet_jet_ratio);
//            input.push_back(njet);
//            input.push_back(ht);
//            input.push_back(vht_pt);
//            input.push_back(dphijet1met);
//            input.push_back(dphijet2met);
//            input.push_back(dphijjmet);
 
            input.push_back(chel);                       // doubleNu_producer[6]
            input.push_back(pdark);                      // doubleNu_producer[8]
            input.push_back(dphi_ttbar);                 // doubleNu_producer[7]
            input.push_back(dphi_met_llb);
            input.push_back(mPhi);

    // Input
    PyObject* pList = PyList_New(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        PyList_SetItem(pList, i, PyFloat_FromDouble((double)input[i]));
    }
    PyObject* pArgs = PyTuple_Pack(1, pList);
    Py_DECREF(pList);
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
            Py_DECREF(pValue);
        } else {
            PyErr_Print();
        }
        Py_DECREF(pArgs);
    } else {
        PyErr_Print();
    }

    PyGILState_Release(gstate);

    return (float)result;
}

#endif
