#!/bin/bash
set -e

root -l -b -q 'bTagEff.cc("Full2022v12","medium")'
root -l -b -q 'bTagEff.cc("Full2022EEv12","medium")'
root -l -b -q 'bTagEff.cc("Full2023v12","medium")'
root -l -b -q 'bTagEff.cc("Full2023BPixv12","medium")'
root -l -b -q 'bTagEff.cc("Full2024v15","medium")'
