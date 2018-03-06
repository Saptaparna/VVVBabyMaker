#!/bin/env python

import os
import ROOT
from QFramework import TQSampleFolder, TQEventlistPrinter, TQTaggable
from rooutil import plottery_wrapper as p

from cuts import bkg_types_SS
from cuts import bkg_types_3L

ROOT.gROOT.SetBatch(True)

#samples = TQSampleFolder.loadSampleFolder("output.root:samples")
samples = TQSampleFolder.loadSampleFolder("validate.root:samples")
#samples.getSampleFolder("/samples/bkg/mm/ttV").printContents("trd")
printer = TQEventlistPrinter(samples)
printer.addCut("SSee")
printer.addCut("SSem")
printer.addCut("SSmm")
printer.addCut("TL0SFOS")
printer.addCut("TL1SFOS")
printer.addCut("TL2SFOS")
#printer.addProcess("/sig/whwww")
#printer.addProcess("/typebkg")
printer.addProcess("/samples/www")
printer.writeEventlists("evtlist", "eventlists", "verbose=true");
