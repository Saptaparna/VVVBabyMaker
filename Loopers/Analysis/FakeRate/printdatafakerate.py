#!/bin/env python

import os
import sys
import ROOT
from QFramework import TQSampleFolder, TQXSecParser, TQCut, TQAnalysisSampleVisitor, TQSampleInitializer, TQCutflowAnalysisJob, TQCutflowPrinter, TQHistoMakerAnalysisJob
from rooutil import plottery_wrapper as p
from plottery import plottery as ply
from rooutil.syncfiles.pyfiles.errors import E

ROOT.gROOT.SetBatch(True)
def compute_fake_factor(th2):
    for ix in xrange(0, th2.GetNbinsX()+2):
        for iy in xrange(0, th2.GetNbinsY()+2):
            frnom = th2.GetBinContent(ix, iy)
            frerr = th2.GetBinError(ix, iy)
            fr = E(frnom, frerr) 
            if fr.val != 0 and fr.val != 1:
                ff = fr / (E(1., 0.) - fr)
            else:
                ff = E(0., 0.)
            th2.SetBinContent(ix, iy, ff.val)
            th2.SetBinError(ix, iy, ff.err)

def compute_fake_factor_1d(th1):
    for ix in xrange(0, th1.GetNbinsX()+2):
        frnom = th1.GetBinContent(ix)
        frerr = th1.GetBinError(ix)
        fr = E(frnom, frerr) 
        if fr.val != 0 and fr.val != 1:
            ff = fr / (E(1., 0.) - fr)
        else:
            ff = E(0., 0.)
        th1.SetBinContent(ix, ff.val)
        th1.SetBinError(ix, ff.err)

samples_nom = TQSampleFolder.loadSampleFolder("output.root:samples")
samples_up = TQSampleFolder.loadSampleFolder("output_up.root:samples")
samples_dn = TQSampleFolder.loadSampleFolder("output_dn.root:samples")

def ewksf(cut, channel, samples):
    data = samples.getHistogram("/data/{}{}".format(channel, channel) , cut+"/MTOneLepFixed").Clone("Data")
    bgs = []
    bgs.append(samples.getHistogram("/top" , cut+"/MTOneLepFixed").Clone("Top"))
    bgs.append(samples.getHistogram("/Zonelep" , cut+"/MTOneLepFixed").Clone("DY"))
    bgs.append(samples.getHistogram("/Wonelep" , cut+"/MTOneLepFixed").Clone("W"))
    totalbkg = p.get_total_hist(bgs)
    data.Rebin(4)
    totalbkg.Rebin(4)
    data.Divide(totalbkg)
    return data.GetBinContent(3)

def move_in_overflow(th2):
    for ix in xrange(0, th2.GetNbinsX()+2):
        overflowbc = th2.GetBinContent(ix, th2.GetNbinsY()+1)
        overflowbe = th2.GetBinError(ix, th2.GetNbinsY()+1)
        overflowb = E(overflowbc, overflowbe)
        highbinbc = th2.GetBinContent(ix, th2.GetNbinsY())
        highbinbe = th2.GetBinError(ix, th2.GetNbinsY())
        highbinb = E(highbinbc, highbinbe)
        newhighbinb = overflowb + highbinb
        th2.SetBinContent(ix, th2.GetNbinsY(), newhighbinb.val)
        th2.SetBinError(ix, th2.GetNbinsY(), newhighbinb.err)
        th2.SetBinContent(ix, th2.GetNbinsY()+1, newhighbinb.val)
        th2.SetBinError(ix, th2.GetNbinsY()+1, newhighbinb.err)


def get_datadriven_fakerate(channel, samples):
    cut = "OneElLoose" if channel == "e" else "OneMuLoose"
    histname = cut+"/lep_ptcorrcoarse_vs_etacoarse"
    dataname = "/data/{}{}".format(channel, channel)
    elsf = ewksf("OneElTightEWKCR", "e", samples)
    musf = ewksf("OneMuTightEWKCR", "m", samples)

    sf = -elsf if channel == "e" else -musf

    bgs  = [
           samples.getHistogram("/top", histname).Clone("Top"),
           samples.getHistogram("/Zonelep", histname).Clone("DY"),
           samples.getHistogram("/Wonelep", histname).Clone("W"),
           ]
    ddqcd = samples.getHistogram(dataname, histname).Clone("Data")
    totalbkg = p.get_total_hist(bgs)
    move_in_overflow(totalbkg)
    move_in_overflow(ddqcd)
    ddqcd.Add(totalbkg, sf)

    bgstight  = [
                samples.getHistogram("/top", str(histname).replace("Loose","Tight")).Clone("Top"),
                samples.getHistogram("/Zonelep", str(histname).replace("Loose","Tight")).Clone("DY"),
                samples.getHistogram("/Wonelep", str(histname).replace("Loose","Tight")).Clone("W"),
                ]
    ddqcdtight = samples.getHistogram(dataname, str(histname).replace("Loose","Tight")).Clone("Data")
    totalbkgtight = p.get_total_hist(bgstight)
    move_in_overflow(totalbkgtight)
    move_in_overflow(ddqcdtight)
    ddqcdtight.Add(totalbkgtight, sf)

    ddqcdtight.Divide(ddqcd)
    compute_fake_factor(ddqcdtight)
    return ddqcdtight

def fullsyst_datadriven_fakerate(channel):
    fr_nom = get_datadriven_fakerate(channel, samples_nom)
    fr_up = get_datadriven_fakerate(channel, samples_up)
    fr_dn = get_datadriven_fakerate(channel, samples_dn)

f = ROOT.TFile("data_fakerates.root", "recreate")
fr_el_nom = get_datadriven_fakerate("e", samples_nom)
fr_mu_nom = get_datadriven_fakerate("m", samples_nom)
fr_el_up = get_datadriven_fakerate("e", samples_up)
fr_mu_up = get_datadriven_fakerate("m", samples_up)
fr_el_dn = get_datadriven_fakerate("e", samples_dn)
fr_mu_dn = get_datadriven_fakerate("m", samples_dn)
fr_el_nom.SetName("fr_el_nom")
fr_mu_nom.SetName("fr_mu_nom")
fr_el_up.SetName("fr_el_up")
fr_mu_up.SetName("fr_mu_up")
fr_el_dn.SetName("fr_el_dn")
fr_mu_dn.SetName("fr_mu_dn")
fr_el_nom.SetTitle("fr_el_nom")
fr_mu_nom.SetTitle("fr_mu_nom")
fr_el_up.SetTitle("fr_el_up")
fr_mu_up.SetTitle("fr_mu_up")
fr_el_dn.SetTitle("fr_el_dn")
fr_mu_dn.SetTitle("fr_mu_dn")
fr_el_nom.Write()
fr_mu_nom.Write()
fr_el_up.Write()
fr_mu_up.Write()
fr_el_dn.Write()
fr_mu_dn.Write()
