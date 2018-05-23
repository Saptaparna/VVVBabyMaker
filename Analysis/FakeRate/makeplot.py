#!/bin/env python

import os
import sys
import ROOT
from QFramework import TQSampleFolder, TQXSecParser, TQCut, TQAnalysisSampleVisitor, TQSampleInitializer, TQCutflowAnalysisJob, TQCutflowPrinter, TQHistoMakerAnalysisJob
from rooutil import plottery_wrapper as p
from plottery import plottery as ply

ROOT.gROOT.SetBatch(True)
samples = TQSampleFolder.loadSampleFolder("output.root:samples")
#samples.printContents("trd")

output_plot_dir = "plots"

doW = False

testsample = "/top"
testsamplename = "tt"
testsamplelegendname = "t#bar{t}"
if doW:
    testsample = "/W"
    testsamplename = "W"
    testsamplelegendname = "W"

#_____________________________________________________________________________________
def plot(histname, output_name, systs=None, options={}, plotfunc=p.plot_hist):
    # Options
    alloptions= {
                "ratio_range":[0.4,1.6],
                #"nbins": 30,
                "autobin": True,
                "legend_scalex": 1.4,
                "legend_scaley": 1.1,
                "output_name": "{}/{}.pdf".format(output_plot_dir, output_name)
                }
    alloptions.update(options)
    sigs = []
    bgs  = [ samples.getHistogram("/qcd/mu", histname).Clone("QCD"),
             ]
    data = None
    colors = [ 2005, 920 ]
    plotfunc(
            sigs = sigs,
            bgs  = bgs,
            data = data,
            colors = colors,
            syst = systs,
            options=alloptions)

#_____________________________________________________________________________________
def fakerate(histname, output_name, systs=None, options={}, plotfunc=p.plot_hist):
    if str(histname).find("Loose") == -1: return
    isqcd = str(histname).find("One") != -1
    ismu = str(histname).find("Mu") != -1
    ispredict = str(histname).find("Predict") != -1

    sample = testsample
    if isqcd and ismu:
        sample = "/qcd/mu"
    elif isqcd and not ismu:
        sample = "/qcd/el/bcToE"
        #sample = "/qcd/el/EM"
        #sample = "/qcd/el"
    samplename = "ttbar estimation" if ispredict else ("QCD Loose" if isqcd else "ttbar Loose")
    color = 920 if isqcd else 2005
    # Options
    alloptions= {
                "ratio_range":[0.0, 2.0] if ispredict else [0.0,0.45],
                #"nbins": 30,
                "autobin": True,
                "legend_scalex": 0.8,
                "legend_scaley": 0.8,
                "legend_datalabel": samplename.replace("estimation", "prediction") if ispredict else samplename.replace("Loose", "Tight"),
                "output_name": "{}/fr_{}.pdf".format(output_plot_dir, output_name),
                "hist_disable_xerrors": True if str(histname).find("varbin") != -1 else False,
                }
    histnum = samples.getHistogram(sample, str(histname).replace("Loose", "Tight")).Clone(samplename)
    histden = samples.getHistogram(sample, histname).Clone(samplename)
    alloptions.update(options)
    sigs = []
    bgs  = [ histden ]
    data = histnum
    colors = [ color ]
    try:
        plotfunc(
                sigs = sigs,
                bgs  = bgs,
                data = data,
                colors = colors,
                syst = systs,
                options=alloptions)
    except:
        print "ERROR: failed plotting {} {}".format(histname, sample)

#_____________________________________________________________________________________
def fakerate2d(histname, output_name, systs=None, options={}, plotfunc=p.plot_hist):
    if str(histname).find("Loose") == -1: return
    isqcd = str(histname).find("One") != -1
    ismu = str(histname).find("Mu") != -1
    ispredict = str(histname).find("Predict") != -1

    sample = testsample
    if isqcd and ismu:
        sample = "/qcd/mu"
    elif isqcd and not ismu:
        #sample = "/qcd/el/bcToE"
        #sample = "/qcd/el/EM"
        sample = "/qcd/el"
    samplename = "ttbar estimation" if ispredict else ("QCD Loose" if isqcd else "ttbar Loose")
    color = 920 if isqcd else 2005
    # Options
    alloptions= {
                "ratio_range":[0.0, 2.0] if ispredict else [0.0,0.3],
                #"nbins": 30,
                "autobin": True,
                "legend_scalex": 0.8,
                "legend_scaley": 0.8,
                "legend_datalabel": samplename.replace("estimation", "prediction") if ispredict else samplename.replace("Loose", "Tight"),
                "output_name": "{}/fr_{}.pdf".format(output_plot_dir, output_name),
                "hist_disable_xerrors": True if str(histname).find("varbin") != -1 else False,
                }
    #samples.getHistogram(sample, str(histname).replace("Loose", "Tight")).Clone(samplename).Print("all")
    histnum = p.flatten_th2(samples.getHistogram(sample, str(histname).replace("Loose", "Tight")).Clone(samplename))
    histden = p.flatten_th2(samples.getHistogram(sample, histname).Clone(samplename))
    alloptions.update(options)
    sigs = []
    bgs  = [ histden ]
    data = histnum
    colors = [ color ]
    plotfunc(
            sigs = sigs,
            bgs  = bgs,
            data = data,
            colors = colors,
            syst = systs,
            options=alloptions)
    histnum.Divide(histden)

#_____________________________________________________________________________________
def fakeratecomp(histname, output_name, systs=None, options={}, plotfunc=p.plot_hist):
    if str(histname).find("Loose") == -1: return
    if str(histname).find("One") == -1: return
    if str(histname).find("varbin") == -1: return
    ismu = str(histname).find("Mu") != -1

    ispredict = str(histname).find("Predict") != -1

    qcdsample = "/qcd/mu" if ismu else "/qcd/el/bcToE"
    ttbarsample = testsample
    qcdsamplename = "QCD" if ismu else "QCD"
    ttbarsamplename = testsamplename

    qcdhistname = histname
    ttbarhistname = histname.replace("One", "Two")
    ttbarhistname = ttbarhistname.replace("lep_", "mu_") if ismu else ttbarhistname.replace("lep_", "el_")

    print qcdhistname, ttbarhistname

    qcdhistnum = p.move_overflow(samples.getHistogram(qcdsample, str(qcdhistname).replace("Loose", "Tight")).Clone(qcdsamplename))
    qcdhistden = p.move_overflow(samples.getHistogram(qcdsample, qcdhistname).Clone(qcdsamplename))
    ttbarhistnum = p.move_overflow(samples.getHistogram(ttbarsample, str(ttbarhistname).replace("Loose", "Tight")).Clone(ttbarsamplename))
    ttbarhistden = p.move_overflow(samples.getHistogram(ttbarsample, ttbarhistname).Clone(ttbarsamplename))

    qcdhistnum.Divide(qcdhistden)
    ttbarhistnum.Divide(ttbarhistden)

    # Options
    alloptions= {
                "ratio_range":[0.0, 2.0] if ispredict else [0.5 if ismu else 0.0, 1.5 if ismu else 3.0],
                "yaxis_range":[0.0,0.25 if ismu else 0.4],
                #"nbins": 30,
                "yaxis_log": False,
                "draw_points": True,
                "autobin": True,
                "legend_scalex": 0.8,
                "legend_scaley": 0.8,
                "legend_datalabel": testsamplelegendname,
                "output_name": "{}/fr_closure_{}.pdf".format(output_plot_dir, output_name)
                }
    #samples.getHistogram(sample, str(histname).replace("Loose", "Tight")).Clone(samplename).Print("all")
    #histnum = p.flatten_th2(samples.getHistogram(sample, str(histname).replace("Loose", "Tight")).Clone(samplename))
    #histden = p.flatten_th2(samples.getHistogram(sample, histname).Clone(samplename))
    alloptions.update(options)
    sigs = []
    bgs  = [ qcdhistnum ]
    data = ttbarhistnum
    colors = [ 2 ]
    plotfunc(
            sigs = sigs,
            bgs  = bgs,
            data = data,
            colors = colors,
            syst = systs,
            options=alloptions)

#_____________________________________________________________________________________
def fakerate2dcomp(histname, output_name, systs=None, options={}, plotfunc=p.plot_hist):
    if str(histname).find("Loose") == -1: return
    if str(histname).find("One") == -1: return
    if str(histname).find("corr") == -1: return
    ismu = str(histname).find("Mu") != -1

    qcdsample = "/qcd/mu" if ismu else "/qcd/el/bcToE"
    ttbarsample = testsample
    qcdsamplename = "QCD" if ismu else "QCD"
    ttbarsamplename = testsamplename

    qcdhistname = histname
    ttbarhistname = histname.replace("One", "Two")
    ttbarhistname = ttbarhistname.replace("lep_", "mu_") if ismu else ttbarhistname.replace("lep_", "el_")

    print qcdhistname, ttbarhistname

    qcdhistnum = p.flatten_th2(samples.getHistogram(qcdsample, str(qcdhistname).replace("Loose", "Tight")).Clone(qcdsamplename))
    qcdhistden = p.flatten_th2(samples.getHistogram(qcdsample, qcdhistname).Clone(qcdsamplename))
    ttbarhistnum = p.flatten_th2(samples.getHistogram(ttbarsample, str(ttbarhistname).replace("Loose", "Tight")).Clone(ttbarsamplename))
    ttbarhistden = p.flatten_th2(samples.getHistogram(ttbarsample, ttbarhistname).Clone(ttbarsamplename))

    qcdhistnum.Divide(qcdhistden)
    ttbarhistnum.Divide(ttbarhistden)

    # Options
    alloptions= {
                "ratio_range":[0.0,2.0 if ismu else 3.0],
                "yaxis_range":[0.0,0.15 if ismu else 0.3],
                #"nbins": 30,
                "draw_points": True,
                "autobin": True,
                "legend_scalex": 0.8,
                "legend_scaley": 0.8,
                "legend_datalabel": testsamplelegendname,
                "output_name": "{}/fr_closure_{}.pdf".format(output_plot_dir, output_name)
                }
    #samples.getHistogram(sample, str(histname).replace("Loose", "Tight")).Clone(samplename).Print("all")
    #histnum = p.flatten_th2(samples.getHistogram(sample, str(histname).replace("Loose", "Tight")).Clone(samplename))
    #histden = p.flatten_th2(samples.getHistogram(sample, histname).Clone(samplename))
    alloptions.update(options)
    sigs = []
    bgs  = [ qcdhistnum ]
    data = ttbarhistnum
    colors = [ 2 ]
    plotfunc(
            sigs = sigs,
            bgs  = bgs,
            data = data,
            colors = colors,
            syst = systs,
            options=alloptions)

if __name__ == "__main__":

    if len(sys.argv) < 2:

        import multiprocessing

        jobs = []
        for histname in samples.getListOfHistogramNames():
            hname = str(histname)
            if hname.find("_vs_") != -1:
                hfilename = hname.replace("/", "_")
                #proc = multiprocessing.Process(target=fakerate2d, args=[hname, hfilename], kwargs={"systs":None, "options":{}, "plotfunc": ply.plot_hist_2d})
                proc = multiprocessing.Process(target=fakerate2d, args=[hname, hfilename], kwargs={"systs":None, "options":{"autobin":False, "nbins":15, "lumi_value":35.9, "yaxis_log":True}, "plotfunc": p.plot_hist})
                jobs.append(proc)
                proc.start()
                proc = multiprocessing.Process(target=fakerate2dcomp, args=[hname, hfilename], kwargs={"systs":None, "options":{"autobin":False, "nbins":15, "lumi_value":35.9, "yaxis_log":False}, "plotfunc": p.plot_hist})
                jobs.append(proc)
                proc.start()
            else:
                hfilename = hname.replace("/", "_")
                proc = multiprocessing.Process(target=fakerate, args=[hname, hfilename], kwargs={"systs":None, "options":{"autobin":False, "nbins":15, "lumi_value":35.9, "yaxis_log":False}, "plotfunc": p.plot_hist})
                jobs.append(proc)
                proc.start()
                proc = multiprocessing.Process(target=fakeratecomp, args=[hname, hfilename], kwargs={"systs":None, "options":{"autobin":False, "nbins":15, "lumi_value":35.9, "yaxis_log":False}, "plotfunc": p.plot_hist})
                jobs.append(proc)
                proc.start()

        for job in jobs:
            job.join()

    else:
        hname = str(sys.argv[1])
        hfilename = hname.replace("/", "_")
        ratio = plot(hname, hfilename, options={"legend_scalex":1.5, "autobin":False, "blind":False, "nbins":15, "signal_scale":7, "yaxis_log":False, "yaxis_range":[float(sys.argv[2]),float(sys.argv[3])]}, plotfunc=p.plot_hist).Clone("nvtx")
        #ratio.Print("all")

        f = ROOT.TFile("nvtx_rewgt_mm.root","recreate")
        ratio.Write()
