#list of sample

samples = ["DYJetsToLL_M-50-madgraphMLM","DYJetsToLL_M-50-madgraphMLM_ext1","DYJetsToLL_M-10to50-madgraphMLM", "DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8", "DY2JetsToLL_M-50-madgraphMLM", "DY3JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8", "DY4JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8", "TTTo2L2Nu", "TTToSemiLeptonic", "TTToHadronic", "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8", "W1JetsToLNu", "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8", "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8", "W4JetsToLNu", "ST_tW_antitop_5f_NoFullyHadronicDecays", "ST_tW_top_5f_NoFullyHadronicDecays", "ST_t-channel_antitop_4f_InclusiveDecays", "ST_t-channel_top_4f_InclusiveDecays", "WW", "WZ", "ZZ"]

path = "/eos/cms/store/group/phys_tau/TauFW/nanoV10/Run2_2018/"

# import required module
import os
import ROOT

for s in samples:

        total = 0.0

        path_dir = path+s+"/"

        for filename in os.listdir(path_dir):

            f = os.path.join(path_dir, filename)
	    # checking if it is a filif os.path.isfile(f):
            if os.path.isfile(f):
                print(f)
                _file0 = ROOT.TFile.Open(f)
                t = _file0.Get("Runs")
                for ientry in t:
                    sumw = ientry.genEventSumw
                    total = total + sumw 

        print("Sample: "+s+" EvtSumw: "+str(total))
        #print("EvtSumw: ", total)
