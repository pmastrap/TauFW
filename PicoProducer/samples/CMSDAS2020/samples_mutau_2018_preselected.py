from TauFW.PicoProducer.storage.Sample import MC as M
from TauFW.PicoProducer.storage.Sample import Data as D
storage  = '/eos/cms/store/group/phys_tau/CMSDAS2020/nano/2018/$DAS'
url      = 'root://eoscms.cern.ch/'
filelist = None #"samples/files/2016/$SAMPLE.txt"
samples  = [
 
  # SINGLE MUON
  D('Data','SingleMuon_Run2018A',"/SingleMuon/Run2018A-Nano25Oct2019-v1/NANOAOD",
   store=storage,url=url,file=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2018B',"/SingleMuon/Run2018B-Nano25Oct2019-v1/NANOAOD",
    store=storage,url=url,file=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2018C',"/SingleMuon/Run2018C-Nano25Oct2019-v1/NANOAOD",
   store=storage,url=url,file=filelist,channels=["skim*",'mutau','mumu','emu'],
  ),
  D('Data','SingleMuon_Run2018D',"/SingleMuon/Run2018D-Nano25Oct2019-v1/NANOAOD",
   store=storage,url=url,file=filelist,channels=["skim*",'mutau','mumu','emu'],
  )
]
