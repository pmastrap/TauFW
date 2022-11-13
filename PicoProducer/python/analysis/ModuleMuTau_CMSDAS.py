# Author: Izaak Neutelings (May 2020)
# Description: Simple module to pre-select mutau events
from ROOT import TFile, TTree, TH1D
from ROOT import TMath
from ROOT import Math, TLorentzVector
import numpy as np
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object


# Inspired by 'Object' class from NanoAODTools.
# Convenient to do so to be able to add MET as 4-momentum to other physics objects using p4()
class Met(Object):
  def __init__(self,event,prefix,index=None):
    self.eta = 0.0
    self.mass = 0.0
    Object.__init__(self,event,prefix,index)


class ModuleMuTau(Module):
  
  def __init__(self,fname,**kwargs):
    self.outfile = TFile(fname,'RECREATE')
    self.default_float = -999.0
    self.default_int = -999
    self.dtype      = kwargs.get('dtype', 'data')
    self.ismc       = self.dtype=='mc'
    self.isdata     = self.dtype=='data'
    self.tes        = kwargs.get('tes',1.0)  #tau energy scale
  
  def beginJob(self):
    """Prepare output analysis tree and cutflow histogram."""
   
    # CUTFLOW HISTOGRAM
    self.cutflow           = TH1D('cutflow','cutflow',25,0,25)
    self.cut_none          = 0
    self.cut_trig          = 1
    self.cut_muon          = 2
    self.cut_muon_veto     = 3
    self.cut_tau           = 4
    self.cut_electron_veto = 5
    self.cut_pair          = 6
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_none,           "no cut"        )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_trig,           "trigger"       )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_muon,           "muon"          )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_muon_veto,      "muon     veto" )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_tau,            "tau"           )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_electron_veto,  "electron veto" )
    self.cutflow.GetXaxis().SetBinLabel(1+self.cut_pair,           "pair"          )
    
    # TREE
    self.tree        = TTree('tree','tree')
    self.pt_1        = np.zeros(1,dtype='f')
    self.eta_1       = np.zeros(1,dtype='f')
    self.m_1         = np.zeros(1,dtype='f')
    self.e_1         = np.zeros(1,dtype='f')
    self.phi_1       = np.zeros(1,dtype='f')
    self.q_1         = np.zeros(1,dtype='i')
    self.id_1        = np.zeros(1,dtype='?')
    self.iso_1       = np.zeros(1,dtype='f')
    self.genmatch_1  = np.zeros(1,dtype='f')
    self.decayMode_1 = np.zeros(1,dtype='i')
    self.pt_2        = np.zeros(1,dtype='f')
    self.eta_2       = np.zeros(1,dtype='f')
    self.m_2         = np.zeros(1,dtype='f')
    self.e_2         = np.zeros(1,dtype='f')
    self.phi_2       = np.zeros(1,dtype='f')
    self.q_2         = np.zeros(1,dtype='i')
    self.id_2        = np.zeros(1,dtype='i')
    self.anti_e_2    = np.zeros(1,dtype='i')
    self.anti_mu_2   = np.zeros(1,dtype='i')
    self.iso_2       = np.zeros(1,dtype='f')
    self.genmatch_2  = np.zeros(1,dtype='f')
    self.decayMode_2 = np.zeros(1,dtype='i')
    self.m_vis       = np.zeros(1,dtype='f')
    self.genWeight   = np.zeros(1,dtype='f')
    self.isoWeight   = np.zeros(1,dtype='f')
    self.idWeight    = np.zeros(1,dtype='f')

    self.met_pt      = np.zeros(1,dtype='f')
    self.met_phi     = np.zeros(1,dtype='f')
    self.met_sumEt   = np.zeros(1,dtype='f')
    self.puppimet_pt      = np.zeros(1,dtype='f')
    self.puppimet_phi     = np.zeros(1,dtype='f')
    self.puppimet_sumEt   = np.zeros(1,dtype='f')
    self.mT_mu_met        = np.zeros(1,dtype='f')
    self.mT_mu_puppimet   = np.zeros(1,dtype='f')
   
    self.pT_vis_Z               = np.zeros(1,dtype='f')
    self.pT_visPLUSmet_Z        = np.zeros(1,dtype='f')
    self.pT_visPLUSpuppimet_Z   = np.zeros(1,dtype='f')
   
    self.m_visPLUSmet_Z        = np.zeros(1,dtype='f')
    self.m_visPLUSpuppimet_Z   = np.zeros(1,dtype='f')

    self.p_miss_zeta            = np.zeros(1,dtype='f')
    self.p_vis_zeta             = np.zeros(1,dtype='f')
    self.Dzeta                  = np.zeros(1,dtype='f')
    self.p_puppimiss_zeta       = np.zeros(1,dtype='f')
    self.Dzeta_puppi            = np.zeros(1,dtype='f') 

    self.deltaR_mu_tau          = np.zeros(1,dtype='f')

    self.pt_j1       = np.zeros(1,dtype='f')
    self.pt_j2       = np.zeros(1,dtype='f')
    self.eta_j1      = np.zeros(1,dtype='f')
    self.eta_j2      = np.zeros(1,dtype='f')
    self.pt_bj1      = np.zeros(1,dtype='f')
    self.pt_bj2      = np.zeros(1,dtype='f')
    self.eta_bj1     = np.zeros(1,dtype='f')
    self.eta_bj2     = np.zeros(1,dtype='f')

    self.pu_density  = np.zeros(1,dtype='f')
    self.npvs        = np.zeros(1,dtype='i')    

    self.pu_True_density = np.zeros(1,dtype='f')

    self.tree.Branch('pt_1',         self.pt_1,        'pt_1/F'       )
    self.tree.Branch('eta_1',        self.eta_1,       'eta_1/F'      )
    self.tree.Branch('m_1',          self.m_1,         'm_1/F'        )
    self.tree.Branch('e_1',          self.e_1,         'e_1/F'        )
    self.tree.Branch('phi_1',        self.phi_1,       'phi_1/F'      )
    self.tree.Branch('q_1',          self.q_1,         'q_1/I'        )
    self.tree.Branch('id_1',         self.id_1,        'id_1/O'       )
    self.tree.Branch('iso_1',        self.iso_1,       'iso_1/F'      )
    self.tree.Branch('genmatch_1',   self.genmatch_1,  'genmatch_1/F' )
    self.tree.Branch('decayMode_1',  self.decayMode_1, 'decayMode_1/I')
    self.tree.Branch('pt_2',         self.pt_2,  'pt_2/F'             )
    self.tree.Branch('eta_2',        self.eta_2, 'eta_2/F'            )
    self.tree.Branch('m_2',          self.m_2,         'm_2/F'        )
    self.tree.Branch('e_2',          self.e_2,         'e_2/F'        )
    self.tree.Branch('phi_2',        self.phi_2,       'phi_2/F'      )
    self.tree.Branch('q_2',          self.q_2,   'q_2/I'              )
    self.tree.Branch('id_2',         self.id_2,  'id_2/I'             )
    self.tree.Branch('anti_e_2',     self.anti_e_2,   'anti_e_2/I'    )
    self.tree.Branch('anti_mu_2',    self.anti_mu_2,  'anti_mu_2/I'   )
    self.tree.Branch('iso_2',        self.iso_2, 'iso_2/F'            )
    self.tree.Branch('genmatch_2',   self.genmatch_2,  'genmatch_2/F' )
    self.tree.Branch('decayMode_2',  self.decayMode_2, 'decayMode_2/I')
    self.tree.Branch('m_vis',        self.m_vis, 'm_vis/F'            )
    self.tree.Branch('genWeight',    self.genWeight,   'genWeight/F'  )

    self.tree.Branch('idWeight',     self.idWeight,    'idWeight/F'  )
    self.tree.Branch('isoWeight',    self.isoWeight,   'isoWeight/F'  )

    self.tree.Branch('met_pt',        self.met_pt,        'met_pt/F'       )
    self.tree.Branch('met_phi',       self.met_phi,       'met_phi/F'      )
    self.tree.Branch('met_sumEt',     self.met_sumEt,     'met_sumEt/F'    )
  
    self.tree.Branch('puppimet_pt',   self.puppimet_pt,   'puppimet_pt/F'       )
    self.tree.Branch('puppimet_phi',  self.puppimet_phi,  'puppimet_phi/F'      )
    self.tree.Branch('puppimet_sumEt',self.puppimet_sumEt,'puppimet_sumEt/F'    )
    self.tree.Branch('mT_mu_met',     self.mT_mu_met,     'mT_mu_met/F'         )
    self.tree.Branch('mT_mu_puppimet',self.mT_mu_puppimet, 'mT_mu_puppimet/F'    )

    self.tree.Branch('pT_vis_Z',            self.pT_vis_Z,             'pT_vis_Z/F'            )
    self.tree.Branch('pT_visPLUSmet_Z ',    self.pT_visPLUSmet_Z ,     'pT_visPLUSmet_Z/F'    )
    self.tree.Branch('pT_visPLUSpuppimet_Z',self.pT_visPLUSpuppimet_Z, 'pT_visPLUSpuppimet_Z/F' )

    self.tree.Branch('m_visPLUSmet_Z ',    self.m_visPLUSmet_Z ,     'm_visPLUSmet_Z/F'    )
    self.tree.Branch('m_visPLUSpuppimet_Z',self.m_visPLUSpuppimet_Z, 'm_visPLUSpuppimet_Z/F' )

    self.tree.Branch('p_miss_zeta',        self.p_miss_zeta,        'p_miss_zeta/F'      )    
    self.tree.Branch('p_vis_zeta',         self.p_vis_zeta,         'p_vis_zeta/F'       )   
    self.tree.Branch('Dzeta',              self.Dzeta,              'Dzeta/F'            ) 
    self.tree.Branch('p_puppimiss_zeta',   self.p_puppimiss_zeta,   'p_puppimiss_zeta/F' )
    self.tree.Branch('Dzeta_puppi',        self.Dzeta_puppi,        'Dzeta_puppi/F'      )    

    self.tree.Branch('deltaR_mu_tau',      self.deltaR_mu_tau,      'deltaR_mu_tau/F'    )

    self.tree.Branch('pt_j1',         self.pt_j1,        'pt_j1/F'       )
    self.tree.Branch('eta_j1',        self.eta_j1,       'eta_j1/F'      )
    self.tree.Branch('pt_j2',         self.pt_j2,        'pt_j2/F'       )
    self.tree.Branch('eta_j2',        self.eta_j2,       'eta_j2/F'      )

    self.tree.Branch('pt_bj1',         self.pt_bj1,        'pt_bj1/F'       )
    self.tree.Branch('eta_bj1',        self.eta_bj1,       'eta_bj1/F'      )
    self.tree.Branch('pt_bj2',         self.pt_bj2,        'pt_bj2/F'       )
    self.tree.Branch('eta_bj2',        self.eta_bj2,       'eta_bj2/F'      )

    self.tree.Branch('pu_density',     self.pu_density,    'pu_density/F'      )
    self.tree.Branch('npvs',     self.npvs,    'npvs/I'      )
    self.tree.Branch('pu_True_density', self.pu_True_density, 'pu_True_density/F')

  def endJob(self):
    """Wrap up after running on all events and files"""
    self.outfile.Write()
    self.outfile.Close()
  
  def analyze(self, event):
    """Process event, return True (pass, go to next module) or False (fail, go to next event)."""
 
    ###################################################
    ## Functions to take SFs from root files #########
    ###################################################

    def get_SF_ID(wp_name,pt,eta):
	    infile=TFile("/eos/user/p/pmastrap/TauLongCMSDAS2020/CMSSW_10_6_20_patch1/src/TauFW/PicoProducer/data/lepton/MuonPOG/Run2018/RunABCD_SF_ID.root","read")
	    #es: NUM_MediumID_DEN_genTracks_pt_abseta
	    hin = infile.Get(wp_name)
	    xaxis = hin.GetXaxis()
	    yaxis = hin.GetYaxis()
	    binx = xaxis.FindBin(pt)
	    biny = yaxis.FindBin(abs(eta))
	    return hin.GetBinContent(binx,biny)

    def get_SF_ISO(wp_name,pt,eta):
            infile=TFile("/eos/user/p/pmastrap/TauLongCMSDAS2020/CMSSW_10_6_20_patch1/src/TauFW/PicoProducer/data/lepton/MuonPOG/Run2018/RunABCD_SF_ISO.root","read")
            #es: NUM_LooseRelIso_DEN_mediumID_pt_abseta
            hin = infile.Get(wp_name)
            xaxis = hin.GetXaxis()
            yaxis = hin.GetYaxis()
            binx = xaxis.FindBin(pt)
            biny = yaxis.FindBin(abs(eta))
            return hin.GetBinContent(binx,biny)
 
    # NO CUT
    self.cutflow.Fill(self.cut_none)
    
    # TRIGGER
    if not event.HLT_IsoMu27: return False
    self.cutflow.Fill(self.cut_trig)
    
    # SELECT MUON
    muons = [ ]
    veto_muons = [ ]

    for muon in Collection(event,'Muon'):
      good_muon = muon.mediumId and muon.pfRelIso04_all < 0.5 and abs(muon.eta) < 2.5
      signal_muon = good_muon and muon.pt > 28.0
      veto_muon   = good_muon and muon.pt > 15.0  
      if signal_muon:
        muons.append(muon)
      if veto_muon: # CAUTION: that's NOT an elif here and intended in that way!
        veto_muons.append(muon)
     
    if len(muons) == 0: return False
    self.cutflow.Fill(self.cut_muon)
    
    if len(veto_muons) > 1: return False
    self.cutflow.Fill(self.cut_muon_veto)

   

    # SELECT TAU
    # TODO section 6: Which decay modes of a tau should be considered for an analysis? Extend tau selection accordingly
    taus = [ ]
    for tau in Collection(event,'Tau'):

      #Tau Energy Scale
      if self.ismc and tau.genPartFlav == 5:
         tau_p4 = TLorentzVector() 
         tau_p4.SetPxPyPzE(tau.p4().Px()*self.tes , tau.p4().Py()*self.tes, tau.p4().Pz()*self.tes, tau.p4().E()*self.tes)
         #print("Old vector")
         #print("Energy: ",tau.p4().E())
         #print("Mass: ",  tau.p4().M())
         #print("Eta: ",   tau.p4().Eta())
         #print("Phi: ",   tau.p4().Phi())
         tau.pt = tau_p4.Pt()
         tau.eta = tau_p4.Eta()
         tau.phi = tau_p4.Phi()
         tau.mass = tau_p4.M()
         #print("New vector")
         #print("Energy: ",tau.p4().E())
         #print("Mass: ",  tau.p4().M())
         #print("Eta: ",   tau.p4().Eta())
         #print("Phi: ",   tau.p4().Phi())


      good_tau = tau.pt > 18.0 and tau.idDeepTau2017v2p1VSe >= 1 and tau.idDeepTau2017v2p1VSmu >= 1 and tau.idDeepTau2017v2p1VSjet >= 1
      if good_tau:
        taus.append(tau)
    if len(taus)<1: return False
    self.cutflow.Fill(self.cut_tau)
    
    print("Tau energy scale: ", self.tes)

    # SELECT ELECTRON
    electrons = []
    for electron in Collection(event,'Electron'):
      good_electron = electron.mvaFall17V2noIso_WPL and electron.pfRelIso03_all < 0.5
      veto_electron = good_electron and electron.pt > 15.0 
      if veto_electron:
         electrons.append(electron)
    if len(electrons) > 0: return False
    self.cutflow.Fill(self.cut_electron_veto)
    
    muon = max(muons,key=lambda p: p.pt)
    tau  = max(taus,key=lambda p: p.pt)
    
    #############################################################
    #pair building algorithm 
    '''
    print("---------------------------------------------")
    print("In this event there are " +str(len(muons))+" muons")
    print("In this event there are " +str(len(taus))+" taus")
    print("Taus pt: ")
    print([t.pt for t in taus] )
    print("Taus iso: ")
    print([t.rawDeepTau2017v2p1VSjet for t in taus] )
    print("Muons pt: ")
    print([m.pt for m in muons] )
    print("Muons iso: ")
    print([m.pfRelIso04_all for m in muons] )

    '''

    mu_tau_pair = []
    indices_pair = []
    for i, mu in enumerate(muons):
        #print("Muon: ", i)
        for k, ta in enumerate(taus):
            #print("Tau: ", k)
            if [k,i] in indices_pair:
               print("done already")
            else:
               indices_pair.append([i,k])
               if mu.DeltaR(ta) > 0.4:
                  mu_tau_pair.append([mu,ta])

    if len(mu_tau_pair) < 1 : return False
    self.cutflow.Fill(self.cut_pair)

    if len(mu_tau_pair) > 1:    
	    mu_tau_pair.sort(key=lambda p: p[0].pfRelIso04_all)
	 
	    if mu_tau_pair[0][0].pfRelIso04_all < mu_tau_pair[1][0].pfRelIso04_all:

		muon = mu_tau_pair[0][0]
		tau  = mu_tau_pair[0][1]

	    elif mu_tau_pair[0][0].pfRelIso04_all == mu_tau_pair[1][0].pfRelIso04_all:

		mu_tau_pair.sort(key=lambda p: p[0].pt, reverse = True)

		if mu_tau_pair[0][0].pt >mu_tau_pair[1][0].pt or len(taus) == 1:

		   muon = mu_tau_pair[0][0]
		   tau  = mu_tau_pair[0][1]

	 
		elif mu_tau_pair[0][0].pt == mu_tau_pair[1][0].pt:
                   #print("Tau comparison")
	 
		   mu_tau_pair.sort(key=lambda p: p[1].rawDeepTau2017v2p1VSjet,reverse = True)
		   
		   if mu_tau_pair[0][1].rawDeepTau2017v2p1VSjet > mu_tau_pair[1][1].rawDeepTau2017v2p1VSjet:
		      
		      muon = mu_tau_pair[0][0]
		      tau  = mu_tau_pair[0][1]

		   elif mu_tau_pair[0][1].rawDeepTau2017v2p1VSjet == mu_tau_pair[1][1].rawDeepTau2017v2p1VSjet:
		      
		      mu_tau_pair.sort(key=lambda p: p[1].pt, reverse = True)
		      #if pair_sort_byPt[0][0].pt > pair_sort_byPt[1][0].pt:
		      muon = mu_tau_pair[0][0]
		      tau  = mu_tau_pair[0][1]

    if len(mu_tau_pair) == 1:
         muon = mu_tau_pair[0][0]
         tau  = mu_tau_pair[0][1]

   # print("chosen Muon Pt : ",muon.pt)
   # print("chosen Muon Iso: ",muon.pfRelIso04_all)
   # print("chosen Tau Pt : ", tau.pt)
   # print("chosen Tau Iso: ", tau.rawDeepTau2017v2p1VSjet)
   # for j, pair in enumerate(mu_tau_pair):
   #     print("Pair: ",j)
   #     print("Muon pt: ", pair[0].pt)
   #     print("Tau pt: ", pair[1].pt)    

    ############################################################################################
    
    # SELECT JET
    hard_jets = []
    btag_jets  = []

    for jet in Collection(event,'Jet'):
        basic_jet = jet.pt > 20.0 and abs(jet.eta) < 4.7 and (jet.puId & (1<<0)) and (jet.jetId & (1<<2)) and jet.DeltaR(muon) > 0.5 and jet.DeltaR(tau) > 0.5 
        hard_jet  = basic_jet and jet.pt > 30
        btagged_jet =  basic_jet and abs(jet.eta) < 2.5 and jet.btagDeepFlavB > 0.2770
        if hard_jet:
           hard_jets.append(jet)
        if btagged_jet:
           btag_jets.append(jet)

    hard_jets.sort(key=lambda p: p.pt, reverse = True)
    btag_jets.sort(key=lambda p: p.pt, reverse = True)


    #MET
    puppimet = Met(event, 'PuppiMET')
    met = Met(event, 'MET')

    #print(met.sumEt)
    #print(met.p4().E())

    #change met according to es
    if self.ismc and tau.genPartFlav == 5:

       met.sumEt = met.sumEt*self.tes
       puppimet.sumEt = puppimet.sumEt*self.tes

       pT_met_vec = Math.Polar2DVectorF()
       pT_met_vec.SetR(met.pt)
       pT_met_vec.SetPhi(met.phi)
       pT_tau_vec = Math.Polar2DVectorF()  
       pT_tau_vec.SetR(abs(1.0-self.tes)*tau.pt)
       pT_tau_vec.SetPhi(tau.phi)
       if self.tes >= 1.00:
          pT_met_tes = pT_met_vec - pT_tau_vec

          met.sumEt = met.sumEt+abs(1.0-self.tes)*tau.p4().Et()
          puppimet.sumEt = puppimet.sumEt+abs(1.0-self.tes)*tau.p4().Et()    

       elif self.tes < 1.00:
          pT_met_tes = pT_met_vec + pT_tau_vec

          met.sumEt = met.sumEt-abs(1.0-self.tes)*tau.p4().Et()
          puppimet.sumEt = puppimet.sumEt-abs(1.0-self.tes)*tau.p4().Et()

       met.pt = pT_met_tes.R()
  
       pT_met_vec.SetR(puppimet.pt)
       pT_met_vec.SetPhi(puppimet.phi)
       if self.tes >= 1.00:
          pT_met_tes = pT_met_vec - pT_tau_vec
       elif self.tes < 1.00:
          pT_met_tes = pT_met_vec + pT_tau_vec
       puppimet.pt = pT_met_tes.R()


    #PILEUP
    pileup = Object(event,'Pileup')

    #PV 
    pv = Object(event, 'PV')

    #TOPOLOGICAL DISCRIMINANT Dzeta
    pT_miss_vec = Math.Polar2DVectorF()
    pT_miss_vec.SetR(met.pt)
    pT_miss_vec.SetPhi(met.phi)

    pT_puppimiss_vec = Math.Polar2DVectorF()
    pT_puppimiss_vec.SetR(puppimet.pt)
    pT_puppimiss_vec.SetPhi(puppimet.phi)

    pT_vis_tau  = Math.Polar2DVectorF()
    pT_vis_tau.SetR(tau.pt)
    pT_vis_tau.SetPhi(tau.phi)
    pT_vis_mu   = Math.Polar2DVectorF()
    pT_vis_mu.SetR(tau.pt)
    pT_vis_mu.SetPhi(tau.phi)
    pT_vis_vec = pT_vis_tau + pT_vis_mu
    
    zeta = (pT_vis_mu.Unit() + pT_vis_tau.Unit()).Unit()

    #check
    #print("pT_vis_mu R: ", pT_vis_mu.R())
    #print("pT_vis_mu Phi: ", pT_vis_mu.Phi())
    #print("X, Y calculated: ", pT_vis_mu.R() * TMath.Cos(pT_vis_mu.Phi()), pT_vis_mu.R() * TMath.Sin(pT_vis_mu.Phi()))
    #print("X, Y function: ", pT_vis_mu.X(), pT_vis_mu.Y() )

    p_miss_zeta = (pT_miss_vec.X()*zeta.X()) + (pT_miss_vec.Y()*zeta.Y()) 
    p_puppimiss_zeta = (pT_puppimiss_vec.X()*zeta.X()) + (pT_puppimiss_vec.Y()*zeta.Y())
    
    p_vis_zeta  = (pT_vis_vec.X()*zeta.X())  + (pT_vis_vec.Y()*zeta.Y())
    
    Dzeta       = p_miss_zeta -0.85 * p_vis_zeta
    Dzeta_puppi = p_puppimiss_zeta - 0.85 * p_vis_zeta

    #########################################################################################
    ##Put variables in the TTree ############################################################
    #########################################################################################

    self.pt_1[0]        = muon.pt
    self.eta_1[0]       = muon.eta
    self.m_1[0]         = muon.mass
    self.e_1[0]         = muon.p4().E()
    self.phi_1[0]       = muon.phi
    self.q_1[0]         = muon.charge
    self.id_1[0]        = muon.mediumId
    self.iso_1[0]       = muon.pfRelIso04_all # keep in mind: the SMALLER the value, the more the muon is isolated
    self.decayMode_1[0] = self.default_int # not needed for a muon
    self.pt_2[0]        = tau.pt
    self.eta_2[0]       = tau.eta
    self.m_2[0]         = tau.mass
    self.e_2[0]         = tau.p4().E()
    self.phi_2[0]       = tau.phi
    self.q_2[0]         = tau.charge
    self.id_2[0]        = tau.idDeepTau2017v2p1VSjet
    self.anti_e_2[0]    = tau.idDeepTau2017v2p1VSe
    self.anti_mu_2[0]   = tau.idDeepTau2017v2p1VSmu
    self.iso_2[0]       = tau.rawDeepTau2017v2p1VSjet # keep in mind: the HIGHER the value of the discriminator, the more the tau is isolated
    self.decayMode_2[0] = tau.decayMode
    self.m_vis[0]       = (muon.p4()+tau.p4()).M()
  
    self.pT_vis_Z[0]             = (muon.p4()+tau.p4()).Pt()
    self.pT_visPLUSmet_Z[0]      = (muon.p4()+tau.p4()+met.p4()).Pt()  
    self.pT_visPLUSpuppimet_Z[0] = (muon.p4()+tau.p4()+puppimet.p4()).Pt()  
    self.m_visPLUSmet_Z[0]       = (muon.p4()+tau.p4()+met.p4()).M()
    self.m_visPLUSpuppimet_Z[0]  = (muon.p4()+tau.p4()+puppimet.p4()).M()
    self.mT_mu_met[0]            = TMath.Sqrt(2* muon.pt * met.sumEt*(1-TMath.Cos(muon.p4().DeltaPhi(met.p4()))))
    self.mT_mu_puppimet[0]       = TMath.Sqrt(2* muon.pt * puppimet.sumEt*(1-TMath.Cos(muon.p4().DeltaPhi(puppimet.p4()))))
    
    self.deltaR_mu_tau[0]  = muon.DeltaR(tau)
    #print(muon.DeltaR(tau))
 
    self.met_phi[0]     = met.phi
    self.met_pt[0]      = met.pt
    self.met_sumEt[0]   = met.sumEt

    self.puppimet_phi[0]     = puppimet.phi
    self.puppimet_pt[0]      = puppimet.pt
    self.puppimet_sumEt[0]   = puppimet.sumEt

    self.p_miss_zeta[0]  = p_miss_zeta
    self.p_vis_zeta[0]   = p_vis_zeta
    self.Dzeta[0]        = Dzeta
    self.p_puppimiss_zeta[0]  = p_puppimiss_zeta
    self.Dzeta_puppi[0]        = Dzeta_puppi 

    if len(hard_jets) > 0:
       self.pt_j1[0]        = hard_jets[0].pt
       self.eta_j1[0]       = hard_jets[0].eta
       if len(hard_jets) > 1: 
          self.pt_j2[0]        = hard_jets[1].pt
          self.eta_j2[0]       = hard_jets[1].eta
       else:
          self.pt_j2[0]        = 0.0
          self.eta_j2[0]       = 0.0
    else:
       self.pt_j1[0]        = 0.0
       self.eta_j1[0]       = 0.0
       self.pt_j2[0]        = 0.0         
       self.eta_j2[0]       = 0.0

    if len(btag_jets) > 0:
       self.pt_bj1[0]        = btag_jets[0].pt
       self.eta_bj1[0]       = btag_jets[0].eta
       if len(btag_jets) > 1:
          self.pt_bj2[0]        = btag_jets[1].pt
          self.eta_bj2[0]       = btag_jets[1].eta
       else:
          self.pt_bj2[0]        = 0.0         
          self.eta_bj2[0]       = 0.0
    else:
       self.pt_bj1[0]        = 0.0
       self.eta_bj1[0]       = 0.0
       self.pt_bj2[0]        = 0.0
       self.eta_bj2[0]       = 0.0

    #event variables
   
    self.npvs[0]       = pv.npvs

    if self.ismc:
      self.genmatch_1[0]  = muon.genPartFlav # in case of muons: 1 == prompt muon, 15 == muon from tau decay, also other values available for muons from jets
      self.genmatch_2[0]  = tau.genPartFlav # in case of taus: 0 == unmatched (corresponds then to jet),
                                            #                  1 == prompt electron, 2 == prompt muon, 3 == electron from tau decay,
                                            #                  4 == muon from tau decay, 5 == hadronic tau decay
      self.genWeight[0] = event.genWeight
      self.idWeight[0]  = get_SF_ID("NUM_MediumID_DEN_genTracks_pt_abseta",muon.pt,muon.eta)
      self.isoWeight[0] = get_SF_ISO("NUM_LooseRelIso_DEN_MediumID_pt_abseta", muon.pt,muon.eta)
      self.pu_True_density[0] = pileup.gpudensity
      self.pu_density[0] = pileup.pudensity

    self.tree.Fill()
    
    return True
