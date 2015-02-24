import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltHM85 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM85.HLTPaths = ['HLT_PixelTracks_Multiplicity70']
process.hltHM85.andOr = cms.bool(True)
process.hltHM85.throw = cms.bool(False)

process.hltHM100 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM100.HLTPaths = ['HLT_PixelTracks_Multiplicity70','HLT_PixelTracks_Multiplicity85']
process.hltHM100.andOr = cms.bool(True)
process.hltHM100.throw = cms.bool(False)

process.hltHM120 = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltHM120.HLTPaths = ['HLT_PixelTracks_Multiplicity70','HLT_PixelTracks_Multiplicity85','HLT_PixelTracks_Multiplicity100']
process.hltHM120.andOr = cms.bool(True)
process.hltHM120.throw = cms.bool(False)

process.hltHM.andOr = cms.bool(True)
process.hltHM.throw = cms.bool(False)

process.ana_pp = cms.EDAnalyzer('V0AnalyzerHisto_pp',
                      
                vertexSrc = cms.string('offlinePrimaryVertices'),
                trackSrc = cms.InputTag('generalTracks'),             
                generalV0_ks = cms.InputTag('generalV0CandidatesNew:Kshort'),
                generalV0_la = cms.InputTag('generalV0CandidatesNew:Lambda'),
        			  generalV0_xi = cms.InputTag('generalV0CandidatesNew:Xi'),
        			  genParticleSrc = cms.InputTag('genParticles'),
                doGenParticle = cms.untracked.bool(True),

        			  multmin = cms.untracked.int32(0),
        			  multmax = cms.untracked.int32(1000000),
        			  mult = cms.untracked.int32(0)
)

### standard includes
process.load("Configuration.StandardSequences.Digi_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")

### conditions
#Global Tag:give access to different data info, needs to be changed for different data.
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#for MC data:
#process.GlobalTag.globaltag = 'STARTHI53_V17::All'

#globaltag for pp data
#process.GlobalTag.globaltag = 'GR_P_V43F::All'


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)
process.options   = cms.untracked.PSet( wantSummary =
cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

	'root://xrootd3.cmsaf.mit.edu/'

)

)

process.TFileService = cms.Service("TFileService",fileName = cms.string("V0reco_pPb_histo_pp.root"))
process.p = cms.Path(  # process.hltHM*
		  
                  			process.ana_pp
                  					)
