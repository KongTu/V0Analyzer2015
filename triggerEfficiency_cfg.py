import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.ana_pp = cms.EDAnalyzer('TriggerEfficiency',
              
                        vertexSrc = cms.string('offlinePrimaryVertices'),
                trackSrc = cms.InputTag('generalTracks'),             
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

process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')

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

	

)

)

process.TFileService = cms.Service("TFileService",fileName = cms.string("TriggerEfficiency_test.root"))
process.p = cms.Path(  
		  					#process.hfCoincFilter*
                  			process.ana_pp
                  					)
