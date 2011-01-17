import FWCore.ParameterSet.Config as cms

#############################################
#############################################
def addVariables(module):

    Variables = cms.PSet(
        mass = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        pt = cms.vstring("Probe p_{T}", "0", "100", "GeV/c"),
        abseta = cms.vstring("Probe |#eta|", "-2.5", "2.5", "")
        )

    module.Variables = Variables
    

#############################################
#############################################
def addEfficiencies(module, absetaBinEdges, ptBinEdges, binToPDFmap):

    Efficiencies = cms.PSet(

        pt_tauAntiEMVA = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("tauAntiEMVA","passed"),
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = cms.PSet(
        #mcTrue = cms.vstring("true"),
        pt = ptBinEdges
        ),
        BinToPDFmap = binToPDFmap
        ),

        abseta_tauAntiEMVA = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("tauAntiEMVA","passed"),
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = cms.PSet(
        #mcTrue = cms.vstring("true"),
        abseta = absetaBinEdges
        ),
        BinToPDFmap = binToPDFmap
        ),

        #pt_tauAntiECrackRem = cms.PSet(
        #EfficiencyCategoryAndState = cms.vstring("tauAntiECrackRem","passed"),
        #UnbinnedVariables = cms.vstring("mass"),
        #BinnedVariables = cms.PSet(
        #pt = ptBinEdges
        #),
        #BinToPDFmap = binToPDFmap
        #),

        #abseta_tauAntiECrackRem = cms.PSet(
        #EfficiencyCategoryAndState = cms.vstring("tauAntiECrackRem","passed"),
        #UnbinnedVariables = cms.vstring("mass"),
        #BinnedVariables = cms.PSet(
        #abseta = absetaBinEdges
        #),
        #BinToPDFmap = binToPDFmap
        #),
        
        pt_tauAntiE2D = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("tauAntiE2D","passed"),
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = cms.PSet(
        #mcTrue = cms.vstring("true"),
        pt = ptBinEdges
        ),
        BinToPDFmap = binToPDFmap
        ),

        abseta_tauAntiE2D = cms.PSet(
        EfficiencyCategoryAndState = cms.vstring("tauAntiE2D","passed"),
        UnbinnedVariables = cms.vstring("mass"),
        BinnedVariables = cms.PSet(
        #mcTrue = cms.vstring("true"),
        abseta = absetaBinEdges
        ),
        BinToPDFmap = binToPDFmap
        ),

        #pt_tauAntiEMVAandCrackRem = cms.PSet(
        #EfficiencyCategoryAndState = cms.vstring("tauAntiEMVAandCrackRem","passed"),
        #UnbinnedVariables = cms.vstring("mass"),
        #BinnedVariables = cms.PSet(
        #pt = ptBinEdges
        #),
        #BinToPDFmap = binToPDFmap
        #),

        #abseta_tauAntiEMVAandCrackRem = cms.PSet(
        #EfficiencyCategoryAndState = cms.vstring("tauAntiEMVAandCrackRem","passed"),
        #UnbinnedVariables = cms.vstring("mass"),
        #BinnedVariables = cms.PSet(
        #abseta = absetaBinEdges
        #),
        #BinToPDFmap = binToPDFmap
        #),

        #pt_tauAntiE2DAndCrackRem = cms.PSet(
        #EfficiencyCategoryAndState = cms.vstring("tauAntiE2DandCrackRem","passed"),
        #UnbinnedVariables = cms.vstring("mass"),
        #BinnedVariables = cms.PSet(
        #pt = ptBinEdges
        #),
        #BinToPDFmap = binToPDFmap
        #),

        #abseta_tauAntiE2DAndCrackRem = cms.PSet(
        #EfficiencyCategoryAndState = cms.vstring("tauAntiE2DandCrackRem","passed"),
        #UnbinnedVariables = cms.vstring("mass"),
        #BinnedVariables = cms.PSet(
        #abseta = absetaBinEdges
        #),
        #BinToPDFmap = binToPDFmap
        #)

        )

    module.Efficiencies = Efficiencies

#############################################
#############################################


#############################################
#############################################
def addCategories(module):
    
    Categories = cms.PSet(
        #mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
        tauAntiEMVA = cms.vstring("tau anti-E mva", "dummy[passed=1,fail=0]"),
        tauAntiE2D = cms.vstring("tau anti-E 2D", "dummy[passed=1,fail=0]"),
        #tauAntiECrackRem = cms.vstring("tau anti-E crack removal", "dummy[passed=1,fail=0]"),
        #tauAntiEMVAandCrackRem = cms.vstring("tau anti-E mva AND crack rem", "dummy[passed=1,fail=0]"),
        #tauAntiE2DandCrackRem = cms.vstring("tau anti-E mva AND crack removal", "dummy[passed=1,fail=0]"),
        )

    module.Categories = Categories

#############################################
#############################################

#############################################
#############################################
    
def addPDFs(module):

    PDFs = cms.PSet(
        gaussPlusLinear = cms.vstring(
        "Gaussian::signalFail(mass, meanFail[91.2, 89.0, 93.0], sigmaFail[2.3, 0.5, 10.0])",
        "Gaussian::signal(mass, mean[85.2, 80.0, 93.0], sigma[5, 0.5, 10.0])",
        "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
        "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
        "efficiency[0.1,0,1]",
        "signalFractionInPassing[0.1]"
        ),
        CBPlusLinear = cms.vstring(
        "RooCBShape::signal(mass,m1[87,85,95],sigma[5.0,3.0,8.5],alfa[1.8,0.,10.],n[1,1.,10.])",
        "RooCBShape::signalFail(mass,m1Fail[91.2,85,95],sigmaFail[2.3,0.5,10.],alfaFail[1.8,0.,10.],nFail[1,0.,10.])",
        "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
        "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
        "efficiency[0.1,0,1]",
        "signalFractionInPassing[0.1]"
        ),
        CBGaussPlusLinear = cms.vstring(
        "Gaussian::signal(mass,m1[87,85,90],sigma[5.0,3.0,8.5])",
        "RooCBShape::signalFail(mass,m1Fail[91.2,85,95],sigmaFail[2.3,0.5,10.],alfaFail[1.8,0.,10.],nFail[1,0.,10.])",
        "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
        "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
        "efficiency[0.1,0,1]",
        "signalFractionInPassing[0.1]"
        ),
        CBPlusExpLinear = cms.vstring(
        "RooCBShape::signal(mass,m1[87,85,90],sigma[5.0,3.0,8.5],alfa[1.0,0.,15.],n[1,0.,15.])",
        "RooCBShape::signalFail(mass,m1Fail[91.2,85,95],sigmaFail[2.3,0.5,10.],alfaFail[1.8,0.,10.],nFail[1,0.,10.])",
        "SUM::backgroundPass(fract[0.5,0,1]*RooExponential::b1(mass, cPass1[0,-10,10]), RooChebychev::b2(mass, cPass2[0,-10,10]))",
        #"RooExponential::backgroundPass(mass, cPass[0,-10,10])",
        "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
        "efficiency[0.1,0,1]",
        "signalFractionInPassing[0.1]"
        ),
        VoigtPlusBifurcPlusExpLinear = cms.vstring(
        "mean[90,80,100]",
        "meanFail[90,80,100]",
        "SUM::signal(fVB[0.5,0,1]*RooBifurGauss::bif(mass,mean,sigmaL[10,0.5,40],sigmaR[1,0.,10]), RooVoigtian::voig(mass, mean, width[2.9], sigmaVoig[5,0.5,10])  )",
        "SUM::signalFail(fVBFail[0.5,0,1]*RooBifurGauss::bifFail(mass,meanFail,sigmaLFail[10,0.5,40],sigmaRFail[1,0.,10]), RooVoigtian::voigFail(mass, meanFail, widthFail[2.9], sigmaVoigFail[5,0.5,10])  )",
        "RooExponential::backgroundPass(mass, cPass[0,-10,10])",
        "RooExponential::backgroundFail(mass, cFail[0,-10,10])",
        "efficiency[0.1,0,1]",
        "signalFractionInPassing[0.1]"
        ),
        
    )

    module.PDFs = PDFs

#############################################
#############################################
    

def cloneNewFitter(process, tnp="",suffix=""):
    setattr(process,"etoTau"+tnp,getattr(process,"etoTau").clone())
    getattr(process,"etoTau"+tnp).InputDirectoryName = "etoTau"+tnp
    getattr(process,"etoTau"+tnp).OutputFileName = "etoTau"+tnp+suffix+".root"
    

#############################################
#############################################
