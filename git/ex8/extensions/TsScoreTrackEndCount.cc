// Scorer for TrackEndCount

#include "TsScoreTrackEndCount.hh"

TsScoreTrackEndCount::TsScoreTrackEndCount(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
										   G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
: TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	SetUnit("");
}


TsScoreTrackEndCount::~TsScoreTrackEndCount() {;}


G4bool TsScoreTrackEndCount::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}

    /*
	if (aStep->GetPostStepPoint()->GetKineticEnergy() == 0) {
		G4double wgt = aStep->GetPreStepPoint()->GetWeight();
		AccumulateHit(aStep, wgt);
		return true;
	}
    */

    //Is following equivalent?
    //are there any possibilities to set track status to kill while energy is non-zero?
    switch(  aStep->GetTrack()->GetTrackStatus()  ){
    case fAlive:
        return false;
    case fPostponeToNextEvent:
        return false;
    case fStopAndKill: 
    case fStopButAlive:
    case fKillTrackAndSecondaries:
    default:
    	G4double wgt = aStep->GetPreStepPoint()->GetWeight();
		AccumulateHit(aStep, wgt);
        return true;
    }

        
	return false;
}
