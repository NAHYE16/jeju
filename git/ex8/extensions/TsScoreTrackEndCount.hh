#ifndef TsScoreTrackEndCount_hh
#define TsScoreTrackEndCount_hh

#include "TsVBinnedScorer.hh"

class TsScoreTrackEndCount : public TsVBinnedScorer
{
public:
	TsScoreTrackEndCount(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
							 G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);

	virtual ~TsScoreTrackEndCount();

	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
};

#endif
