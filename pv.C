#include "CommonDataFormat/RangeReference.h"
#include "CommonDataFormat/TimeStamp.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsCommonDataFormats/NameConf.h"
#include "Field/MagneticField.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/MatchInfoTOF.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "gsl/span"
#include <TChain.h>
#include <numeric>
#include <string>

#include "TH1F.h"
#include "TH2F.h"

#include "Vertexer.h"

namespace o2d = o2::dataformats;

using V2TRef = o2d::RangeReference<int, int>;
using TimeEst = o2d::TimeStampWithError<float, float>;
using TrackVF = Vertexer::Track;
using Track = o2::track::TrackParCov;

struct TimeVal {
  TimeEst timeEst;
  int idStart = -1;
  int idEnd = -1;
  
  void addTrack(int i, const TimeEst& trcT) {
    auto trcTErr2 = trcT.getTimeStampError()*trcT.getTimeStampError();
    auto trcTErr2Inv = 1./trcTErr2;
    if (idStart<0) {
      idStart = idEnd = i;
      timeEst.setTimeStamp(trcT.getTimeStamp());
      timeEst.setTimeStampError(trcTErr2); // we store the error^2
    }
    else {
      auto vtxTErr2Inv = 1./timeEst.getTimeStampError();
      auto vtxTErr2UpdInv = trcTErr2Inv + vtxTErr2Inv;
      auto vtxTErr2Upd = 1./vtxTErr2UpdInv;
      timeEst.setTimeStamp( (timeEst.getTimeStamp()*vtxTErr2Inv +  trcT.getTimeStamp()*trcTErr2Inv)*vtxTErr2Upd );
      timeEst.setTimeStampError( vtxTErr2Upd );
      if (i>idEnd) {
        idEnd = i;
      }
    }
  }

  bool isCompatible(const TimeEst& trcT, double margin2, double cut) {
    if (idStart<0) {
      return true;
    }
    auto dt = timeEst.getTimeStamp() - trcT.getTimeStamp();
    auto trcTErr2 = trcT.getTimeStampError()*trcT.getTimeStampError();
    auto err = trcTErr2 + timeEst.getTimeStampError() + margin2;
    return dt*dt/err<cut;
  }

  void print() const {
    printf("**** %e/%e / %d : %d -> %d\n", timeEst.getTimeStamp(), timeEst.getTimeStampError(), idStart, idEnd, idEnd - idStart + 1);
  }

};

TH1F *ht = nullptr;
TH1F *he = nullptr;
TH1F *htv = nullptr;
TH2F *htz = nullptr;

//------------------------------------------
// will be taken from DPL input
std::vector<o2d::TrackTPCITS> *mITSTPCtracksPtr = nullptr;
std::vector<o2d::MatchInfoTOF> *mTOFtracksPtr = nullptr;
//------------------------------------------

gsl::span<o2d::TrackTPCITS> mITSTPCTracks;
gsl::span<o2d::MatchInfoTOF> mTOFTracks;
o2d::VertexBase mMeanVertex;
std::vector<TrackVF> mTracksPool;
std::vector<int> mSortID;

std::vector<TimeVal> intTimes;

std::vector<Vertex> mVertices;
std::vector<int> mVertexTrackIDs;
std::vector<V2TRef> mV2TRefs;

float mBz = 0;

void createTracksPool();
void clusterizeTime(gsl::span<TrackVF> pool, gsl::span<int> idxsort, double margin2=1., double cut=25);
void finalizeVertex(Vertex vtx, gsl::span<TrackVF> pool, gsl::span<int> idxsort);

void pv(std::string dataDir = "./data",
        const std::string &itstpcTrFName = "o2match_itstpc.root",
        const std::string &tofTrFName = "o2match_tof.root") {
  if (dataDir.empty()) {
    dataDir = "./";
  }
  if (dataDir.back() != '/') {
    dataDir += '/';
  }

  TChain itstpcChain("matchTPCITS");
  itstpcChain.AddFile(o2::utils::concat_string(dataDir, itstpcTrFName).c_str());
  itstpcChain.SetBranchAddress("TPCITS", &mITSTPCtracksPtr);
  itstpcChain.GetEntry(0);
  mITSTPCTracks = gsl::span<o2d::TrackTPCITS>(mITSTPCtracksPtr->data(),
                                              mITSTPCtracksPtr->size());

  //
  TChain tofChain("matchTOF");
  tofChain.AddFile(o2::utils::concat_string(dataDir, tofTrFName).c_str());
  tofChain.SetBranchAddress("TOFMatchInfo", &mTOFtracksPtr);
  tofChain.GetEntry(0);
  mTOFTracks = gsl::span<o2d::MatchInfoTOF>(mTOFtracksPtr->data(),
                                            mTOFtracksPtr->size());
  //
  mMeanVertex.setSigmaX2(0.1 * 01);
  mMeanVertex.setSigmaY2(0.1 * 01);
  mMeanVertex.setSigmaZ2(6. * 6.);

  o2::base::GeometryManager::loadGeometry(
      o2::utils::concat_string(dataDir, o2::base::NameConf::getGeomFileName()));
  o2::base::Propagator::initFieldFromGRP(
      o2::utils::concat_string(dataDir, o2::base::NameConf::getGRPFileName()));
  mBz = o2::base::Propagator::Instance()->getNominalBz();

  createTracksPool();
}

void createTracksPool() {
  // create pull of all candidate tracks in a global array ordered in time
  mTracksPool.clear();
  mSortID.clear();

  mTracksPool.reserve(mITSTPCTracks.size());
  // check all containers
  float vtxErr2 = 0.5 * (mMeanVertex.getSigmaX2() + mMeanVertex.getSigmaY2());
  float dcaToler = 1.0 + 3. * std::sqrt(vtxErr2);
  float pullIniCut = 9.;
  o2d::DCA dca;

  auto ntGlo = mITSTPCTracks.size();
  for (uint32_t i = 0; i < ntGlo; i++) {
    Track trc = mITSTPCTracks[i];
    if (!trc.propagateToDCA(mMeanVertex, mBz, &dca, dcaToler) ||
        dca.getY() * dca.getY() / (dca.getSigmaY2() + vtxErr2) > pullIniCut) {
      continue;
    }
    mTracksPool.emplace_back(trc, mITSTPCTracks[i].getTimeMUS(), i, 0);
  }
  // TODO: try to narrow timestamps using tof times

  if (mTracksPool.empty()) {
    return;
  }

  // estimate mean error
  float erTMean = 0.;
  for (const auto &trc : mTracksPool) {
    erTMean += trc.timeEst.getTimeStampError();
  }
  erTMean /= mTracksPool.size();

  //
  mSortID.resize(mTracksPool.size());
  std::iota(mSortID.begin(), mSortID.end(), 0);

  std::sort(mSortID.begin(), mSortID.end(), [](int i, int j) {
    return mTracksPool[i].timeEst.getTimeStamp() <
           mTracksPool[j].timeEst.getTimeStamp();
  });

  auto tMin = mTracksPool[mSortID.front()].timeEst.getTimeStamp();
  auto tMax = mTracksPool[mSortID.back()].timeEst.getTimeStamp();

  int itmin = int(tMin), itmax = std::ceil(tMax);
  ht = new TH1F("ht", "time", 10*(itmax - itmin + 1), itmin, itmax);
  htv = new TH1F("htv", "timeV", 10*(itmax - itmin + 1), itmin, itmax);
  he = new TH1F("et", "time er", 1000, 0, 2.);
  htz = new TH2F("tz","t z", 3*(itmax - itmin + 1+2), itmin-1, itmax+1, 2000, -20,20);
  for (auto i : mSortID) {
    ht->Fill(mTracksPool[i].timeEst.getTimeStamp());
    he->Fill(mTracksPool[i].timeEst.getTimeStampError());
    htz->Fill(mTracksPool[i].timeEst.getTimeStamp(), mTracksPool[i].z);
  }
}

void finalizeVertex(Vertex vtx, gsl::span<TrackVF> pool, gsl::span<int> idxsort) {
  int lastID = mVertices.size();
  mVertices.emplace_back(vtx);
  auto &ref = mV2TRefs.emplace_back(mVertexTrackIDs.size(), 0);
  for (int i : idxsort) {
    if (pool[i].canAssign()) {
      mVertexTrackIDs.push_back(pool[i].entry);
      pool[i].vtxID = lastID;
      // remove track from ZSeeds histo
      auto bin = getZSeedBin( pool[i].getZForXY(mXYConstraint[0],mXYConstraint[1]) );
      //      mZSeedHisto[bin]--;
      //      mZSeedsFilled--;
    }
  }
  ref.setEntries(vtx.getNContributors());
}

void clusterizeTime(gsl::span<TrackVF> pool, gsl::span<int> idxsort, double margin2, double cut)
{
  intTimes.clear();
  auto* vtxT =  &intTimes.emplace_back();
  int ntr = idxsort.size();
  int cntTmp = 0;
  for (int i=0; i<ntr;i++) {
    int icur = idxsort[i];
    const auto& trc = pool[icur];
    if (!vtxT->isCompatible(trc.timeEst, margin2, cut)) {
      htv->Fill(vtxT->timeEst.getTimeStamp(),vtxT->idEnd-vtxT->idStart+1); // old one
      vtxT->print();
      vtxT =  &intTimes.emplace_back();
      cntTmp = 0;
    }
    vtxT->addTrack(i,trc.timeEst);    
    printf("%4d/%3d %5d %e/%e -> %e/%e\n", int(intTimes.size()-1), cntTmp++, i,
           trc.timeEst.getTimeStamp(), trc.timeEst.getTimeStampError(),
           vtxT->timeEst.getTimeStamp(), std::sqrt(vtxT->timeEst.getTimeStampError()));
  }
}


void fitVertex(Vertexer &vtf, gsl::span<TrackVF> pool, gsl::span<int> idxsort) {
  std::vector<TrackVF> seeds(idxsort.size());
  Vertex vtx;
  float scaleSigma2 = 100;
  while (vtf.fitVertex(pool, idxsort, vtx, scaleSigma2, false, false)) {
    finalizeVertex(vtx, pool, idxsort);
  }
}
