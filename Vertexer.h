#ifndef O2_VERTEXER_H
#define O2_VERTEXER_H

#include <array>

#include "CommonDataFormat/TimeStamp.h"
#include "CommonDataFormat/RangeReference.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "gsl/span"

using TimeEst = o2::dataformats::TimeStampWithError<float, float>; // FIXME
using Vertex = o2::dataformats::Vertex<TimeEst>;
using V2TRef = o2::dataformats::RangeReference<int, int>;


class Vertexer {

public:
  enum class FitStatus : int { Failure, PoolEmpty, NotEnoughTracks, OK };

  struct StatAccumulator {
    // mean / RMS accumulator
    double sum  = 0.;
    double sum2 = 0.;
    double wsum = 0.;
    int n = 0;
    void add(float v, float w=1.) {
      auto c = v*w;
      sum += c;
      sum2 += c*v;
      wsum += w;
      n++;
    }
    double getMean() const { return wsum > 0. ? sum/wsum : 0.;}
    bool getMeanRMS2(double &mean, double &rms2) const
    {
      if (!wsum) {
        mean = rms2 = 0;
        return false;
      }
      auto wi = 1./wsum;
      mean = sum*wi;
      rms2 = sum2*wi - mean*mean;
      return true;
    }
    bool getMeanRMS2(float &mean, float &rms2) const
    {
      if (!wsum) {
        mean = rms2 = 0;
        return false;
      }
      auto wi = 1./wsum;
      mean = sum*wi;
      rms2 = sum2*wi - mean*mean;
      return true;      
    }
    StatAccumulator& operator+=(const StatAccumulator& other) {
      sum += other.sum;
      sum2 += other.sum2;
      wsum += other.wsum;
      return *this;
    }
    
    StatAccumulator operator+(const StatAccumulator& other) const {
      StatAccumulator res = *this;
      res += other;
      return res;
    }
    
    void clear() {
      sum = sum2 = wsum = 0.;
      n = 0;
    }
  };
  
  ///< weights and scaling params for current vertex
  struct VertexSeed : public Vertex {
    double wghSum = 0.;  // sum of tracks weights
    double wghChi2 = 0.; // sum of tracks weighted chi2's
    double tMeanAcc = 0.; // sum of track times * inv.err^2
    double tMeanAccErr = 0.; // some of tracks times inv.err^2
    double cxx = 0., cyy = 0., czz = 0., cxy = 0., cxz = 0., cyz = 0., cx0 = 0., cy0 = 0., cz0 = 0.; // elements of lin.equation matrix
    float scaleSigma2 = 1.e9;  // Tukey scaling parameter
    float scaleSigma2Prev = 1.;
    float scaleSig2ITuk2I = 0; // inverse squared Tukey parameter scaled by scaleSigma2
    bool useConstraint = true;
    bool fillErrors = true;
    int nStuck = 0;
    StatAccumulator statDZNeg, statDZPos;
    StatAccumulator statDTNeg, statDTPos;
    
    void setScale(float scale2, float tukey2I) {
      scaleSigma2Prev = scaleSigma2;
      scaleSigma2 = scale2;
      scaleSig2ITuk2I = tukey2I / scale2;
    }

    void resetForNewIteration() {
      setNContributors(0);
      //setTimeStamp({0., 0.});
      wghSum = 0;
      wghChi2 = 0;
      tMeanAcc = 0;
      tMeanAccErr = 0;
      cxx = cyy = czz = cxy = cxz = cyz = cx0 = cy0 = cz0 = 0.;
      statDZNeg.clear();
      statDZPos.clear();
      statDTNeg.clear();
      statDTPos.clear();
    }

    VertexSeed() = default;
    VertexSeed(const Vertex &vtx, bool _constraint, bool _errors)
        : Vertex(vtx), useConstraint(_constraint), fillErrors(_errors) {}

    void print() const
    {
      auto terr2 = tMeanAccErr > 0 ? 1./tMeanAccErr : 0.;
      printf("VtxSeed: Scale: %+e ScalePrev: %+e |NStuck: %d| WChi2: %e WSum: %e | TMean: %e TMeanE: %e\n",
             scaleSigma2, scaleSigma2Prev, nStuck, wghChi2, wghSum, tMeanAcc*terr2, std::sqrt(terr2));
      double dZP, rmsZP, dZN, rmsZN, dTP, rmsTP, dTN, rmsTN;
      double dZ, rmsZ, dT, rmsT;
      statDZNeg.getMeanRMS2(dZN, rmsZN);
      statDZPos.getMeanRMS2(dZP, rmsZP);
      statDTNeg.getMeanRMS2(dTN, rmsTN);
      statDTPos.getMeanRMS2(dTP, rmsTP);
      auto statDZ = statDZNeg+statDZPos;
      auto statDT = statDTNeg+statDTPos;
      statDZ.getMeanRMS2(dZ, rmsZ);
      statDT.getMeanRMS2(dT, rmsT);
      printf("DZN: %+e / %e w: %.1e/%4d | DZP: %+e / %e w: %.1e/%4d || DZ: %+e / %e\n"
             "DTN: %+e / %e w: %.1e/%4d | DTP: %+e / %e w: %.1e/%4d || DT: %+e / %e\n",
             dZN, std::sqrt(rmsZN), statDZNeg.wsum, statDZNeg.n, dZP, std::sqrt(rmsZP), statDZPos.wsum, statDZPos.n, dZ, std::sqrt(rmsZ),
             dTN, std::sqrt(rmsTN), statDTNeg.wsum, statDTNeg.n, dTP, std::sqrt(rmsTP), statDTPos.wsum, statDTPos.n, dT, std::sqrt(rmsT) );
      Vertex::print();
    }
  };

  struct Track {
    /** Straight track parameterization in the frame defined by alpha angle.
        Assumed to be defined in the proximity to vertex, so that the
       straight-line extrapolation Y=mY+mTgP*(x-mX) and Z=mZ+mTgL*(x-mX) is
       precise
    */
    enum { kUsed, kNoVtx = -1, kDiscarded = kNoVtx - 1 };
    float x;      ///< reference X
    float y;      ///< Y at X
    float z;      ///< Z at X
    float sig2YI; ///< YY component of inverse cov.matrix
    float sig2ZI; ///< ZZ component of inverse cov.matrix
    float sigYZI; ///< YZ component of inverse cov.matrix
    float tgP;    ///< tangent(phi) in tracking frame
    float tgL;    ///< tangent(lambda)
    float cosAlp; ///< cos of alpha frame
    float sinAlp; ///< sin of alpha frame

    TimeEst timeEst;
    float wgh = 0.; ///< track weight wrt current vertex seed
    uint32_t entry = 0;
    int16_t bin = -1; // seeds histo bin
    uint8_t srcID = 0;
    uint8_t flags = 0;
    int vtxID = kNoVtx; ///< assigned vertex

    //
    bool canAssign() const { return wgh > 0. && vtxID == kNoVtx; }
    bool canUse() const { return vtxID == kNoVtx; }
    bool canUse(float zmin, float zmax) const {
      return canUse() && (z > zmin && z < zmax);
    }
    bool operator<(const Track &trc) const { return z < trc.z; }

    float getZForXY(float vx, float vy) const {
      return z + tgL * (vx * cosAlp + vy * sinAlp - x);
    }
    
    float getResiduals(const Vertex &vtx, float &dy, float &dz) const {
      // get residuals (Y and Z DCA in track frame) and calculate chi2
      float dx = vtx.getX() * cosAlp + vtx.getY() * sinAlp - x; // VX rotated to track frame - trackX
      dy = y + tgP * dx - (-vtx.getX() * sinAlp + vtx.getY() * cosAlp);
      dz = z + tgL * dx - vtx.getZ();
      return (dy * dy * sig2YI + dz * dz * sig2ZI) + 2.* dy * dz * sigYZI;
    }

    Track() = default;
    Track(const o2::track::TrackParCov &src, const TimeEst &t_est, uint32_t _entry, uint8_t _srcID)
    : x(src.getX()), y(src.getY()), z(src.getZ()), tgL(src.getTgl()), tgP(src.getSnp() / std::sqrt(1. - src.getSnp()) * (1. + src.getSnp())),
      timeEst(t_est), entry(_entry), srcID(_srcID) {
      o2::utils::sincosf(src.getAlpha(), sinAlp, cosAlp);
      auto det = src.getSigmaY2() * src.getSigmaZ2() -
                 src.getSigmaZY() * src.getSigmaZY();
      auto detI = 1. / det;
      sig2YI = src.getSigmaZ2() * detI;
      sig2ZI = src.getSigmaY2() * detI;
      sigYZI = -src.getSigmaZY() * detI;
    }
  };

  struct VertexingInput {
    gsl::span<Vertexer::Track> tracks;
    gsl::span<int> idxsort;
    float scaleSigma2 = 10;
    bool useConstraint = false;
    bool fillErrors = true;

    ClassDefNV(VertexingInput,1);
  };

  struct SeedHisto {
    float range = 20;
    float binSize = 0.5;
    float binSizeInv = 0.;
    int nFilled = 0;
    std::vector<int> data;

    int size() const { return data.size(); }
    
    void init(float _range = 20., float _binsize = 0.5) {
      range = _range;
      binSize = _binsize;
      auto zr = 2*range;
      int nzb = zr/binSize;
      if (nzb*binSize<zr-Vertexer::kAlmost0F) {
        nzb++;
      }
      binSizeInv = 1./binSize;
      range = nzb*binSize/2.;
      data.resize(nzb);
    }

    void fill(float z) {
      incrementBin(findBin(z));
    }

    void incrementBin(int bin) {
      data[bin]++;
      nFilled++;
    }

    void decrementBin(int bin) {
      data[bin]--;
      nFilled--;
    }
    
    int findBin(float z) {
      auto d = z+range;
      if (d<0.) {
        return 0;
      }
      uint32_t n = d*binSizeInv;
      return n<data.size() ? n : data.size()-1;
    }

    int findHighestPeakBin() const {
      if (nFilled<2) {
        return -1;
      }
      int n = data.size(), maxBin = -1, maxv = 0;
      for (int i=0;i<n;i++) {
        if (data[i]>maxv) {
          maxv = data[(maxBin=i)];
        }
      }
      return maxBin;
    }

    bool isValidBin(int ib) const {
      return static_cast<uint32_t>(ib)<data.size();
    }
    
    float getBinCenter(int ib) const {
      return (ib+0.5)*binSize - range; // no check for being in the range!!!
    }

    void discardBin(int ib) { // no check for being in the range!!!
      nFilled -= data[ib];
      data[ib] = 0;
    }
  };
  
  void init();
  bool findVertices(const VertexingInput& input, std::vector<Vertex> &vertices, std::vector<int>& vertexTrackIDs, std::vector<V2TRef>& v2tRefs);
  bool findVertex(const VertexingInput& input, Vertex &vtx);
  FitStatus fitIteration(const VertexingInput& input, Vertexer::VertexSeed &vtxSeed) const;
  void setConstraint(float x, float y, float sigyy, float sigyz, float sigzz);

  void setTukey(float t) {
    mTukey2I = t > 0.f ? 1.f / (t * t) : 1.f / (kDefTukey * kDefTukey);
  }
  float getTukey() const;

  float getMinScale2() const { return mMinScale2; }
  void setMinScale2(float v) { mMinScale2 = v < 0.01 ? 1. : v; }

  float getMaxScale2() const { return mMaxScale2; }
  void setMaxScale2(float v) { mMaxScale2 = v < 1.e6 ? 1.e6 : v; }

  void setMaxIterations(int n) { mMaxIterations = n > 2 ? n : 2; }
  int getMaxIterations() const { return mMaxIterations; }

  void setMaxChi2Change(float v) { mMaxChi2Change = v > 1.e-3 ? v : 1.e-3; }
  float getMaxChi2Change() const { return mMaxChi2Change; }

  void finalizeVertex(const VertexingInput& input, const Vertex& vtx, std::vector<Vertex>& vertices, std::vector<V2TRef>& v2tRefs, std::vector<int> &vertexTrackIDs, SeedHisto& histo);
  
 private:
  void accountTrack(Track &trc, VertexSeed &vtxSeed) const;
  bool solveVertex(VertexSeed &vtxSeed) const;
  bool stopIterations(VertexSeed &vtxSeed, Vertex &vtx) const;
  TimeEst timeEstimate(const VertexingInput& input) const;
  float findZSeedHistoPeak() const;

  
  //___________________________________________________________________
  void applyConstraint(VertexSeed &vtxSeed) const {
    // impose meanVertex constraint, i.e. account terms
    // (V_i-Constrain_i)^2/sig2constr_i for i=X,Y in the fit chi2 definition
    vtxSeed.cxx += mXYConstraintInvErr[0];
    vtxSeed.cyy += mXYConstraintInvErr[1];
    vtxSeed.cx0 += mXYConstraintInvErr[0] * mXYConstraint[0];
    vtxSeed.cy0 += mXYConstraintInvErr[1] * mXYConstraint[1];
  }

  //___________________________________________________________________
  bool upscaleSigma(VertexSeed &vtxSeed) const {
    // scale upward the scaleSigma2 if needes
    if (vtxSeed.scaleSigma2 < mMaxScale2) {
      auto s = vtxSeed.scaleSigma2 * 9.;
      vtxSeed.setScale(s > mMaxScale2 ? mMaxScale2 : s, mTukey2I);
      return true;
    }
    return false;
  }

  std::array<float, 2> mXYConstraint = {0.f, 0.f}; ///< nominal vertex constraint
  std::array<float, 3> mXYConstraintInvErr = {1.0f, 0.f, 1.0f}; ///< nominal vertex constraint inverted errors^2
  //
  float mTukey2I = 1. / (kDefTukey * kDefTukey); ///< 1./[Tukey parameter]^2
  float mMinScale2 = 1.;                         ///< min slaling factor^2
  float mMaxScale2 = 1.e6;                       ///< max slaling factor^2
  float mMaxChi2Change = 0.001; ///< max chi2 change to continue iterations
  
  int mMinTracksPerVtx = 2;                  ///< min N tracks per vertex
  int mMaxIterations = 100;                  ///< max iterations per vertex fit
  int mMaxNStuck = 1;                        ///< max number of scaling-non-decreasing iterations
  static constexpr float kDefTukey = 5.0f;   ///< def.value for tukey constant
  static constexpr float kHugeF = 1.e12;     ///< very large float
  static constexpr float kAlmost0F = 1e-12;  ///< tiny float
  static constexpr double kAlmost0D = 1e-16; ///< tiny double

  ClassDefNV(Vertexer, 1);
};

#endif
