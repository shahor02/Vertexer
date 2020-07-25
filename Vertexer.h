#ifndef O2_VERTEXER_H
#define O2_VERTEXER_H

#include <array>

#include "CommonDataFormat/TimeStamp.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "MathUtils/Utils.h"
#include "gsl/span"

using TimeEst = o2::dataformats::TimeStampWithError<float, float>; // FIXME
using Vertex = o2::dataformats::Vertex<TimeEst>;

class Vertexer {

 public:
  
  ///< weights and scaling params for current vertex
  struct Weights {
    int nTracks = 0;
    double wghSum = 0.;  // sum of tracks weights
    double wghChi2 = 0.; // sum of tracks weighted chi2's
    double cxx=0.,cyy=0.,czz=0.,cxy=0.,cxz=0.,cyz=0., cx0=0.,cy0=0.,cz0=0.; // elements of lin.equation matrix
    float tstamp = 0;
    float scaleSigma2 = 1.; // Tukey scaling parameter
    float scaleSig2ITuk2I = 0; // inverse squared Tukey parameter scaled by scaleSigma2    

    void setScale(float scale2, float tukey2I) {
      scaleSigma2 = scale2;
      scaleSig2ITuk2I = tukey2I/scale2;
    }
    void resetForNewIteration(float tukey2I) {
      scaleSig2ITuk2I = tukey2I/scaleSigma2;
      nTracks = 0;
      wghSum = 0;
      wghChi2 = 0;
      tstamp = 0;
      cxx = cyy = czz = cxy = cxz = cyz = cx0 = cy0 = cz0 = 0.;
    }
  };
  
  struct Track {
    /** Straight track parameterization in the frame defined by alpha angle.
	Assumed to be defined in the proximity to vertex, so that the straight-line
	extrapolation Y=mY+mTgP*(x-mX) and Z=mZ+mTgL*(x-mX) is precise
    */
    enum {kUsed, kNoVtx=-1,kDiscarded=kNoVtx-1};
    unsigned int tstamp;///< track stamp (time, event id, etc)
    float x;           ///< reference X
    float y;           ///< Y at X
    float z;           ///< Z at X
    float sig2YI;      ///< YY component of inverse cov.matrix
    float sig2ZI;      ///< ZZ component of inverse cov.matrix
    float sigYZI;      ///< YZ component of inverse cov.matrix
    float tgP;         ///< tangent(phi) in tracking frame
    float tgL;         ///< tangent(lambda)
    float cosAlp;      ///< cos of alpha frame
    float sinAlp;      ///< sin of alpha frame

    TimeEst timeEst;
    float wgh = 0.;         ///< track weight wrt current vertex seed
    uint32_t entry = 0;
    uint8_t srcID = 0;  
    uint8_t flags = 0;
    short vtxID = kNoVtx;   ///< assigned vertex
    
    //
    bool canUse()                      const {return vtxID==kNoVtx;}
    bool canUse(float zmin,float zmax) const {return vtxID==kNoVtx && z>zmin && z<zmax;}
    bool operator < (const Track& trc) const {return z<trc.z;}

    float getResiduals(const Vertex& vtx, float& dy, float& dz) const
    {
      // get residuals (Y and Z DCA in track frame) and calculate chi2
      float dx = vtx.getX()*cosAlp + vtx.getY()*sinAlp - x;     // VX rotated to track frame - trackX
      dy = y + tgP*dx - (-vtx.getX()*sinAlp + vtx.getY()*cosAlp); 
      dz = z + tgP*dx - vtx.getZ();
      return 0.5f*(dy*dy*sig2YI  + dz*dz*sig2ZI) + dy*dz*sigYZI;
    }

    Track() = default;
    Track(const o2::track::TrackParCov& src, const TimeEst& t_est, uint32_t _entry, uint8_t srcID)
    {
      x = src.getX();
      y = src.getY();
      z = src.getZ();
      o2::utils::sincosf(src.getAlpha(), sinAlp, cosAlp);
      tgL = src.getTgl();
      tgP = src.getSnp() / std::sqrt(1.-src.getSnp())*(1.+src.getSnp());
      auto det = src.getSigmaY2()*src.getSigmaZ2() - src.getSigmaZY()*src.getSigmaZY();
      auto detI = 1./det;
      sig2YI = src.getSigmaZ2()*detI;
      sig2ZI = src.getSigmaY2()*detI;
      sigYZI =-src.getSigmaZY()*detI;
      timeEst = t_est;
    }

    //    ClassDefNV(Track,1);
  };

  bool fitVertex(gsl::span<Track> pool,gsl::span<int> idxsort, Vertex &vtx, float &scaleSigma2, bool useConstraint, bool fillError);
  void accountTrack(Track &trc, Vertex &vtx, Weights &cumw) const;
  void setConstraint(float x, float y, float sigyy, float sigyz, float sigzz);
  void  setTukey(float t)  {mTukey2I = t>0.f ? 1.f/(t*t) : 1.f/(kDefTukey*kDefTukey);}
  float getTukey()   const;

 private:
  
  bool solveVertex(Vertex &vtx, Weights &cumw, bool fillErrors) const;
    
  //___________________________________________________________________
  void applyConstraint(Weights &cumw) const
  {
    // impose meanVertex constraint, i.e. account terms (V_i-Constrain_i)^2/sig2constr_i for i=X,Y in the fit chi2 definition
    cumw.cxx += mXYConstraintInvErr[0];
    cumw.cyy += mXYConstraintInvErr[1];
    cumw.cx0 += mXYConstraintInvErr[0]*mXYConstraint[0];
    cumw.cy0 += mXYConstraintInvErr[1]*mXYConstraint[1];
  }

  //___________________________________________________________________
  bool upscale(Weights &cumw) const
  {
    // scale upward the scaleSigma2 if needes
    cumw.scaleSigma2 *=2;
    return true;
  }
  
  std::array<float,2> mXYConstraint = {0.f,0.f};               ///< nominal vertex constraint
  std::array<float,3> mXYConstraintInvErr = {1.0f, 0.f, 1.0f}; ///< nominal vertex constraint inverted errors^2
  //
  float mTukey2I = 1./(kDefTukey*kDefTukey);                ///< 1./[Tukey parameter]^2
  int mMinTracksPerVtx = 2;                                 ///< min N tracks per vertex
  static constexpr float kDefTukey = 5.0f;                  ///< def.value for tukey constant
  static constexpr float kHugeF = 1.e12;                    ///< very large float
  static constexpr float kAlmost0F = 1e-12;                 ///< tiny float
  static constexpr double kAlmost0D = 1e-16;                ///< tiny double

  ClassDefNV(Vertexer,1);
};


#endif
