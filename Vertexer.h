#ifndef O2_VERTEXER_H
#define O2_VERTEXER_H

#include <array>

class Vertexer {

 public:
  
  struct Vertex {
    float xyz[3];      ///< vertex position
    float cov[6];      ///< sym.cov. matrix in lower triangle representation
    float chi2;        ///< total chi^2
    int   nTracks;     ///< number of tracks associated
    unsigned int tstamp;///< vertex stamp (time, event id, etc)
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
    float wgh;         ///< track weight wrt current vertex seed
    short vtxID;       ///< assigned vertex
    //
    bool CanUse()                      const {return vtxID==kNoVtx;}
    bool CanUse(float zmin,float zmax) const {return vtxID==kNoVtx && z>zmin && z<zmax;}
    bool operator < (const Track& trc) const {return z<trc.z;}

    float getResiduals(const Vertex& vtx, float& dy, float& dz) const
    {
      // get residuals (Y and Z DCA in track frame) and calculate chi2
      float dx = vtx.xyz[0]*cosAlp + vtx.xyz[1]*sinAlp - x;     // VX rotated to track frame - trackX
      dy = y + tgP*dx - (-vtx.xyz[0]*sinAlp + vtx.xyz[1]*cosAlp); 
      dz = z + tgP*dx - vtx.xyz[2];
      return 0.5f*(dy*dy*sig2YI  + dz*dz*sig2ZI) + dy*dz*sigYZI;
    }
  };

  bool FitVertex(std::vector<Track> &tracks, Vertex &vtx, float &scaleSigma2, bool useConstraint, bool fillError);

  void setConstraint(float x, float y, float sigyy, float sigyz, float sigzz);
  void  setTukey(float t)  {mTukey2I = t>0.f ? 1.f/(t*t) : 1.f/(kDefTukey*kDefTukey);}
  float getTukey()   const;

 private:
  
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
