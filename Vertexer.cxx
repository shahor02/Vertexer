#include "Vertexer.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"

constexpr float  Vertexer::kAlmost0F;
constexpr double Vertexer::kAlmost0D;
constexpr float  Vertexer::kHugeF;
constexpr float  Vertexer::kDefTukey;

//______________________________________________
bool Vertexer::fitVertex(gsl::span<Vertexer::Track> tracks, gsl::span<int> idxsort, Vertex &vtx, 
			 float &scaleSigma2, bool useConstraint, bool fillErrors)
{
  // fit vertex taking provided vertex as a seed
  // tracks pool may contain arbitrary number of tracks, only those which are in the idxsort (indices of tracks sorted in time)
  // will be used.
  
  int ntAcc=0,ntr = idxsort.size();
  if (ntr<mMinTracksPerVtx) return false;
  //
  Weights cumw;
  cumw.setScale(scaleSigma2, mTukey2I);
  //
  int it = 100;
  while(it-- && cumw.scaleSigma2>0.5) {
    for (int i : idxsort) {
      if (tracks[i].canUse()) {
        accountTrack(tracks[i], vtx, cumw);
      }
    }
    //
    if (cumw.nTracks<mMinTracksPerVtx) {
      if (upscale(cumw)) {
        continue; // redo with stronger rescaling
      }
    }
    if (useConstraint) {
      applyConstraint(cumw);
    }
    if (!solveVertex(vtx, cumw, fillErrors)) {
      return false;
    }
    vtx.print();
    printf("Sigma: %e T: %e\n",cumw.scaleSigma2, cumw.tstamp);
    cumw.resetForNewIteration(mTukey2I);
  }
  return true;
}

//___________________________________________________________________
void Vertexer::accountTrack(Vertexer::Track &trc, Vertex &vtx, Vertexer::Weights &cumw) const
{
  float dy,dz, chi2T = trc.getResiduals(vtx, dy, dz); // track-vertex residuals and chi2
  float wghT = (1.f-chi2T*cumw.scaleSig2ITuk2I);           // weighted distance to vertex
  if (wghT<kAlmost0F)  {
    trc.wgh = 0.f;
    return;
  }
  float syyI(trc.sig2YI),szzI(trc.sig2ZI),syzI(trc.sigYZI);
    
  cumw.wghSum  += wghT;
  cumw.wghChi2 += wghT*chi2T;
  //
  syyI *= wghT;
  syzI *= wghT;
  szzI *= wghT;
  trc.wgh = wghT;
  //
  // aux variables
  double tmpSP = trc.sinAlp*trc.tgP, tmpCP = trc.cosAlp*trc.tgP
    ,tmpSC = trc.sinAlp+tmpCP, tmpCS = -trc.cosAlp+tmpSP
    ,tmpCL = trc.cosAlp*trc.tgL, tmpSL = trc.sinAlp*trc.tgL
    ,tmpYXP = trc.y-trc.tgP*trc.x, tmpZXL = trc.z-trc.tgL*trc.x
    ,tmpCLzz = tmpCL*szzI, tmpSLzz = tmpSL*szzI, tmpSCyz = tmpSC*syzI, tmpCSyz = tmpCS*syzI
    ,tmpCSyy = tmpCS*syyI, tmpSCyy = tmpSC*syyI, tmpSLyz = tmpSL*syzI, tmpCLyz = tmpCL*syzI;
  //
  // symmetric matrix equation 
  cumw.cxx += tmpCL*(tmpCLzz+tmpSCyz+tmpSCyz)+tmpSC*tmpSCyy;          // dchi^2/dx/dx
  cumw.cxy += tmpCL*(tmpSLzz+tmpCSyz)+tmpSL*tmpSCyz+tmpSC*tmpCSyy;    // dchi^2/dx/dy
  cumw.cxz += -trc.sinAlp*syzI-tmpCLzz-tmpCP*syzI;                   // dchi^2/dx/dz
  cumw.cx0 += -(tmpCLyz+tmpSCyy)*tmpYXP-(tmpCLzz+tmpSCyz)*tmpZXL;     // RHS 
  //
  cumw.cyy += tmpSL*(tmpSLzz+tmpCSyz+tmpCSyz)+tmpCS*tmpCSyy;          // dchi^2/dy/dy
  cumw.cyz += -(tmpCSyz+tmpSLzz);                                     // dchi^2/dy/dz
  cumw.cy0 += -tmpYXP*(tmpCSyy+tmpSLyz)-tmpZXL*(tmpCSyz+tmpSLzz);     // RHS
  //
  cumw.czz += szzI;                                                    // dchi^2/dz/dz
  cumw.cz0 += tmpZXL*szzI+tmpYXP*syzI;                                 // RHS
  //
  cumw.tstamp = float(cumw.tstamp*cumw.nTracks + trc.tstamp)/(1+cumw.nTracks);
  cumw.nTracks++;
}

//___________________________________________________________________
bool Vertexer::solveVertex(Vertex &vtx, Vertexer::Weights &cumw, bool fillErrors) const
{
  ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> mat;
  mat(0,0) = cumw.cxx;
  mat(0,1) = cumw.cxy;
  mat(0,2) = cumw.cxz;
  mat(1,1) = cumw.cyy;
  mat(1,2) = cumw.cyz;
  mat(2,2) = cumw.czz;
  if (!mat.InvertFast()) {
    printf("Failed to invert matrix\n");
    std::cout << mat << "\n";
    return false;
  }
  ROOT::Math::SVector<double, 3> rhs(cumw.cx0, cumw.cy0, cumw.cz0);
  auto sol = mat*rhs;
  vtx.setXYZ(sol(0),sol(1),sol(2));
  if (fillErrors) {
    vtx.setCov(mat(0,0), mat(1,0), mat(1,1), mat(2,0), mat(2,1), mat(2,2));
  }
  vtx.setNContributors(cumw.nTracks);
  vtx.setChi2( 2.f*(cumw.nTracks-cumw.wghSum)/cumw.scaleSig2ITuk2I );         // calculate chi^2
  cumw.scaleSigma2 = cumw.wghChi2/cumw.wghSum;
  return true;
}
  
//___________________________________________________________________
void Vertexer::setConstraint(float x, float y, float sigyy, float sigyz, float sigzz)
{
  // set mean vertex constraint and its errors 
  mXYConstraint[0] = x;
  mXYConstraint[1] = y;
  double det = sigyy*sigzz - sigyz*sigyz;
  if (det<=kAlmost0D || sigyy<kAlmost0D || sigzz<kAlmost0D) {
    printf("Error: wrong settings for vertex constraint: %e %e | %e %e %e\n",x,y,sigyy,sigyz,sigzz);
    exit(1);
  }
  mXYConstraintInvErr[0] = sigzz / det;
  mXYConstraintInvErr[2] = sigyy / det;
  mXYConstraintInvErr[1] =-sigyz / det;
}

//______________________________________________
float Vertexer::getTukey() const
{
  // convert 1/tukey^2 to tukey
  return sqrtf(1./mTukey2I);
}
