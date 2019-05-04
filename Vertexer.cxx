#include "Vertexer.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"

constexpr float  Vertexer::kAlmost0F;
constexpr double Vertexer::kAlmost0D;
constexpr float  Vertexer::kHugeF;
constexpr float  Vertexer::kDefTukey;

//______________________________________________
bool Vertexer::FitVertex(std::vector<Vertexer::Track> &tracks, Vertexer::Vertex &vtx,
			 float &scaleSigma2, bool useConstraint, bool fillErrors)
{
  // fit vertex taking provided vertex as a seed
  
  int ntAcc=0,ntr = tracks.size();
  if (ntr<2) return false;
  //
  double wghSum=0,wghChi2=0; 
  double cxx=0,cxy=0,cxz=0,cx0=0,cyy=0,cyz=0,cy0=0,czz=0,cz0=0;
  float scaleSig2ITuk2I = mTukey2I/scaleSigma2;
  //
  for (int itr=ntr;itr--;) {
    Track &trc = tracks[itr];
    if (!trc.CanUse()) continue;
    /* check stamps compatibility
    if (mCheckStamps) {
      if (!CheckVertexTrackStamps(vtx,trc)) continue; // the track is invalidated or out of range, skip
    }
    */
    float dy,dz, chi2T = trc.getResiduals(vtx, dy, dz); // track-vertex residuals and chi2
    float wghT = (1.f-chi2T*scaleSig2ITuk2I);           // weighted distance to vertex
    if (wghT<kAlmost0F)  {
      trc.wgh = 0.f;
      continue;
    }
    float syyI(trc.sig2YI),szzI(trc.sig2ZI),syzI(trc.sigYZI);
    
    wghSum  += wghT;
    wghChi2 += wghT*chi2T;
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
    cxx += tmpCL*(tmpCLzz+tmpSCyz+tmpSCyz)+tmpSC*tmpSCyy;          // dchi^2/dx/dx
    cxy += tmpCL*(tmpSLzz+tmpCSyz)+tmpSL*tmpSCyz+tmpSC*tmpCSyy;    // dchi^2/dx/dy
    cxz += -trc.sinAlp*syzI-tmpCLzz-tmpCP*syzI;                   // dchi^2/dx/dz
    cx0 += -(tmpCLyz+tmpSCyy)*tmpYXP-(tmpCLzz+tmpSCyz)*tmpZXL;     // RHS 
    //
    cyy += tmpSL*(tmpSLzz+tmpCSyz+tmpCSyz)+tmpCS*tmpCSyy;          // dchi^2/dy/dy
    cyz += -(tmpCSyz+tmpSLzz);                                     // dchi^2/dy/dz
    cy0 += -tmpYXP*(tmpCSyy+tmpSLyz)-tmpZXL*(tmpCSyz+tmpSLzz);     // RHS
    //
    czz += szzI;                                                    // dchi^2/dz/dz
    cz0 += tmpZXL*szzI+tmpYXP*syzI;                                 // RHS
    //
    vtx.tstamp = float(vtx.tstamp*ntAcc + trc.tstamp*100)/(1+ntAcc);
    ntAcc++;
  }
  //
  vtx.nTracks = ntAcc;
  if (ntAcc<mMinTracksPerVtx) return false;
  //
  if (useConstraint) {
    // impose meanVertex constraint, i.e. account terms (V_i-Constrain_i)^2/sig2constr_i for i=X,Y 
    // in the fit chi2 definition
    cxx += mXYConstraintInvErr[0];
    cx0 += mXYConstraintInvErr[0]*mXYConstraint[0];
    cyy += mXYConstraintInvErr[1];
    cy0 += mXYConstraintInvErr[1]*mXYConstraint[1];
  }
  ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> mat;  
  mat(0,0) = cxx;
  mat(0,1) = cxy;
  mat(0,2) = cxz;
  mat(1,1) = cyy;
  mat(1,2) = cyz;
  mat(2,2) = czz;
  if (!mat.InvertFast()) {
    printf("Failed to invert matrix\n");
    std::cout << mat << "\n";
    return false;
  }
  ROOT::Math::SVector<double, 3> rhs(cx0,cy0,cz0);
  auto sol = mat*rhs;
  for (int i=3;i--;) {
    vtx.xyz[i] = sol(i);
  }    
  if (fillErrors) {
    vtx.cov[0] = mat(0,0);
    vtx.cov[1] = mat(1,0);
    vtx.cov[2] = mat(1,1);
    vtx.cov[3] = mat(2,0);
    vtx.cov[4] = mat(2,1);
    vtx.cov[5] = mat(2,2);
    vtx.nTracks = ntAcc;
    vtx.chi2 = 2.f*(ntAcc-wghSum)/scaleSig2ITuk2I;         // calculate chi^2
  }
  scaleSigma2 = wghChi2/wghSum;
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
