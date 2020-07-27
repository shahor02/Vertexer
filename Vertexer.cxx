#include "Vertexer.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"

constexpr float Vertexer::kAlmost0F;
constexpr double Vertexer::kAlmost0D;
constexpr float Vertexer::kHugeF;
constexpr float Vertexer::kDefTukey;

//______________________________________________
bool Vertexer::fitVertex(gsl::span<Vertexer::Track> tracks,
                         gsl::span<int> idxsort, Vertex &vtx, float scaleSigma2,
                         bool useConstraint, bool fillErrors) {
  // fit vertex taking provided vertex as a seed
  // tracks pool may contain arbitrary number of tracks, only those which are in
  // the idxsort (indices of tracks sorted in time) will be used.

  int ntAcc = 0, ntr = idxsort.size();
  if (ntr < mMinTracksPerVtx)
    return false;
  //
  VertexSeed vtxseed(vtx, useConstraint, fillErrors);
  vtxseed.setScale(scaleSigma2, mTukey2I);
  vtx.setChi2(1.e30);
  //
  FitStatus result;
  int it = mMaxIterations;
  bool found = false;
  while (it-- && vtxseed.scaleSigma2 > mMinScale2) {
    vtxseed.resetForNewIteration();
    result = fitIteration(tracks, idxsort, vtxseed);

    if (result == FitStatus::OK) {
      found = stopIterations(vtxseed, vtx);
      vtx.print();
      LOG(INFO) << "#" << it << " Sigma: " << vtxseed.scaleSigma2
                << " T: " << vtx.getTimeStamp();
      if (found) {
        break;
      }
    } else if (result == FitStatus::NotEnoughTracks) {
      if (upscaleSigma(vtxseed)) {
        LOG(INFO) << "Upscaling scale to " << vtxseed.scaleSigma2;
        continue; // redo with stronger rescaling
      } else {
        break;
      }
    } else if (result == FitStatus::PoolEmpty || result == FitStatus::Failure) {
      break;
    } else {
      LOG(FATAL) << "Unknown fit status " << int(result);
    }
  }
  return found;
}

//___________________________________________________________________
Vertexer::FitStatus
Vertexer::fitIteration(gsl::span<Vertexer::Track> tracks,
                       gsl::span<int> idxsort,
                       Vertexer::VertexSeed &vtxseed) const {
  int nTested = 0;
  for (int i : idxsort) {
    if (tracks[i].canUse()) {
      nTested++;
      accountTrack(tracks[i], vtxseed);
    }
  }
  if (vtxseed.getNContributors() < mMinTracksPerVtx) {
    return nTested < mMinTracksPerVtx ? FitStatus::PoolEmpty
                                      : FitStatus::NotEnoughTracks;
  }
  if (vtxseed.useConstraint) {
    applyConstraint(vtxseed);
  }
  if (!solveVertex(vtxseed)) {
    return FitStatus::Failure;
  }
  return FitStatus::OK;
}

//___________________________________________________________________
void Vertexer::accountTrack(Vertexer::Track &trc,
                            Vertexer::VertexSeed &vtxseed) const {
  float dy, dz,
      chi2T =
          trc.getResiduals(vtxseed, dy, dz); // track-vertex residuals and chi2
  float wghT =
      (1.f - chi2T * vtxseed.scaleSig2ITuk2I); // weighted distance to vertex
  if (wghT < kAlmost0F) {
    trc.wgh = 0.f;
    return;
  }
  float syyI(trc.sig2YI), szzI(trc.sig2ZI), syzI(trc.sigYZI);

  vtxseed.wghSum += wghT;
  vtxseed.wghChi2 += wghT * chi2T;
  //
  syyI *= wghT;
  syzI *= wghT;
  szzI *= wghT;
  trc.wgh = wghT;
  //
  // aux variables
  double tmpSP = trc.sinAlp * trc.tgP, tmpCP = trc.cosAlp * trc.tgP,
         tmpSC = trc.sinAlp + tmpCP, tmpCS = -trc.cosAlp + tmpSP,
         tmpCL = trc.cosAlp * trc.tgL, tmpSL = trc.sinAlp * trc.tgL,
         tmpYXP = trc.y - trc.tgP * trc.x, tmpZXL = trc.z - trc.tgL * trc.x,
         tmpCLzz = tmpCL * szzI, tmpSLzz = tmpSL * szzI, tmpSCyz = tmpSC * syzI,
         tmpCSyz = tmpCS * syzI, tmpCSyy = tmpCS * syyI, tmpSCyy = tmpSC * syyI,
         tmpSLyz = tmpSL * syzI, tmpCLyz = tmpCL * syzI;
  //
  // symmetric matrix equation
  vtxseed.cxx +=
      tmpCL * (tmpCLzz + tmpSCyz + tmpSCyz) + tmpSC * tmpSCyy; // dchi^2/dx/dx
  vtxseed.cxy += tmpCL * (tmpSLzz + tmpCSyz) + tmpSL * tmpSCyz +
                 tmpSC * tmpCSyy;                             // dchi^2/dx/dy
  vtxseed.cxz += -trc.sinAlp * syzI - tmpCLzz - tmpCP * syzI; // dchi^2/dx/dz
  vtxseed.cx0 +=
      -(tmpCLyz + tmpSCyy) * tmpYXP - (tmpCLzz + tmpSCyz) * tmpZXL; // RHS
  //
  vtxseed.cyy +=
      tmpSL * (tmpSLzz + tmpCSyz + tmpCSyz) + tmpCS * tmpCSyy; // dchi^2/dy/dy
  vtxseed.cyz += -(tmpCSyz + tmpSLzz);                         // dchi^2/dy/dz
  vtxseed.cy0 +=
      -tmpYXP * (tmpCSyy + tmpSLyz) - tmpZXL * (tmpCSyz + tmpSLzz); // RHS
  //
  vtxseed.czz += szzI;                          // dchi^2/dz/dz
  vtxseed.cz0 += tmpZXL * szzI + tmpYXP * syzI; // RHS
  //
  auto &timeV = vtxseed.getTimeStamp();
  auto &timeT = trc.timeEst;
  auto trErr2I = 1. / (timeT.getTimeStampError() * timeT.getTimeStampError());
  auto tAv = timeV.getTimeStamp() * timeV.getTimeStampError() +
             wghT * timeT.getTimeStamp() *
                 trErr2I; // <t_prev>*(1/err_prev)^2 [until end of the fit we
                          // keep inverse error^2]
  auto tAvE = timeV.getTimeStampError() + wghT * trErr2I;
  timeV.setTimeStamp(tAv / tAvE);
  timeV.setTimeStampError(tAvE);
  vtxseed.addContributor();
}

//___________________________________________________________________
bool Vertexer::solveVertex(Vertexer::VertexSeed &vtxseed) const {
  ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> mat;
  mat(0, 0) = vtxseed.cxx;
  mat(0, 1) = vtxseed.cxy;
  mat(0, 2) = vtxseed.cxz;
  mat(1, 1) = vtxseed.cyy;
  mat(1, 2) = vtxseed.cyz;
  mat(2, 2) = vtxseed.czz;
  if (!mat.InvertFast()) {
    printf("Failed to invert matrix\n");
    std::cout << mat << "\n";
    return false;
  }
  ROOT::Math::SVector<double, 3> rhs(vtxseed.cx0, vtxseed.cy0, vtxseed.cz0);
  auto sol = mat * rhs;
  vtxseed.setXYZ(sol(0), sol(1), sol(2));
  if (vtxseed.fillErrors) {
    vtxseed.setCov(mat(0, 0), mat(1, 0), mat(1, 1), mat(2, 0), mat(2, 1),
                   mat(2, 2));
  }
  vtxseed.setChi2(2.f * (vtxseed.getNContributors() - vtxseed.wghSum) /
                  vtxseed.scaleSig2ITuk2I); // calculate chi^2
  auto newScale = vtxseed.wghChi2 / vtxseed.wghSum;
  vtxseed.setScale(newScale < mMinScale2 ? mMinScale2 : newScale, mTukey2I);
  return true;
}

//___________________________________________________________________
bool Vertexer::stopIterations(VertexSeed &vtxSeed, Vertex &vtx) const {
  // decide if new iteration should be done, prepare next one if needed
  // if scaleSigma2 reached its lower limit stop
  bool stop = false;
  while (1) {
    if (vtxSeed.scaleSigma2 <= mMinScale2 + kAlmost0F) {
      stop = true;
      LOG(INFO) << "stop on simga :" << vtxSeed.scaleSigma2
                << " prev: " << vtxSeed.scaleSigma2Prev;
      //      break;
    }
    auto dchi = (vtx.getChi2() - vtxSeed.getChi2()) / vtxSeed.getChi2();
    if (dchi > 0 && dchi < mMaxChi2Change) {
      stop = true;
      LOG(INFO) << "stop on chi2 :";
      //      break;
    }
    auto dx = vtxSeed.getX() - vtx.getX(), dy = vtxSeed.getY() - vtx.getY(),
         dz = vtxSeed.getZ() - vtx.getZ();
    auto dst = std::sqrt(dx * dx + dy * dy + dz * dz);

    LOG(INFO) << "dChi:" << vtx.getChi2() << "->" << vtxSeed.getChi2()
              << " :-> " << dchi;
    LOG(INFO) << "dx: " << dx << " dy: " << dy << " dz: " << dz << " -> "
              << dst;

    break;
  }

  vtx = vtxSeed;
  auto &vTStamp = vtx.getTimeStamp();
  if (vTStamp.getTimeStampError() > 0) {
    vTStamp.setTimeStampError(
        1. / std::sqrt(vTStamp.getTimeStampError())); // seed accumulates
                                                      // inverse squared error!
  }
  if (!stop) {
    auto scforce = vtxSeed.scaleSigma2Prev * 0.5;
    if (vtxSeed.scaleSigma2 > scforce) {
      auto sav = vtxSeed.scaleSigma2Prev;
      vtxSeed.setScale(scforce < mMinScale2 ? mMinScale2 : scforce, mTukey2I);
      vtxSeed.scaleSigma2Prev = sav;
      LOG(INFO) << "Forcing scale from " << vtxSeed.scaleSigma2Prev << " to "
                << vtxSeed.scaleSigma2;
    }
  }

  return stop;
}

//___________________________________________________________________
void Vertexer::setConstraint(float x, float y, float sigyy, float sigyz,
                             float sigzz) {
  // set mean vertex constraint and its errors
  mXYConstraint[0] = x;
  mXYConstraint[1] = y;
  double det = sigyy * sigzz - sigyz * sigyz;
  if (det <= kAlmost0D || sigyy < kAlmost0D || sigzz < kAlmost0D) {
    printf("Error: wrong settings for vertex constraint: %e %e | %e %e %e\n", x,
           y, sigyy, sigyz, sigzz);
    exit(1);
  }
  mXYConstraintInvErr[0] = sigzz / det;
  mXYConstraintInvErr[2] = sigyy / det;
  mXYConstraintInvErr[1] = -sigyz / det;
}

//______________________________________________
float Vertexer::getTukey() const {
  // convert 1/tukey^2 to tukey
  return sqrtf(1. / mTukey2I);
}
