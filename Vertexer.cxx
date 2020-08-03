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
                         bool useConstraint, bool fillErrors)
{
  // fit vertex taking provided vertex as a seed
  // tracks pool may contain arbitrary number of tracks, only those which are in
  // the idxsort (indices of tracks sorted in time) will be used.

  int ntAcc = 0, ntr = idxsort.size();
  if (ntr < mMinTracksPerVtx) {
    return false;
  }
  //
  fillZSeedHisto(tracks, idxsort);
  
  VertexSeed vtxSeed(vtx, useConstraint, fillErrors);
  vtxSeed.setScale(scaleSigma2, mTukey2I);
  vtxSeed.scaleSigma2Prev = scaleSigma2;
  vtxSeed.setTimeStamp( timeEstimate(tracks, idxsort) );
  LOG(INFO) << "Start time guess: " << vtxSeed.getTimeStamp();
  vtx.setChi2(1.e30);
  //
  FitStatus result;
  int it = mMaxIterations;
  bool found = false;
  while (it-- && (vtxSeed.scaleSigma2 > mMinScale2 || vtxSeed.scaleSigma2!=vtxSeed.scaleSigma2Prev)) {
    LOG(INFO) << "iter " << it << " with scale=" << vtxSeed.scaleSigma2 << " prevScale=" << vtxSeed.scaleSigma2Prev;
    vtxSeed.resetForNewIteration();
    result = fitIteration(tracks, idxsort, vtxSeed);

    if (result == FitStatus::OK) {
      found = stopIterations(vtxSeed, vtx);
      //      vtx.print();
      LOG(INFO) << "#" << it << " Sigma: " << vtxSeed.scaleSigma2 << " T: " << vtx.getTimeStamp();
      if (found) {
        break;
      }
    } else if (result == FitStatus::NotEnoughTracks) {
      if (upscaleSigma(vtxSeed)) {
        LOG(INFO) << "Upscaling scale to " << vtxSeed.scaleSigma2;
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
  LOG(INFO) << "Stopped with scale=" << vtxSeed.scaleSigma2 << " prevScale=" << vtxSeed.scaleSigma2Prev;
  return found;
}

//___________________________________________________________________
Vertexer::FitStatus Vertexer::fitIteration(gsl::span<Vertexer::Track> tracks, gsl::span<int> idxsort, Vertexer::VertexSeed &vtxSeed) const
{
  int nTested = 0;
  for (int i : idxsort) {
    if (tracks[i].canUse()) {
      nTested++;
      accountTrack(tracks[i], vtxSeed);
    }
  }
  if (vtxSeed.getNContributors() < mMinTracksPerVtx) {
    return nTested < mMinTracksPerVtx ? FitStatus::PoolEmpty
                                      : FitStatus::NotEnoughTracks;
  }
  if (vtxSeed.useConstraint) {
    applyConstraint(vtxSeed);
  }
  if (!solveVertex(vtxSeed)) {
    return FitStatus::Failure;
  }
  return FitStatus::OK;
}

//___________________________________________________________________
void Vertexer::accountTrack(Vertexer::Track &trc, Vertexer::VertexSeed &vtxSeed) const
{
  // deltas defined as track - vertex
  float dy, dz, chi2T = trc.getResiduals(vtxSeed, dy, dz); // track-vertex residuals and chi2
  auto &timeV = vtxSeed.getTimeStamp();
  auto &timeT = trc.timeEst;

  auto dt = timeT.getTimeStamp() - timeV.getTimeStamp();
  auto trErr2I = 1. / (timeT.getTimeStampError() * timeT.getTimeStampError());
  //  chi2T += dt*dt*trErr2I;
  //  chi2T *= 1./3;
  chi2T *= 0.5;
  float wghT = (1.f - chi2T * vtxSeed.scaleSig2ITuk2I); // weighted distance to vertex
  if (wghT < kAlmost0F) {
    trc.wgh = 0.f;
    return;
  }
  float syyI(trc.sig2YI), szzI(trc.sig2ZI), syzI(trc.sigYZI);

  // special statistics
  if (dz<0) {
    vtxSeed.statDZNeg.add(dz,wghT);
    vtxSeed.statDTNeg.add(dt,wghT);    
  }
  else {
    vtxSeed.statDZPos.add(dz,wghT);
    vtxSeed.statDTPos.add(dt,wghT);
  }
  
  //
  vtxSeed.wghSum += wghT;
  vtxSeed.wghChi2 += wghT * chi2T;
  //
  syyI *= wghT;
  syzI *= wghT;
  szzI *= wghT;
  trErr2I *= wghT;
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
  vtxSeed.cxx += tmpCL * (tmpCLzz + tmpSCyz + tmpSCyz) + tmpSC * tmpSCyy;                 // dchi^2/dx/dx
  vtxSeed.cxy += tmpCL * (tmpSLzz + tmpCSyz) + tmpSL * tmpSCyz + tmpSC * tmpCSyy;         // dchi^2/dx/dy
  vtxSeed.cxz += -trc.sinAlp * syzI - tmpCLzz - tmpCP * syzI;                             // dchi^2/dx/dz
  vtxSeed.cx0 += -(tmpCLyz + tmpSCyy) * tmpYXP - (tmpCLzz + tmpSCyz) * tmpZXL;            // RHS
  //
  vtxSeed.cyy += tmpSL * (tmpSLzz + tmpCSyz + tmpCSyz) + tmpCS * tmpCSyy;                 // dchi^2/dy/dy
  vtxSeed.cyz += -(tmpCSyz + tmpSLzz);                                                    // dchi^2/dy/dz
  vtxSeed.cy0 += -tmpYXP * (tmpCSyy + tmpSLyz) - tmpZXL * (tmpCSyz + tmpSLzz);            // RHS
  //
  vtxSeed.czz += szzI;                                                                    // dchi^2/dz/dz
  vtxSeed.cz0 += tmpZXL * szzI + tmpYXP * syzI;                                           // RHS
  //
  vtxSeed.tMeanAcc += timeT.getTimeStamp() * trErr2I;
  vtxSeed.tMeanAccErr += trErr2I;
  vtxSeed.addContributor();
}

//___________________________________________________________________
bool Vertexer::solveVertex(Vertexer::VertexSeed &vtxSeed) const
{
  ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> mat;
  mat(0, 0) = vtxSeed.cxx;
  mat(0, 1) = vtxSeed.cxy;
  mat(0, 2) = vtxSeed.cxz;
  mat(1, 1) = vtxSeed.cyy;
  mat(1, 2) = vtxSeed.cyz;
  mat(2, 2) = vtxSeed.czz;
  if (!mat.InvertFast()) {
    printf("Failed to invert matrix\n");
    std::cout << mat << "\n";
    return false;
  }
  ROOT::Math::SVector<double, 3> rhs(vtxSeed.cx0, vtxSeed.cy0, vtxSeed.cz0);
  auto sol = mat * rhs;
  vtxSeed.setXYZ(sol(0), sol(1), sol(2));
  if (vtxSeed.fillErrors) {
    vtxSeed.setCov(mat(0, 0), mat(1, 0), mat(1, 1), mat(2, 0), mat(2, 1), mat(2, 2));
  }
  if (vtxSeed.tMeanAccErr>0.) {
    auto err2 = 1./vtxSeed.tMeanAccErr;
    vtxSeed.setTimeStamp({ float(vtxSeed.tMeanAcc*err2), float(std::sqrt(err2)) });
  }
  
  vtxSeed.setChi2(2.f * (vtxSeed.getNContributors() - vtxSeed.wghSum) / vtxSeed.scaleSig2ITuk2I); // calculate chi^2
  auto newScale = vtxSeed.wghChi2 / vtxSeed.wghSum;
  LOG(INFO) << "Solve: wghChi2=" << vtxSeed.wghChi2 << " wghSum=" << vtxSeed.wghSum << " -> scale= " << newScale << " old scale " << vtxSeed.scaleSigma2Prev;
  vtxSeed.print();
  vtxSeed.setScale(newScale < mMinScale2 ? mMinScale2 : newScale, mTukey2I);
  return true;
}

//___________________________________________________________________
bool Vertexer::stopIterations(VertexSeed &vtxSeed, Vertex &vtx) const
{
  // decide if new iteration should be done, prepare next one if needed
  // if scaleSigma2 reached its lower limit stop
  bool stop = false;
  while (1) {
    if (vtxSeed.scaleSigma2 <= mMinScale2 + kAlmost0F) {
      stop = true;
      LOG(INFO) << "stop on simga :" << vtxSeed.scaleSigma2 << " prev: " << vtxSeed.scaleSigma2Prev;
      //      break;
    }
    auto dchi = (vtx.getChi2() - vtxSeed.getChi2()) / vtxSeed.getChi2();
    if (dchi > 0 && dchi < mMaxChi2Change) {
      stop = true;
      LOG(INFO) << "stop on chi2 :";
      //      break;
    }
    auto dx = vtxSeed.getX() - vtx.getX(), dy = vtxSeed.getY() - vtx.getY(), dz = vtxSeed.getZ() - vtx.getZ();
    auto dst = std::sqrt(dx * dx + dy * dy + dz * dz);

    LOG(INFO) << "dChi:" << vtx.getChi2() << "->" << vtxSeed.getChi2() << " :-> " << dchi;
    LOG(INFO) << "dx: " << dx << " dy: " << dy << " dz: " << dz << " -> " << dst;

    break;
  }

  vtx = vtxSeed;
  /*
  auto &vTStamp = vtx.getTimeStamp();
  if (vTStamp.getTimeStampError() > 0) {
    vTStamp.setTimeStampError( 1. / std::sqrt(vTStamp.getTimeStampError())); // seed accumulates inverse squared error!
  }
  */
  if (!stop) {
    auto scforce = vtxSeed.scaleSigma2Prev * 0.9;
    if (vtxSeed.scaleSigma2 > scforce) {
      if (vtxSeed.nStuck>=mMaxNStuck) {
        auto &timeV = vtxSeed.getTimeStamp();
        float tmean, trms2;
        // seed is stuck between 2 real vertx attractors, push it towards better one (delta was defined as T - V)
        if (vtxSeed.statDZNeg.wsum > vtxSeed.statDZPos.wsum) {
          vtxSeed.setZ(vtxSeed.getZ() + vtxSeed.statDZNeg.getMean());
          vtxSeed.statDTNeg.getMeanRMS2(tmean, trms2);        
          vtxSeed.setTimeStamp( {timeV.getTimeStamp()+tmean, trms2}  );
        }
        else {
          vtxSeed.setZ(vtxSeed.getZ() + vtxSeed.statDZPos.getMean());
          vtxSeed.statDTPos.getMeanRMS2(tmean, trms2);        
          vtxSeed.setTimeStamp( {timeV.getTimeStamp()+tmean, trms2}  );
        }
        vtxSeed.nStuck = 0;
        auto sav = vtxSeed.scaleSigma2Prev;
        vtxSeed.setScale(scforce < mMinScale2 ? mMinScale2 : scforce, mTukey2I);
        vtxSeed.scaleSigma2Prev = sav;
        LOG(INFO) << "Forcing scale from " << vtxSeed.scaleSigma2Prev << " to " << vtxSeed.scaleSigma2 << " Pushed Z from " << vtx.getZ() << " to " << vtxSeed.getZ();
        printf("\n\n");
      }
      else {
        vtxSeed.nStuck++;
      }
    }
    else {
      vtxSeed.nStuck = 0;
    }
  }

  return stop;
}

//___________________________________________________________________
void Vertexer::setConstraint(float x, float y, float sigyy, float sigyz, float sigzz)
{
  // set mean vertex constraint and its errors
  mXYConstraint[0] = x;
  mXYConstraint[1] = y;
  double det = sigyy * sigzz - sigyz * sigyz;
  if (det <= kAlmost0D || sigyy < kAlmost0D || sigzz < kAlmost0D) {
    printf("Error: wrong settings for vertex constraint: %e %e | %e %e %e\n", x, y, sigyy, sigyz, sigzz);
    exit(1);
  }
  mXYConstraintInvErr[0] = sigzz / det;
  mXYConstraintInvErr[2] = sigyy / det;
  mXYConstraintInvErr[1] = -sigyz / det;
}

//______________________________________________
float Vertexer::getTukey() const
{
  // convert 1/tukey^2 to tukey
  return sqrtf(1. / mTukey2I);
}

//___________________________________________________________________
TimeEst Vertexer::timeEstimate(gsl::span<Vertexer::Track> tracks, gsl::span<int> idxsort) const
{
  StatAccumulator test;
  for (int i : idxsort) {
    if (tracks[i].canUse()) {
      const auto &timeT = tracks[i].timeEst;
      auto trErr2I = 1. / (timeT.getTimeStampError() * timeT.getTimeStampError());
      test.add(timeT.getTimeStamp(), trErr2I);
    }
  }
  if (test.wsum>0) {
    float t = 0., te2 = 0.;
    test.getMeanRMS2(t, te2);
    return {t, te2};
  }
  else {
    return {0.,0.};
  }
  /*
  double t = 0., te = 0.;
  for (int i : idxsort) {
    if (tracks[i].canUse()) {
      const auto &timeT = tracks[i].timeEst;
      auto trErr2I = 1. / (timeT.getTimeStampError() * timeT.getTimeStampError());
      t += timeT.getTimeStamp() * trErr2I;
      te += trErr2I;
    }
  }
  if (te>0) {
    te = 1./te;
    return {float(t*te), float(std::sqrt(te))};
  }
  else {
    return {0.,0.};
  }
  */
}

//___________________________________________________________________
void Vertexer::init()
{
  auto zr = 2*mSeedZRange;
  int nzb = zr/mSeedZBinSize;
  if (nzb*zr<zr-1e-6) {
    nzb++;
  }
  mSeedZBinSizeInv = 1./mSeedZBinSize;
  mSeedZRange = nzb*zr/2.;
  mSeedZHisto.resize(nzb);
}

//___________________________________________________________________
void Vertexer::fillZSeedHisto(gsl::span<Vertexer::Track> tracks, gsl::span<int> idxsort)
{
  memset(mSeedZHisto.data(),0,mSeedZHisto.size()*sizeof(int));
  mZSeedsFilled = 0;
  for (int i : idxsort) {
    if (tracks[i].canUse()) {
      auto bin = getZSeedBin( tracks[i].getZForXY(mXYConstraint[0],mXYConstraint[1]) );
      mSeedZHisto[bin]++;
      mZSeedsFilled++;
    }
  }
}
