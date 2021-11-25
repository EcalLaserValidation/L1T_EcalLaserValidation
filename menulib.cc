/* automatically generated from L1Menu_Collisions2018_v2_1_0 with menu2lib.py */
/* https://gitlab.cern.ch/cms-l1t-utm/scripts */

#include <algorithm>
#include <map>
#include <string>
#include <sstream>
#include <cmath>

#include "menulib.hh"

//
// common functions for algorithm implementations
//
std::pair<double, double>
get_missing_et(L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower,
               const int max_eta,
               const double threshold)
{
  // https://github.com/cms-sw/cmssw/blob/CMSSW_9_0_X/L1Trigger/L1TCalorimeter/src/CaloTools.cc#L13=L15
  const int64_t cos_coeff[72] = {1023, 1019, 1007, 988, 961, 927, 886, 838, 784, 723, 658, 587, 512, 432, 350, 265, 178, 89, 0, -89, -178, -265, -350, -432, -512, -587, -658, -723, -784, -838, -886, -927, -961, -988, -1007, -1019, -1023, -1019, -1007, -988, -961, -927, -886, -838, -784, -723, -658, -587, -512, -432, -350, -265, -178, -89, 0, 89, 178, 265, 350, 432, 511, 587, 658, 723, 784, 838, 886, 927, 961, 988, 1007, 1019};

  const int64_t sin_coeff[72] = {0, 89, 178, 265, 350, 432, 512, 587, 658, 723, 784, 838, 886, 927, 961, 988, 1007, 1019, 1023, 1019, 1007, 988, 961, 927, 886, 838, 784, 723, 658, 587, 512, 432, 350, 265, 178, 89, 0, -89, -178, -265, -350, -432, -512, -587, -658, -723, -784, -838, -886, -927, -961, -988, -1007, -1019, -1023, -1019, -1007, -988, -961, -927, -886, -838, -784, -723, -658, -587, -512, -432, -350, -265, -178, -89};

  if (not calo_tower) return std::make_pair(-1., -9999.);

  double ex = 0.;
  double ey = 0.;

  for (int ii = 0; ii < calo_tower->nTower; ii++)
  {
    if (abs(calo_tower->ieta.at(ii)) <= max_eta)
    {
      const double et = calo_tower->iet.at(ii) * 0.5;
      if (et > threshold)
      {
        const int index = calo_tower->iphi.at(ii) - 1;
        ex += (et*cos_coeff[index]/1024.);
        ey += (et*sin_coeff[index]/1024.);
      }
    }
  }

  return std::make_pair(sqrt(ex*ex + ey*ey), atan2(-ey, -ex));
}


double
get_total_ht(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade,
             const int max_eta,
             const double threshold)
{
  double sum = 0.;

  for (int ii = 0; ii < upgrade->nJets; ii++)
  {
    if (upgrade->jetBx.at(ii) != 0) continue;

    if (abs(upgrade->jetIEta.at(ii)) <= 2*max_eta)
    {
      const double et = upgrade->jetEt.at(ii);
      if (et > threshold)
      {
        sum += et;
      }
    }
  }

  return sum;
} 


double
get_transverse_mass(L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade,
                    const double threshold_eg,
                    const double threshold_met)
{
  double mt = -1.;

  if (upgrade->nEGs == 0) return mt;

  // leading-eg
  int id_leading_eg = -1;
  for (int ii = 0; ii < upgrade->nEGs; ii++)
  {
    if (upgrade->egBx.at(ii) != 0) continue;
    if (id_leading_eg < 0)
    {
      id_leading_eg = ii;
      break;
    }
  }

  if (id_leading_eg < 0) return mt;

  const double eg_et = upgrade->egEt.at(id_leading_eg);
  const double eg_phi = upgrade->egPhi.at(id_leading_eg);

  if (eg_et < threshold_eg) return mt;


  // missing-Et
  int id_missing_et = -1;
  for (int ii = 0; ii < upgrade->nSums; ii++)
  {
    if (upgrade->sumBx.at(ii) != 0) continue;
    if (upgrade->sumType.at(ii) == L1Analysis::kMissingEt)
    {
      id_missing_et = ii;
      break;
    }
  }

  if (id_missing_et < 0) return mt;

  const double met_et = upgrade->sumEt.at(id_missing_et);
  const double met_phi = upgrade->sumPhi.at(id_missing_et);

  if (met_et < threshold_met) return mt;


  // mt
  double delta_phi = eg_phi - met_phi;
  while (delta_phi >= M_PI) delta_phi -= 2.*M_PI;
  while (delta_phi < -M_PI) delta_phi += 2.*M_PI;

  mt = sqrt(2.*eg_et*met_et*(1. - cos(delta_phi)));
  return mt;
}


// utility methods
void
getCombination(int N,
               int K,
               std::vector<std::vector<int> >& combination)
{
  std::string bitmask(K, 1);
  bitmask.resize(N, 0);

  do
  {
    std::vector<int> set;
    for (int ii = 0; ii < N; ++ii)
    {
      if (bitmask[ii]) set.push_back(ii);
    }
    combination.push_back(set);
  }
  while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}


void
getPermutation(int N,
               std::vector<std::vector<int> >& permutation)
{
  std::vector<int> indicies(N);
  for (int ii = 0; ii < N; ii++) indicies.at(ii) = ii;

  do
  {
    std::vector<int> set;
    for (int ii = 0; ii < N; ++ii)
    {
      set.push_back(indicies.at(ii));
    }
    permutation.push_back(set);
  }
  while (std::next_permutation(indicies.begin(), indicies.end()));
}




//
// NB: tmEventSetup.XxxWithOverlapRemoval was removed between utm-overlapRemoval-xsd330 and utm_0.6.5
//
// generate conditions
    





bool
CaloCaloCorrelation_i119
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(ii) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
              
                                      // JET34: ET >= 68 at BX = 0
      if (not (data->jetIEt.at(jj) >= 68)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(jj)) and (data->jetIEta.at(jj) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_JET[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_JET[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_i120
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(ii) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
              
                                      // JET34: ET >= 68 at BX = 0
      if (not (data->jetIEt.at(jj) >= 68)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(jj)) and (data->jetIEta.at(jj) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_JET[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_JET[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_i121
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(ii) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
              
                                      // JET34: ET >= 68 at BX = 0
      if (not (data->jetIEt.at(jj) >= 68)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(jj)) and (data->jetIEta.at(jj) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_JET[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_JET[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_i132
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 200)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 52));
            
          if (not etaWindow1) continue;
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 200)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 52));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloCaloCorrelation_i133
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET112: ET >= 224 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 224)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 52));
            
          if (not etaWindow1) continue;
                                // JET112: ET >= 224 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 224)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 52));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloCaloCorrelation_i135
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 64)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 52));
            
          if (not etaWindow1) continue;
                                // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 64)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 52));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloCaloCorrelation_i137
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 52));
            
          if (not etaWindow1) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 52));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloCaloCorrelation_i152
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(ii) >= 44)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
              
                                      // TAU70: ET >= 140 at BX = 0
      if (not (data->tauIEt.at(jj) >= 140)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_i153
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(ii) >= 44)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
              
                                      // TAU26: ET >= 52 at BX = 0
      if (not (data->tauIEt.at(jj) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(jj)) & 1)) continue;

          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_i154
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(ii) >= 48)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
              
                                      // TAU27: ET >= 54 at BX = 0
      if (not (data->tauIEt.at(jj) >= 54)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(jj)) and (data->tauIEta.at(jj) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(jj)) & 1)) continue;

          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];
  
    int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_EG_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.09 * POW10[6]);
  maximum = (long long)(139.24 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
CaloCaloCorrelation_i255
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 48));
            
          if (not etaWindow1) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.6 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloMuonCorrelation_i134
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET32: ET >= 64 at BX = 0
      if (not (data->jetIEt.at(ii) >= 64)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 52));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(jj) >= 21)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000624999999997
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(jj)) and (data->muonIEtaAtVtx.at(jj) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i136
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(ii) >= 80)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 52));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(jj) >= 25)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000624999999997
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(jj)) and (data->muonIEtaAtVtx.at(jj) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i254
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(ii) >= 80)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 48));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(jj) >= 25)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000624999999997
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(jj)) and (data->muonIEtaAtVtx.at(jj) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i272
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET16: ET >= 32 at BX = 0
      if (not (data->jetIEt.at(ii) >= 32)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i273
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(ii) >= 120)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i274
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(ii) >= 240)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i275
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(ii) >= 70)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i276
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(ii) >= 160)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.40
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.161 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i277
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(ii) >= 240)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.641 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i303
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(ii) >= 180)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(jj) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.641 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i305
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
          
        bool etaWindow1;
                              // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(ii) >= 180)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(ii)) and (data->jetIEta.at(ii) <= 57));
            
          if (not etaWindow1) continue;
    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
                                      // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(jj) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(jj)) & 1)) continue;

          
          long long minimum;
  long long maximum;
    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->jetIEta.at(ii);
      if (iEta < 0) iEta += 256;
    iEta = LUT_ETA_JET2MU[iEta];
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(jj));
      unsigned int deltaEta = LUT_DETA_JET_MU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
      iPhi = LUT_PHI_JET2MU[iPhi];
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(jj));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(0.641 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i144
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i145
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i173
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i206
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i214
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i246
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i247
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i268
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i269
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG25: ET >= 50 at BX = 0
      if (not (data->egIEt.at(idx) >= 50)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i278
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i279
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i280
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG27: ET >= 54 at BX = 0
      if (not (data->egIEt.at(idx) >= 54)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG14: ET >= 28 at BX = 0
      if (not (data->egIEt.at(idx) >= 28)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i61
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i62
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG25: ET >= 50 at BX = 0
      if (not (data->egIEt.at(idx) >= 50)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i63
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG25: ET >= 50 at BX = 0
      if (not (data->egIEt.at(idx) >= 50)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG14: ET >= 28 at BX = 0
      if (not (data->egIEt.at(idx) >= 28)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i64
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG22: ET >= 44 at BX = 0
      if (not (data->egIEt.at(idx) >= 44)) continue;

                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i65
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i146
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx) >= 60)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i148
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET110: ET >= 220 at BX = 0
      if (not (data->jetIEt.at(idx) >= 220)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i150
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET115: ET >= 230 at BX = 0
      if (not (data->jetIEt.at(idx) >= 230)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i258
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx) >= 60)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i260
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET70: ET >= 140 at BX = 0
      if (not (data->jetIEt.at(idx) >= 140)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i262
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET65: ET >= 130 at BX = 0
      if (not (data->jetIEt.at(idx) >= 130)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i309
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx) >= 90)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i312
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i313
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET150: ET >= 300 at BX = 0
      if (not (data->jetIEt.at(idx) >= 300)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET150: ET >= 300 at BX = 0
      if (not (data->jetIEt.at(idx) >= 300)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i87
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i88
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i90
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET85: ET >= 170 at BX = 0
      if (not (data->jetIEt.at(idx) >= 170)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i165
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i169
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx) >= 31)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i170
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx) >= 31)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i176
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i177
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i192
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU9: ET >= 19 at BX = 0
      if (not (data->muonIEt.at(idx) >= 19)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU9: ET >= 19 at BX = 0
      if (not (data->muonIEt.at(idx) >= 19)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i193
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i194
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i195
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i196
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 184));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 184));
            
          if (not etaWindow1) continue;
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i197
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i198
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      // charge correlation
      bool equal = true;
      bool invalid = false;
      for (size_t mm = 0; mm < 2 -1; mm++)
      {
        int idx0 = candidates.at(set.at(indicies.at(mm)));
        int idx1 = candidates.at(set.at(indicies.at(mm+1)));
        if ((data->muonChg.at(idx0) == 0) or (data->muonChg.at(idx1) == 0))
        {
          invalid = true;
          break;
        }
        if (data->muonChg.at(idx0) != data->muonChg.at(idx1))
        {
          equal = false;
          break;
        }
      }
      if (invalid) continue;

      // charge correlation: "os"
      if (equal) continue;

      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i243
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU8: ET >= 17 at BX = 0
      if (not (data->muonIEt.at(idx) >= 17)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU8: ET >= 17 at BX = 0
      if (not (data->muonIEt.at(idx) >= 17)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i256
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i28
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i29
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i292
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i30
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx) >= 31)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i31
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125000000003 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125000000003 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i356
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU30: ET >= 60 at BX = 0
      if (not (data->tauIEt.at(idx) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU30: ET >= 60 at BX = 0
      if (not (data->tauIEt.at(idx) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i357
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU28: ET >= 56 at BX = 0
      if (not (data->tauIEt.at(idx) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU28: ET >= 56 at BX = 0
      if (not (data->tauIEt.at(idx) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i69
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU70: ET >= 140 at BX = 0
      if (not (data->tauIEt.at(idx) >= 140)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU70: ET >= 140 at BX = 0
      if (not (data->tauIEt.at(idx) >= 140)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i70
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU32: ET >= 64 at BX = 0
      if (not (data->tauIEt.at(idx) >= 64)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU32: ET >= 64 at BX = 0
      if (not (data->tauIEt.at(idx) >= 64)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i71
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU34: ET >= 68 at BX = 0
      if (not (data->tauIEt.at(idx) >= 68)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU34: ET >= 68 at BX = 0
      if (not (data->tauIEt.at(idx) >= 68)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i72
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU36: ET >= 72 at BX = 0
      if (not (data->tauIEt.at(idx) >= 72)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // TAU36: ET >= 72 at BX = 0
      if (not (data->tauIEt.at(idx) >= 72)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

                          

  



bool
InvariantMassOvRm_i264
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
        // remove overlap -- reference: TAU40
  std::vector<int> reference;
    size_t nref = 0;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nref++;
          if (nref > 12) break;
          
                              // TAU40: ET >= 80 at BX = 0
      if (not (data->tauIEt.at(ii) >= 80)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(ii)) & 1)) continue;

          
    reference.push_back(ii);
  }
  if (not reference.size()) return false;

    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                     int iEta = -9999999; unsigned int deltaIEta = 9999999;
                          
  // remove overlap -- target: JET80
       // 0.00 <= DeltaR <= 0.20
  long long minDeltaR2 = std::numeric_limits<long long>::max();
  const long long cutDeltaR2Min = (long long)(0.0 * POW10[6]);
  const long long cutDeltaR2Max = (long long)(0.041 * POW10[6]);
         
  // compute minimum distance to reference objects
  for (size_t _jj = 0; _jj < reference.size(); _jj++)
  {
    const int index = reference.at(_jj);
                iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_TAU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
      if (deltaR2 < minDeltaR2) minDeltaR2 = deltaR2;
                  }

  // skip if needed
      if ((cutDeltaR2Min <= minDeltaR2) and (minDeltaR2 <= cutDeltaR2Max)) continue;

        
        candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 160)) continue;

          
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
            // 420.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(88200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                          

  



bool
InvariantMassOvRm_i316
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
        // remove overlap -- reference: TAU45
  std::vector<int> reference;
    size_t nref = 0;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nref++;
          if (nref > 12) break;
          
                              // TAU45: ET >= 90 at BX = 0
      if (not (data->tauIEt.at(ii) >= 90)) continue;

                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(ii)) & 1)) continue;

          
    reference.push_back(ii);
  }
  if (not reference.size()) return false;

    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                     int iEta = -9999999; unsigned int deltaIEta = 9999999;
                          
  // remove overlap -- target: JET35
       // 0.00 <= DeltaR <= 0.20
  long long minDeltaR2 = std::numeric_limits<long long>::max();
  const long long cutDeltaR2Min = (long long)(0.0 * POW10[6]);
  const long long cutDeltaR2Max = (long long)(0.041 * POW10[6]);
         
  // compute minimum distance to reference objects
  for (size_t _jj = 0; _jj < reference.size(); _jj++)
  {
    const int index = reference.at(_jj);
                iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_TAU[deltaIEta];
  
    int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_JET_TAU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
      if (deltaR2 < minDeltaR2) minDeltaR2 = deltaR2;
                  }

  // skip if needed
      if ((cutDeltaR2Min <= minDeltaR2) and (minDeltaR2 <= cutDeltaR2Max)) continue;

        
        candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 70)) continue;

          
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 70)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
            // 450.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(101250.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i113
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

          // 300.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(45000.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i114
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 250.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(31250.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i115
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 330.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(54450.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i116
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 360.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(64800.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i147
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

          
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i149
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 70)) continue;

          
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 70)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i151
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

          
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i171
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 31)) continue;

          
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 15)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 1.0 <= mass <= 151982.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.5 * POW10[6]);
  maximum = (long long)(11549264162.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i175
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

          
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 1.0 <= mass <= 151982.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.5 * POW10[6]);
  maximum = (long long)(11549264162.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i182
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 184));
            
          if (not etaWindow1) continue;
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 184));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 7.0 <= mass <= 151982.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(24.5 * POW10[6]);
  maximum = (long long)(11549264162.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i187
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 6)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 5.0 <= mass <= 17.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(12.5 * POW10[6]);
  maximum = (long long)(144.5 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i190
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.0 <= mass <= 9.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(40.5 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i200
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000624999999997
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 211));
            
          if (not etaWindow1) continue;
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000624999999997
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 8.0 <= mass <= 14.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(32.0 * POW10[6]);
  maximum = (long long)(98.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i201
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // EG3: ET >= 6 at BX = 0
      if (not (data->egIEt.at(idx0) >= 6)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx0)) and (data->egIEta.at(idx0) <= 48));
            
          if (not etaWindow1) continue;
                                // EG3: ET >= 6 at BX = 0
      if (not (data->egIEt.at(idx1) >= 6)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx1)) and (data->egIEta.at(idx1) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 20.0
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
  
    int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_EG_EG[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_EG_EG[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(idx0)];
  long long pt1 = LUT_EG_ET[data->egIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[5]);
  maximum = (long long)(200.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i202
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000624999999997
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 211));
            
          if (not etaWindow1) continue;
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.3000625 <= eta <= 2.3000624999999997
              etaWindow1 = ((-211 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 211));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 14.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(98.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

          const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i203
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // EG7p5: ET >= 15 at BX = 0
      if (not (data->egIEt.at(idx0) >= 15)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx0)) and (data->egIEta.at(idx0) <= 48));
            
          if (not etaWindow1) continue;
                                // EG7p5: ET >= 15 at BX = 0
      if (not (data->egIEt.at(idx1) >= 15)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx1)) and (data->egIEta.at(idx1) <= 48));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 20.0
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
  
    int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_EG_EG[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_EG_EG[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(idx0)];
  long long pt1 = LUT_EG_ET[data->egIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[5]);
  maximum = (long long)(200.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i205
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 11)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 6)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 5.0 <= mass <= 17.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(12.5 * POW10[6]);
  maximum = (long long)(144.5 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i265
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

          // 150.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(11250.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i266
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 57));
            
          if (not etaWindow1) continue;
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 57));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    minimum = (long long)(0.0 * POW10[3]);
  maximum = (long long)(1.5 * POW10[3]);
  if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;

          // 200.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(20000.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i270
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 160)) continue;

          
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 60)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 420.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(88200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i310
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 90)) continue;

          
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i326
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -2.6969999999999996 <= eta <= 2.6969999999999996
              etaWindow1 = ((-62 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 61));
            
          if (not etaWindow1) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.6969999999999996 <= eta <= 2.6969999999999996
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i327
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 120)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i328
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= -70));
            
                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i329
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.6969999999999996 <= eta <= 2.6969999999999996
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i330
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= -70));
            
                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -2.6969999999999996 <= eta <= 2.6969999999999996
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i331
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 80)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= -70));
            
                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 80)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= -70));
            
                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i332
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 90)) continue;

                        // -2.6969999999999996 <= eta <= 2.6969999999999996
              etaWindow1 = ((-62 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 61));
            
          if (not etaWindow1) continue;
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -2.6969999999999996 <= eta <= 2.6969999999999996
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i333
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= -70));
            
                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i334
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 120)) continue;

          
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -2.6969999999999996 <= eta <= 2.6969999999999996
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i335
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 90)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= -70));
            
                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -2.6969999999999996 <= eta <= 2.6969999999999996
              etaWindow1 = ((-62 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 61));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i336
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1, etaWindow2;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx0) >= 90)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= -70));
            
                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx0)) and (data->jetIEta.at(idx0) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx1) >= 90)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= -70));
            
                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow2 = ((69 <= data->jetIEta.at(idx1)) and (data->jetIEta.at(idx1) <= 114));
            
          if (not (etaWindow1 or etaWindow2)) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
    int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
  long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
  long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(192200.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i338
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 184));
            
          if (not etaWindow1) continue;
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 184));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 7.0 <= mass <= 18.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
  long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
  long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(24.5 * POW10[6]);
  maximum = (long long)(162.0 * POW10[6]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i354
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
        if (nobj0 > 12) break;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // TAU28: ET >= 56 at BX = 0
      if (not (data->tauIEt.at(idx0) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx0)) and (data->tauIEta.at(idx0) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx0)) & 1)) continue;

          if (not etaWindow1) continue;
                                // TAU28: ET >= 56 at BX = 0
      if (not (data->tauIEt.at(idx1) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx1)) and (data->tauIEta.at(idx1) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx1)) & 1)) continue;

          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 80.0
      iEta = data->tauIEta.at(idx0);
    deltaIEta = abs(iEta - data->tauIEta.at(idx1));
  
    int iPhi = data->tauIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_TAU_TAU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_TAU_TAU[deltaIPhi];
  long long pt0 = LUT_TAU_ET[data->tauIEt.at(idx0)];
  long long pt1 = LUT_TAU_ET[data->tauIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[5]);
  maximum = (long long)(3200.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i355
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
        if (nobj0 > 12) break;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // TAU28: ET >= 56 at BX = 0
      if (not (data->tauIEt.at(idx0) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx0)) and (data->tauIEta.at(idx0) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx0)) & 1)) continue;

          if (not etaWindow1) continue;
                                // TAU28: ET >= 56 at BX = 0
      if (not (data->tauIEt.at(idx1) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx1)) and (data->tauIEta.at(idx1) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx1)) & 1)) continue;

          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 90.0
      iEta = data->tauIEta.at(idx0);
    deltaIEta = abs(iEta - data->tauIEta.at(idx1));
  
    int iPhi = data->tauIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_TAU_TAU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_TAU_TAU[deltaIPhi];
  long long pt0 = LUT_TAU_ET[data->tauIEt.at(idx0)];
  long long pt1 = LUT_TAU_ET[data->tauIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[5]);
  maximum = (long long)(4050.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i358
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
        if (nobj0 > 12) break;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // TAU30: ET >= 60 at BX = 0
      if (not (data->tauIEt.at(idx0) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx0)) and (data->tauIEta.at(idx0) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx0)) & 1)) continue;

          if (not etaWindow1) continue;
                                // TAU30: ET >= 60 at BX = 0
      if (not (data->tauIEt.at(idx1) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx1)) and (data->tauIEta.at(idx1) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx1)) & 1)) continue;

          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 90.0
      iEta = data->tauIEta.at(idx0);
    deltaIEta = abs(iEta - data->tauIEta.at(idx1));
  
    int iPhi = data->tauIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_TAU_TAU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_TAU_TAU[deltaIPhi];
  long long pt0 = LUT_TAU_ET[data->tauIEt.at(idx0)];
  long long pt1 = LUT_TAU_ET[data->tauIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[5]);
  maximum = (long long)(4050.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i359
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
        if (nobj0 > 12) break;
                   candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // TAU30: ET >= 60 at BX = 0
      if (not (data->tauIEt.at(idx0) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx0)) and (data->tauIEta.at(idx0) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx0)) & 1)) continue;

          if (not etaWindow1) continue;
                                // TAU30: ET >= 60 at BX = 0
      if (not (data->tauIEt.at(idx1) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx1)) and (data->tauIEta.at(idx1) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx1)) & 1)) continue;

          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 80.0
      iEta = data->tauIEta.at(idx0);
    deltaIEta = abs(iEta - data->tauIEta.at(idx1));
  
    int iPhi = data->tauIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long coshDeltaEta = LUT_COSH_DETA_TAU_TAU[deltaIEta];
  long long cosDeltaPhi = LUT_COS_DPHI_TAU_TAU[deltaIPhi];
  long long pt0 = LUT_TAU_ET[data->tauIEt.at(idx0)];
  long long pt1 = LUT_TAU_ET[data->tauIEt.at(idx1)];
  long long mass2 = pt0*pt1*(coshDeltaEta - cosDeltaPhi);
    minimum = (long long)(0.0 * POW10[5]);
  maximum = (long long)(3200.0 * POW10[5]);
  if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
      





bool
MuonMuonCorrelation_i117
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 184));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 184));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i118
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 138));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 138));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i181
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 138));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 138));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i183
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 9)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.20
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.441 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i184
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -1.4083124999999999 <= eta <= 1.4083124999999999
              etaWindow1 = ((-129 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 129));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -1.4083124999999999 <= eta <= 1.4083124999999999
              etaWindow1 = ((-129 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 129));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i185
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU4p5: ET >= 10 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 10)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.20
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.441 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i235
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == -1)) continue;
    nobj0++;
      
        const int idx0 = ii;
    bool etaWindow1;bool phiWindow1;
                              // MU3-1: ET >= 7 at BX = -1
      if (not (data->muonIEt.at(idx0) >= 7)) continue;

                        // -1.2016874999999998 <= eta <= 1.2016875
              etaWindow1 = ((-110 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 110));
            
                        // 0.5235987755982988 <= phi <= 2.6179938779914944
              phiWindow1 = ((48 <= data->muonIPhiAtVtx.at(idx0)) and (data->muonIPhiAtVtx.at(idx0) <= 239));
            
                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

          if (not etaWindow1) continue;
	  if (not phiWindow1) continue;

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
        
            const int idx1 = jj;
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 7)) continue;

                        // -1.2016874999999998 <= eta <= 1.2016875
              etaWindow1 = ((-110 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 110));
            
                        // 3.665191429188092 <= phi <= 5.759586531581287
              phiWindow1 = ((336 <= data->muonIPhiAtVtx.at(idx1)) and (data->muonIPhiAtVtx.at(idx1) <= 527));
            
                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

          if (not etaWindow1) continue;
	  if (not phiWindow1) continue;
          long long minimum;
  long long maximum;

  
        // 2.618 <= DeltaPhi <= 3.142
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    minimum = (long long)(2.618 * POW10[3]);
  maximum = (long long)(3.142 * POW10[3]);
  if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;

    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
MuonMuonCorrelation_i267
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      bool etaWindow1;
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx0)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx0)) and (data->muonIEtaAtVtx.at(idx0) <= 184));
            
          if (not etaWindow1) continue;
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx1)) & 1)) continue;

                        // -2.0064375 <= eta <= 2.0064374999999997
              etaWindow1 = ((-184 <= data->muonIEtaAtVtx.at(idx1)) and (data->muonIEtaAtVtx.at(idx1) <= 184));
            
          if (not etaWindow1) continue;
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        const std::string OS = "os";
    const std::string SS = "ls";
    if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if ("os" == OS)
    {
      if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    }
    else if ("os" == SS)
    {
      if (not (data->muonChg.at(idx0) == data->muonChg.at(idx1))) continue;
    }
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(1.961 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i302
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.60
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(2.561 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i304
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 2, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(2, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx0) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx0)) & 1)) continue;

          
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx1) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx1)) & 1)) continue;

          
          long long minimum;
  long long maximum;

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.60
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
      unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
  long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    minimum = (long long)(0.0 * POW10[6]);
  maximum = (long long)(2.561 * POW10[6]);
  if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      


bool
QuadJET_i129
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET70: ET >= 140 at BX = 0
      if (not (data->jetIEt.at(idx) >= 140)) continue;

                        // -2.3924999999999996 <= eta <= 2.3924999999999996
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET55: ET >= 110 at BX = 0
      if (not (data->jetIEt.at(idx) >= 110)) continue;

                        // -2.3924999999999996 <= eta <= 2.3924999999999996
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.3924999999999996 <= eta <= 2.3924999999999996
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -2.3924999999999996 <= eta <= 2.3924999999999996
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i130
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET70: ET >= 140 at BX = 0
      if (not (data->jetIEt.at(idx) >= 140)) continue;

                        // -2.3924999999999996 <= eta <= 2.3924999999999996
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET55: ET >= 110 at BX = 0
      if (not (data->jetIEt.at(idx) >= 110)) continue;

                        // -2.3924999999999996 <= eta <= 2.3924999999999996
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.3924999999999996 <= eta <= 2.3924999999999996
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.3924999999999996 <= eta <= 2.3924999999999996
              etaWindow1 = ((-55 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 54));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i131
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx) >= 90)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 52));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET40: ET >= 80 at BX = 0
      if (not (data->jetIEt.at(idx) >= 80)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 52));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i142
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET36: ET >= 72 at BX = 0
      if (not (data->jetIEt.at(idx) >= 72)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i239
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET50: ET >= 100 at BX = 0
      if (not (data->jetIEt.at(idx) >= 100)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 52));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET45: ET >= 90 at BX = 0
      if (not (data->jetIEt.at(idx) >= 90)) continue;

                        // -2.3055 <= eta <= 2.3055
              etaWindow1 = ((-53 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 52));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i306
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET95: ET >= 190 at BX = 0
      if (not (data->jetIEt.at(idx) >= 190)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET65: ET >= 130 at BX = 0
      if (not (data->jetIEt.at(idx) >= 130)) continue;

          
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i91
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(3)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadMU_i284
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(3)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadMU_i285
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(3)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadMU_i36
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 4) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 4, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(4, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(3)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i122
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i124
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i125
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i126
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG23: ET >= 46 at BX = 0
      if (not (data->egIEt.at(idx) >= 46)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i127
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i128
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG23: ET >= 46 at BX = 0
      if (not (data->egIEt.at(idx) >= 46)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i178
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG9: ET >= 18 at BX = 0
      if (not (data->egIEt.at(idx) >= 18)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i248
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i317
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i318
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i319
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i320
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i321
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i337
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG60: ET >= 120 at BX = 0
      if (not (data->egIEt.at(idx) >= 120)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i339
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i340
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i341
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i342
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i343
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -1.5225 <= eta <= 1.5225
              etaWindow1 = ((-35 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 34));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i344
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i345
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i346
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // 2.5229999999999997 <= eta <= 5.0
              etaWindow1 = ((58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i347
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -5.0 <= eta <= -2.5229999999999997
              etaWindow1 = ((-115 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= -59));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i348
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // 2.5229999999999997 <= eta <= 5.0
              etaWindow1 = ((58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 114));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i349
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -5.0 <= eta <= -2.5229999999999997
              etaWindow1 = ((-115 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= -59));
            
                        // isolation : 0xc
      if (not ((12 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i350
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // 2.5229999999999997 <= eta <= 5.0
              etaWindow1 = ((58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 114));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i351
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -5.0 <= eta <= -2.5229999999999997
              etaWindow1 = ((-115 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= -59));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i352
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG20: ET >= 40 at BX = 0
      if (not (data->egIEt.at(idx) >= 40)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i37
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i38
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG10: ET >= 20 at BX = 0
      if (not (data->egIEt.at(idx) >= 20)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i39
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i40
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i41
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG45: ET >= 90 at BX = 0
      if (not (data->egIEt.at(idx) >= 90)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i42
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG50: ET >= 100 at BX = 0
      if (not (data->egIEt.at(idx) >= 100)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i43
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG34: ET >= 68 at BX = 0
      if (not (data->egIEt.at(idx) >= 68)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i44
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG36: ET >= 72 at BX = 0
      if (not (data->egIEt.at(idx) >= 72)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i45
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG38: ET >= 76 at BX = 0
      if (not (data->egIEt.at(idx) >= 76)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i46
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG40: ET >= 80 at BX = 0
      if (not (data->egIEt.at(idx) >= 80)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i47
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG42: ET >= 84 at BX = 0
      if (not (data->egIEt.at(idx) >= 84)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i48
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG24: ET >= 48 at BX = 0
      if (not (data->egIEt.at(idx) >= 48)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i49
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i50
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i51
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i52
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(idx) >= 64)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 48));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i53
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG26: ET >= 52 at BX = 0
      if (not (data->egIEt.at(idx) >= 52)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i54
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG28: ET >= 56 at BX = 0
      if (not (data->egIEt.at(idx) >= 56)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i55
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG30: ET >= 60 at BX = 0
      if (not (data->egIEt.at(idx) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i56
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(idx) >= 64)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i57
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG34: ET >= 68 at BX = 0
      if (not (data->egIEt.at(idx) >= 68)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      
bool
SingleETMHF_i101
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF100: ET >= 200 at BX = 0
      if (not (data->sumIEt.at(ii) >= 200)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i102
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF110: ET >= 220 at BX = 0
      if (not (data->sumIEt.at(ii) >= 220)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i103
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF120: ET >= 240 at BX = 0
      if (not (data->sumIEt.at(ii) >= 240)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i104
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF130: ET >= 260 at BX = 0
      if (not (data->sumIEt.at(ii) >= 260)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i139
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF90: ET >= 180 at BX = 0
      if (not (data->sumIEt.at(ii) >= 180)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i140
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF80: ET >= 160 at BX = 0
      if (not (data->sumIEt.at(ii) >= 160)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i166
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF50: ET >= 100 at BX = 0
      if (not (data->sumIEt.at(ii) >= 100)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i240
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF60: ET >= 120 at BX = 0
      if (not (data->sumIEt.at(ii) >= 120)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i251
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF40: ET >= 80 at BX = 0
      if (not (data->sumIEt.at(ii) >= 80)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i263
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF140: ET >= 280 at BX = 0
      if (not (data->sumIEt.at(ii) >= 280)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i311
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF150: ET >= 300 at BX = 0
      if (not (data->sumIEt.at(ii) >= 300)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i353
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETMHF70: ET >= 140 at BX = 0
      if (not (data->sumIEt.at(ii) >= 140)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETM_i314
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM120: ET >= 240 at BX = 0
      if (not (data->sumIEt.at(ii) >= 240)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETM_i315
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETM150: ET >= 300 at BX = 0
      if (not (data->sumIEt.at(ii) >= 300)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETT_i322
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT1200: ET >= 2400 at BX = 0
      if (not (data->sumIEt.at(ii) >= 2400)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETT_i323
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT1600: ET >= 3200 at BX = 0
      if (not (data->sumIEt.at(ii) >= 3200)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETT_i324
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // ETT2000: ET >= 4000 at BX = 0
      if (not (data->sumIEt.at(ii) >= 4000)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleEXT_i0
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i1
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i105
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i107
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i108
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i2
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i207
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i208
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i209
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i210
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i211
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i212
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i213
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i217
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i218
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i219
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i220
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i221
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i222
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i223
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i224
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i225
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i226
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i227
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i228
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i229
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i230
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i231
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i232
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i233
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i234
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i3
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i4
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleEXT_i5
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // for now return false always
  // could check decision available in ugt data
  bool pass = false;
  return pass;
}

      
bool
SingleHTT_i100
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT450: ET >= 900 at BX = 0
      if (not (data->sumIEt.at(ii) >= 900)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i123
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT100: ET >= 200 at BX = 0
      if (not (data->sumIEt.at(ii) >= 200)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i163
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT240: ET >= 480 at BX = 0
      if (not (data->sumIEt.at(ii) >= 480)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i164
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT250: ET >= 500 at BX = 0
      if (not (data->sumIEt.at(ii) >= 500)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i168
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT220: ET >= 440 at BX = 0
      if (not (data->sumIEt.at(ii) >= 440)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i174
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT260: ET >= 520 at BX = 0
      if (not (data->sumIEt.at(ii) >= 520)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i179
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT60: ET >= 120 at BX = 0
      if (not (data->sumIEt.at(ii) >= 120)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i241
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT300: ET >= 600 at BX = 0
      if (not (data->sumIEt.at(ii) >= 600)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i242
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT340: ET >= 680 at BX = 0
      if (not (data->sumIEt.at(ii) >= 680)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i92
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT120: ET >= 240 at BX = 0
      if (not (data->sumIEt.at(ii) >= 240)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i93
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT160: ET >= 320 at BX = 0
      if (not (data->sumIEt.at(ii) >= 320)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i94
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT200: ET >= 400 at BX = 0
      if (not (data->sumIEt.at(ii) >= 400)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i95
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT255: ET >= 510 at BX = 0
      if (not (data->sumIEt.at(ii) >= 510)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i96
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT280: ET >= 560 at BX = 0
      if (not (data->sumIEt.at(ii) >= 560)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i97
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT320: ET >= 640 at BX = 0
      if (not (data->sumIEt.at(ii) >= 640)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i98
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT360: ET >= 720 at BX = 0
      if (not (data->sumIEt.at(ii) >= 720)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i99
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // HTT400: ET >= 800 at BX = 0
      if (not (data->sumIEt.at(ii) >= 800)) continue;
      

    pass = true;
    break;
  }

  return pass;
}

      


bool
SingleJET_i110
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i111
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET43: ET >= 86 at BX = 0
      if (not (data->jetIEt.at(idx) >= 86)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i112
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET46: ET >= 92 at BX = 0
      if (not (data->jetIEt.at(idx) >= 92)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i167
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i180
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET30: ET >= 60 at BX = 0
      if (not (data->jetIEt.at(idx) >= 60)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i236
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET8: ET >= 16 at BX = 0
      if (not (data->jetIEt.at(idx) >= 16)) continue;

                        // -2.9579999999999997 <= eta <= -1.392
              etaWindow1 = ((-68 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -33));
            
                        // 1.392 <= eta <= 2.9579999999999997
              etaWindow2 = ((32 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 67));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i237
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET10: ET >= 20 at BX = 0
      if (not (data->jetIEt.at(idx) >= 20)) continue;

                        // -2.9579999999999997 <= eta <= -1.392
              etaWindow1 = ((-68 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -33));
            
                        // 1.392 <= eta <= 2.9579999999999997
              etaWindow2 = ((32 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 67));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i238
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET12: ET >= 24 at BX = 0
      if (not (data->jetIEt.at(idx) >= 24)) continue;

                        // -2.9579999999999997 <= eta <= -1.392
              etaWindow1 = ((-68 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -33));
            
                        // 1.392 <= eta <= 2.9579999999999997
              etaWindow2 = ((32 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 67));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i250
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i253
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET140: ET >= 280 at BX = 0
      if (not (data->jetIEt.at(idx) >= 280)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i287
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i288
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i289
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i290
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET160: ET >= 320 at BX = 0
      if (not (data->jetIEt.at(idx) >= 320)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i291
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET180: ET >= 360 at BX = 0
      if (not (data->jetIEt.at(idx) >= 360)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i307
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i308
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET20: ET >= 40 at BX = 0
      if (not (data->jetIEt.at(idx) >= 40)) continue;

                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i325
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET115: ET >= 230 at BX = 0
      if (not (data->jetIEt.at(idx) >= 230)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i73
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i74
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i75
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i76
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i77
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET180: ET >= 360 at BX = 0
      if (not (data->jetIEt.at(idx) >= 360)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i78
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET200: ET >= 400 at BX = 0
      if (not (data->jetIEt.at(idx) >= 400)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i79
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i80
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET35: ET >= 70 at BX = 0
      if (not (data->jetIEt.at(idx) >= 70)) continue;

                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i81
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i82
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET60: ET >= 120 at BX = 0
      if (not (data->jetIEt.at(idx) >= 120)) continue;

                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i83
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i84
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET90: ET >= 180 at BX = 0
      if (not (data->jetIEt.at(idx) >= 180)) continue;

                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i85
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // -5.0 <= eta <= -3.0014999999999996
              etaWindow1 = ((-115 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= -70));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i86
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET120: ET >= 240 at BX = 0
      if (not (data->jetIEt.at(idx) >= 240)) continue;

                        // 3.0014999999999996 <= eta <= 5.0
              etaWindow1 = ((69 <= data->jetIEta.at(idx)) and (data->jetIEta.at(idx) <= 114));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      
bool
SingleMBT0HFM_i216
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFM0)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // MBT0HFM1: Count >= 1 at BX = 0
      if (not (data->sumIEt.at(ii) >= 1)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleMBT0HFP_i215
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  
  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFP0)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
                        // MBT0HFP1: Count >= 1 at BX = 0
      if (not (data->sumIEt.at(ii) >= 1)) continue;
          

    pass = true;
    break;
  }

  return pass;
}

      


bool
SingleMU_i10
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // 1.2451875 <= eta <= 2.45
              etaWindow1 = ((115 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 225));
            
                        // -2.45 <= eta <= -1.2451875
              etaWindow2 = ((-225 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -115));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i106
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.4083124999999999 <= eta <= 1.4083124999999999
              etaWindow1 = ((-129 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 129));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i109
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.1038124999999999 <= eta <= 1.1038124999999999
              etaWindow1 = ((-101 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 101));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i11
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -0.7993125 <= eta <= 0.7993125
              etaWindow1 = ((-73 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 73));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i12
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // 0.7993125 <= eta <= 1.2451875
              etaWindow1 = ((74 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 114));
            
                        // -1.2451875 <= eta <= -0.7993125
              etaWindow2 = ((-114 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -74));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i13
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // 1.2451875 <= eta <= 2.45
              etaWindow1 = ((115 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 225));
            
                        // -2.45 <= eta <= -1.2451875
              etaWindow2 = ((-225 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -115));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i14
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i141
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125000000003 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i15
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i155
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -2.1043125000000003 <= eta <= 2.1043125
              etaWindow1 = ((-193 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 193));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i16
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i162
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU6: ET >= 13 at BX = 0
      if (not (data->muonIEt.at(idx) >= 13)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i17
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i18
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

                        // -0.7993125 <= eta <= 0.7993125
              etaWindow1 = ((-73 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 73));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i19
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

                        // 0.7993125 <= eta <= 1.2451875
              etaWindow1 = ((74 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 114));
            
                        // -1.2451875 <= eta <= -0.7993125
              etaWindow2 = ((-114 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -74));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i20
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

                        // 1.2451875 <= eta <= 2.45
              etaWindow1 = ((115 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 225));
            
                        // -2.45 <= eta <= -1.2451875
              etaWindow2 = ((-225 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -115));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i21
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i22
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU20: ET >= 41 at BX = 0
      if (not (data->muonIEt.at(idx) >= 41)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i23
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i24
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -0.7993125 <= eta <= 0.7993125
              etaWindow1 = ((-73 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 73));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i249
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i25
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // 0.7993125 <= eta <= 1.2451875
              etaWindow1 = ((74 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 114));
            
                        // -1.2451875 <= eta <= -0.7993125
              etaWindow2 = ((-114 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -74));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i252
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i257
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU15: ET >= 31 at BX = 0
      if (not (data->muonIEt.at(idx) >= 31)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i26
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU22: ET >= 45 at BX = 0
      if (not (data->muonIEt.at(idx) >= 45)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // 1.2451875 <= eta <= 2.45
              etaWindow1 = ((115 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 225));
            
                        // -2.45 <= eta <= -1.2451875
              etaWindow2 = ((-225 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -115));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i27
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU25: ET >= 51 at BX = 0
      if (not (data->muonIEt.at(idx) >= 51)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i271
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU8: ET >= 17 at BX = 0
      if (not (data->muonIEt.at(idx) >= 17)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i293
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU6: ET >= 13 at BX = 0
      if (not (data->muonIEt.at(idx) >= 13)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i294
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU8: ET >= 17 at BX = 0
      if (not (data->muonIEt.at(idx) >= 17)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i295
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU7: ET >= 15 at BX = 0
      if (not (data->muonIEt.at(idx) >= 15)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i296
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU9: ET >= 19 at BX = 0
      if (not (data->muonIEt.at(idx) >= 19)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i297
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU10: ET >= 21 at BX = 0
      if (not (data->muonIEt.at(idx) >= 21)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i298
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU12: ET >= 25 at BX = 0
      if (not (data->muonIEt.at(idx) >= 25)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i299
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU14: ET >= 29 at BX = 0
      if (not (data->muonIEt.at(idx) >= 29)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i300
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU16: ET >= 33 at BX = 0
      if (not (data->muonIEt.at(idx) >= 33)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i301
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU18: ET >= 37 at BX = 0
      if (not (data->muonIEt.at(idx) >= 37)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

                        // -1.5061874999999998 <= eta <= 1.5061875
              etaWindow1 = ((-138 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 138));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i6
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i7
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i8
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // -0.7993125 <= eta <= 0.7993125
              etaWindow1 = ((-73 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 73));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i9
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1, etaWindow2;
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // 0.7993125 <= eta <= 1.2451875
              etaWindow1 = ((74 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= 114));
            
                        // -1.2451875 <= eta <= -0.7993125
              etaWindow2 = ((-114 <= data->muonIEtaAtVtx.at(idx)) and (data->muonIEtaAtVtx.at(idx) <= -74));
            
          if (not (etaWindow1 or etaWindow2)) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i138
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU40: ET >= 80 at BX = 0
      if (not (data->tauIEt.at(idx) >= 80)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i143
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU52: ET >= 104 at BX = 0
      if (not (data->tauIEt.at(idx) >= 104)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i156
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU24: ET >= 48 at BX = 0
      if (not (data->tauIEt.at(idx) >= 48)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i157
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU26: ET >= 52 at BX = 0
      if (not (data->tauIEt.at(idx) >= 52)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i158
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU70: ET >= 140 at BX = 0
      if (not (data->tauIEt.at(idx) >= 140)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i159
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU32: ET >= 64 at BX = 0
      if (not (data->tauIEt.at(idx) >= 64)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i160
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU34: ET >= 68 at BX = 0
      if (not (data->tauIEt.at(idx) >= 68)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i161
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU36: ET >= 72 at BX = 0
      if (not (data->tauIEt.at(idx) >= 72)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i360
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU30: ET >= 60 at BX = 0
      if (not (data->tauIEt.at(idx) >= 60)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i361
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU28: ET >= 56 at BX = 0
      if (not (data->tauIEt.at(idx) >= 56)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
                        // isolation : 0xe
      if (not ((14 >> data->tauIso.at(idx)) & 1)) continue;

          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i67
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU120: ET >= 240 at BX = 0
      if (not (data->tauIEt.at(idx) >= 240)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i68
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
          if (nobj > 12) break;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 1) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 1, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(1, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // TAU130: ET >= 260 at BX = 0
      if (not (data->tauIEt.at(idx) >= 260)) continue;

                        // -2.1315 <= eta <= 2.1315
              etaWindow1 = ((-49 <= data->tauIEta.at(idx)) and (data->tauIEta.at(idx) <= 48));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

                    




bool
TransverseMass_i58
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  bool etaWindow1;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
          
                                  // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(ii) >= 64)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEt)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
                          // ETM10: ET >= 20 at BX = 0
      if (not (data->sumIEt.at(jj) >= 20)) continue;
      
          long long minimum;
  long long maximum;
        // 40.0 <= Mt <= 151982.0
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->sumIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long cosDeltaPhi = LUT_COS_DPHI_EG_ETM[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(ii)];
  long long pt1 = LUT_ETM_ET[data->sumIEt.at(jj)];
  long long mt2 = pt0*pt1*(1*POW10[3] - cosDeltaPhi);
    minimum = (long long)(800.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mt2) and (mt2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}
    
                    




bool
TransverseMass_i59
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  bool etaWindow1;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
          
                                  // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(ii) >= 64)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEt)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
                          // ETM10: ET >= 20 at BX = 0
      if (not (data->sumIEt.at(jj) >= 20)) continue;
      
          long long minimum;
  long long maximum;
        // 44.0 <= Mt <= 151982.0
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->sumIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long cosDeltaPhi = LUT_COS_DPHI_EG_ETM[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(ii)];
  long long pt1 = LUT_ETM_ET[data->sumIEt.at(jj)];
  long long mt2 = pt0*pt1*(1*POW10[3] - cosDeltaPhi);
    minimum = (long long)(968.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mt2) and (mt2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}
    
                    




bool
TransverseMass_i60
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  bool etaWindow1;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
          
                                  // EG32: ET >= 64 at BX = 0
      if (not (data->egIEt.at(ii) >= 64)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(ii)) and (data->egIEta.at(ii) <= 57));
            
                        // isolation : 0xa
      if (not ((10 >> data->egIso.at(ii)) & 1)) continue;

          if (not etaWindow1) continue;

    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEt)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
                          // ETM10: ET >= 20 at BX = 0
      if (not (data->sumIEt.at(jj) >= 20)) continue;
      
          long long minimum;
  long long maximum;
        // 48.0 <= Mt <= 151982.0
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->sumIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
  long long cosDeltaPhi = LUT_COS_DPHI_EG_ETM[deltaIPhi];
  long long pt0 = LUT_EG_ET[data->egIEt.at(ii)];
  long long pt1 = LUT_ETM_ET[data->sumIEt.at(jj)];
  long long mt2 = pt0*pt1*(1*POW10[3] - cosDeltaPhi);
    minimum = (long long)(1152.0 * POW10[5]);
  maximum = (long long)(11549264162.0 * POW10[5]);
  if (not ((minimum <= mt2) and (mt2 <= maximum))) continue;

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}
    
      


bool
TripleEG_i244
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_i281
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_i282
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG16: ET >= 32 at BX = 0
      if (not (data->egIEt.at(idx) >= 32)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG15: ET >= 30 at BX = 0
      if (not (data->egIEt.at(idx) >= 30)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_i283
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG12: ET >= 24 at BX = 0
      if (not (data->egIEt.at(idx) >= 24)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_i66
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      bool etaWindow1;
            idx = candidates.at(set.at(indicies.at(0)));
                                // EG18: ET >= 36 at BX = 0
      if (not (data->egIEt.at(idx) >= 36)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(1)));
                                // EG17: ET >= 34 at BX = 0
      if (not (data->egIEt.at(idx) >= 34)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            idx = candidates.at(set.at(indicies.at(2)));
                                // EG8: ET >= 16 at BX = 0
      if (not (data->egIEt.at(idx) >= 16)) continue;

                        // -2.5229999999999997 <= eta <= 2.5229999999999997
              etaWindow1 = ((-58 <= data->egIEta.at(idx)) and (data->egIEta.at(idx) <= 57));
            
          if (not etaWindow1) continue;
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleJET_i259
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET100: ET >= 200 at BX = 0
      if (not (data->jetIEt.at(idx) >= 200)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET80: ET >= 160 at BX = 0
      if (not (data->jetIEt.at(idx) >= 160)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET70: ET >= 140 at BX = 0
      if (not (data->jetIEt.at(idx) >= 140)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleJET_i261
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET95: ET >= 190 at BX = 0
      if (not (data->jetIEt.at(idx) >= 190)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET65: ET >= 130 at BX = 0
      if (not (data->jetIEt.at(idx) >= 130)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleJET_i89
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // JET105: ET >= 210 at BX = 0
      if (not (data->jetIEt.at(idx) >= 210)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // JET85: ET >= 170 at BX = 0
      if (not (data->jetIEt.at(idx) >= 170)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // JET75: ET >= 150 at BX = 0
      if (not (data->jetIEt.at(idx) >= 150)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i172
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i186
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3p5: ET >= 8 at BX = 0
      if (not (data->muonIEt.at(idx) >= 8)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx) >= 6)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i188
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU4: ET >= 9 at BX = 0
      if (not (data->muonIEt.at(idx) >= 9)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx) >= 6)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i189
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i191
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i199
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i204
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3p5: ET >= 8 at BX = 0
      if (not (data->muonIEt.at(idx) >= 8)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU2p5: ET >= 6 at BX = 0
      if (not (data->muonIEt.at(idx) >= 6)) continue;

                        // quality : 0xfff0
      if (not ((65520 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i245
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i286
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xf000
      if (not ((61440 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i32
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU0: ET >= 1 at BX = 0
      if (not (data->muonIEt.at(idx) >= 1)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i33
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i34
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i35
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
               candidates.push_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  std::vector<std::vector<int> > combination;
  getCombination(candidates.size(), 3, combination);
  std::vector<std::vector<int> > permutation;
  getPermutation(3, permutation);

  for (size_t ii = 0; ii < combination.size(); ii++)
  {
    const std::vector<int>& set = combination.at(ii);
    for (size_t jj = 0; jj < permutation.size(); jj++)
    {
      const std::vector<int>& indicies = permutation.at(jj);
      int idx = -1;
      
            idx = candidates.at(set.at(indicies.at(0)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(1)));
                                // MU5: ET >= 11 at BX = 0
      if (not (data->muonIEt.at(idx) >= 11)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            idx = candidates.at(set.at(indicies.at(2)));
                                // MU3: ET >= 7 at BX = 0
      if (not (data->muonIEt.at(idx) >= 7)) continue;

                        // quality : 0xff00
      if (not ((65280 >> data->muonQual.at(idx)) & 1)) continue;

          
            
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

  

// generate algorithms
bool
L1_AlwaysTrue(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i208(data) or ( not SingleEXT_i208(data));
}
bool
L1_BPTX_AND_Ref1_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i223(data);
}
bool
L1_BPTX_AND_Ref3_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i227(data);
}
bool
L1_BPTX_AND_Ref4_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i231(data);
}
bool
L1_BPTX_BeamGas_B1_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i221(data);
}
bool
L1_BPTX_BeamGas_B2_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i222(data);
}
bool
L1_BPTX_BeamGas_Ref1_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i219(data);
}
bool
L1_BPTX_BeamGas_Ref2_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i220(data);
}
bool
L1_BPTX_NotOR_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i226(data);
}
bool
L1_BPTX_OR_Ref3_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i228(data);
}
bool
L1_BPTX_OR_Ref4_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i232(data);
}
bool
L1_BPTX_RefAND_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i229(data);
}
bool
L1_BptxMinus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i218(data);
}
bool
L1_BptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i105(data);
}
bool
L1_BptxPlus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i217(data);
}
bool
L1_BptxXOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return (SingleEXT_i217(data) and ( not SingleEXT_i218(data))) or (SingleEXT_i218(data) and ( not SingleEXT_i217(data)));
}
bool
L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i235(data);
}
bool
L1_DoubleEG8er2p5_HTT260er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i173(data) and SingleHTT_i174(data);
}
bool
L1_DoubleEG8er2p5_HTT280er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i173(data) and SingleHTT_i96(data);
}
bool
L1_DoubleEG8er2p5_HTT300er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i173(data) and SingleHTT_i241(data);
}
bool
L1_DoubleEG8er2p5_HTT320er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i173(data) and SingleHTT_i97(data);
}
bool
L1_DoubleEG8er2p5_HTT340er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i173(data) and SingleHTT_i242(data);
}
bool
L1_DoubleEG_15_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i278(data);
}
bool
L1_DoubleEG_20_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i279(data);
}
bool
L1_DoubleEG_22_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i61(data);
}
bool
L1_DoubleEG_25_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i62(data);
}
bool
L1_DoubleEG_25_14_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i63(data);
}
bool
L1_DoubleEG_27_14_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i280(data);
}
bool
L1_DoubleEG_LooseIso20_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i145(data);
}
bool
L1_DoubleEG_LooseIso22_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i144(data);
}
bool
L1_DoubleEG_LooseIso22_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i268(data);
}
bool
L1_DoubleEG_LooseIso25_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i269(data);
}
bool
L1_DoubleIsoTau28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i357(data);
}
bool
L1_DoubleIsoTau28er2p1_Mass_Max80(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i354(data);
}
bool
L1_DoubleIsoTau28er2p1_Mass_Max90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i355(data);
}
bool
L1_DoubleIsoTau30er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i356(data);
}
bool
L1_DoubleIsoTau30er2p1_Mass_Max80(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i359(data);
}
bool
L1_DoubleIsoTau30er2p1_Mass_Max90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i358(data);
}
bool
L1_DoubleIsoTau32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i70(data);
}
bool
L1_DoubleIsoTau34er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i71(data);
}
bool
L1_DoubleIsoTau36er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i72(data);
}
bool
L1_DoubleJet100er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i132(data);
}
bool
L1_DoubleJet100er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i88(data);
}
bool
L1_DoubleJet112er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i133(data);
}
bool
L1_DoubleJet120er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i312(data);
}
bool
L1_DoubleJet150er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i313(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i265(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i266(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i114(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i113(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i115(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i116(data);
}
bool
L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMassOvRm_i316(data);
}
bool
L1_DoubleJet40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i87(data);
}
bool
L1_DoubleJet_100_30_DoubleJet30_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i146(data) and InvariantMass_i147(data);
}
bool
L1_DoubleJet_110_35_DoubleJet35_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i148(data) and InvariantMass_i149(data);
}
bool
L1_DoubleJet_115_40_DoubleJet40_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i150(data) and InvariantMass_i151(data);
}
bool
L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i325(data) and (InvariantMass_i326(data) or InvariantMass_i327(data) or InvariantMass_i328(data) or InvariantMass_i329(data) or InvariantMass_i330(data) or InvariantMass_i331(data));
}
bool
L1_DoubleJet_120_45_DoubleJet45_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i309(data) and InvariantMass_i310(data);
}
bool
L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i76(data) and (InvariantMass_i332(data) or InvariantMass_i327(data) or InvariantMass_i333(data) or InvariantMass_i334(data) or InvariantMass_i335(data) or InvariantMass_i336(data));
}
bool
L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i270(data) and DoubleMU_i193(data);
}
bool
L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMassOvRm_i264(data);
}
bool
L1_DoubleJet_80_30_Mass_Min420_Mu8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i270(data) and SingleMU_i271(data);
}
bool
L1_DoubleJet_90_30_DoubleJet30_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i258(data) and InvariantMass_i147(data);
}
bool
L1_DoubleLooseIsoEG22er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i64(data);
}
bool
L1_DoubleLooseIsoEG24er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i65(data);
}
bool
L1_DoubleMu0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i28(data);
}
bool
L1_DoubleMu0_Mass_Min1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i175(data);
}
bool
L1_DoubleMu0_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i176(data);
}
bool
L1_DoubleMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i193(data);
}
bool
L1_DoubleMu0_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i194(data);
}
bool
L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i304(data) and CaloMuonCorrelation_i305(data);
}
bool
L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i184(data);
}
bool
L1_DoubleMu0er1p5_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i292(data);
}
bool
L1_DoubleMu0er1p5_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i197(data);
}
bool
L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i181(data);
}
bool
L1_DoubleMu0er1p5_SQ_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i118(data);
}
bool
L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i267(data);
}
bool
L1_DoubleMu0er2p0_SQ_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i117(data);
}
bool
L1_DoubleMu18er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i31(data);
}
bool
L1_DoubleMu3_OS_DoubleEG7p5Upsilon(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i202(data) and InvariantMass_i203(data);
}
bool
L1_DoubleMu3_SQ_ETMHF50_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i165(data) and SingleETMHF_i166(data) and SingleHTT_i179(data);
}
bool
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i165(data) and SingleETMHF_i166(data) and SingleJET_i167(data);
}
bool
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i165(data) and SingleETMHF_i166(data) and (SingleJET_i167(data) or DoubleJET_i87(data));
}
bool
L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i165(data) and SingleETMHF_i240(data) and SingleJET_i167(data);
}
bool
L1_DoubleMu3_SQ_HTT220er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i165(data) and SingleHTT_i168(data);
}
bool
L1_DoubleMu3_SQ_HTT240er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i165(data) and SingleHTT_i163(data);
}
bool
L1_DoubleMu3_SQ_HTT260er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i165(data) and SingleHTT_i174(data);
}
bool
L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i302(data) and CaloMuonCorrelation_i303(data);
}
bool
L1_DoubleMu4_SQ_EG9er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i256(data) and SingleEG_i178(data);
}
bool
L1_DoubleMu4_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i198(data);
}
bool
L1_DoubleMu4_SQ_OS_dR_Max1p2(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i183(data);
}
bool
L1_DoubleMu4p5_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i195(data);
}
bool
L1_DoubleMu4p5_SQ_OS_dR_Max1p2(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i185(data);
}
bool
L1_DoubleMu4p5er2p0_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i196(data);
}
bool
L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i338(data);
}
bool
L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i182(data);
}
bool
L1_DoubleMu5Upsilon_OS_DoubleEG3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i200(data) and InvariantMass_i201(data);
}
bool
L1_DoubleMu5_SQ_EG9er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i177(data) and SingleEG_i178(data);
}
bool
L1_DoubleMu8_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i243(data);
}
bool
L1_DoubleMu9_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i192(data);
}
bool
L1_DoubleMu_12_5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i29(data);
}
bool
L1_DoubleMu_15_5_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i169(data);
}
bool
L1_DoubleMu_15_7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i30(data);
}
bool
L1_DoubleMu_15_7_Mass_Min1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i171(data);
}
bool
L1_DoubleMu_15_7_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i170(data);
}
bool
L1_DoubleTau70er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i69(data);
}
bool
L1_ETM120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETM_i314(data);
}
bool
L1_ETM150(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETM_i315(data);
}
bool
L1_ETMHF100(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i101(data);
}
bool
L1_ETMHF100_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i101(data) and SingleHTT_i179(data);
}
bool
L1_ETMHF110(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i102(data);
}
bool
L1_ETMHF110_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i102(data) and SingleHTT_i179(data);
}
bool
L1_ETMHF110_HTT60er_NotSecondBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i102(data) and SingleHTT_i179(data) and ((SingleEXT_i207(data)) or ( not SingleEXT_i213(data)) or ( not SingleEXT_i208(data)) or ( not SingleEXT_i210(data)) or ( not SingleEXT_i211(data)));
}
bool
L1_ETMHF120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i103(data);
}
bool
L1_ETMHF120_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i103(data) and SingleHTT_i179(data);
}
bool
L1_ETMHF120_NotSecondBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i103(data) and ((SingleEXT_i207(data)) or ( not SingleEXT_i213(data)) or ( not SingleEXT_i208(data)) or ( not SingleEXT_i210(data)) or ( not SingleEXT_i211(data)));
}
bool
L1_ETMHF130(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i104(data);
}
bool
L1_ETMHF130_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i104(data) and SingleHTT_i179(data);
}
bool
L1_ETMHF140(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i263(data);
}
bool
L1_ETMHF150(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i311(data);
}
bool
L1_ETMHF90_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i139(data) and SingleHTT_i179(data);
}
bool
L1_ETT1200(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETT_i322(data);
}
bool
L1_ETT1600(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETT_i323(data);
}
bool
L1_ETT2000(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETT_i324(data);
}
bool
L1_FirstBunchAfterTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i212(data) and SingleEXT_i213(data) and ( not SingleEXT_i105(data)) and ( not SingleEXT_i108(data)) and ( not SingleEXT_i209(data));
}
bool
L1_FirstBunchBeforeTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_i207(data)) and ( not SingleEXT_i107(data)) and ( not SingleEXT_i105(data)) and SingleEXT_i210(data) and SingleEXT_i211(data);
}
bool
L1_FirstBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_i207(data)) and ( not SingleEXT_i107(data)) and SingleEXT_i208(data) and SingleEXT_i210(data) and SingleEXT_i211(data);
}
bool
L1_FirstCollisionInOrbit(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i230(data);
}
bool
L1_FirstCollisionInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i234(data);
}
bool
L1_HCAL_LaserMon_Trig(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i0(data);
}
bool
L1_HCAL_LaserMon_Veto(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i1(data);
}
bool
L1_HTT120er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i92(data);
}
bool
L1_HTT160er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i93(data);
}
bool
L1_HTT200er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i94(data);
}
bool
L1_HTT255er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i95(data);
}
bool
L1_HTT280er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i96(data);
}
bool
L1_HTT280er_QuadJet_70_55_40_35_er2p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i96(data) and QuadJET_i129(data);
}
bool
L1_HTT320er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i97(data);
}
bool
L1_HTT320er_QuadJet_70_55_40_40_er2p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i97(data) and QuadJET_i130(data);
}
bool
L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i97(data) and QuadJET_i131(data);
}
bool
L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i97(data) and QuadJET_i239(data);
}
bool
L1_HTT360er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i98(data);
}
bool
L1_HTT400er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i99(data);
}
bool
L1_HTT450er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i100(data);
}
bool
L1_IsoEG32er2p5_Mt40(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TransverseMass_i58(data);
}
bool
L1_IsoEG32er2p5_Mt44(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TransverseMass_i59(data);
}
bool
L1_IsoEG32er2p5_Mt48(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TransverseMass_i60(data);
}
bool
L1_IsoTau40er2p1_ETMHF100(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i138(data) and SingleETMHF_i101(data);
}
bool
L1_IsoTau40er2p1_ETMHF110(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i138(data) and SingleETMHF_i102(data);
}
bool
L1_IsoTau40er2p1_ETMHF80(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i138(data) and SingleETMHF_i140(data);
}
bool
L1_IsoTau40er2p1_ETMHF90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i138(data) and SingleETMHF_i139(data);
}
bool
L1_IsolatedBunch(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_i207(data)) and ( not SingleEXT_i107(data)) and SingleEXT_i208(data) and ( not SingleEXT_i108(data)) and ( not SingleEXT_i209(data));
}
bool
L1_LastBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i212(data) and SingleEXT_i213(data) and SingleEXT_i208(data) and ( not SingleEXT_i108(data)) and ( not SingleEXT_i209(data));
}
bool
L1_LastCollisionInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i233(data);
}
bool
L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i153(data);
}
bool
L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i152(data);
}
bool
L1_LooseIsoEG24er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i122(data) and SingleHTT_i123(data);
}
bool
L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i154(data);
}
bool
L1_LooseIsoEG26er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i124(data) and SingleHTT_i123(data);
}
bool
L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i119(data);
}
bool
L1_LooseIsoEG28er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i125(data) and SingleHTT_i123(data);
}
bool
L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i120(data);
}
bool
L1_LooseIsoEG30er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i248(data) and SingleHTT_i123(data);
}
bool
L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i121(data);
}
bool
L1_MinimumBiasHF0_AND_BptxAND(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return (SingleMBT0HFP_i215(data) and SingleMBT0HFM_i216(data)) and SingleEXT_i208(data);
}
bool
L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i134(data) and CaloCaloCorrelation_i135(data);
}
bool
L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i254(data) and CaloCaloCorrelation_i255(data);
}
bool
L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i136(data) and CaloCaloCorrelation_i137(data);
}
bool
L1_Mu18er2p1_Tau24er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i155(data) and SingleTAU_i156(data);
}
bool
L1_Mu18er2p1_Tau26er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i155(data) and SingleTAU_i157(data);
}
bool
L1_Mu20_EG10er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i22(data) and SingleEG_i38(data);
}
bool
L1_Mu22er2p1_IsoTau28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i141(data) and SingleTAU_i361(data);
}
bool
L1_Mu22er2p1_IsoTau30er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i141(data) and SingleTAU_i360(data);
}
bool
L1_Mu22er2p1_IsoTau32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i141(data) and SingleTAU_i159(data);
}
bool
L1_Mu22er2p1_IsoTau34er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i141(data) and SingleTAU_i160(data);
}
bool
L1_Mu22er2p1_IsoTau36er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i141(data) and SingleTAU_i161(data);
}
bool
L1_Mu22er2p1_IsoTau40er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i141(data) and SingleTAU_i138(data);
}
bool
L1_Mu22er2p1_Tau70er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i141(data) and SingleTAU_i158(data);
}
bool
L1_Mu3_Jet120er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i274(data);
}
bool
L1_Mu3_Jet120er2p5_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i277(data);
}
bool
L1_Mu3_Jet16er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i272(data);
}
bool
L1_Mu3_Jet30er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i14(data) and SingleJET_i180(data);
}
bool
L1_Mu3_Jet35er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i275(data);
}
bool
L1_Mu3_Jet60er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i273(data);
}
bool
L1_Mu3_Jet80er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i276(data);
}
bool
L1_Mu3er1p5_Jet100er2p5_ETMHF40(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i249(data) and SingleJET_i250(data) and SingleETMHF_i251(data);
}
bool
L1_Mu3er1p5_Jet100er2p5_ETMHF50(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i249(data) and SingleJET_i250(data) and SingleETMHF_i166(data);
}
bool
L1_Mu5_EG23er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i15(data) and SingleEG_i126(data);
}
bool
L1_Mu5_LooseIsoEG20er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i15(data) and SingleEG_i127(data);
}
bool
L1_Mu6_DoubleEG10er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i162(data) and DoubleEG_i214(data);
}
bool
L1_Mu6_DoubleEG12er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i162(data) and DoubleEG_i246(data);
}
bool
L1_Mu6_DoubleEG15er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i162(data) and DoubleEG_i247(data);
}
bool
L1_Mu6_DoubleEG17er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i162(data) and DoubleEG_i206(data);
}
bool
L1_Mu6_HTT240er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i162(data) and SingleHTT_i163(data);
}
bool
L1_Mu6_HTT250er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i162(data) and SingleHTT_i164(data);
}
bool
L1_Mu7_EG20er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i16(data) and SingleEG_i352(data);
}
bool
L1_Mu7_EG23er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i16(data) and SingleEG_i126(data);
}
bool
L1_Mu7_LooseIsoEG20er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i16(data) and SingleEG_i127(data);
}
bool
L1_Mu7_LooseIsoEG23er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i16(data) and SingleEG_i128(data);
}
bool
L1_NotBptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return not SingleEXT_i105(data);
}
bool
L1_QuadJet36er2p5_IsoTau52er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadJET_i142(data) and SingleTAU_i143(data);
}
bool
L1_QuadJet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadJET_i91(data);
}
bool
L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadJET_i306(data) and DoubleJET_i262(data) and (SingleJET_i307(data) or SingleJET_i308(data));
}
bool
L1_QuadMu0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadMU_i36(data);
}
bool
L1_QuadMu0_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadMU_i284(data);
}
bool
L1_QuadMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadMU_i285(data);
}
bool
L1_SecondBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_i207(data)) and SingleEXT_i213(data) and SingleEXT_i208(data) and SingleEXT_i210(data) and SingleEXT_i211(data);
}
bool
L1_SecondLastBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i212(data) and SingleEXT_i213(data) and SingleEXT_i208(data) and SingleEXT_i210(data) and ( not SingleEXT_i209(data));
}
bool
L1_SingleEG10er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i38(data);
}
bool
L1_SingleEG15er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i39(data);
}
bool
L1_SingleEG26er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i40(data);
}
bool
L1_SingleEG28_FWD2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i346(data) or SingleEG_i347(data);
}
bool
L1_SingleEG28er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i343(data);
}
bool
L1_SingleEG28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i344(data);
}
bool
L1_SingleEG28er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i345(data);
}
bool
L1_SingleEG34er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i43(data);
}
bool
L1_SingleEG36er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i44(data);
}
bool
L1_SingleEG38er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i45(data);
}
bool
L1_SingleEG40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i46(data);
}
bool
L1_SingleEG42er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i47(data);
}
bool
L1_SingleEG45er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i41(data);
}
bool
L1_SingleEG50(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i42(data);
}
bool
L1_SingleEG60(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i337(data);
}
bool
L1_SingleEG8er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i37(data);
}
bool
L1_SingleIsoEG24er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i317(data);
}
bool
L1_SingleIsoEG24er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i48(data);
}
bool
L1_SingleIsoEG26er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i318(data);
}
bool
L1_SingleIsoEG26er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i49(data);
}
bool
L1_SingleIsoEG26er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i53(data);
}
bool
L1_SingleIsoEG28_FWD2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i350(data) or SingleEG_i351(data);
}
bool
L1_SingleIsoEG28er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i319(data);
}
bool
L1_SingleIsoEG28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i50(data);
}
bool
L1_SingleIsoEG28er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i54(data);
}
bool
L1_SingleIsoEG30er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i51(data);
}
bool
L1_SingleIsoEG30er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i55(data);
}
bool
L1_SingleIsoEG32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i52(data);
}
bool
L1_SingleIsoEG32er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i56(data);
}
bool
L1_SingleIsoEG34er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i57(data);
}
bool
L1_SingleJet10erHE(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i237(data);
}
bool
L1_SingleJet120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i76(data);
}
bool
L1_SingleJet120_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i85(data) or SingleJET_i86(data);
}
bool
L1_SingleJet120er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i287(data);
}
bool
L1_SingleJet12erHE(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i238(data);
}
bool
L1_SingleJet140er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i253(data);
}
bool
L1_SingleJet140er2p5_ETMHF70(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i253(data) and SingleETMHF_i353(data);
}
bool
L1_SingleJet140er2p5_ETMHF80(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i253(data) and SingleETMHF_i140(data);
}
bool
L1_SingleJet140er2p5_ETMHF90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i253(data) and SingleETMHF_i139(data);
}
bool
L1_SingleJet160er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i290(data);
}
bool
L1_SingleJet180(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i77(data);
}
bool
L1_SingleJet180er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i291(data);
}
bool
L1_SingleJet200(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i78(data);
}
bool
L1_SingleJet20er2p5_NotBptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i110(data) and ( not SingleEXT_i105(data));
}
bool
L1_SingleJet20er2p5_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i110(data) and ( not SingleEXT_i107(data)) and ( not SingleEXT_i105(data)) and ( not SingleEXT_i108(data));
}
bool
L1_SingleJet35(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i73(data);
}
bool
L1_SingleJet35_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i79(data) or SingleJET_i80(data);
}
bool
L1_SingleJet35er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i288(data);
}
bool
L1_SingleJet43er2p5_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i111(data) and ( not SingleEXT_i107(data)) and ( not SingleEXT_i105(data)) and ( not SingleEXT_i108(data));
}
bool
L1_SingleJet46er2p5_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i112(data) and ( not SingleEXT_i107(data)) and ( not SingleEXT_i105(data)) and ( not SingleEXT_i108(data));
}
bool
L1_SingleJet60(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i74(data);
}
bool
L1_SingleJet60_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i81(data) or SingleJET_i82(data);
}
bool
L1_SingleJet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i167(data);
}
bool
L1_SingleJet8erHE(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i236(data);
}
bool
L1_SingleJet90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i75(data);
}
bool
L1_SingleJet90_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i83(data) or SingleJET_i84(data);
}
bool
L1_SingleJet90er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i289(data);
}
bool
L1_SingleLooseIsoEG26er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i339(data);
}
bool
L1_SingleLooseIsoEG26er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i340(data);
}
bool
L1_SingleLooseIsoEG28_FWD2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i348(data) or SingleEG_i349(data);
}
bool
L1_SingleLooseIsoEG28er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i320(data);
}
bool
L1_SingleLooseIsoEG28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i125(data);
}
bool
L1_SingleLooseIsoEG28er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i341(data);
}
bool
L1_SingleLooseIsoEG30er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i321(data);
}
bool
L1_SingleLooseIsoEG30er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i342(data);
}
bool
L1_SingleMu0_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i11(data);
}
bool
L1_SingleMu0_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i252(data);
}
bool
L1_SingleMu0_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i13(data);
}
bool
L1_SingleMu0_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i12(data);
}
bool
L1_SingleMu10er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i297(data);
}
bool
L1_SingleMu12_DQ_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i18(data);
}
bool
L1_SingleMu12_DQ_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i20(data);
}
bool
L1_SingleMu12_DQ_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i19(data);
}
bool
L1_SingleMu12er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i298(data);
}
bool
L1_SingleMu14er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i299(data);
}
bool
L1_SingleMu15_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i257(data);
}
bool
L1_SingleMu16er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i300(data);
}
bool
L1_SingleMu18(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i21(data);
}
bool
L1_SingleMu18er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i301(data);
}
bool
L1_SingleMu20(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i22(data);
}
bool
L1_SingleMu22(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i23(data);
}
bool
L1_SingleMu22_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i24(data);
}
bool
L1_SingleMu22_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i26(data);
}
bool
L1_SingleMu22_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i25(data);
}
bool
L1_SingleMu25(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i27(data);
}
bool
L1_SingleMu3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i14(data);
}
bool
L1_SingleMu5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i15(data);
}
bool
L1_SingleMu6er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i293(data);
}
bool
L1_SingleMu7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i16(data);
}
bool
L1_SingleMu7_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i17(data);
}
bool
L1_SingleMu7er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i295(data);
}
bool
L1_SingleMu8er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i294(data);
}
bool
L1_SingleMu9er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i296(data);
}
bool
L1_SingleMuCosmics(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i6(data);
}
bool
L1_SingleMuCosmics_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i8(data);
}
bool
L1_SingleMuCosmics_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i10(data);
}
bool
L1_SingleMuCosmics_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i9(data);
}
bool
L1_SingleMuOpen(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i7(data);
}
bool
L1_SingleMuOpen_NotBptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i7(data) and ( not SingleEXT_i105(data));
}
bool
L1_SingleMuOpen_er1p1_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i109(data) and ( not SingleEXT_i107(data)) and ( not SingleEXT_i105(data)) and ( not SingleEXT_i108(data));
}
bool
L1_SingleMuOpen_er1p4_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i106(data) and ( not SingleEXT_i107(data)) and ( not SingleEXT_i105(data)) and ( not SingleEXT_i108(data));
}
bool
L1_SingleTau120er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i67(data);
}
bool
L1_SingleTau130er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i68(data);
}
bool
L1_TOTEM_1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i2(data);
}
bool
L1_TOTEM_2(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i3(data);
}
bool
L1_TOTEM_3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i4(data);
}
bool
L1_TOTEM_4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i5(data);
}
bool
L1_TripleEG16er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i244(data);
}
bool
L1_TripleEG_16_12_8_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i281(data);
}
bool
L1_TripleEG_16_15_8_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i282(data);
}
bool
L1_TripleEG_18_17_8_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i66(data);
}
bool
L1_TripleEG_18_18_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i283(data);
}
bool
L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleJET_i259(data) and DoubleJET_i260(data);
}
bool
L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleJET_i89(data) and DoubleJET_i90(data);
}
bool
L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleJET_i261(data) and DoubleJET_i262(data);
}
bool
L1_TripleMu0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i32(data);
}
bool
L1_TripleMu0_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i199(data);
}
bool
L1_TripleMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i245(data);
}
bool
L1_TripleMu3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i33(data);
}
bool
L1_TripleMu3_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i172(data);
}
bool
L1_TripleMu_5SQ_3SQ_0OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i189(data);
}
bool
L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i189(data) and InvariantMass_i190(data);
}
bool
L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i191(data) and InvariantMass_i190(data);
}
bool
L1_TripleMu_5_3_3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i34(data);
}
bool
L1_TripleMu_5_3_3_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i286(data);
}
bool
L1_TripleMu_5_3p5_2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i186(data);
}
bool
L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i186(data) and InvariantMass_i187(data);
}
bool
L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i204(data) and InvariantMass_i205(data);
}
bool
L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i188(data) and InvariantMass_i187(data);
}
bool
L1_TripleMu_5_5_3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i35(data);
}
bool
L1_UnpairedBunchBptxMinus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i225(data);
}
bool
L1_UnpairedBunchBptxPlus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i224(data);
}
bool
L1_ZeroBias(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i208(data);
}
bool
L1_ZeroBias_copy(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i208(data);
}


std::string getNameFromId(const int index)
{
  static const std::pair<int, std::string> id2name[] = {
          std::make_pair(458, "L1_AlwaysTrue"),          std::make_pair(486, "L1_BPTX_AND_Ref1_VME"),          std::make_pair(487, "L1_BPTX_AND_Ref3_VME"),          std::make_pair(488, "L1_BPTX_AND_Ref4_VME"),          std::make_pair(491, "L1_BPTX_BeamGas_B1_VME"),          std::make_pair(492, "L1_BPTX_BeamGas_B2_VME"),          std::make_pair(489, "L1_BPTX_BeamGas_Ref1_VME"),          std::make_pair(490, "L1_BPTX_BeamGas_Ref2_VME"),          std::make_pair(482, "L1_BPTX_NotOR_VME"),          std::make_pair(483, "L1_BPTX_OR_Ref3_VME"),          std::make_pair(484, "L1_BPTX_OR_Ref4_VME"),          std::make_pair(485, "L1_BPTX_RefAND_VME"),          std::make_pair(467, "L1_BptxMinus"),          std::make_pair(464, "L1_BptxOR"),          std::make_pair(466, "L1_BptxPlus"),          std::make_pair(465, "L1_BptxXOR"),          std::make_pair(494, "L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142"),          std::make_pair(247, "L1_DoubleEG8er2p5_HTT260er"),          std::make_pair(248, "L1_DoubleEG8er2p5_HTT280er"),          std::make_pair(249, "L1_DoubleEG8er2p5_HTT300er"),          std::make_pair(250, "L1_DoubleEG8er2p5_HTT320er"),          std::make_pair(251, "L1_DoubleEG8er2p5_HTT340er"),          std::make_pair(205, "L1_DoubleEG_15_10_er2p5"),          std::make_pair(206, "L1_DoubleEG_20_10_er2p5"),          std::make_pair(207, "L1_DoubleEG_22_10_er2p5"),          std::make_pair(208, "L1_DoubleEG_25_12_er2p5"),          std::make_pair(209, "L1_DoubleEG_25_14_er2p5"),          std::make_pair(210, "L1_DoubleEG_27_14_er2p5"),          std::make_pair(212, "L1_DoubleEG_LooseIso20_10_er2p5"),          std::make_pair(213, "L1_DoubleEG_LooseIso22_10_er2p5"),          std::make_pair(214, "L1_DoubleEG_LooseIso22_12_er2p5"),          std::make_pair(215, "L1_DoubleEG_LooseIso25_12_er2p5"),          std::make_pair(269, "L1_DoubleIsoTau28er2p1"),          std::make_pair(275, "L1_DoubleIsoTau28er2p1_Mass_Max80"),          std::make_pair(274, "L1_DoubleIsoTau28er2p1_Mass_Max90"),          std::make_pair(270, "L1_DoubleIsoTau30er2p1"),          std::make_pair(277, "L1_DoubleIsoTau30er2p1_Mass_Max80"),          std::make_pair(276, "L1_DoubleIsoTau30er2p1_Mass_Max90"),          std::make_pair(271, "L1_DoubleIsoTau32er2p1"),          std::make_pair(272, "L1_DoubleIsoTau34er2p1"),          std::make_pair(273, "L1_DoubleIsoTau36er2p1"),          std::make_pair(345, "L1_DoubleJet100er2p3_dEta_Max1p6"),          std::make_pair(341, "L1_DoubleJet100er2p5"),          std::make_pair(346, "L1_DoubleJet112er2p3_dEta_Max1p6"),          std::make_pair(342, "L1_DoubleJet120er2p5"),          std::make_pair(343, "L1_DoubleJet150er2p5"),          std::make_pair(348, "L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5"),          std::make_pair(349, "L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5"),          std::make_pair(350, "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5"),          std::make_pair(351, "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5"),          std::make_pair(352, "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5"),          std::make_pair(353, "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5"),          std::make_pair(363, "L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp"),          std::make_pair(340, "L1_DoubleJet40er2p5"),          std::make_pair(356, "L1_DoubleJet_100_30_DoubleJet30_Mass_Min620"),          std::make_pair(357, "L1_DoubleJet_110_35_DoubleJet35_Mass_Min620"),          std::make_pair(358, "L1_DoubleJet_115_40_DoubleJet40_Mass_Min620"),          std::make_pair(360, "L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28"),          std::make_pair(359, "L1_DoubleJet_120_45_DoubleJet45_Mass_Min620"),          std::make_pair(361, "L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28"),          std::make_pair(366, "L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ"),          std::make_pair(364, "L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp"),          std::make_pair(365, "L1_DoubleJet_80_30_Mass_Min420_Mu8"),          std::make_pair(355, "L1_DoubleJet_90_30_DoubleJet30_Mass_Min620"),          std::make_pair(217, "L1_DoubleLooseIsoEG22er2p1"),          std::make_pair(218, "L1_DoubleLooseIsoEG24er2p1"),          std::make_pair(40, "L1_DoubleMu0"),          std::make_pair(43, "L1_DoubleMu0_Mass_Min1"),          std::make_pair(39, "L1_DoubleMu0_OQ"),          std::make_pair(41, "L1_DoubleMu0_SQ"),          std::make_pair(42, "L1_DoubleMu0_SQ_OS"),          std::make_pair(142, "L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8"),          std::make_pair(59, "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"),          std::make_pair(55, "L1_DoubleMu0er1p5_SQ"),          std::make_pair(56, "L1_DoubleMu0er1p5_SQ_OS"),          std::make_pair(58, "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"),          std::make_pair(57, "L1_DoubleMu0er1p5_SQ_dR_Max1p4"),          std::make_pair(54, "L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4"),          std::make_pair(53, "L1_DoubleMu0er2p0_SQ_dR_Max1p4"),          std::make_pair(51, "L1_DoubleMu18er2p1"),          std::make_pair(112, "L1_DoubleMu3_OS_DoubleEG7p5Upsilon"),          std::make_pair(145, "L1_DoubleMu3_SQ_ETMHF50_HTT60er"),          std::make_pair(147, "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5"),          std::make_pair(146, "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5"),          std::make_pair(148, "L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5"),          std::make_pair(150, "L1_DoubleMu3_SQ_HTT220er"),          std::make_pair(151, "L1_DoubleMu3_SQ_HTT240er"),          std::make_pair(152, "L1_DoubleMu3_SQ_HTT260er"),          std::make_pair(143, "L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8"),          std::make_pair(109, "L1_DoubleMu4_SQ_EG9er2p5"),          std::make_pair(60, "L1_DoubleMu4_SQ_OS"),          std::make_pair(61, "L1_DoubleMu4_SQ_OS_dR_Max1p2"),          std::make_pair(62, "L1_DoubleMu4p5_SQ_OS"),          std::make_pair(63, "L1_DoubleMu4p5_SQ_OS_dR_Max1p2"),          std::make_pair(64, "L1_DoubleMu4p5er2p0_SQ_OS"),          std::make_pair(66, "L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18"),          std::make_pair(65, "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7"),          std::make_pair(113, "L1_DoubleMu5Upsilon_OS_DoubleEG3"),          std::make_pair(110, "L1_DoubleMu5_SQ_EG9er2p5"),          std::make_pair(44, "L1_DoubleMu8_SQ"),          std::make_pair(45, "L1_DoubleMu9_SQ"),          std::make_pair(46, "L1_DoubleMu_12_5"),          std::make_pair(47, "L1_DoubleMu_15_5_SQ"),          std::make_pair(48, "L1_DoubleMu_15_7"),          std::make_pair(50, "L1_DoubleMu_15_7_Mass_Min1"),          std::make_pair(49, "L1_DoubleMu_15_7_SQ"),          std::make_pair(267, "L1_DoubleTau70er2p1"),          std::make_pair(416, "L1_ETM120"),          std::make_pair(417, "L1_ETM150"),          std::make_pair(421, "L1_ETMHF100"),          std::make_pair(429, "L1_ETMHF100_HTT60er"),          std::make_pair(422, "L1_ETMHF110"),          std::make_pair(430, "L1_ETMHF110_HTT60er"),          std::make_pair(444, "L1_ETMHF110_HTT60er_NotSecondBunchInTrain"),          std::make_pair(423, "L1_ETMHF120"),          std::make_pair(431, "L1_ETMHF120_HTT60er"),          std::make_pair(443, "L1_ETMHF120_NotSecondBunchInTrain"),          std::make_pair(424, "L1_ETMHF130"),          std::make_pair(432, "L1_ETMHF130_HTT60er"),          std::make_pair(425, "L1_ETMHF140"),          std::make_pair(426, "L1_ETMHF150"),          std::make_pair(428, "L1_ETMHF90_HTT60er"),          std::make_pair(410, "L1_ETT1200"),          std::make_pair(411, "L1_ETT1600"),          std::make_pair(412, "L1_ETT2000"),          std::make_pair(477, "L1_FirstBunchAfterTrain"),          std::make_pair(472, "L1_FirstBunchBeforeTrain"),          std::make_pair(473, "L1_FirstBunchInTrain"),          std::make_pair(480, "L1_FirstCollisionInOrbit"),          std::make_pair(479, "L1_FirstCollisionInTrain"),          std::make_pair(500, "L1_HCAL_LaserMon_Trig"),          std::make_pair(501, "L1_HCAL_LaserMon_Veto"),          std::make_pair(398, "L1_HTT120er"),          std::make_pair(399, "L1_HTT160er"),          std::make_pair(400, "L1_HTT200er"),          std::make_pair(401, "L1_HTT255er"),          std::make_pair(402, "L1_HTT280er"),          std::make_pair(384, "L1_HTT280er_QuadJet_70_55_40_35_er2p4"),          std::make_pair(403, "L1_HTT320er"),          std::make_pair(385, "L1_HTT320er_QuadJet_70_55_40_40_er2p4"),          std::make_pair(386, "L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3"),          std::make_pair(387, "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3"),          std::make_pair(404, "L1_HTT360er"),          std::make_pair(405, "L1_HTT400er"),          std::make_pair(406, "L1_HTT450er"),          std::make_pair(197, "L1_IsoEG32er2p5_Mt40"),          std::make_pair(198, "L1_IsoEG32er2p5_Mt44"),          std::make_pair(199, "L1_IsoEG32er2p5_Mt48"),          std::make_pair(293, "L1_IsoTau40er2p1_ETMHF100"),          std::make_pair(294, "L1_IsoTau40er2p1_ETMHF110"),          std::make_pair(291, "L1_IsoTau40er2p1_ETMHF80"),          std::make_pair(292, "L1_IsoTau40er2p1_ETMHF90"),          std::make_pair(471, "L1_IsolatedBunch"),          std::make_pair(476, "L1_LastBunchInTrain"),          std::make_pair(478, "L1_LastCollisionInTrain"),          std::make_pair(257, "L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3"),          std::make_pair(259, "L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3"),          std::make_pair(238, "L1_LooseIsoEG24er2p1_HTT100er"),          std::make_pair(258, "L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3"),          std::make_pair(239, "L1_LooseIsoEG26er2p1_HTT100er"),          std::make_pair(234, "L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3"),          std::make_pair(240, "L1_LooseIsoEG28er2p1_HTT100er"),          std::make_pair(235, "L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3"),          std::make_pair(241, "L1_LooseIsoEG30er2p1_HTT100er"),          std::make_pair(236, "L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3"),          std::make_pair(461, "L1_MinimumBiasHF0_AND_BptxAND"),          std::make_pair(134, "L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6"),          std::make_pair(136, "L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6"),          std::make_pair(135, "L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6"),          std::make_pair(279, "L1_Mu18er2p1_Tau24er2p1"),          std::make_pair(280, "L1_Mu18er2p1_Tau26er2p1"),          std::make_pair(99, "L1_Mu20_EG10er2p5"),          std::make_pair(282, "L1_Mu22er2p1_IsoTau28er2p1"),          std::make_pair(283, "L1_Mu22er2p1_IsoTau30er2p1"),          std::make_pair(284, "L1_Mu22er2p1_IsoTau32er2p1"),          std::make_pair(285, "L1_Mu22er2p1_IsoTau34er2p1"),          std::make_pair(286, "L1_Mu22er2p1_IsoTau36er2p1"),          std::make_pair(287, "L1_Mu22er2p1_IsoTau40er2p1"),          std::make_pair(289, "L1_Mu22er2p1_Tau70er2p1"),          std::make_pair(126, "L1_Mu3_Jet120er2p5_dR_Max0p4"),          std::make_pair(125, "L1_Mu3_Jet120er2p5_dR_Max0p8"),          std::make_pair(121, "L1_Mu3_Jet16er2p5_dR_Max0p4"),          std::make_pair(119, "L1_Mu3_Jet30er2p5"),          std::make_pair(122, "L1_Mu3_Jet35er2p5_dR_Max0p4"),          std::make_pair(123, "L1_Mu3_Jet60er2p5_dR_Max0p4"),          std::make_pair(124, "L1_Mu3_Jet80er2p5_dR_Max0p4"),          std::make_pair(128, "L1_Mu3er1p5_Jet100er2p5_ETMHF40"),          std::make_pair(129, "L1_Mu3er1p5_Jet100er2p5_ETMHF50"),          std::make_pair(96, "L1_Mu5_EG23er2p5"),          std::make_pair(100, "L1_Mu5_LooseIsoEG20er2p5"),          std::make_pair(104, "L1_Mu6_DoubleEG10er2p5"),          std::make_pair(105, "L1_Mu6_DoubleEG12er2p5"),          std::make_pair(106, "L1_Mu6_DoubleEG15er2p5"),          std::make_pair(107, "L1_Mu6_DoubleEG17er2p5"),          std::make_pair(131, "L1_Mu6_HTT240er"),          std::make_pair(132, "L1_Mu6_HTT250er"),          std::make_pair(97, "L1_Mu7_EG20er2p5"),          std::make_pair(98, "L1_Mu7_EG23er2p5"),          std::make_pair(101, "L1_Mu7_LooseIsoEG20er2p5"),          std::make_pair(102, "L1_Mu7_LooseIsoEG23er2p5"),          std::make_pair(463, "L1_NotBptxOR"),          std::make_pair(298, "L1_QuadJet36er2p5_IsoTau52er2p1"),          std::make_pair(382, "L1_QuadJet60er2p5"),          std::make_pair(376, "L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0"),          std::make_pair(89, "L1_QuadMu0"),          std::make_pair(88, "L1_QuadMu0_OQ"),          std::make_pair(90, "L1_QuadMu0_SQ"),          std::make_pair(474, "L1_SecondBunchInTrain"),          std::make_pair(475, "L1_SecondLastBunchInTrain"),          std::make_pair(160, "L1_SingleEG10er2p5"),          std::make_pair(161, "L1_SingleEG15er2p5"),          std::make_pair(162, "L1_SingleEG26er2p5"),          std::make_pair(163, "L1_SingleEG28_FWD2p5"),          std::make_pair(166, "L1_SingleEG28er1p5"),          std::make_pair(165, "L1_SingleEG28er2p1"),          std::make_pair(164, "L1_SingleEG28er2p5"),          std::make_pair(167, "L1_SingleEG34er2p5"),          std::make_pair(168, "L1_SingleEG36er2p5"),          std::make_pair(169, "L1_SingleEG38er2p5"),          std::make_pair(170, "L1_SingleEG40er2p5"),          std::make_pair(171, "L1_SingleEG42er2p5"),          std::make_pair(172, "L1_SingleEG45er2p5"),          std::make_pair(173, "L1_SingleEG50"),          std::make_pair(174, "L1_SingleEG60"),          std::make_pair(159, "L1_SingleEG8er2p5"),          std::make_pair(184, "L1_SingleIsoEG24er1p5"),          std::make_pair(183, "L1_SingleIsoEG24er2p1"),          std::make_pair(187, "L1_SingleIsoEG26er1p5"),          std::make_pair(186, "L1_SingleIsoEG26er2p1"),          std::make_pair(185, "L1_SingleIsoEG26er2p5"),          std::make_pair(188, "L1_SingleIsoEG28_FWD2p5"),          std::make_pair(191, "L1_SingleIsoEG28er1p5"),          std::make_pair(190, "L1_SingleIsoEG28er2p1"),          std::make_pair(189, "L1_SingleIsoEG28er2p5"),          std::make_pair(193, "L1_SingleIsoEG30er2p1"),          std::make_pair(192, "L1_SingleIsoEG30er2p5"),          std::make_pair(195, "L1_SingleIsoEG32er2p1"),          std::make_pair(194, "L1_SingleIsoEG32er2p5"),          std::make_pair(196, "L1_SingleIsoEG34er2p5"),          std::make_pair(330, "L1_SingleJet10erHE"),          std::make_pair(312, "L1_SingleJet120"),          std::make_pair(327, "L1_SingleJet120_FWD3p0"),          std::make_pair(319, "L1_SingleJet120er2p5"),          std::make_pair(331, "L1_SingleJet12erHE"),          std::make_pair(320, "L1_SingleJet140er2p5"),          std::make_pair(332, "L1_SingleJet140er2p5_ETMHF70"),          std::make_pair(333, "L1_SingleJet140er2p5_ETMHF80"),          std::make_pair(334, "L1_SingleJet140er2p5_ETMHF90"),          std::make_pair(321, "L1_SingleJet160er2p5"),          std::make_pair(313, "L1_SingleJet180"),          std::make_pair(322, "L1_SingleJet180er2p5"),          std::make_pair(314, "L1_SingleJet200"),          std::make_pair(450, "L1_SingleJet20er2p5_NotBptxOR"),          std::make_pair(451, "L1_SingleJet20er2p5_NotBptxOR_3BX"),          std::make_pair(309, "L1_SingleJet35"),          std::make_pair(324, "L1_SingleJet35_FWD3p0"),          std::make_pair(316, "L1_SingleJet35er2p5"),          std::make_pair(452, "L1_SingleJet43er2p5_NotBptxOR_3BX"),          std::make_pair(453, "L1_SingleJet46er2p5_NotBptxOR_3BX"),          std::make_pair(310, "L1_SingleJet60"),          std::make_pair(325, "L1_SingleJet60_FWD3p0"),          std::make_pair(317, "L1_SingleJet60er2p5"),          std::make_pair(329, "L1_SingleJet8erHE"),          std::make_pair(311, "L1_SingleJet90"),          std::make_pair(326, "L1_SingleJet90_FWD3p0"),          std::make_pair(318, "L1_SingleJet90er2p5"),          std::make_pair(176, "L1_SingleLooseIsoEG26er1p5"),          std::make_pair(175, "L1_SingleLooseIsoEG26er2p5"),          std::make_pair(177, "L1_SingleLooseIsoEG28_FWD2p5"),          std::make_pair(180, "L1_SingleLooseIsoEG28er1p5"),          std::make_pair(179, "L1_SingleLooseIsoEG28er2p1"),          std::make_pair(178, "L1_SingleLooseIsoEG28er2p5"),          std::make_pair(182, "L1_SingleLooseIsoEG30er1p5"),          std::make_pair(181, "L1_SingleLooseIsoEG30er2p5"),          std::make_pair(6, "L1_SingleMu0_BMTF"),          std::make_pair(5, "L1_SingleMu0_DQ"),          std::make_pair(8, "L1_SingleMu0_EMTF"),          std::make_pair(7, "L1_SingleMu0_OMTF"),          std::make_pair(29, "L1_SingleMu10er1p5"),          std::make_pair(13, "L1_SingleMu12_DQ_BMTF"),          std::make_pair(15, "L1_SingleMu12_DQ_EMTF"),          std::make_pair(14, "L1_SingleMu12_DQ_OMTF"),          std::make_pair(30, "L1_SingleMu12er1p5"),          std::make_pair(31, "L1_SingleMu14er1p5"),          std::make_pair(16, "L1_SingleMu15_DQ"),          std::make_pair(32, "L1_SingleMu16er1p5"),          std::make_pair(17, "L1_SingleMu18"),          std::make_pair(33, "L1_SingleMu18er1p5"),          std::make_pair(18, "L1_SingleMu20"),          std::make_pair(19, "L1_SingleMu22"),          std::make_pair(20, "L1_SingleMu22_BMTF"),          std::make_pair(22, "L1_SingleMu22_EMTF"),          std::make_pair(21, "L1_SingleMu22_OMTF"),          std::make_pair(23, "L1_SingleMu25"),          std::make_pair(9, "L1_SingleMu3"),          std::make_pair(10, "L1_SingleMu5"),          std::make_pair(25, "L1_SingleMu6er1p5"),          std::make_pair(12, "L1_SingleMu7"),          std::make_pair(11, "L1_SingleMu7_DQ"),          std::make_pair(26, "L1_SingleMu7er1p5"),          std::make_pair(27, "L1_SingleMu8er1p5"),          std::make_pair(28, "L1_SingleMu9er1p5"),          std::make_pair(0, "L1_SingleMuCosmics"),          std::make_pair(1, "L1_SingleMuCosmics_BMTF"),          std::make_pair(3, "L1_SingleMuCosmics_EMTF"),          std::make_pair(2, "L1_SingleMuCosmics_OMTF"),          std::make_pair(4, "L1_SingleMuOpen"),          std::make_pair(446, "L1_SingleMuOpen_NotBptxOR"),          std::make_pair(448, "L1_SingleMuOpen_er1p1_NotBptxOR_3BX"),          std::make_pair(447, "L1_SingleMuOpen_er1p4_NotBptxOR_3BX"),          std::make_pair(264, "L1_SingleTau120er2p1"),          std::make_pair(265, "L1_SingleTau130er2p1"),          std::make_pair(503, "L1_TOTEM_1"),          std::make_pair(504, "L1_TOTEM_2"),          std::make_pair(505, "L1_TOTEM_3"),          std::make_pair(506, "L1_TOTEM_4"),          std::make_pair(228, "L1_TripleEG16er2p5"),          std::make_pair(224, "L1_TripleEG_16_12_8_er2p5"),          std::make_pair(225, "L1_TripleEG_16_15_8_er2p5"),          std::make_pair(226, "L1_TripleEG_18_17_8_er2p5"),          std::make_pair(227, "L1_TripleEG_18_18_12_er2p5"),          std::make_pair(373, "L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5"),          std::make_pair(374, "L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5"),          std::make_pair(372, "L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5"),          std::make_pair(72, "L1_TripleMu0"),          std::make_pair(71, "L1_TripleMu0_OQ"),          std::make_pair(73, "L1_TripleMu0_SQ"),          std::make_pair(74, "L1_TripleMu3"),          std::make_pair(75, "L1_TripleMu3_SQ"),          std::make_pair(76, "L1_TripleMu_5SQ_3SQ_0OQ"),          std::make_pair(85, "L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9"),          std::make_pair(86, "L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"),          std::make_pair(78, "L1_TripleMu_5_3_3"),          std::make_pair(79, "L1_TripleMu_5_3_3_SQ"),          std::make_pair(77, "L1_TripleMu_5_3p5_2p5"),          std::make_pair(83, "L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17"),          std::make_pair(82, "L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17"),          std::make_pair(84, "L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17"),          std::make_pair(80, "L1_TripleMu_5_5_3"),          std::make_pair(469, "L1_UnpairedBunchBptxMinus"),          std::make_pair(468, "L1_UnpairedBunchBptxPlus"),          std::make_pair(459, "L1_ZeroBias"),          std::make_pair(460, "L1_ZeroBias_copy")      };

  static const std::map<int, std::string> Id2Name(id2name, id2name + sizeof(id2name) / sizeof(id2name[0]));
  const std::map<int, std::string>::const_iterator rc = Id2Name.find(index);
  std::string name;
  if (rc != Id2Name.end()) name = rc->second;
  return name;
}


int getIdFromName(const std::string& name)
{
  static const std::pair<std::string, int> name2id[] = {
          std::make_pair("L1_AlwaysTrue", 458),          std::make_pair("L1_BPTX_AND_Ref1_VME", 486),          std::make_pair("L1_BPTX_AND_Ref3_VME", 487),          std::make_pair("L1_BPTX_AND_Ref4_VME", 488),          std::make_pair("L1_BPTX_BeamGas_B1_VME", 491),          std::make_pair("L1_BPTX_BeamGas_B2_VME", 492),          std::make_pair("L1_BPTX_BeamGas_Ref1_VME", 489),          std::make_pair("L1_BPTX_BeamGas_Ref2_VME", 490),          std::make_pair("L1_BPTX_NotOR_VME", 482),          std::make_pair("L1_BPTX_OR_Ref3_VME", 483),          std::make_pair("L1_BPTX_OR_Ref4_VME", 484),          std::make_pair("L1_BPTX_RefAND_VME", 485),          std::make_pair("L1_BptxMinus", 467),          std::make_pair("L1_BptxOR", 464),          std::make_pair("L1_BptxPlus", 466),          std::make_pair("L1_BptxXOR", 465),          std::make_pair("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", 494),          std::make_pair("L1_DoubleEG8er2p5_HTT260er", 247),          std::make_pair("L1_DoubleEG8er2p5_HTT280er", 248),          std::make_pair("L1_DoubleEG8er2p5_HTT300er", 249),          std::make_pair("L1_DoubleEG8er2p5_HTT320er", 250),          std::make_pair("L1_DoubleEG8er2p5_HTT340er", 251),          std::make_pair("L1_DoubleEG_15_10_er2p5", 205),          std::make_pair("L1_DoubleEG_20_10_er2p5", 206),          std::make_pair("L1_DoubleEG_22_10_er2p5", 207),          std::make_pair("L1_DoubleEG_25_12_er2p5", 208),          std::make_pair("L1_DoubleEG_25_14_er2p5", 209),          std::make_pair("L1_DoubleEG_27_14_er2p5", 210),          std::make_pair("L1_DoubleEG_LooseIso20_10_er2p5", 212),          std::make_pair("L1_DoubleEG_LooseIso22_10_er2p5", 213),          std::make_pair("L1_DoubleEG_LooseIso22_12_er2p5", 214),          std::make_pair("L1_DoubleEG_LooseIso25_12_er2p5", 215),          std::make_pair("L1_DoubleIsoTau28er2p1", 269),          std::make_pair("L1_DoubleIsoTau28er2p1_Mass_Max80", 275),          std::make_pair("L1_DoubleIsoTau28er2p1_Mass_Max90", 274),          std::make_pair("L1_DoubleIsoTau30er2p1", 270),          std::make_pair("L1_DoubleIsoTau30er2p1_Mass_Max80", 277),          std::make_pair("L1_DoubleIsoTau30er2p1_Mass_Max90", 276),          std::make_pair("L1_DoubleIsoTau32er2p1", 271),          std::make_pair("L1_DoubleIsoTau34er2p1", 272),          std::make_pair("L1_DoubleIsoTau36er2p1", 273),          std::make_pair("L1_DoubleJet100er2p3_dEta_Max1p6", 345),          std::make_pair("L1_DoubleJet100er2p5", 341),          std::make_pair("L1_DoubleJet112er2p3_dEta_Max1p6", 346),          std::make_pair("L1_DoubleJet120er2p5", 342),          std::make_pair("L1_DoubleJet150er2p5", 343),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", 348),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", 349),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", 350),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", 351),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", 352),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", 353),          std::make_pair("L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", 363),          std::make_pair("L1_DoubleJet40er2p5", 340),          std::make_pair("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", 356),          std::make_pair("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", 357),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", 358),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", 360),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", 359),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", 361),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", 366),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", 364),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_Mu8", 365),          std::make_pair("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", 355),          std::make_pair("L1_DoubleLooseIsoEG22er2p1", 217),          std::make_pair("L1_DoubleLooseIsoEG24er2p1", 218),          std::make_pair("L1_DoubleMu0", 40),          std::make_pair("L1_DoubleMu0_Mass_Min1", 43),          std::make_pair("L1_DoubleMu0_OQ", 39),          std::make_pair("L1_DoubleMu0_SQ", 41),          std::make_pair("L1_DoubleMu0_SQ_OS", 42),          std::make_pair("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", 142),          std::make_pair("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", 59),          std::make_pair("L1_DoubleMu0er1p5_SQ", 55),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS", 56),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", 58),          std::make_pair("L1_DoubleMu0er1p5_SQ_dR_Max1p4", 57),          std::make_pair("L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", 54),          std::make_pair("L1_DoubleMu0er2p0_SQ_dR_Max1p4", 53),          std::make_pair("L1_DoubleMu18er2p1", 51),          std::make_pair("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", 112),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_HTT60er", 145),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", 147),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", 146),          std::make_pair("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", 148),          std::make_pair("L1_DoubleMu3_SQ_HTT220er", 150),          std::make_pair("L1_DoubleMu3_SQ_HTT240er", 151),          std::make_pair("L1_DoubleMu3_SQ_HTT260er", 152),          std::make_pair("L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", 143),          std::make_pair("L1_DoubleMu4_SQ_EG9er2p5", 109),          std::make_pair("L1_DoubleMu4_SQ_OS", 60),          std::make_pair("L1_DoubleMu4_SQ_OS_dR_Max1p2", 61),          std::make_pair("L1_DoubleMu4p5_SQ_OS", 62),          std::make_pair("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", 63),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS", 64),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", 66),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", 65),          std::make_pair("L1_DoubleMu5Upsilon_OS_DoubleEG3", 113),          std::make_pair("L1_DoubleMu5_SQ_EG9er2p5", 110),          std::make_pair("L1_DoubleMu8_SQ", 44),          std::make_pair("L1_DoubleMu9_SQ", 45),          std::make_pair("L1_DoubleMu_12_5", 46),          std::make_pair("L1_DoubleMu_15_5_SQ", 47),          std::make_pair("L1_DoubleMu_15_7", 48),          std::make_pair("L1_DoubleMu_15_7_Mass_Min1", 50),          std::make_pair("L1_DoubleMu_15_7_SQ", 49),          std::make_pair("L1_DoubleTau70er2p1", 267),          std::make_pair("L1_ETM120", 416),          std::make_pair("L1_ETM150", 417),          std::make_pair("L1_ETMHF100", 421),          std::make_pair("L1_ETMHF100_HTT60er", 429),          std::make_pair("L1_ETMHF110", 422),          std::make_pair("L1_ETMHF110_HTT60er", 430),          std::make_pair("L1_ETMHF110_HTT60er_NotSecondBunchInTrain", 444),          std::make_pair("L1_ETMHF120", 423),          std::make_pair("L1_ETMHF120_HTT60er", 431),          std::make_pair("L1_ETMHF120_NotSecondBunchInTrain", 443),          std::make_pair("L1_ETMHF130", 424),          std::make_pair("L1_ETMHF130_HTT60er", 432),          std::make_pair("L1_ETMHF140", 425),          std::make_pair("L1_ETMHF150", 426),          std::make_pair("L1_ETMHF90_HTT60er", 428),          std::make_pair("L1_ETT1200", 410),          std::make_pair("L1_ETT1600", 411),          std::make_pair("L1_ETT2000", 412),          std::make_pair("L1_FirstBunchAfterTrain", 477),          std::make_pair("L1_FirstBunchBeforeTrain", 472),          std::make_pair("L1_FirstBunchInTrain", 473),          std::make_pair("L1_FirstCollisionInOrbit", 480),          std::make_pair("L1_FirstCollisionInTrain", 479),          std::make_pair("L1_HCAL_LaserMon_Trig", 500),          std::make_pair("L1_HCAL_LaserMon_Veto", 501),          std::make_pair("L1_HTT120er", 398),          std::make_pair("L1_HTT160er", 399),          std::make_pair("L1_HTT200er", 400),          std::make_pair("L1_HTT255er", 401),          std::make_pair("L1_HTT280er", 402),          std::make_pair("L1_HTT280er_QuadJet_70_55_40_35_er2p4", 384),          std::make_pair("L1_HTT320er", 403),          std::make_pair("L1_HTT320er_QuadJet_70_55_40_40_er2p4", 385),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", 386),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", 387),          std::make_pair("L1_HTT360er", 404),          std::make_pair("L1_HTT400er", 405),          std::make_pair("L1_HTT450er", 406),          std::make_pair("L1_IsoEG32er2p5_Mt40", 197),          std::make_pair("L1_IsoEG32er2p5_Mt44", 198),          std::make_pair("L1_IsoEG32er2p5_Mt48", 199),          std::make_pair("L1_IsoTau40er2p1_ETMHF100", 293),          std::make_pair("L1_IsoTau40er2p1_ETMHF110", 294),          std::make_pair("L1_IsoTau40er2p1_ETMHF80", 291),          std::make_pair("L1_IsoTau40er2p1_ETMHF90", 292),          std::make_pair("L1_IsolatedBunch", 471),          std::make_pair("L1_LastBunchInTrain", 476),          std::make_pair("L1_LastCollisionInTrain", 478),          std::make_pair("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", 257),          std::make_pair("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", 259),          std::make_pair("L1_LooseIsoEG24er2p1_HTT100er", 238),          std::make_pair("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", 258),          std::make_pair("L1_LooseIsoEG26er2p1_HTT100er", 239),          std::make_pair("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", 234),          std::make_pair("L1_LooseIsoEG28er2p1_HTT100er", 240),          std::make_pair("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", 235),          std::make_pair("L1_LooseIsoEG30er2p1_HTT100er", 241),          std::make_pair("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", 236),          std::make_pair("L1_MinimumBiasHF0_AND_BptxAND", 461),          std::make_pair("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", 134),          std::make_pair("L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", 136),          std::make_pair("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", 135),          std::make_pair("L1_Mu18er2p1_Tau24er2p1", 279),          std::make_pair("L1_Mu18er2p1_Tau26er2p1", 280),          std::make_pair("L1_Mu20_EG10er2p5", 99),          std::make_pair("L1_Mu22er2p1_IsoTau28er2p1", 282),          std::make_pair("L1_Mu22er2p1_IsoTau30er2p1", 283),          std::make_pair("L1_Mu22er2p1_IsoTau32er2p1", 284),          std::make_pair("L1_Mu22er2p1_IsoTau34er2p1", 285),          std::make_pair("L1_Mu22er2p1_IsoTau36er2p1", 286),          std::make_pair("L1_Mu22er2p1_IsoTau40er2p1", 287),          std::make_pair("L1_Mu22er2p1_Tau70er2p1", 289),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p4", 126),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p8", 125),          std::make_pair("L1_Mu3_Jet16er2p5_dR_Max0p4", 121),          std::make_pair("L1_Mu3_Jet30er2p5", 119),          std::make_pair("L1_Mu3_Jet35er2p5_dR_Max0p4", 122),          std::make_pair("L1_Mu3_Jet60er2p5_dR_Max0p4", 123),          std::make_pair("L1_Mu3_Jet80er2p5_dR_Max0p4", 124),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF40", 128),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF50", 129),          std::make_pair("L1_Mu5_EG23er2p5", 96),          std::make_pair("L1_Mu5_LooseIsoEG20er2p5", 100),          std::make_pair("L1_Mu6_DoubleEG10er2p5", 104),          std::make_pair("L1_Mu6_DoubleEG12er2p5", 105),          std::make_pair("L1_Mu6_DoubleEG15er2p5", 106),          std::make_pair("L1_Mu6_DoubleEG17er2p5", 107),          std::make_pair("L1_Mu6_HTT240er", 131),          std::make_pair("L1_Mu6_HTT250er", 132),          std::make_pair("L1_Mu7_EG20er2p5", 97),          std::make_pair("L1_Mu7_EG23er2p5", 98),          std::make_pair("L1_Mu7_LooseIsoEG20er2p5", 101),          std::make_pair("L1_Mu7_LooseIsoEG23er2p5", 102),          std::make_pair("L1_NotBptxOR", 463),          std::make_pair("L1_QuadJet36er2p5_IsoTau52er2p1", 298),          std::make_pair("L1_QuadJet60er2p5", 382),          std::make_pair("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", 376),          std::make_pair("L1_QuadMu0", 89),          std::make_pair("L1_QuadMu0_OQ", 88),          std::make_pair("L1_QuadMu0_SQ", 90),          std::make_pair("L1_SecondBunchInTrain", 474),          std::make_pair("L1_SecondLastBunchInTrain", 475),          std::make_pair("L1_SingleEG10er2p5", 160),          std::make_pair("L1_SingleEG15er2p5", 161),          std::make_pair("L1_SingleEG26er2p5", 162),          std::make_pair("L1_SingleEG28_FWD2p5", 163),          std::make_pair("L1_SingleEG28er1p5", 166),          std::make_pair("L1_SingleEG28er2p1", 165),          std::make_pair("L1_SingleEG28er2p5", 164),          std::make_pair("L1_SingleEG34er2p5", 167),          std::make_pair("L1_SingleEG36er2p5", 168),          std::make_pair("L1_SingleEG38er2p5", 169),          std::make_pair("L1_SingleEG40er2p5", 170),          std::make_pair("L1_SingleEG42er2p5", 171),          std::make_pair("L1_SingleEG45er2p5", 172),          std::make_pair("L1_SingleEG50", 173),          std::make_pair("L1_SingleEG60", 174),          std::make_pair("L1_SingleEG8er2p5", 159),          std::make_pair("L1_SingleIsoEG24er1p5", 184),          std::make_pair("L1_SingleIsoEG24er2p1", 183),          std::make_pair("L1_SingleIsoEG26er1p5", 187),          std::make_pair("L1_SingleIsoEG26er2p1", 186),          std::make_pair("L1_SingleIsoEG26er2p5", 185),          std::make_pair("L1_SingleIsoEG28_FWD2p5", 188),          std::make_pair("L1_SingleIsoEG28er1p5", 191),          std::make_pair("L1_SingleIsoEG28er2p1", 190),          std::make_pair("L1_SingleIsoEG28er2p5", 189),          std::make_pair("L1_SingleIsoEG30er2p1", 193),          std::make_pair("L1_SingleIsoEG30er2p5", 192),          std::make_pair("L1_SingleIsoEG32er2p1", 195),          std::make_pair("L1_SingleIsoEG32er2p5", 194),          std::make_pair("L1_SingleIsoEG34er2p5", 196),          std::make_pair("L1_SingleJet10erHE", 330),          std::make_pair("L1_SingleJet120", 312),          std::make_pair("L1_SingleJet120_FWD3p0", 327),          std::make_pair("L1_SingleJet120er2p5", 319),          std::make_pair("L1_SingleJet12erHE", 331),          std::make_pair("L1_SingleJet140er2p5", 320),          std::make_pair("L1_SingleJet140er2p5_ETMHF70", 332),          std::make_pair("L1_SingleJet140er2p5_ETMHF80", 333),          std::make_pair("L1_SingleJet140er2p5_ETMHF90", 334),          std::make_pair("L1_SingleJet160er2p5", 321),          std::make_pair("L1_SingleJet180", 313),          std::make_pair("L1_SingleJet180er2p5", 322),          std::make_pair("L1_SingleJet200", 314),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR", 450),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR_3BX", 451),          std::make_pair("L1_SingleJet35", 309),          std::make_pair("L1_SingleJet35_FWD3p0", 324),          std::make_pair("L1_SingleJet35er2p5", 316),          std::make_pair("L1_SingleJet43er2p5_NotBptxOR_3BX", 452),          std::make_pair("L1_SingleJet46er2p5_NotBptxOR_3BX", 453),          std::make_pair("L1_SingleJet60", 310),          std::make_pair("L1_SingleJet60_FWD3p0", 325),          std::make_pair("L1_SingleJet60er2p5", 317),          std::make_pair("L1_SingleJet8erHE", 329),          std::make_pair("L1_SingleJet90", 311),          std::make_pair("L1_SingleJet90_FWD3p0", 326),          std::make_pair("L1_SingleJet90er2p5", 318),          std::make_pair("L1_SingleLooseIsoEG26er1p5", 176),          std::make_pair("L1_SingleLooseIsoEG26er2p5", 175),          std::make_pair("L1_SingleLooseIsoEG28_FWD2p5", 177),          std::make_pair("L1_SingleLooseIsoEG28er1p5", 180),          std::make_pair("L1_SingleLooseIsoEG28er2p1", 179),          std::make_pair("L1_SingleLooseIsoEG28er2p5", 178),          std::make_pair("L1_SingleLooseIsoEG30er1p5", 182),          std::make_pair("L1_SingleLooseIsoEG30er2p5", 181),          std::make_pair("L1_SingleMu0_BMTF", 6),          std::make_pair("L1_SingleMu0_DQ", 5),          std::make_pair("L1_SingleMu0_EMTF", 8),          std::make_pair("L1_SingleMu0_OMTF", 7),          std::make_pair("L1_SingleMu10er1p5", 29),          std::make_pair("L1_SingleMu12_DQ_BMTF", 13),          std::make_pair("L1_SingleMu12_DQ_EMTF", 15),          std::make_pair("L1_SingleMu12_DQ_OMTF", 14),          std::make_pair("L1_SingleMu12er1p5", 30),          std::make_pair("L1_SingleMu14er1p5", 31),          std::make_pair("L1_SingleMu15_DQ", 16),          std::make_pair("L1_SingleMu16er1p5", 32),          std::make_pair("L1_SingleMu18", 17),          std::make_pair("L1_SingleMu18er1p5", 33),          std::make_pair("L1_SingleMu20", 18),          std::make_pair("L1_SingleMu22", 19),          std::make_pair("L1_SingleMu22_BMTF", 20),          std::make_pair("L1_SingleMu22_EMTF", 22),          std::make_pair("L1_SingleMu22_OMTF", 21),          std::make_pair("L1_SingleMu25", 23),          std::make_pair("L1_SingleMu3", 9),          std::make_pair("L1_SingleMu5", 10),          std::make_pair("L1_SingleMu6er1p5", 25),          std::make_pair("L1_SingleMu7", 12),          std::make_pair("L1_SingleMu7_DQ", 11),          std::make_pair("L1_SingleMu7er1p5", 26),          std::make_pair("L1_SingleMu8er1p5", 27),          std::make_pair("L1_SingleMu9er1p5", 28),          std::make_pair("L1_SingleMuCosmics", 0),          std::make_pair("L1_SingleMuCosmics_BMTF", 1),          std::make_pair("L1_SingleMuCosmics_EMTF", 3),          std::make_pair("L1_SingleMuCosmics_OMTF", 2),          std::make_pair("L1_SingleMuOpen", 4),          std::make_pair("L1_SingleMuOpen_NotBptxOR", 446),          std::make_pair("L1_SingleMuOpen_er1p1_NotBptxOR_3BX", 448),          std::make_pair("L1_SingleMuOpen_er1p4_NotBptxOR_3BX", 447),          std::make_pair("L1_SingleTau120er2p1", 264),          std::make_pair("L1_SingleTau130er2p1", 265),          std::make_pair("L1_TOTEM_1", 503),          std::make_pair("L1_TOTEM_2", 504),          std::make_pair("L1_TOTEM_3", 505),          std::make_pair("L1_TOTEM_4", 506),          std::make_pair("L1_TripleEG16er2p5", 228),          std::make_pair("L1_TripleEG_16_12_8_er2p5", 224),          std::make_pair("L1_TripleEG_16_15_8_er2p5", 225),          std::make_pair("L1_TripleEG_18_17_8_er2p5", 226),          std::make_pair("L1_TripleEG_18_18_12_er2p5", 227),          std::make_pair("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", 373),          std::make_pair("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", 374),          std::make_pair("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", 372),          std::make_pair("L1_TripleMu0", 72),          std::make_pair("L1_TripleMu0_OQ", 71),          std::make_pair("L1_TripleMu0_SQ", 73),          std::make_pair("L1_TripleMu3", 74),          std::make_pair("L1_TripleMu3_SQ", 75),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ", 76),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", 85),          std::make_pair("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", 86),          std::make_pair("L1_TripleMu_5_3_3", 78),          std::make_pair("L1_TripleMu_5_3_3_SQ", 79),          std::make_pair("L1_TripleMu_5_3p5_2p5", 77),          std::make_pair("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", 83),          std::make_pair("L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", 82),          std::make_pair("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", 84),          std::make_pair("L1_TripleMu_5_5_3", 80),          std::make_pair("L1_UnpairedBunchBptxMinus", 469),          std::make_pair("L1_UnpairedBunchBptxPlus", 468),          std::make_pair("L1_ZeroBias", 459),          std::make_pair("L1_ZeroBias_copy", 460)      };

  static const std::map<std::string, int> Name2Id(name2id, name2id + sizeof(name2id) / sizeof(name2id[0]));
  const std::map<std::string, int>::const_iterator rc = Name2Id.find(name);
  int id = -1;
  if (rc != Name2Id.end()) id = rc->second;
  return id;
}


AlgorithmFunction getFuncFromId(const int index)
{
  static const std::pair<int, AlgorithmFunction> id2func[] = {
          std::make_pair(458, &L1_AlwaysTrue),          std::make_pair(486, &L1_BPTX_AND_Ref1_VME),          std::make_pair(487, &L1_BPTX_AND_Ref3_VME),          std::make_pair(488, &L1_BPTX_AND_Ref4_VME),          std::make_pair(491, &L1_BPTX_BeamGas_B1_VME),          std::make_pair(492, &L1_BPTX_BeamGas_B2_VME),          std::make_pair(489, &L1_BPTX_BeamGas_Ref1_VME),          std::make_pair(490, &L1_BPTX_BeamGas_Ref2_VME),          std::make_pair(482, &L1_BPTX_NotOR_VME),          std::make_pair(483, &L1_BPTX_OR_Ref3_VME),          std::make_pair(484, &L1_BPTX_OR_Ref4_VME),          std::make_pair(485, &L1_BPTX_RefAND_VME),          std::make_pair(467, &L1_BptxMinus),          std::make_pair(464, &L1_BptxOR),          std::make_pair(466, &L1_BptxPlus),          std::make_pair(465, &L1_BptxXOR),          std::make_pair(494, &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142),          std::make_pair(247, &L1_DoubleEG8er2p5_HTT260er),          std::make_pair(248, &L1_DoubleEG8er2p5_HTT280er),          std::make_pair(249, &L1_DoubleEG8er2p5_HTT300er),          std::make_pair(250, &L1_DoubleEG8er2p5_HTT320er),          std::make_pair(251, &L1_DoubleEG8er2p5_HTT340er),          std::make_pair(205, &L1_DoubleEG_15_10_er2p5),          std::make_pair(206, &L1_DoubleEG_20_10_er2p5),          std::make_pair(207, &L1_DoubleEG_22_10_er2p5),          std::make_pair(208, &L1_DoubleEG_25_12_er2p5),          std::make_pair(209, &L1_DoubleEG_25_14_er2p5),          std::make_pair(210, &L1_DoubleEG_27_14_er2p5),          std::make_pair(212, &L1_DoubleEG_LooseIso20_10_er2p5),          std::make_pair(213, &L1_DoubleEG_LooseIso22_10_er2p5),          std::make_pair(214, &L1_DoubleEG_LooseIso22_12_er2p5),          std::make_pair(215, &L1_DoubleEG_LooseIso25_12_er2p5),          std::make_pair(269, &L1_DoubleIsoTau28er2p1),          std::make_pair(275, &L1_DoubleIsoTau28er2p1_Mass_Max80),          std::make_pair(274, &L1_DoubleIsoTau28er2p1_Mass_Max90),          std::make_pair(270, &L1_DoubleIsoTau30er2p1),          std::make_pair(277, &L1_DoubleIsoTau30er2p1_Mass_Max80),          std::make_pair(276, &L1_DoubleIsoTau30er2p1_Mass_Max90),          std::make_pair(271, &L1_DoubleIsoTau32er2p1),          std::make_pair(272, &L1_DoubleIsoTau34er2p1),          std::make_pair(273, &L1_DoubleIsoTau36er2p1),          std::make_pair(345, &L1_DoubleJet100er2p3_dEta_Max1p6),          std::make_pair(341, &L1_DoubleJet100er2p5),          std::make_pair(346, &L1_DoubleJet112er2p3_dEta_Max1p6),          std::make_pair(342, &L1_DoubleJet120er2p5),          std::make_pair(343, &L1_DoubleJet150er2p5),          std::make_pair(348, &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5),          std::make_pair(349, &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5),          std::make_pair(350, &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5),          std::make_pair(351, &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5),          std::make_pair(352, &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5),          std::make_pair(353, &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5),          std::make_pair(363, &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp),          std::make_pair(340, &L1_DoubleJet40er2p5),          std::make_pair(356, &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620),          std::make_pair(357, &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620),          std::make_pair(358, &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620),          std::make_pair(360, &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28),          std::make_pair(359, &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620),          std::make_pair(361, &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28),          std::make_pair(366, &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ),          std::make_pair(364, &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp),          std::make_pair(365, &L1_DoubleJet_80_30_Mass_Min420_Mu8),          std::make_pair(355, &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620),          std::make_pair(217, &L1_DoubleLooseIsoEG22er2p1),          std::make_pair(218, &L1_DoubleLooseIsoEG24er2p1),          std::make_pair(40, &L1_DoubleMu0),          std::make_pair(43, &L1_DoubleMu0_Mass_Min1),          std::make_pair(39, &L1_DoubleMu0_OQ),          std::make_pair(41, &L1_DoubleMu0_SQ),          std::make_pair(42, &L1_DoubleMu0_SQ_OS),          std::make_pair(142, &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair(59, &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4),          std::make_pair(55, &L1_DoubleMu0er1p5_SQ),          std::make_pair(56, &L1_DoubleMu0er1p5_SQ_OS),          std::make_pair(58, &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4),          std::make_pair(57, &L1_DoubleMu0er1p5_SQ_dR_Max1p4),          std::make_pair(54, &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4),          std::make_pair(53, &L1_DoubleMu0er2p0_SQ_dR_Max1p4),          std::make_pair(51, &L1_DoubleMu18er2p1),          std::make_pair(112, &L1_DoubleMu3_OS_DoubleEG7p5Upsilon),          std::make_pair(145, &L1_DoubleMu3_SQ_ETMHF50_HTT60er),          std::make_pair(147, &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5),          std::make_pair(146, &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5),          std::make_pair(148, &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5),          std::make_pair(150, &L1_DoubleMu3_SQ_HTT220er),          std::make_pair(151, &L1_DoubleMu3_SQ_HTT240er),          std::make_pair(152, &L1_DoubleMu3_SQ_HTT260er),          std::make_pair(143, &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair(109, &L1_DoubleMu4_SQ_EG9er2p5),          std::make_pair(60, &L1_DoubleMu4_SQ_OS),          std::make_pair(61, &L1_DoubleMu4_SQ_OS_dR_Max1p2),          std::make_pair(62, &L1_DoubleMu4p5_SQ_OS),          std::make_pair(63, &L1_DoubleMu4p5_SQ_OS_dR_Max1p2),          std::make_pair(64, &L1_DoubleMu4p5er2p0_SQ_OS),          std::make_pair(66, &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18),          std::make_pair(65, &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7),          std::make_pair(113, &L1_DoubleMu5Upsilon_OS_DoubleEG3),          std::make_pair(110, &L1_DoubleMu5_SQ_EG9er2p5),          std::make_pair(44, &L1_DoubleMu8_SQ),          std::make_pair(45, &L1_DoubleMu9_SQ),          std::make_pair(46, &L1_DoubleMu_12_5),          std::make_pair(47, &L1_DoubleMu_15_5_SQ),          std::make_pair(48, &L1_DoubleMu_15_7),          std::make_pair(50, &L1_DoubleMu_15_7_Mass_Min1),          std::make_pair(49, &L1_DoubleMu_15_7_SQ),          std::make_pair(267, &L1_DoubleTau70er2p1),          std::make_pair(416, &L1_ETM120),          std::make_pair(417, &L1_ETM150),          std::make_pair(421, &L1_ETMHF100),          std::make_pair(429, &L1_ETMHF100_HTT60er),          std::make_pair(422, &L1_ETMHF110),          std::make_pair(430, &L1_ETMHF110_HTT60er),          std::make_pair(444, &L1_ETMHF110_HTT60er_NotSecondBunchInTrain),          std::make_pair(423, &L1_ETMHF120),          std::make_pair(431, &L1_ETMHF120_HTT60er),          std::make_pair(443, &L1_ETMHF120_NotSecondBunchInTrain),          std::make_pair(424, &L1_ETMHF130),          std::make_pair(432, &L1_ETMHF130_HTT60er),          std::make_pair(425, &L1_ETMHF140),          std::make_pair(426, &L1_ETMHF150),          std::make_pair(428, &L1_ETMHF90_HTT60er),          std::make_pair(410, &L1_ETT1200),          std::make_pair(411, &L1_ETT1600),          std::make_pair(412, &L1_ETT2000),          std::make_pair(477, &L1_FirstBunchAfterTrain),          std::make_pair(472, &L1_FirstBunchBeforeTrain),          std::make_pair(473, &L1_FirstBunchInTrain),          std::make_pair(480, &L1_FirstCollisionInOrbit),          std::make_pair(479, &L1_FirstCollisionInTrain),          std::make_pair(500, &L1_HCAL_LaserMon_Trig),          std::make_pair(501, &L1_HCAL_LaserMon_Veto),          std::make_pair(398, &L1_HTT120er),          std::make_pair(399, &L1_HTT160er),          std::make_pair(400, &L1_HTT200er),          std::make_pair(401, &L1_HTT255er),          std::make_pair(402, &L1_HTT280er),          std::make_pair(384, &L1_HTT280er_QuadJet_70_55_40_35_er2p4),          std::make_pair(403, &L1_HTT320er),          std::make_pair(385, &L1_HTT320er_QuadJet_70_55_40_40_er2p4),          std::make_pair(386, &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3),          std::make_pair(387, &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3),          std::make_pair(404, &L1_HTT360er),          std::make_pair(405, &L1_HTT400er),          std::make_pair(406, &L1_HTT450er),          std::make_pair(197, &L1_IsoEG32er2p5_Mt40),          std::make_pair(198, &L1_IsoEG32er2p5_Mt44),          std::make_pair(199, &L1_IsoEG32er2p5_Mt48),          std::make_pair(293, &L1_IsoTau40er2p1_ETMHF100),          std::make_pair(294, &L1_IsoTau40er2p1_ETMHF110),          std::make_pair(291, &L1_IsoTau40er2p1_ETMHF80),          std::make_pair(292, &L1_IsoTau40er2p1_ETMHF90),          std::make_pair(471, &L1_IsolatedBunch),          std::make_pair(476, &L1_LastBunchInTrain),          std::make_pair(478, &L1_LastCollisionInTrain),          std::make_pair(257, &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3),          std::make_pair(259, &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3),          std::make_pair(238, &L1_LooseIsoEG24er2p1_HTT100er),          std::make_pair(258, &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3),          std::make_pair(239, &L1_LooseIsoEG26er2p1_HTT100er),          std::make_pair(234, &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair(240, &L1_LooseIsoEG28er2p1_HTT100er),          std::make_pair(235, &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair(241, &L1_LooseIsoEG30er2p1_HTT100er),          std::make_pair(236, &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair(461, &L1_MinimumBiasHF0_AND_BptxAND),          std::make_pair(134, &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6),          std::make_pair(136, &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6),          std::make_pair(135, &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6),          std::make_pair(279, &L1_Mu18er2p1_Tau24er2p1),          std::make_pair(280, &L1_Mu18er2p1_Tau26er2p1),          std::make_pair(99, &L1_Mu20_EG10er2p5),          std::make_pair(282, &L1_Mu22er2p1_IsoTau28er2p1),          std::make_pair(283, &L1_Mu22er2p1_IsoTau30er2p1),          std::make_pair(284, &L1_Mu22er2p1_IsoTau32er2p1),          std::make_pair(285, &L1_Mu22er2p1_IsoTau34er2p1),          std::make_pair(286, &L1_Mu22er2p1_IsoTau36er2p1),          std::make_pair(287, &L1_Mu22er2p1_IsoTau40er2p1),          std::make_pair(289, &L1_Mu22er2p1_Tau70er2p1),          std::make_pair(126, &L1_Mu3_Jet120er2p5_dR_Max0p4),          std::make_pair(125, &L1_Mu3_Jet120er2p5_dR_Max0p8),          std::make_pair(121, &L1_Mu3_Jet16er2p5_dR_Max0p4),          std::make_pair(119, &L1_Mu3_Jet30er2p5),          std::make_pair(122, &L1_Mu3_Jet35er2p5_dR_Max0p4),          std::make_pair(123, &L1_Mu3_Jet60er2p5_dR_Max0p4),          std::make_pair(124, &L1_Mu3_Jet80er2p5_dR_Max0p4),          std::make_pair(128, &L1_Mu3er1p5_Jet100er2p5_ETMHF40),          std::make_pair(129, &L1_Mu3er1p5_Jet100er2p5_ETMHF50),          std::make_pair(96, &L1_Mu5_EG23er2p5),          std::make_pair(100, &L1_Mu5_LooseIsoEG20er2p5),          std::make_pair(104, &L1_Mu6_DoubleEG10er2p5),          std::make_pair(105, &L1_Mu6_DoubleEG12er2p5),          std::make_pair(106, &L1_Mu6_DoubleEG15er2p5),          std::make_pair(107, &L1_Mu6_DoubleEG17er2p5),          std::make_pair(131, &L1_Mu6_HTT240er),          std::make_pair(132, &L1_Mu6_HTT250er),          std::make_pair(97, &L1_Mu7_EG20er2p5),          std::make_pair(98, &L1_Mu7_EG23er2p5),          std::make_pair(101, &L1_Mu7_LooseIsoEG20er2p5),          std::make_pair(102, &L1_Mu7_LooseIsoEG23er2p5),          std::make_pair(463, &L1_NotBptxOR),          std::make_pair(298, &L1_QuadJet36er2p5_IsoTau52er2p1),          std::make_pair(382, &L1_QuadJet60er2p5),          std::make_pair(376, &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0),          std::make_pair(89, &L1_QuadMu0),          std::make_pair(88, &L1_QuadMu0_OQ),          std::make_pair(90, &L1_QuadMu0_SQ),          std::make_pair(474, &L1_SecondBunchInTrain),          std::make_pair(475, &L1_SecondLastBunchInTrain),          std::make_pair(160, &L1_SingleEG10er2p5),          std::make_pair(161, &L1_SingleEG15er2p5),          std::make_pair(162, &L1_SingleEG26er2p5),          std::make_pair(163, &L1_SingleEG28_FWD2p5),          std::make_pair(166, &L1_SingleEG28er1p5),          std::make_pair(165, &L1_SingleEG28er2p1),          std::make_pair(164, &L1_SingleEG28er2p5),          std::make_pair(167, &L1_SingleEG34er2p5),          std::make_pair(168, &L1_SingleEG36er2p5),          std::make_pair(169, &L1_SingleEG38er2p5),          std::make_pair(170, &L1_SingleEG40er2p5),          std::make_pair(171, &L1_SingleEG42er2p5),          std::make_pair(172, &L1_SingleEG45er2p5),          std::make_pair(173, &L1_SingleEG50),          std::make_pair(174, &L1_SingleEG60),          std::make_pair(159, &L1_SingleEG8er2p5),          std::make_pair(184, &L1_SingleIsoEG24er1p5),          std::make_pair(183, &L1_SingleIsoEG24er2p1),          std::make_pair(187, &L1_SingleIsoEG26er1p5),          std::make_pair(186, &L1_SingleIsoEG26er2p1),          std::make_pair(185, &L1_SingleIsoEG26er2p5),          std::make_pair(188, &L1_SingleIsoEG28_FWD2p5),          std::make_pair(191, &L1_SingleIsoEG28er1p5),          std::make_pair(190, &L1_SingleIsoEG28er2p1),          std::make_pair(189, &L1_SingleIsoEG28er2p5),          std::make_pair(193, &L1_SingleIsoEG30er2p1),          std::make_pair(192, &L1_SingleIsoEG30er2p5),          std::make_pair(195, &L1_SingleIsoEG32er2p1),          std::make_pair(194, &L1_SingleIsoEG32er2p5),          std::make_pair(196, &L1_SingleIsoEG34er2p5),          std::make_pair(330, &L1_SingleJet10erHE),          std::make_pair(312, &L1_SingleJet120),          std::make_pair(327, &L1_SingleJet120_FWD3p0),          std::make_pair(319, &L1_SingleJet120er2p5),          std::make_pair(331, &L1_SingleJet12erHE),          std::make_pair(320, &L1_SingleJet140er2p5),          std::make_pair(332, &L1_SingleJet140er2p5_ETMHF70),          std::make_pair(333, &L1_SingleJet140er2p5_ETMHF80),          std::make_pair(334, &L1_SingleJet140er2p5_ETMHF90),          std::make_pair(321, &L1_SingleJet160er2p5),          std::make_pair(313, &L1_SingleJet180),          std::make_pair(322, &L1_SingleJet180er2p5),          std::make_pair(314, &L1_SingleJet200),          std::make_pair(450, &L1_SingleJet20er2p5_NotBptxOR),          std::make_pair(451, &L1_SingleJet20er2p5_NotBptxOR_3BX),          std::make_pair(309, &L1_SingleJet35),          std::make_pair(324, &L1_SingleJet35_FWD3p0),          std::make_pair(316, &L1_SingleJet35er2p5),          std::make_pair(452, &L1_SingleJet43er2p5_NotBptxOR_3BX),          std::make_pair(453, &L1_SingleJet46er2p5_NotBptxOR_3BX),          std::make_pair(310, &L1_SingleJet60),          std::make_pair(325, &L1_SingleJet60_FWD3p0),          std::make_pair(317, &L1_SingleJet60er2p5),          std::make_pair(329, &L1_SingleJet8erHE),          std::make_pair(311, &L1_SingleJet90),          std::make_pair(326, &L1_SingleJet90_FWD3p0),          std::make_pair(318, &L1_SingleJet90er2p5),          std::make_pair(176, &L1_SingleLooseIsoEG26er1p5),          std::make_pair(175, &L1_SingleLooseIsoEG26er2p5),          std::make_pair(177, &L1_SingleLooseIsoEG28_FWD2p5),          std::make_pair(180, &L1_SingleLooseIsoEG28er1p5),          std::make_pair(179, &L1_SingleLooseIsoEG28er2p1),          std::make_pair(178, &L1_SingleLooseIsoEG28er2p5),          std::make_pair(182, &L1_SingleLooseIsoEG30er1p5),          std::make_pair(181, &L1_SingleLooseIsoEG30er2p5),          std::make_pair(6, &L1_SingleMu0_BMTF),          std::make_pair(5, &L1_SingleMu0_DQ),          std::make_pair(8, &L1_SingleMu0_EMTF),          std::make_pair(7, &L1_SingleMu0_OMTF),          std::make_pair(29, &L1_SingleMu10er1p5),          std::make_pair(13, &L1_SingleMu12_DQ_BMTF),          std::make_pair(15, &L1_SingleMu12_DQ_EMTF),          std::make_pair(14, &L1_SingleMu12_DQ_OMTF),          std::make_pair(30, &L1_SingleMu12er1p5),          std::make_pair(31, &L1_SingleMu14er1p5),          std::make_pair(16, &L1_SingleMu15_DQ),          std::make_pair(32, &L1_SingleMu16er1p5),          std::make_pair(17, &L1_SingleMu18),          std::make_pair(33, &L1_SingleMu18er1p5),          std::make_pair(18, &L1_SingleMu20),          std::make_pair(19, &L1_SingleMu22),          std::make_pair(20, &L1_SingleMu22_BMTF),          std::make_pair(22, &L1_SingleMu22_EMTF),          std::make_pair(21, &L1_SingleMu22_OMTF),          std::make_pair(23, &L1_SingleMu25),          std::make_pair(9, &L1_SingleMu3),          std::make_pair(10, &L1_SingleMu5),          std::make_pair(25, &L1_SingleMu6er1p5),          std::make_pair(12, &L1_SingleMu7),          std::make_pair(11, &L1_SingleMu7_DQ),          std::make_pair(26, &L1_SingleMu7er1p5),          std::make_pair(27, &L1_SingleMu8er1p5),          std::make_pair(28, &L1_SingleMu9er1p5),          std::make_pair(0, &L1_SingleMuCosmics),          std::make_pair(1, &L1_SingleMuCosmics_BMTF),          std::make_pair(3, &L1_SingleMuCosmics_EMTF),          std::make_pair(2, &L1_SingleMuCosmics_OMTF),          std::make_pair(4, &L1_SingleMuOpen),          std::make_pair(446, &L1_SingleMuOpen_NotBptxOR),          std::make_pair(448, &L1_SingleMuOpen_er1p1_NotBptxOR_3BX),          std::make_pair(447, &L1_SingleMuOpen_er1p4_NotBptxOR_3BX),          std::make_pair(264, &L1_SingleTau120er2p1),          std::make_pair(265, &L1_SingleTau130er2p1),          std::make_pair(503, &L1_TOTEM_1),          std::make_pair(504, &L1_TOTEM_2),          std::make_pair(505, &L1_TOTEM_3),          std::make_pair(506, &L1_TOTEM_4),          std::make_pair(228, &L1_TripleEG16er2p5),          std::make_pair(224, &L1_TripleEG_16_12_8_er2p5),          std::make_pair(225, &L1_TripleEG_16_15_8_er2p5),          std::make_pair(226, &L1_TripleEG_18_17_8_er2p5),          std::make_pair(227, &L1_TripleEG_18_18_12_er2p5),          std::make_pair(373, &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5),          std::make_pair(374, &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5),          std::make_pair(372, &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5),          std::make_pair(72, &L1_TripleMu0),          std::make_pair(71, &L1_TripleMu0_OQ),          std::make_pair(73, &L1_TripleMu0_SQ),          std::make_pair(74, &L1_TripleMu3),          std::make_pair(75, &L1_TripleMu3_SQ),          std::make_pair(76, &L1_TripleMu_5SQ_3SQ_0OQ),          std::make_pair(85, &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair(86, &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair(78, &L1_TripleMu_5_3_3),          std::make_pair(79, &L1_TripleMu_5_3_3_SQ),          std::make_pair(77, &L1_TripleMu_5_3p5_2p5),          std::make_pair(83, &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair(82, &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17),          std::make_pair(84, &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair(80, &L1_TripleMu_5_5_3),          std::make_pair(469, &L1_UnpairedBunchBptxMinus),          std::make_pair(468, &L1_UnpairedBunchBptxPlus),          std::make_pair(459, &L1_ZeroBias),          std::make_pair(460, &L1_ZeroBias_copy)      };

  static const std::map<int, AlgorithmFunction> Id2Func(id2func, id2func + sizeof(id2func) / sizeof(id2func[0]));
  const std::map<int, AlgorithmFunction>::const_iterator rc = Id2Func.find(index);
  AlgorithmFunction fp = 0;
  if (rc != Id2Func.end()) fp = rc->second;
  return fp;
}


AlgorithmFunction getFuncFromName(const std::string& name)
{
  static const std::pair<std::string, AlgorithmFunction> name2func[] = {
          std::make_pair("L1_AlwaysTrue", &L1_AlwaysTrue),          std::make_pair("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME),          std::make_pair("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME),          std::make_pair("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME),          std::make_pair("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME),          std::make_pair("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME),          std::make_pair("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME),          std::make_pair("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME),          std::make_pair("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME),          std::make_pair("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME),          std::make_pair("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME),          std::make_pair("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME),          std::make_pair("L1_BptxMinus", &L1_BptxMinus),          std::make_pair("L1_BptxOR", &L1_BptxOR),          std::make_pair("L1_BptxPlus", &L1_BptxPlus),          std::make_pair("L1_BptxXOR", &L1_BptxXOR),          std::make_pair("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142),          std::make_pair("L1_DoubleEG8er2p5_HTT260er", &L1_DoubleEG8er2p5_HTT260er),          std::make_pair("L1_DoubleEG8er2p5_HTT280er", &L1_DoubleEG8er2p5_HTT280er),          std::make_pair("L1_DoubleEG8er2p5_HTT300er", &L1_DoubleEG8er2p5_HTT300er),          std::make_pair("L1_DoubleEG8er2p5_HTT320er", &L1_DoubleEG8er2p5_HTT320er),          std::make_pair("L1_DoubleEG8er2p5_HTT340er", &L1_DoubleEG8er2p5_HTT340er),          std::make_pair("L1_DoubleEG_15_10_er2p5", &L1_DoubleEG_15_10_er2p5),          std::make_pair("L1_DoubleEG_20_10_er2p5", &L1_DoubleEG_20_10_er2p5),          std::make_pair("L1_DoubleEG_22_10_er2p5", &L1_DoubleEG_22_10_er2p5),          std::make_pair("L1_DoubleEG_25_12_er2p5", &L1_DoubleEG_25_12_er2p5),          std::make_pair("L1_DoubleEG_25_14_er2p5", &L1_DoubleEG_25_14_er2p5),          std::make_pair("L1_DoubleEG_27_14_er2p5", &L1_DoubleEG_27_14_er2p5),          std::make_pair("L1_DoubleEG_LooseIso20_10_er2p5", &L1_DoubleEG_LooseIso20_10_er2p5),          std::make_pair("L1_DoubleEG_LooseIso22_10_er2p5", &L1_DoubleEG_LooseIso22_10_er2p5),          std::make_pair("L1_DoubleEG_LooseIso22_12_er2p5", &L1_DoubleEG_LooseIso22_12_er2p5),          std::make_pair("L1_DoubleEG_LooseIso25_12_er2p5", &L1_DoubleEG_LooseIso25_12_er2p5),          std::make_pair("L1_DoubleIsoTau28er2p1", &L1_DoubleIsoTau28er2p1),          std::make_pair("L1_DoubleIsoTau28er2p1_Mass_Max80", &L1_DoubleIsoTau28er2p1_Mass_Max80),          std::make_pair("L1_DoubleIsoTau28er2p1_Mass_Max90", &L1_DoubleIsoTau28er2p1_Mass_Max90),          std::make_pair("L1_DoubleIsoTau30er2p1", &L1_DoubleIsoTau30er2p1),          std::make_pair("L1_DoubleIsoTau30er2p1_Mass_Max80", &L1_DoubleIsoTau30er2p1_Mass_Max80),          std::make_pair("L1_DoubleIsoTau30er2p1_Mass_Max90", &L1_DoubleIsoTau30er2p1_Mass_Max90),          std::make_pair("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1),          std::make_pair("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1),          std::make_pair("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1),          std::make_pair("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6),          std::make_pair("L1_DoubleJet100er2p5", &L1_DoubleJet100er2p5),          std::make_pair("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6),          std::make_pair("L1_DoubleJet120er2p5", &L1_DoubleJet120er2p5),          std::make_pair("L1_DoubleJet150er2p5", &L1_DoubleJet150er2p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5),          std::make_pair("L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp),          std::make_pair("L1_DoubleJet40er2p5", &L1_DoubleJet40er2p5),          std::make_pair("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620),          std::make_pair("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_Mu8", &L1_DoubleJet_80_30_Mass_Min420_Mu8),          std::make_pair("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620),          std::make_pair("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1),          std::make_pair("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1),          std::make_pair("L1_DoubleMu0", &L1_DoubleMu0),          std::make_pair("L1_DoubleMu0_Mass_Min1", &L1_DoubleMu0_Mass_Min1),          std::make_pair("L1_DoubleMu0_OQ", &L1_DoubleMu0_OQ),          std::make_pair("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ),          std::make_pair("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS),          std::make_pair("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er1p5_SQ", &L1_DoubleMu0er1p5_SQ),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4),          std::make_pair("L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er2p0_SQ_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_dR_Max1p4),          std::make_pair("L1_DoubleMu18er2p1", &L1_DoubleMu18er2p1),          std::make_pair("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", &L1_DoubleMu3_OS_DoubleEG7p5Upsilon),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_HTT60er", &L1_DoubleMu3_SQ_ETMHF50_HTT60er),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5),          std::make_pair("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5),          std::make_pair("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er),          std::make_pair("L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er),          std::make_pair("L1_DoubleMu3_SQ_HTT260er", &L1_DoubleMu3_SQ_HTT260er),          std::make_pair("L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair("L1_DoubleMu4_SQ_EG9er2p5", &L1_DoubleMu4_SQ_EG9er2p5),          std::make_pair("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS),          std::make_pair("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2),          std::make_pair("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS),          std::make_pair("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7),          std::make_pair("L1_DoubleMu5Upsilon_OS_DoubleEG3", &L1_DoubleMu5Upsilon_OS_DoubleEG3),          std::make_pair("L1_DoubleMu5_SQ_EG9er2p5", &L1_DoubleMu5_SQ_EG9er2p5),          std::make_pair("L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ),          std::make_pair("L1_DoubleMu9_SQ", &L1_DoubleMu9_SQ),          std::make_pair("L1_DoubleMu_12_5", &L1_DoubleMu_12_5),          std::make_pair("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ),          std::make_pair("L1_DoubleMu_15_7", &L1_DoubleMu_15_7),          std::make_pair("L1_DoubleMu_15_7_Mass_Min1", &L1_DoubleMu_15_7_Mass_Min1),          std::make_pair("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ),          std::make_pair("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1),          std::make_pair("L1_ETM120", &L1_ETM120),          std::make_pair("L1_ETM150", &L1_ETM150),          std::make_pair("L1_ETMHF100", &L1_ETMHF100),          std::make_pair("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er),          std::make_pair("L1_ETMHF110", &L1_ETMHF110),          std::make_pair("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er),          std::make_pair("L1_ETMHF110_HTT60er_NotSecondBunchInTrain", &L1_ETMHF110_HTT60er_NotSecondBunchInTrain),          std::make_pair("L1_ETMHF120", &L1_ETMHF120),          std::make_pair("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er),          std::make_pair("L1_ETMHF120_NotSecondBunchInTrain", &L1_ETMHF120_NotSecondBunchInTrain),          std::make_pair("L1_ETMHF130", &L1_ETMHF130),          std::make_pair("L1_ETMHF130_HTT60er", &L1_ETMHF130_HTT60er),          std::make_pair("L1_ETMHF140", &L1_ETMHF140),          std::make_pair("L1_ETMHF150", &L1_ETMHF150),          std::make_pair("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er),          std::make_pair("L1_ETT1200", &L1_ETT1200),          std::make_pair("L1_ETT1600", &L1_ETT1600),          std::make_pair("L1_ETT2000", &L1_ETT2000),          std::make_pair("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain),          std::make_pair("L1_FirstBunchBeforeTrain", &L1_FirstBunchBeforeTrain),          std::make_pair("L1_FirstBunchInTrain", &L1_FirstBunchInTrain),          std::make_pair("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit),          std::make_pair("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain),          std::make_pair("L1_HCAL_LaserMon_Trig", &L1_HCAL_LaserMon_Trig),          std::make_pair("L1_HCAL_LaserMon_Veto", &L1_HCAL_LaserMon_Veto),          std::make_pair("L1_HTT120er", &L1_HTT120er),          std::make_pair("L1_HTT160er", &L1_HTT160er),          std::make_pair("L1_HTT200er", &L1_HTT200er),          std::make_pair("L1_HTT255er", &L1_HTT255er),          std::make_pair("L1_HTT280er", &L1_HTT280er),          std::make_pair("L1_HTT280er_QuadJet_70_55_40_35_er2p4", &L1_HTT280er_QuadJet_70_55_40_35_er2p4),          std::make_pair("L1_HTT320er", &L1_HTT320er),          std::make_pair("L1_HTT320er_QuadJet_70_55_40_40_er2p4", &L1_HTT320er_QuadJet_70_55_40_40_er2p4),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3),          std::make_pair("L1_HTT360er", &L1_HTT360er),          std::make_pair("L1_HTT400er", &L1_HTT400er),          std::make_pair("L1_HTT450er", &L1_HTT450er),          std::make_pair("L1_IsoEG32er2p5_Mt40", &L1_IsoEG32er2p5_Mt40),          std::make_pair("L1_IsoEG32er2p5_Mt44", &L1_IsoEG32er2p5_Mt44),          std::make_pair("L1_IsoEG32er2p5_Mt48", &L1_IsoEG32er2p5_Mt48),          std::make_pair("L1_IsoTau40er2p1_ETMHF100", &L1_IsoTau40er2p1_ETMHF100),          std::make_pair("L1_IsoTau40er2p1_ETMHF110", &L1_IsoTau40er2p1_ETMHF110),          std::make_pair("L1_IsoTau40er2p1_ETMHF80", &L1_IsoTau40er2p1_ETMHF80),          std::make_pair("L1_IsoTau40er2p1_ETMHF90", &L1_IsoTau40er2p1_ETMHF90),          std::make_pair("L1_IsolatedBunch", &L1_IsolatedBunch),          std::make_pair("L1_LastBunchInTrain", &L1_LastBunchInTrain),          std::make_pair("L1_LastCollisionInTrain", &L1_LastCollisionInTrain),          std::make_pair("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_LooseIsoEG30er2p1_HTT100er", &L1_LooseIsoEG30er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND),          std::make_pair("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6),          std::make_pair("L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6),          std::make_pair("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6),          std::make_pair("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1),          std::make_pair("L1_Mu18er2p1_Tau26er2p1", &L1_Mu18er2p1_Tau26er2p1),          std::make_pair("L1_Mu20_EG10er2p5", &L1_Mu20_EG10er2p5),          std::make_pair("L1_Mu22er2p1_IsoTau28er2p1", &L1_Mu22er2p1_IsoTau28er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau30er2p1", &L1_Mu22er2p1_IsoTau30er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1),          std::make_pair("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p4", &L1_Mu3_Jet120er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p8", &L1_Mu3_Jet120er2p5_dR_Max0p8),          std::make_pair("L1_Mu3_Jet16er2p5_dR_Max0p4", &L1_Mu3_Jet16er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5),          std::make_pair("L1_Mu3_Jet35er2p5_dR_Max0p4", &L1_Mu3_Jet35er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet60er2p5_dR_Max0p4", &L1_Mu3_Jet60er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet80er2p5_dR_Max0p4", &L1_Mu3_Jet80er2p5_dR_Max0p4),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF40", &L1_Mu3er1p5_Jet100er2p5_ETMHF40),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF50", &L1_Mu3er1p5_Jet100er2p5_ETMHF50),          std::make_pair("L1_Mu5_EG23er2p5", &L1_Mu5_EG23er2p5),          std::make_pair("L1_Mu5_LooseIsoEG20er2p5", &L1_Mu5_LooseIsoEG20er2p5),          std::make_pair("L1_Mu6_DoubleEG10er2p5", &L1_Mu6_DoubleEG10er2p5),          std::make_pair("L1_Mu6_DoubleEG12er2p5", &L1_Mu6_DoubleEG12er2p5),          std::make_pair("L1_Mu6_DoubleEG15er2p5", &L1_Mu6_DoubleEG15er2p5),          std::make_pair("L1_Mu6_DoubleEG17er2p5", &L1_Mu6_DoubleEG17er2p5),          std::make_pair("L1_Mu6_HTT240er", &L1_Mu6_HTT240er),          std::make_pair("L1_Mu6_HTT250er", &L1_Mu6_HTT250er),          std::make_pair("L1_Mu7_EG20er2p5", &L1_Mu7_EG20er2p5),          std::make_pair("L1_Mu7_EG23er2p5", &L1_Mu7_EG23er2p5),          std::make_pair("L1_Mu7_LooseIsoEG20er2p5", &L1_Mu7_LooseIsoEG20er2p5),          std::make_pair("L1_Mu7_LooseIsoEG23er2p5", &L1_Mu7_LooseIsoEG23er2p5),          std::make_pair("L1_NotBptxOR", &L1_NotBptxOR),          std::make_pair("L1_QuadJet36er2p5_IsoTau52er2p1", &L1_QuadJet36er2p5_IsoTau52er2p1),          std::make_pair("L1_QuadJet60er2p5", &L1_QuadJet60er2p5),          std::make_pair("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0),          std::make_pair("L1_QuadMu0", &L1_QuadMu0),          std::make_pair("L1_QuadMu0_OQ", &L1_QuadMu0_OQ),          std::make_pair("L1_QuadMu0_SQ", &L1_QuadMu0_SQ),          std::make_pair("L1_SecondBunchInTrain", &L1_SecondBunchInTrain),          std::make_pair("L1_SecondLastBunchInTrain", &L1_SecondLastBunchInTrain),          std::make_pair("L1_SingleEG10er2p5", &L1_SingleEG10er2p5),          std::make_pair("L1_SingleEG15er2p5", &L1_SingleEG15er2p5),          std::make_pair("L1_SingleEG26er2p5", &L1_SingleEG26er2p5),          std::make_pair("L1_SingleEG28_FWD2p5", &L1_SingleEG28_FWD2p5),          std::make_pair("L1_SingleEG28er1p5", &L1_SingleEG28er1p5),          std::make_pair("L1_SingleEG28er2p1", &L1_SingleEG28er2p1),          std::make_pair("L1_SingleEG28er2p5", &L1_SingleEG28er2p5),          std::make_pair("L1_SingleEG34er2p5", &L1_SingleEG34er2p5),          std::make_pair("L1_SingleEG36er2p5", &L1_SingleEG36er2p5),          std::make_pair("L1_SingleEG38er2p5", &L1_SingleEG38er2p5),          std::make_pair("L1_SingleEG40er2p5", &L1_SingleEG40er2p5),          std::make_pair("L1_SingleEG42er2p5", &L1_SingleEG42er2p5),          std::make_pair("L1_SingleEG45er2p5", &L1_SingleEG45er2p5),          std::make_pair("L1_SingleEG50", &L1_SingleEG50),          std::make_pair("L1_SingleEG60", &L1_SingleEG60),          std::make_pair("L1_SingleEG8er2p5", &L1_SingleEG8er2p5),          std::make_pair("L1_SingleIsoEG24er1p5", &L1_SingleIsoEG24er1p5),          std::make_pair("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1),          std::make_pair("L1_SingleIsoEG26er1p5", &L1_SingleIsoEG26er1p5),          std::make_pair("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1),          std::make_pair("L1_SingleIsoEG26er2p5", &L1_SingleIsoEG26er2p5),          std::make_pair("L1_SingleIsoEG28_FWD2p5", &L1_SingleIsoEG28_FWD2p5),          std::make_pair("L1_SingleIsoEG28er1p5", &L1_SingleIsoEG28er1p5),          std::make_pair("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1),          std::make_pair("L1_SingleIsoEG28er2p5", &L1_SingleIsoEG28er2p5),          std::make_pair("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1),          std::make_pair("L1_SingleIsoEG30er2p5", &L1_SingleIsoEG30er2p5),          std::make_pair("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1),          std::make_pair("L1_SingleIsoEG32er2p5", &L1_SingleIsoEG32er2p5),          std::make_pair("L1_SingleIsoEG34er2p5", &L1_SingleIsoEG34er2p5),          std::make_pair("L1_SingleJet10erHE", &L1_SingleJet10erHE),          std::make_pair("L1_SingleJet120", &L1_SingleJet120),          std::make_pair("L1_SingleJet120_FWD3p0", &L1_SingleJet120_FWD3p0),          std::make_pair("L1_SingleJet120er2p5", &L1_SingleJet120er2p5),          std::make_pair("L1_SingleJet12erHE", &L1_SingleJet12erHE),          std::make_pair("L1_SingleJet140er2p5", &L1_SingleJet140er2p5),          std::make_pair("L1_SingleJet140er2p5_ETMHF70", &L1_SingleJet140er2p5_ETMHF70),          std::make_pair("L1_SingleJet140er2p5_ETMHF80", &L1_SingleJet140er2p5_ETMHF80),          std::make_pair("L1_SingleJet140er2p5_ETMHF90", &L1_SingleJet140er2p5_ETMHF90),          std::make_pair("L1_SingleJet160er2p5", &L1_SingleJet160er2p5),          std::make_pair("L1_SingleJet180", &L1_SingleJet180),          std::make_pair("L1_SingleJet180er2p5", &L1_SingleJet180er2p5),          std::make_pair("L1_SingleJet200", &L1_SingleJet200),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR", &L1_SingleJet20er2p5_NotBptxOR),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR_3BX", &L1_SingleJet20er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet35", &L1_SingleJet35),          std::make_pair("L1_SingleJet35_FWD3p0", &L1_SingleJet35_FWD3p0),          std::make_pair("L1_SingleJet35er2p5", &L1_SingleJet35er2p5),          std::make_pair("L1_SingleJet43er2p5_NotBptxOR_3BX", &L1_SingleJet43er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet46er2p5_NotBptxOR_3BX", &L1_SingleJet46er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet60", &L1_SingleJet60),          std::make_pair("L1_SingleJet60_FWD3p0", &L1_SingleJet60_FWD3p0),          std::make_pair("L1_SingleJet60er2p5", &L1_SingleJet60er2p5),          std::make_pair("L1_SingleJet8erHE", &L1_SingleJet8erHE),          std::make_pair("L1_SingleJet90", &L1_SingleJet90),          std::make_pair("L1_SingleJet90_FWD3p0", &L1_SingleJet90_FWD3p0),          std::make_pair("L1_SingleJet90er2p5", &L1_SingleJet90er2p5),          std::make_pair("L1_SingleLooseIsoEG26er1p5", &L1_SingleLooseIsoEG26er1p5),          std::make_pair("L1_SingleLooseIsoEG26er2p5", &L1_SingleLooseIsoEG26er2p5),          std::make_pair("L1_SingleLooseIsoEG28_FWD2p5", &L1_SingleLooseIsoEG28_FWD2p5),          std::make_pair("L1_SingleLooseIsoEG28er1p5", &L1_SingleLooseIsoEG28er1p5),          std::make_pair("L1_SingleLooseIsoEG28er2p1", &L1_SingleLooseIsoEG28er2p1),          std::make_pair("L1_SingleLooseIsoEG28er2p5", &L1_SingleLooseIsoEG28er2p5),          std::make_pair("L1_SingleLooseIsoEG30er1p5", &L1_SingleLooseIsoEG30er1p5),          std::make_pair("L1_SingleLooseIsoEG30er2p5", &L1_SingleLooseIsoEG30er2p5),          std::make_pair("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF),          std::make_pair("L1_SingleMu0_DQ", &L1_SingleMu0_DQ),          std::make_pair("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF),          std::make_pair("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF),          std::make_pair("L1_SingleMu10er1p5", &L1_SingleMu10er1p5),          std::make_pair("L1_SingleMu12_DQ_BMTF", &L1_SingleMu12_DQ_BMTF),          std::make_pair("L1_SingleMu12_DQ_EMTF", &L1_SingleMu12_DQ_EMTF),          std::make_pair("L1_SingleMu12_DQ_OMTF", &L1_SingleMu12_DQ_OMTF),          std::make_pair("L1_SingleMu12er1p5", &L1_SingleMu12er1p5),          std::make_pair("L1_SingleMu14er1p5", &L1_SingleMu14er1p5),          std::make_pair("L1_SingleMu15_DQ", &L1_SingleMu15_DQ),          std::make_pair("L1_SingleMu16er1p5", &L1_SingleMu16er1p5),          std::make_pair("L1_SingleMu18", &L1_SingleMu18),          std::make_pair("L1_SingleMu18er1p5", &L1_SingleMu18er1p5),          std::make_pair("L1_SingleMu20", &L1_SingleMu20),          std::make_pair("L1_SingleMu22", &L1_SingleMu22),          std::make_pair("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF),          std::make_pair("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF),          std::make_pair("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF),          std::make_pair("L1_SingleMu25", &L1_SingleMu25),          std::make_pair("L1_SingleMu3", &L1_SingleMu3),          std::make_pair("L1_SingleMu5", &L1_SingleMu5),          std::make_pair("L1_SingleMu6er1p5", &L1_SingleMu6er1p5),          std::make_pair("L1_SingleMu7", &L1_SingleMu7),          std::make_pair("L1_SingleMu7_DQ", &L1_SingleMu7_DQ),          std::make_pair("L1_SingleMu7er1p5", &L1_SingleMu7er1p5),          std::make_pair("L1_SingleMu8er1p5", &L1_SingleMu8er1p5),          std::make_pair("L1_SingleMu9er1p5", &L1_SingleMu9er1p5),          std::make_pair("L1_SingleMuCosmics", &L1_SingleMuCosmics),          std::make_pair("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF),          std::make_pair("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF),          std::make_pair("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF),          std::make_pair("L1_SingleMuOpen", &L1_SingleMuOpen),          std::make_pair("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR),          std::make_pair("L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX),          std::make_pair("L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX),          std::make_pair("L1_SingleTau120er2p1", &L1_SingleTau120er2p1),          std::make_pair("L1_SingleTau130er2p1", &L1_SingleTau130er2p1),          std::make_pair("L1_TOTEM_1", &L1_TOTEM_1),          std::make_pair("L1_TOTEM_2", &L1_TOTEM_2),          std::make_pair("L1_TOTEM_3", &L1_TOTEM_3),          std::make_pair("L1_TOTEM_4", &L1_TOTEM_4),          std::make_pair("L1_TripleEG16er2p5", &L1_TripleEG16er2p5),          std::make_pair("L1_TripleEG_16_12_8_er2p5", &L1_TripleEG_16_12_8_er2p5),          std::make_pair("L1_TripleEG_16_15_8_er2p5", &L1_TripleEG_16_15_8_er2p5),          std::make_pair("L1_TripleEG_18_17_8_er2p5", &L1_TripleEG_18_17_8_er2p5),          std::make_pair("L1_TripleEG_18_18_12_er2p5", &L1_TripleEG_18_18_12_er2p5),          std::make_pair("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5),          std::make_pair("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5),          std::make_pair("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5),          std::make_pair("L1_TripleMu0", &L1_TripleMu0),          std::make_pair("L1_TripleMu0_OQ", &L1_TripleMu0_OQ),          std::make_pair("L1_TripleMu0_SQ", &L1_TripleMu0_SQ),          std::make_pair("L1_TripleMu3", &L1_TripleMu3),          std::make_pair("L1_TripleMu3_SQ", &L1_TripleMu3_SQ),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3),          std::make_pair("L1_TripleMu_5_3_3_SQ", &L1_TripleMu_5_3_3_SQ),          std::make_pair("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5),          std::make_pair("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3),          std::make_pair("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus),          std::make_pair("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus),          std::make_pair("L1_ZeroBias", &L1_ZeroBias),          std::make_pair("L1_ZeroBias_copy", &L1_ZeroBias_copy)      };

  static const std::map<std::string, AlgorithmFunction> Name2Func(name2func, name2func + sizeof(name2func) / sizeof(name2func[0]));
  const std::map<std::string, AlgorithmFunction>::const_iterator rc = Name2Func.find(name);
  AlgorithmFunction fp = 0;
  if (rc != Name2Func.end()) fp = rc->second;
  if (fp == 0)
  {
    std::stringstream ss;
    ss << "fat> algorithm '" << name << "' is not defined in L1Menu_Collisions2018_v2_1_0\n";
    throw std::runtime_error(ss.str());
  }
  return fp;
}


bool addFuncFromName(std::map<std::string, std::function<bool()>> &L1SeedFun,
                     L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade,
                     L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  static const std::pair<std::string, AlgorithmFunction> name2func[] = {
          std::make_pair("L1_AlwaysTrue", &L1_AlwaysTrue),          std::make_pair("L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME),          std::make_pair("L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME),          std::make_pair("L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME),          std::make_pair("L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME),          std::make_pair("L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME),          std::make_pair("L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME),          std::make_pair("L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME),          std::make_pair("L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME),          std::make_pair("L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME),          std::make_pair("L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME),          std::make_pair("L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME),          std::make_pair("L1_BptxMinus", &L1_BptxMinus),          std::make_pair("L1_BptxOR", &L1_BptxOR),          std::make_pair("L1_BptxPlus", &L1_BptxPlus),          std::make_pair("L1_BptxXOR", &L1_BptxXOR),          std::make_pair("L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142),          std::make_pair("L1_DoubleEG8er2p5_HTT260er", &L1_DoubleEG8er2p5_HTT260er),          std::make_pair("L1_DoubleEG8er2p5_HTT280er", &L1_DoubleEG8er2p5_HTT280er),          std::make_pair("L1_DoubleEG8er2p5_HTT300er", &L1_DoubleEG8er2p5_HTT300er),          std::make_pair("L1_DoubleEG8er2p5_HTT320er", &L1_DoubleEG8er2p5_HTT320er),          std::make_pair("L1_DoubleEG8er2p5_HTT340er", &L1_DoubleEG8er2p5_HTT340er),          std::make_pair("L1_DoubleEG_15_10_er2p5", &L1_DoubleEG_15_10_er2p5),          std::make_pair("L1_DoubleEG_20_10_er2p5", &L1_DoubleEG_20_10_er2p5),          std::make_pair("L1_DoubleEG_22_10_er2p5", &L1_DoubleEG_22_10_er2p5),          std::make_pair("L1_DoubleEG_25_12_er2p5", &L1_DoubleEG_25_12_er2p5),          std::make_pair("L1_DoubleEG_25_14_er2p5", &L1_DoubleEG_25_14_er2p5),          std::make_pair("L1_DoubleEG_27_14_er2p5", &L1_DoubleEG_27_14_er2p5),          std::make_pair("L1_DoubleEG_LooseIso20_10_er2p5", &L1_DoubleEG_LooseIso20_10_er2p5),          std::make_pair("L1_DoubleEG_LooseIso22_10_er2p5", &L1_DoubleEG_LooseIso22_10_er2p5),          std::make_pair("L1_DoubleEG_LooseIso22_12_er2p5", &L1_DoubleEG_LooseIso22_12_er2p5),          std::make_pair("L1_DoubleEG_LooseIso25_12_er2p5", &L1_DoubleEG_LooseIso25_12_er2p5),          std::make_pair("L1_DoubleIsoTau28er2p1", &L1_DoubleIsoTau28er2p1),          std::make_pair("L1_DoubleIsoTau28er2p1_Mass_Max80", &L1_DoubleIsoTau28er2p1_Mass_Max80),          std::make_pair("L1_DoubleIsoTau28er2p1_Mass_Max90", &L1_DoubleIsoTau28er2p1_Mass_Max90),          std::make_pair("L1_DoubleIsoTau30er2p1", &L1_DoubleIsoTau30er2p1),          std::make_pair("L1_DoubleIsoTau30er2p1_Mass_Max80", &L1_DoubleIsoTau30er2p1_Mass_Max80),          std::make_pair("L1_DoubleIsoTau30er2p1_Mass_Max90", &L1_DoubleIsoTau30er2p1_Mass_Max90),          std::make_pair("L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1),          std::make_pair("L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1),          std::make_pair("L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1),          std::make_pair("L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6),          std::make_pair("L1_DoubleJet100er2p5", &L1_DoubleJet100er2p5),          std::make_pair("L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6),          std::make_pair("L1_DoubleJet120er2p5", &L1_DoubleJet120er2p5),          std::make_pair("L1_DoubleJet150er2p5", &L1_DoubleJet150er2p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5),          std::make_pair("L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5),          std::make_pair("L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp),          std::make_pair("L1_DoubleJet40er2p5", &L1_DoubleJet40er2p5),          std::make_pair("L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620),          std::make_pair("L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620),          std::make_pair("L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620),          std::make_pair("L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp),          std::make_pair("L1_DoubleJet_80_30_Mass_Min420_Mu8", &L1_DoubleJet_80_30_Mass_Min420_Mu8),          std::make_pair("L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620),          std::make_pair("L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1),          std::make_pair("L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1),          std::make_pair("L1_DoubleMu0", &L1_DoubleMu0),          std::make_pair("L1_DoubleMu0_Mass_Min1", &L1_DoubleMu0_Mass_Min1),          std::make_pair("L1_DoubleMu0_OQ", &L1_DoubleMu0_OQ),          std::make_pair("L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ),          std::make_pair("L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS),          std::make_pair("L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair("L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er1p5_SQ", &L1_DoubleMu0er1p5_SQ),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS),          std::make_pair("L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4),          std::make_pair("L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4),          std::make_pair("L1_DoubleMu0er2p0_SQ_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_dR_Max1p4),          std::make_pair("L1_DoubleMu18er2p1", &L1_DoubleMu18er2p1),          std::make_pair("L1_DoubleMu3_OS_DoubleEG7p5Upsilon", &L1_DoubleMu3_OS_DoubleEG7p5Upsilon),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_HTT60er", &L1_DoubleMu3_SQ_ETMHF50_HTT60er),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5),          std::make_pair("L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5),          std::make_pair("L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5),          std::make_pair("L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er),          std::make_pair("L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er),          std::make_pair("L1_DoubleMu3_SQ_HTT260er", &L1_DoubleMu3_SQ_HTT260er),          std::make_pair("L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8),          std::make_pair("L1_DoubleMu4_SQ_EG9er2p5", &L1_DoubleMu4_SQ_EG9er2p5),          std::make_pair("L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS),          std::make_pair("L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2),          std::make_pair("L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS),          std::make_pair("L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass7to18),          std::make_pair("L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7),          std::make_pair("L1_DoubleMu5Upsilon_OS_DoubleEG3", &L1_DoubleMu5Upsilon_OS_DoubleEG3),          std::make_pair("L1_DoubleMu5_SQ_EG9er2p5", &L1_DoubleMu5_SQ_EG9er2p5),          std::make_pair("L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ),          std::make_pair("L1_DoubleMu9_SQ", &L1_DoubleMu9_SQ),          std::make_pair("L1_DoubleMu_12_5", &L1_DoubleMu_12_5),          std::make_pair("L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ),          std::make_pair("L1_DoubleMu_15_7", &L1_DoubleMu_15_7),          std::make_pair("L1_DoubleMu_15_7_Mass_Min1", &L1_DoubleMu_15_7_Mass_Min1),          std::make_pair("L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ),          std::make_pair("L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1),          std::make_pair("L1_ETM120", &L1_ETM120),          std::make_pair("L1_ETM150", &L1_ETM150),          std::make_pair("L1_ETMHF100", &L1_ETMHF100),          std::make_pair("L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er),          std::make_pair("L1_ETMHF110", &L1_ETMHF110),          std::make_pair("L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er),          std::make_pair("L1_ETMHF110_HTT60er_NotSecondBunchInTrain", &L1_ETMHF110_HTT60er_NotSecondBunchInTrain),          std::make_pair("L1_ETMHF120", &L1_ETMHF120),          std::make_pair("L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er),          std::make_pair("L1_ETMHF120_NotSecondBunchInTrain", &L1_ETMHF120_NotSecondBunchInTrain),          std::make_pair("L1_ETMHF130", &L1_ETMHF130),          std::make_pair("L1_ETMHF130_HTT60er", &L1_ETMHF130_HTT60er),          std::make_pair("L1_ETMHF140", &L1_ETMHF140),          std::make_pair("L1_ETMHF150", &L1_ETMHF150),          std::make_pair("L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er),          std::make_pair("L1_ETT1200", &L1_ETT1200),          std::make_pair("L1_ETT1600", &L1_ETT1600),          std::make_pair("L1_ETT2000", &L1_ETT2000),          std::make_pair("L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain),          std::make_pair("L1_FirstBunchBeforeTrain", &L1_FirstBunchBeforeTrain),          std::make_pair("L1_FirstBunchInTrain", &L1_FirstBunchInTrain),          std::make_pair("L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit),          std::make_pair("L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain),          std::make_pair("L1_HCAL_LaserMon_Trig", &L1_HCAL_LaserMon_Trig),          std::make_pair("L1_HCAL_LaserMon_Veto", &L1_HCAL_LaserMon_Veto),          std::make_pair("L1_HTT120er", &L1_HTT120er),          std::make_pair("L1_HTT160er", &L1_HTT160er),          std::make_pair("L1_HTT200er", &L1_HTT200er),          std::make_pair("L1_HTT255er", &L1_HTT255er),          std::make_pair("L1_HTT280er", &L1_HTT280er),          std::make_pair("L1_HTT280er_QuadJet_70_55_40_35_er2p4", &L1_HTT280er_QuadJet_70_55_40_35_er2p4),          std::make_pair("L1_HTT320er", &L1_HTT320er),          std::make_pair("L1_HTT320er_QuadJet_70_55_40_40_er2p4", &L1_HTT320er_QuadJet_70_55_40_40_er2p4),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3),          std::make_pair("L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3),          std::make_pair("L1_HTT360er", &L1_HTT360er),          std::make_pair("L1_HTT400er", &L1_HTT400er),          std::make_pair("L1_HTT450er", &L1_HTT450er),          std::make_pair("L1_IsoEG32er2p5_Mt40", &L1_IsoEG32er2p5_Mt40),          std::make_pair("L1_IsoEG32er2p5_Mt44", &L1_IsoEG32er2p5_Mt44),          std::make_pair("L1_IsoEG32er2p5_Mt48", &L1_IsoEG32er2p5_Mt48),          std::make_pair("L1_IsoTau40er2p1_ETMHF100", &L1_IsoTau40er2p1_ETMHF100),          std::make_pair("L1_IsoTau40er2p1_ETMHF110", &L1_IsoTau40er2p1_ETMHF110),          std::make_pair("L1_IsoTau40er2p1_ETMHF80", &L1_IsoTau40er2p1_ETMHF80),          std::make_pair("L1_IsoTau40er2p1_ETMHF90", &L1_IsoTau40er2p1_ETMHF90),          std::make_pair("L1_IsolatedBunch", &L1_IsolatedBunch),          std::make_pair("L1_LastBunchInTrain", &L1_LastBunchInTrain),          std::make_pair("L1_LastCollisionInTrain", &L1_LastCollisionInTrain),          std::make_pair("L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3),          std::make_pair("L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_LooseIsoEG30er2p1_HTT100er", &L1_LooseIsoEG30er2p1_HTT100er),          std::make_pair("L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3),          std::make_pair("L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND),          std::make_pair("L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6),          std::make_pair("L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6),          std::make_pair("L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6),          std::make_pair("L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1),          std::make_pair("L1_Mu18er2p1_Tau26er2p1", &L1_Mu18er2p1_Tau26er2p1),          std::make_pair("L1_Mu20_EG10er2p5", &L1_Mu20_EG10er2p5),          std::make_pair("L1_Mu22er2p1_IsoTau28er2p1", &L1_Mu22er2p1_IsoTau28er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau30er2p1", &L1_Mu22er2p1_IsoTau30er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1),          std::make_pair("L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1),          std::make_pair("L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p4", &L1_Mu3_Jet120er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet120er2p5_dR_Max0p8", &L1_Mu3_Jet120er2p5_dR_Max0p8),          std::make_pair("L1_Mu3_Jet16er2p5_dR_Max0p4", &L1_Mu3_Jet16er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5),          std::make_pair("L1_Mu3_Jet35er2p5_dR_Max0p4", &L1_Mu3_Jet35er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet60er2p5_dR_Max0p4", &L1_Mu3_Jet60er2p5_dR_Max0p4),          std::make_pair("L1_Mu3_Jet80er2p5_dR_Max0p4", &L1_Mu3_Jet80er2p5_dR_Max0p4),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF40", &L1_Mu3er1p5_Jet100er2p5_ETMHF40),          std::make_pair("L1_Mu3er1p5_Jet100er2p5_ETMHF50", &L1_Mu3er1p5_Jet100er2p5_ETMHF50),          std::make_pair("L1_Mu5_EG23er2p5", &L1_Mu5_EG23er2p5),          std::make_pair("L1_Mu5_LooseIsoEG20er2p5", &L1_Mu5_LooseIsoEG20er2p5),          std::make_pair("L1_Mu6_DoubleEG10er2p5", &L1_Mu6_DoubleEG10er2p5),          std::make_pair("L1_Mu6_DoubleEG12er2p5", &L1_Mu6_DoubleEG12er2p5),          std::make_pair("L1_Mu6_DoubleEG15er2p5", &L1_Mu6_DoubleEG15er2p5),          std::make_pair("L1_Mu6_DoubleEG17er2p5", &L1_Mu6_DoubleEG17er2p5),          std::make_pair("L1_Mu6_HTT240er", &L1_Mu6_HTT240er),          std::make_pair("L1_Mu6_HTT250er", &L1_Mu6_HTT250er),          std::make_pair("L1_Mu7_EG20er2p5", &L1_Mu7_EG20er2p5),          std::make_pair("L1_Mu7_EG23er2p5", &L1_Mu7_EG23er2p5),          std::make_pair("L1_Mu7_LooseIsoEG20er2p5", &L1_Mu7_LooseIsoEG20er2p5),          std::make_pair("L1_Mu7_LooseIsoEG23er2p5", &L1_Mu7_LooseIsoEG23er2p5),          std::make_pair("L1_NotBptxOR", &L1_NotBptxOR),          std::make_pair("L1_QuadJet36er2p5_IsoTau52er2p1", &L1_QuadJet36er2p5_IsoTau52er2p1),          std::make_pair("L1_QuadJet60er2p5", &L1_QuadJet60er2p5),          std::make_pair("L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0),          std::make_pair("L1_QuadMu0", &L1_QuadMu0),          std::make_pair("L1_QuadMu0_OQ", &L1_QuadMu0_OQ),          std::make_pair("L1_QuadMu0_SQ", &L1_QuadMu0_SQ),          std::make_pair("L1_SecondBunchInTrain", &L1_SecondBunchInTrain),          std::make_pair("L1_SecondLastBunchInTrain", &L1_SecondLastBunchInTrain),          std::make_pair("L1_SingleEG10er2p5", &L1_SingleEG10er2p5),          std::make_pair("L1_SingleEG15er2p5", &L1_SingleEG15er2p5),          std::make_pair("L1_SingleEG26er2p5", &L1_SingleEG26er2p5),          std::make_pair("L1_SingleEG28_FWD2p5", &L1_SingleEG28_FWD2p5),          std::make_pair("L1_SingleEG28er1p5", &L1_SingleEG28er1p5),          std::make_pair("L1_SingleEG28er2p1", &L1_SingleEG28er2p1),          std::make_pair("L1_SingleEG28er2p5", &L1_SingleEG28er2p5),          std::make_pair("L1_SingleEG34er2p5", &L1_SingleEG34er2p5),          std::make_pair("L1_SingleEG36er2p5", &L1_SingleEG36er2p5),          std::make_pair("L1_SingleEG38er2p5", &L1_SingleEG38er2p5),          std::make_pair("L1_SingleEG40er2p5", &L1_SingleEG40er2p5),          std::make_pair("L1_SingleEG42er2p5", &L1_SingleEG42er2p5),          std::make_pair("L1_SingleEG45er2p5", &L1_SingleEG45er2p5),          std::make_pair("L1_SingleEG50", &L1_SingleEG50),          std::make_pair("L1_SingleEG60", &L1_SingleEG60),          std::make_pair("L1_SingleEG8er2p5", &L1_SingleEG8er2p5),          std::make_pair("L1_SingleIsoEG24er1p5", &L1_SingleIsoEG24er1p5),          std::make_pair("L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1),          std::make_pair("L1_SingleIsoEG26er1p5", &L1_SingleIsoEG26er1p5),          std::make_pair("L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1),          std::make_pair("L1_SingleIsoEG26er2p5", &L1_SingleIsoEG26er2p5),          std::make_pair("L1_SingleIsoEG28_FWD2p5", &L1_SingleIsoEG28_FWD2p5),          std::make_pair("L1_SingleIsoEG28er1p5", &L1_SingleIsoEG28er1p5),          std::make_pair("L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1),          std::make_pair("L1_SingleIsoEG28er2p5", &L1_SingleIsoEG28er2p5),          std::make_pair("L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1),          std::make_pair("L1_SingleIsoEG30er2p5", &L1_SingleIsoEG30er2p5),          std::make_pair("L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1),          std::make_pair("L1_SingleIsoEG32er2p5", &L1_SingleIsoEG32er2p5),          std::make_pair("L1_SingleIsoEG34er2p5", &L1_SingleIsoEG34er2p5),          std::make_pair("L1_SingleJet10erHE", &L1_SingleJet10erHE),          std::make_pair("L1_SingleJet120", &L1_SingleJet120),          std::make_pair("L1_SingleJet120_FWD3p0", &L1_SingleJet120_FWD3p0),          std::make_pair("L1_SingleJet120er2p5", &L1_SingleJet120er2p5),          std::make_pair("L1_SingleJet12erHE", &L1_SingleJet12erHE),          std::make_pair("L1_SingleJet140er2p5", &L1_SingleJet140er2p5),          std::make_pair("L1_SingleJet140er2p5_ETMHF70", &L1_SingleJet140er2p5_ETMHF70),          std::make_pair("L1_SingleJet140er2p5_ETMHF80", &L1_SingleJet140er2p5_ETMHF80),          std::make_pair("L1_SingleJet140er2p5_ETMHF90", &L1_SingleJet140er2p5_ETMHF90),          std::make_pair("L1_SingleJet160er2p5", &L1_SingleJet160er2p5),          std::make_pair("L1_SingleJet180", &L1_SingleJet180),          std::make_pair("L1_SingleJet180er2p5", &L1_SingleJet180er2p5),          std::make_pair("L1_SingleJet200", &L1_SingleJet200),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR", &L1_SingleJet20er2p5_NotBptxOR),          std::make_pair("L1_SingleJet20er2p5_NotBptxOR_3BX", &L1_SingleJet20er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet35", &L1_SingleJet35),          std::make_pair("L1_SingleJet35_FWD3p0", &L1_SingleJet35_FWD3p0),          std::make_pair("L1_SingleJet35er2p5", &L1_SingleJet35er2p5),          std::make_pair("L1_SingleJet43er2p5_NotBptxOR_3BX", &L1_SingleJet43er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet46er2p5_NotBptxOR_3BX", &L1_SingleJet46er2p5_NotBptxOR_3BX),          std::make_pair("L1_SingleJet60", &L1_SingleJet60),          std::make_pair("L1_SingleJet60_FWD3p0", &L1_SingleJet60_FWD3p0),          std::make_pair("L1_SingleJet60er2p5", &L1_SingleJet60er2p5),          std::make_pair("L1_SingleJet8erHE", &L1_SingleJet8erHE),          std::make_pair("L1_SingleJet90", &L1_SingleJet90),          std::make_pair("L1_SingleJet90_FWD3p0", &L1_SingleJet90_FWD3p0),          std::make_pair("L1_SingleJet90er2p5", &L1_SingleJet90er2p5),          std::make_pair("L1_SingleLooseIsoEG26er1p5", &L1_SingleLooseIsoEG26er1p5),          std::make_pair("L1_SingleLooseIsoEG26er2p5", &L1_SingleLooseIsoEG26er2p5),          std::make_pair("L1_SingleLooseIsoEG28_FWD2p5", &L1_SingleLooseIsoEG28_FWD2p5),          std::make_pair("L1_SingleLooseIsoEG28er1p5", &L1_SingleLooseIsoEG28er1p5),          std::make_pair("L1_SingleLooseIsoEG28er2p1", &L1_SingleLooseIsoEG28er2p1),          std::make_pair("L1_SingleLooseIsoEG28er2p5", &L1_SingleLooseIsoEG28er2p5),          std::make_pair("L1_SingleLooseIsoEG30er1p5", &L1_SingleLooseIsoEG30er1p5),          std::make_pair("L1_SingleLooseIsoEG30er2p5", &L1_SingleLooseIsoEG30er2p5),          std::make_pair("L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF),          std::make_pair("L1_SingleMu0_DQ", &L1_SingleMu0_DQ),          std::make_pair("L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF),          std::make_pair("L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF),          std::make_pair("L1_SingleMu10er1p5", &L1_SingleMu10er1p5),          std::make_pair("L1_SingleMu12_DQ_BMTF", &L1_SingleMu12_DQ_BMTF),          std::make_pair("L1_SingleMu12_DQ_EMTF", &L1_SingleMu12_DQ_EMTF),          std::make_pair("L1_SingleMu12_DQ_OMTF", &L1_SingleMu12_DQ_OMTF),          std::make_pair("L1_SingleMu12er1p5", &L1_SingleMu12er1p5),          std::make_pair("L1_SingleMu14er1p5", &L1_SingleMu14er1p5),          std::make_pair("L1_SingleMu15_DQ", &L1_SingleMu15_DQ),          std::make_pair("L1_SingleMu16er1p5", &L1_SingleMu16er1p5),          std::make_pair("L1_SingleMu18", &L1_SingleMu18),          std::make_pair("L1_SingleMu18er1p5", &L1_SingleMu18er1p5),          std::make_pair("L1_SingleMu20", &L1_SingleMu20),          std::make_pair("L1_SingleMu22", &L1_SingleMu22),          std::make_pair("L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF),          std::make_pair("L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF),          std::make_pair("L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF),          std::make_pair("L1_SingleMu25", &L1_SingleMu25),          std::make_pair("L1_SingleMu3", &L1_SingleMu3),          std::make_pair("L1_SingleMu5", &L1_SingleMu5),          std::make_pair("L1_SingleMu6er1p5", &L1_SingleMu6er1p5),          std::make_pair("L1_SingleMu7", &L1_SingleMu7),          std::make_pair("L1_SingleMu7_DQ", &L1_SingleMu7_DQ),          std::make_pair("L1_SingleMu7er1p5", &L1_SingleMu7er1p5),          std::make_pair("L1_SingleMu8er1p5", &L1_SingleMu8er1p5),          std::make_pair("L1_SingleMu9er1p5", &L1_SingleMu9er1p5),          std::make_pair("L1_SingleMuCosmics", &L1_SingleMuCosmics),          std::make_pair("L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF),          std::make_pair("L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF),          std::make_pair("L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF),          std::make_pair("L1_SingleMuOpen", &L1_SingleMuOpen),          std::make_pair("L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR),          std::make_pair("L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX),          std::make_pair("L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX),          std::make_pair("L1_SingleTau120er2p1", &L1_SingleTau120er2p1),          std::make_pair("L1_SingleTau130er2p1", &L1_SingleTau130er2p1),          std::make_pair("L1_TOTEM_1", &L1_TOTEM_1),          std::make_pair("L1_TOTEM_2", &L1_TOTEM_2),          std::make_pair("L1_TOTEM_3", &L1_TOTEM_3),          std::make_pair("L1_TOTEM_4", &L1_TOTEM_4),          std::make_pair("L1_TripleEG16er2p5", &L1_TripleEG16er2p5),          std::make_pair("L1_TripleEG_16_12_8_er2p5", &L1_TripleEG_16_12_8_er2p5),          std::make_pair("L1_TripleEG_16_15_8_er2p5", &L1_TripleEG_16_15_8_er2p5),          std::make_pair("L1_TripleEG_18_17_8_er2p5", &L1_TripleEG_18_17_8_er2p5),          std::make_pair("L1_TripleEG_18_18_12_er2p5", &L1_TripleEG_18_18_12_er2p5),          std::make_pair("L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5),          std::make_pair("L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5),          std::make_pair("L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5),          std::make_pair("L1_TripleMu0", &L1_TripleMu0),          std::make_pair("L1_TripleMu0_OQ", &L1_TripleMu0_OQ),          std::make_pair("L1_TripleMu0_SQ", &L1_TripleMu0_SQ),          std::make_pair("L1_TripleMu3", &L1_TripleMu3),          std::make_pair("L1_TripleMu3_SQ", &L1_TripleMu3_SQ),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ),          std::make_pair("L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair("L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9),          std::make_pair("L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3),          std::make_pair("L1_TripleMu_5_3_3_SQ", &L1_TripleMu_5_3_3_SQ),          std::make_pair("L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5),          std::make_pair("L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17),          std::make_pair("L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3),          std::make_pair("L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus),          std::make_pair("L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus),          std::make_pair("L1_ZeroBias", &L1_ZeroBias),          std::make_pair("L1_ZeroBias_copy", &L1_ZeroBias_copy)      };

  for (auto pair : name2func)
  {
    L1SeedFun[pair.first] = std::bind(pair.second, upgrade, calo_tower);
  }

  return true;
}
// eof
