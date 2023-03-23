/* automatically generated from L1Menu_Collisions2022_v1_4_0 with menu2lib.py */
/* https://github.com/cms-l1-dpg/L1MenuTools */

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
  static const int64_t cos_coeff[72] = {1023, 1019, 1007, 988, 961, 927, 886, 838, 784, 723, 658, 587, 512, 432, 350, 265, 178, 89, 0, -89, -178, -265, -350, -432, -512, -587, -658, -723, -784, -838, -886, -927, -961, -988, -1007, -1019, -1023, -1019, -1007, -988, -961, -927, -886, -838, -784, -723, -658, -587, -512, -432, -350, -265, -178, -89, 0, 89, 178, 265, 350, 432, 511, 587, 658, 723, 784, 838, 886, 927, 961, 988, 1007, 1019};
  static const int64_t sin_coeff[72] = {0, 89, 178, 265, 350, 432, 512, 587, 658, 723, 784, 838, 886, 927, 961, 988, 1007, 1019, 1023, 1019, 1007, 988, 961, 927, 886, 838, 784, 723, 658, 587, 512, 432, 350, 265, 178, 89, 0, -89, -178, -265, -350, -432, -512, -587, -658, -723, -784, -838, -886, -927, -961, -988, -1007, -1019, -1023, -1019, -1007, -988, -961, -927, -886, -838, -784, -723, -658, -587, -512, -432, -350, -265, -178, -89};

  if (not calo_tower) return {-1., -9999.};

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
        ex += (et * cos_coeff[index] / 1024.);
        ey += (et * sin_coeff[index] / 1024.);
      }
    }
  }

  return {sqrt(ex*ex + ey*ey), atan2(-ey, -ex)};
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

  mt = sqrt(2. * eg_et * met_et * (1. - cos(delta_phi)));
  return mt;
}


// utility factories

const CombinationFactory::data_t& CombinationFactory::get(const size_t n, const size_t k)
{
  const auto& rc = cache_.find(std::make_pair(n, k));
  if (rc != cache_.end())
    return rc->second;
  return insert(n, k);
}


void CombinationFactory::clear()
{
  cache_.clear();
}


const CombinationFactory::data_t& CombinationFactory::insert(const size_t n, const size_t k)
{
  data_t v;

  std::string bitmask(k, 1);
  bitmask.resize(n, 0);

  do
  {
    std::vector<size_t> set;
    set.reserve(n);
    for (size_t ii = 0; ii < n; ++ii)
    {
      if (bitmask[ii])
      {
        set.emplace_back(ii);
      }
    }
    v.emplace_back(set);
  }
  while (std::prev_permutation(bitmask.begin(), bitmask.end()));

  const auto key = std::make_pair(n, k);
  cache_.emplace(key, v);
  return cache_.at(key);
}


CombinationFactory::cache_t CombinationFactory::cache_ = {};


const PermutationFactory::data_t& PermutationFactory::get(const size_t n)
{
  const auto& rc = cache_.find(n);
  if (rc != cache_.end())
    return rc->second;
  return insert(n);
}


void PermutationFactory::clear()
{
  cache_.clear();
}


const PermutationFactory::data_t& PermutationFactory::insert(const size_t n)
{
  data_t v;

  std::vector<size_t> indicies(n);
  for (size_t ii = 0; ii < n; ++ii)
  {
    indicies[ii] = ii;
  }

  do
  {
    std::vector<size_t> set;
    set.reserve(n);
    for (size_t ii = 0; ii < n; ++ii)
    {
      set.emplace_back(indicies.at(ii));
    }
    v.emplace_back(set);
  }
  while (std::next_permutation(indicies.begin(), indicies.end()));

  cache_.emplace(n, v);
  return cache_.at(n);
}


PermutationFactory::cache_t PermutationFactory::cache_ = {};




//
// NB: tmEventSetup.XxxWithOverlapRemoval was removed between utm-overlapRemoval-xsd330 and utm_0.6.5
//
/////////////////////////
// Generate conditions //
/////////////////////////


    






bool
CaloCaloCorrelation_i105
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET32: ET >= 64 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 64)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

      {
        // JET32: ET >= 64 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 64)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1600; // 1.6 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i107
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1600; // 1.6 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i109
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1600; // 1.6 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i179
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;

      {
        // EG26: ET >= 52 at BX = 0
        if (not (data->egIEt.at(ii) >= 52)) continue;

        const auto eta = data->egIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }


    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
            
            {
        // JET34: ET >= 68 at BX = 0
        if (not (data->jetIEt.at(jj) >= 68)) continue;

        const auto eta = data->jetIEta.at(jj);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_JET[deltaIEta];
  
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_JET[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 90000; // 0.09 * 10^6
      const long long maximum = 139240000; // 139.24 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      






bool
CaloCaloCorrelation_i180
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;

      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(ii) >= 56)) continue;

        const auto eta = data->egIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }


    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
            
            {
        // JET34: ET >= 68 at BX = 0
        if (not (data->jetIEt.at(jj) >= 68)) continue;

        const auto eta = data->jetIEta.at(jj);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_JET[deltaIEta];
  
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_JET[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 90000; // 0.09 * 10^6
      const long long maximum = 139240000; // 139.24 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      






bool
CaloCaloCorrelation_i181
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;

      {
        // EG30: ET >= 60 at BX = 0
        if (not (data->egIEt.at(ii) >= 60)) continue;

        const auto eta = data->egIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }


    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
            
            {
        // JET34: ET >= 68 at BX = 0
        if (not (data->jetIEt.at(jj) >= 68)) continue;

        const auto eta = data->jetIEta.at(jj);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_JET[deltaIEta];
  
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_JET[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 90000; // 0.09 * 10^6
      const long long maximum = 139240000; // 139.24 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      






bool
CaloCaloCorrelation_i191
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;

      {
        // EG22: ET >= 44 at BX = 0
        if (not (data->egIEt.at(ii) >= 44)) continue;

        const auto eta = data->egIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }


    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
            
            {
        // TAU26: ET >= 52 at BX = 0
        if (not (data->tauIEt.at(jj) >= 52)) continue;

        const auto eta = data->tauIEta.at(jj);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(jj)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];
  
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_TAU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 90000; // 0.09 * 10^6
      const long long maximum = 139240000; // 139.24 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      






bool
CaloCaloCorrelation_i192
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;

      {
        // EG24: ET >= 48 at BX = 0
        if (not (data->egIEt.at(ii) >= 48)) continue;

        const auto eta = data->egIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }


    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
            
            {
        // TAU27: ET >= 54 at BX = 0
        if (not (data->tauIEt.at(jj) >= 54)) continue;

        const auto eta = data->tauIEta.at(jj);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(jj)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];
  
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_TAU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 90000; // 0.09 * 10^6
      const long long maximum = 139240000; // 139.24 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      






bool
CaloCaloCorrelation_i193
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;

      {
        // EG22: ET >= 44 at BX = 0
        if (not (data->egIEt.at(ii) >= 44)) continue;

        const auto eta = data->egIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }


    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
            
            {
        // TAU70: ET >= 140 at BX = 0
        if (not (data->tauIEt.at(jj) >= 140)) continue;

        const auto eta = data->tauIEta.at(jj);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

    int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.30 <= DeltaR <= 11.80
      iEta = data->egIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(jj));
      unsigned int deltaEta = LUT_DETA_EG_TAU[deltaIEta];
  
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_TAU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 90000; // 0.09 * 10^6
      const long long maximum = 139240000; // 139.24 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      






bool
CaloCaloCorrelation_i246
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET100: ET >= 200 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 200)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

      {
        // JET100: ET >= 200 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 200)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1600; // 1.6 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i247
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET112: ET >= 224 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 224)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

      {
        // JET112: ET >= 224 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 224)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1600; // 1.6 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i390
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG5: ET >= 10 at BX = 0
        if (not (data->egIEt.at(idx0) >= 10)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG5: ET >= 10 at BX = 0
        if (not (data->egIEt.at(idx1) >= 10)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.90
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 810000; // 0.81 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i391
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG5p5: ET >= 11 at BX = 0
        if (not (data->egIEt.at(idx0) >= 11)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG5p5: ET >= 11 at BX = 0
        if (not (data->egIEt.at(idx1) >= 11)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 641000; // 0.641 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i392
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG6: ET >= 12 at BX = 0
        if (not (data->egIEt.at(idx0) >= 12)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG6: ET >= 12 at BX = 0
        if (not (data->egIEt.at(idx1) >= 12)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 641000; // 0.641 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i393
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG6p5: ET >= 13 at BX = 0
        if (not (data->egIEt.at(idx0) >= 13)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG6p5: ET >= 13 at BX = 0
        if (not (data->egIEt.at(idx1) >= 13)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 641000; // 0.641 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i394
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG7: ET >= 14 at BX = 0
        if (not (data->egIEt.at(idx0) >= 14)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG7: ET >= 14 at BX = 0
        if (not (data->egIEt.at(idx1) >= 14)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.80
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 641000; // 0.641 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i395
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG7p5: ET >= 15 at BX = 0
        if (not (data->egIEt.at(idx0) >= 15)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG7p5: ET >= 15 at BX = 0
        if (not (data->egIEt.at(idx1) >= 15)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.70
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 491000; // 0.491 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i396
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG8: ET >= 16 at BX = 0
        if (not (data->egIEt.at(idx0) >= 16)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG8: ET >= 16 at BX = 0
        if (not (data->egIEt.at(idx1) >= 16)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.70
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 491000; // 0.491 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i397
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG8p5: ET >= 17 at BX = 0
        if (not (data->egIEt.at(idx0) >= 17)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG8p5: ET >= 17 at BX = 0
        if (not (data->egIEt.at(idx1) >= 17)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.70
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 491000; // 0.491 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i398
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG9: ET >= 18 at BX = 0
        if (not (data->egIEt.at(idx0) >= 18)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG9: ET >= 18 at BX = 0
        if (not (data->egIEt.at(idx1) >= 18)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.70
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 491000; // 0.491 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i399
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG9p5: ET >= 19 at BX = 0
        if (not (data->egIEt.at(idx0) >= 19)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG9p5: ET >= 19 at BX = 0
        if (not (data->egIEt.at(idx1) >= 19)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.60
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 361000; // 0.361 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i400
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx0) >= 20)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx1) >= 20)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.60
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 361000; // 0.361 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i401
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG10p5: ET >= 21 at BX = 0
        if (not (data->egIEt.at(idx0) >= 21)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG10p5: ET >= 21 at BX = 0
        if (not (data->egIEt.at(idx1) >= 21)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.60
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 361000; // 0.361 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i402
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG11: ET >= 22 at BX = 0
        if (not (data->egIEt.at(idx0) >= 22)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG11: ET >= 22 at BX = 0
        if (not (data->egIEt.at(idx1) >= 22)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.60
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 361000; // 0.361 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i403
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG4p5: ET >= 9 at BX = 0
        if (not (data->egIEt.at(idx0) >= 9)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG4p5: ET >= 9 at BX = 0
        if (not (data->egIEt.at(idx1) >= 9)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.90
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 810000; // 0.81 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      






bool
CaloCaloCorrelation_i404
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG4: ET >= 8 at BX = 0
        if (not (data->egIEt.at(idx0) >= 8)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

      {
        // EG4: ET >= 8 at BX = 0
        if (not (data->egIEt.at(idx1) >= 8)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -1.218 <= eta <= 1.218
        const bool etaWindow0 = ((-28 <= eta) and (eta <= 27));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 0.90
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_EG_EG[deltaIEta];
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_EG_EG[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 810000; // 0.81 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
CaloEsumCorrelation_i382
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
        
              {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(ii) >= 120)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }


    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEtHF)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
          {
                    // ETMHF90: ET >= 180 at BX = 0
      if (not (data->sumIEt.at(jj) >= 180)) continue;
      
  
    }


              // 2.094 <= DeltaPhi <= 3.142
      int iPhi = data->sumIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(ii));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_ETMHF[deltaIPhi];
  
    {
      const long long minimum = 2094; // 2.094 * 10^3
      const long long maximum = 3142; // 3.142 * 10^3
      if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;
    }

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      





bool
CaloEsumCorrelation_i383
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
        
              {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(ii) >= 120)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }


    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEtHF)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
          {
                    // ETMHF90: ET >= 180 at BX = 0
      if (not (data->sumIEt.at(jj) >= 180)) continue;
      
  
    }


              // 2.618 <= DeltaPhi <= 3.142
      int iPhi = data->sumIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(ii));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_ETMHF[deltaIPhi];
  
    {
      const long long minimum = 2618; // 2.618 * 10^3
      const long long maximum = 3142; // 3.142 * 10^3
      if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;
    }

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      





bool
CaloEsumCorrelation_i384
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
        
              {
        // JET80: ET >= 160 at BX = 0
        if (not (data->jetIEt.at(ii) >= 160)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }


    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEtHF)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
          {
                    // ETMHF90: ET >= 180 at BX = 0
      if (not (data->sumIEt.at(jj) >= 180)) continue;
      
  
    }


              // 2.094 <= DeltaPhi <= 3.142
      int iPhi = data->sumIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(ii));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_ETMHF[deltaIPhi];
  
    {
      const long long minimum = 2094; // 2.094 * 10^3
      const long long maximum = 3142; // 3.142 * 10^3
      if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;
    }

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      





bool
CaloEsumCorrelation_i385
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
        
              {
        // JET80: ET >= 160 at BX = 0
        if (not (data->jetIEt.at(ii) >= 160)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }


    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEtHF)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
          {
                    // ETMHF90: ET >= 180 at BX = 0
      if (not (data->sumIEt.at(jj) >= 180)) continue;
      
  
    }


              // 2.618 <= DeltaPhi <= 3.142
      int iPhi = data->sumIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(ii));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_ETMHF[deltaIPhi];
  
    {
      const long long minimum = 2618; // 2.618 * 10^3
      const long long maximum = 3142; // 3.142 * 10^3
      if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;
    }

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      





bool
CaloEsumCorrelation_i421
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
        
              {
        // JET55: ET >= 110 at BX = 0
        if (not (data->jetIEt.at(ii) >= 110)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }


    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEtHF)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
          {
                    // ETMHF80: ET >= 160 at BX = 0
      if (not (data->sumIEt.at(jj) >= 160)) continue;
      
  
    }


              // 2.094 <= DeltaPhi <= 3.142
      int iPhi = data->sumIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(ii));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_ETMHF[deltaIPhi];
  
    {
      const long long minimum = 2094; // 2.094 * 10^3
      const long long maximum = 3142; // 3.142 * 10^3
      if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;
    }

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      





bool
CaloEsumCorrelation_i422
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;
        
              {
        // JET55: ET >= 110 at BX = 0
        if (not (data->jetIEt.at(ii) >= 110)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }


    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEtHF)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
          {
                    // ETMHF80: ET >= 160 at BX = 0
      if (not (data->sumIEt.at(jj) >= 160)) continue;
      
  
    }


              // 2.618 <= DeltaPhi <= 3.142
      int iPhi = data->sumIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(ii));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_ETMHF[deltaIPhi];
  
    {
      const long long minimum = 2618; // 2.618 * 10^3
      const long long maximum = 3142; // 3.142 * 10^3
      if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;
    }

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i104
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET32: ET >= 64 at BX = 0
        if (not (data->jetIEt.at(ii) >= 64)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU10: ET >= 21 at BX = 0
        if (not (data->muonIEt.at(jj) >= 21)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(jj)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(jj);
        // -2.3000625 <= eta <= 2.3000624999999997
        const bool etaWindow0 = ((-211 <= eta) and (eta <= 211));

        if (not (etaWindow0)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 161000; // 0.161 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i106
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(ii) >= 80)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU12: ET >= 25 at BX = 0
        if (not (data->muonIEt.at(jj) >= 25)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(jj)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(jj);
        // -2.3000625 <= eta <= 2.3000624999999997
        const bool etaWindow0 = ((-211 <= eta) and (eta <= 211));

        if (not (etaWindow0)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 161000; // 0.161 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i108
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(ii) >= 80)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU12: ET >= 25 at BX = 0
        if (not (data->muonIEt.at(jj) >= 25)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(jj)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(jj);
        // -2.3000625 <= eta <= 2.3000624999999997
        const bool etaWindow0 = ((-211 <= eta) and (eta <= 211));

        if (not (etaWindow0)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 161000; // 0.161 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i111
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET90: ET >= 180 at BX = 0
        if (not (data->jetIEt.at(ii) >= 180)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(jj) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(jj)) & 1)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 641000; // 0.641 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i113
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET90: ET >= 180 at BX = 0
        if (not (data->jetIEt.at(ii) >= 180)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(jj) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(jj)) & 1)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 641000; // 0.641 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i92
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET16: ET >= 32 at BX = 0
        if (not (data->jetIEt.at(ii) >= 32)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(jj) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(jj)) & 1)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 161000; // 0.161 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i93
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(ii) >= 70)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(jj) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(jj)) & 1)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 161000; // 0.161 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i94
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(ii) >= 120)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(jj) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(jj)) & 1)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 161000; // 0.161 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i95
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET80: ET >= 160 at BX = 0
        if (not (data->jetIEt.at(ii) >= 160)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(jj) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(jj)) & 1)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 161000; // 0.161 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i96
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET120: ET >= 240 at BX = 0
        if (not (data->jetIEt.at(ii) >= 240)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(jj) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(jj)) & 1)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 641000; // 0.641 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      






bool
CaloMuonCorrelation_i97
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
    
              {
        // JET120: ET >= 240 at BX = 0
        if (not (data->jetIEt.at(ii) >= 240)) continue;

        const auto eta = data->jetIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
  
                  {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(jj) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(jj)) & 1)) continue;

      }

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
        const unsigned int deltaPhi = LUT_DPHI_JET_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 161000; // 0.161 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i162
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG15: ET >= 30 at BX = 0
        if (not (data->egIEt.at(idx) >= 30)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx) >= 20)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i163
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG20: ET >= 40 at BX = 0
        if (not (data->egIEt.at(idx) >= 40)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx) >= 20)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i164
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG22: ET >= 44 at BX = 0
        if (not (data->egIEt.at(idx) >= 44)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx) >= 20)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i165
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG25: ET >= 50 at BX = 0
        if (not (data->egIEt.at(idx) >= 50)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i166
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG25: ET >= 50 at BX = 0
        if (not (data->egIEt.at(idx) >= 50)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG14: ET >= 28 at BX = 0
        if (not (data->egIEt.at(idx) >= 28)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i167
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG27: ET >= 54 at BX = 0
        if (not (data->egIEt.at(idx) >= 54)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG14: ET >= 28 at BX = 0
        if (not (data->egIEt.at(idx) >= 28)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i168
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG20: ET >= 40 at BX = 0
        if (not (data->egIEt.at(idx) >= 40)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx) >= 20)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i169
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG22: ET >= 44 at BX = 0
        if (not (data->egIEt.at(idx) >= 44)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx) >= 20)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i170
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG22: ET >= 44 at BX = 0
        if (not (data->egIEt.at(idx) >= 44)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i171
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG25: ET >= 50 at BX = 0
        if (not (data->egIEt.at(idx) >= 50)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i172
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG22: ET >= 44 at BX = 0
        if (not (data->egIEt.at(idx) >= 44)) continue;

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG22: ET >= 44 at BX = 0
        if (not (data->egIEt.at(idx) >= 44)) continue;

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG24: ET >= 48 at BX = 0
        if (not (data->egIEt.at(idx) >= 48)) continue;

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG24: ET >= 48 at BX = 0
        if (not (data->egIEt.at(idx) >= 48)) continue;

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i186
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG8: ET >= 16 at BX = 0
        if (not (data->egIEt.at(idx) >= 16)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG8: ET >= 16 at BX = 0
        if (not (data->egIEt.at(idx) >= 16)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i354
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG16: ET >= 32 at BX = 0
        if (not (data->egIEt.at(idx) >= 32)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i355
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG18: ET >= 36 at BX = 0
        if (not (data->egIEt.at(idx) >= 36)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i356
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG20: ET >= 40 at BX = 0
        if (not (data->egIEt.at(idx) >= 40)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i357
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG22: ET >= 44 at BX = 0
        if (not (data->egIEt.at(idx) >= 44)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i358
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG25: ET >= 50 at BX = 0
        if (not (data->egIEt.at(idx) >= 50)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i84
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx) >= 20)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx) >= 20)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i85
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i86
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG15: ET >= 30 at BX = 0
        if (not (data->egIEt.at(idx) >= 30)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG15: ET >= 30 at BX = 0
        if (not (data->egIEt.at(idx) >= 30)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleEG_i87
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG17: ET >= 34 at BX = 0
        if (not (data->egIEt.at(idx) >= 34)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG17: ET >= 34 at BX = 0
        if (not (data->egIEt.at(idx) >= 34)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i117
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i243
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET100: ET >= 200 at BX = 0
        if (not (data->jetIEt.at(idx) >= 200)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET100: ET >= 200 at BX = 0
        if (not (data->jetIEt.at(idx) >= 200)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i244
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET120: ET >= 240 at BX = 0
        if (not (data->jetIEt.at(idx) >= 240)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET120: ET >= 240 at BX = 0
        if (not (data->jetIEt.at(idx) >= 240)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i245
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET150: ET >= 300 at BX = 0
        if (not (data->jetIEt.at(idx) >= 300)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET150: ET >= 300 at BX = 0
        if (not (data->jetIEt.at(idx) >= 300)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i254
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET90: ET >= 180 at BX = 0
        if (not (data->jetIEt.at(idx) >= 180)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx) >= 60)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i256
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET100: ET >= 200 at BX = 0
        if (not (data->jetIEt.at(idx) >= 200)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx) >= 60)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i257
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET110: ET >= 220 at BX = 0
        if (not (data->jetIEt.at(idx) >= 220)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(idx) >= 70)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i259
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET115: ET >= 230 at BX = 0
        if (not (data->jetIEt.at(idx) >= 230)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i261
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET120: ET >= 240 at BX = 0
        if (not (data->jetIEt.at(idx) >= 240)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx) >= 90)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i280
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET75: ET >= 150 at BX = 0
        if (not (data->jetIEt.at(idx) >= 150)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET65: ET >= 130 at BX = 0
        if (not (data->jetIEt.at(idx) >= 130)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i282
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET80: ET >= 160 at BX = 0
        if (not (data->jetIEt.at(idx) >= 160)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET70: ET >= 140 at BX = 0
        if (not (data->jetIEt.at(idx) >= 140)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i284
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET85: ET >= 170 at BX = 0
        if (not (data->jetIEt.at(idx) >= 170)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET75: ET >= 150 at BX = 0
        if (not (data->jetIEt.at(idx) >= 150)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleJET_i413
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

        // displaced jet bit : 0x1
        if (not ((0x1 & data->jetHwQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

        // displaced jet bit : 0x1
        if (not ((0x1 & data->jetHwQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i114
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i33
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i34
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i35
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i359
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 16 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 16)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 8 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 8)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i36
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
DoubleMU_i360
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 16 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 16)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 8 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 8)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i361
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 16 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 16)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 8 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 8)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i362
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 16 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 16)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 8 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 8)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i363
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // impact parameter : 0xe
        if (not ((0xe >> data->muonDxy.at(idx)) & 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 7 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 7)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 5 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 5)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i364
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // impact parameter : 0xe
        if (not ((0xe >> data->muonDxy.at(idx)) & 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 7 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 7)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 5 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 5)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i365
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // impact parameter : 0xe
        if (not ((0xe >> data->muonDxy.at(idx)) & 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 7 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 7)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 5 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 5)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i366
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // impact parameter : 0xe
        if (not ((0xe >> data->muonDxy.at(idx)) & 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 7 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 7)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 5 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 5)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i38
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU8: ET >= 17 at BX = 0
        if (not (data->muonIEt.at(idx) >= 17)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU8: ET >= 17 at BX = 0
        if (not (data->muonIEt.at(idx) >= 17)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i39
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU9: ET >= 19 at BX = 0
        if (not (data->muonIEt.at(idx) >= 19)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU9: ET >= 19 at BX = 0
        if (not (data->muonIEt.at(idx) >= 19)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i40
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU12: ET >= 25 at BX = 0
        if (not (data->muonIEt.at(idx) >= 25)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i41
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU15: ET >= 31 at BX = 0
        if (not (data->muonIEt.at(idx) >= 31)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i416
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 6 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 6)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 6 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 6)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i417
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 6 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 6)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 6 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 6)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i418
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 6 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 6)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 6 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 6)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i419
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.45 <= eta <= -1.2451875
        const bool etaWindow0 = ((-225 <= eta) and (eta <= -115));

        // 1.2451875 <= eta <= 2.45
        const bool etaWindow1 = ((115 <= eta) and (eta <= 225));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 6 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 6)) continue;

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        // MU0: UnconstrainedPt >= 6 at BX = 0
        if (not (data->muonIEtUnconstrained.at(idx) >= 6)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i42
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU15: ET >= 31 at BX = 0
        if (not (data->muonIEt.at(idx) >= 31)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU7: ET >= 15 at BX = 0
        if (not (data->muonIEt.at(idx) >= 15)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i43
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU15: ET >= 31 at BX = 0
        if (not (data->muonIEt.at(idx) >= 31)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU7: ET >= 15 at BX = 0
        if (not (data->muonIEt.at(idx) >= 15)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i45
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU18: ET >= 37 at BX = 0
        if (not (data->muonIEt.at(idx) >= 37)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.1043125000000003 <= eta <= 2.1043125
        const bool etaWindow0 = ((-193 <= eta) and (eta <= 193));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU18: ET >= 37 at BX = 0
        if (not (data->muonIEt.at(idx) >= 37)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.1043125000000003 <= eta <= 2.1043125
        const bool etaWindow0 = ((-193 <= eta) and (eta <= 193));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i48
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i49
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
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
DoubleMU_i53
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU4: ET >= 9 at BX = 0
        if (not (data->muonIEt.at(idx) >= 9)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU4: ET >= 9 at BX = 0
        if (not (data->muonIEt.at(idx) >= 9)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
DoubleMU_i55
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
DoubleMU_i57
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

      
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
DoubleMU_i88
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU4: ET >= 9 at BX = 0
        if (not (data->muonIEt.at(idx) >= 9)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU4: ET >= 9 at BX = 0
        if (not (data->muonIEt.at(idx) >= 9)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleMU_i90
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i196
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU70: ET >= 140 at BX = 0
        if (not (data->tauIEt.at(idx) >= 140)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // TAU70: ET >= 140 at BX = 0
        if (not (data->tauIEt.at(idx) >= 140)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i197
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU28: ET >= 56 at BX = 0
        if (not (data->tauIEt.at(idx) >= 56)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // TAU28: ET >= 56 at BX = 0
        if (not (data->tauIEt.at(idx) >= 56)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i198
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU30: ET >= 60 at BX = 0
        if (not (data->tauIEt.at(idx) >= 60)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // TAU30: ET >= 60 at BX = 0
        if (not (data->tauIEt.at(idx) >= 60)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i199
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU32: ET >= 64 at BX = 0
        if (not (data->tauIEt.at(idx) >= 64)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // TAU32: ET >= 64 at BX = 0
        if (not (data->tauIEt.at(idx) >= 64)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i200
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU34: ET >= 68 at BX = 0
        if (not (data->tauIEt.at(idx) >= 68)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // TAU34: ET >= 68 at BX = 0
        if (not (data->tauIEt.at(idx) >= 68)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i201
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU36: ET >= 72 at BX = 0
        if (not (data->tauIEt.at(idx) >= 72)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // TAU36: ET >= 72 at BX = 0
        if (not (data->tauIEt.at(idx) >= 72)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
DoubleTAU_i367
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU35: ET >= 70 at BX = 0
        if (not (data->tauIEt.at(idx) >= 70)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // TAU35: ET >= 70 at BX = 0
        if (not (data->tauIEt.at(idx) >= 70)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

            

  





bool
DoubleTauOvRm_i381
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // remove overlap -- reference: JET55
  std::vector<int> reference;
    size_t nref = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nref++;
        
          {
        // JET55: ET >= 110 at BX = 0
        if (not (data->jetIEt.at(ii) >= 110)) continue;

      }

    reference.emplace_back(ii);
  }
  if (not reference.size()) return false;

  bool pass = false;
  size_t nobj0 = 0;

  //Loop over leg1
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
      if (nobj0 > 12) break;

      {
        // TAU26: ET >= 52 at BX = 0
        if (not (data->tauIEt.at(ii) >= 52)) continue;

        const auto eta = data->tauIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }


    size_t nobj1 = 0;
    //Loop over leg2, starting from index of leg1 + 1 to avoid double counting pairs.
    for (size_t jj = ii+1 ; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
            

      {
        // TAU26: ET >= 52 at BX = 0
        if (not (data->tauIEt.at(jj) >= 52)) continue;

        const auto eta = data->tauIEta.at(jj);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(jj)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

//Check for correlation requirements between leg1 and leg2
  int iEta = -9999999; unsigned int deltaIEta = 9999999;
      

             // 0.00 <= DeltaR <= 0.50
      const long long cutDeltaR2Min = 0; // 0.0 * 10^6
      const long long cutDeltaR2Max = 250000; // 0.25 * 10^6
       
	//Loop over saved reference objects (leg3)
	for (size_t kk = 0; kk < reference.size(); kk++)
	  {
	    const int index = reference.at(kk);

	    //Check for correlation conditions between leg1-leg3 and leg2-leg3
	    	    		{
	        iEta = data->tauIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(index));
      unsigned int deltaEta = LUT_DETA_TAU_TAU[deltaIEta];
  
      int iPhi = data->tauIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_TAU_TAU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

		{
	        iEta = data->tauIEta.at(jj);
    deltaIEta = abs(iEta - data->jetIEta.at(index));
      unsigned int deltaEta = LUT_DETA_TAU_TAU[deltaIEta];
  
      int iPhi = data->tauIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_TAU_TAU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

	    	    
	      pass = true;
	      break;

	  }

      if (pass) break;
    }
    if (pass) break;
  }

  return pass;

}


            

  





bool
DoubleTauOvRm_i408
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // remove overlap -- reference: JET70
  std::vector<int> reference;
    size_t nref = 0;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nref++;
        
          {
        // JET70: ET >= 140 at BX = 0
        if (not (data->jetIEt.at(ii) >= 140)) continue;

      }

    reference.emplace_back(ii);
  }
  if (not reference.size()) return false;

  bool pass = false;
  size_t nobj0 = 0;

  //Loop over leg1
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
      if (nobj0 > 12) break;

      {
        // TAU26: ET >= 52 at BX = 0
        if (not (data->tauIEt.at(ii) >= 52)) continue;

        const auto eta = data->tauIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }


    size_t nobj1 = 0;
    //Loop over leg2, starting from index of leg1 + 1 to avoid double counting pairs.
    for (size_t jj = ii+1 ; jj < data->tauBx.size(); jj++)
    {
      if (not (data->tauBx.at(jj) == 0)) continue;
      nobj1++;
              if (nobj1 > 12) break;
            

      {
        // TAU26: ET >= 52 at BX = 0
        if (not (data->tauIEt.at(jj) >= 52)) continue;

        const auto eta = data->tauIEta.at(jj);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(jj)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

//Check for correlation requirements between leg1 and leg2
  int iEta = -9999999; unsigned int deltaIEta = 9999999;
      

             // 0.00 <= DeltaR <= 0.50
      const long long cutDeltaR2Min = 0; // 0.0 * 10^6
      const long long cutDeltaR2Max = 250000; // 0.25 * 10^6
       
	//Loop over saved reference objects (leg3)
	for (size_t kk = 0; kk < reference.size(); kk++)
	  {
	    const int index = reference.at(kk);

	    //Check for correlation conditions between leg1-leg3 and leg2-leg3
	    	    		{
	        iEta = data->tauIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(index));
      unsigned int deltaEta = LUT_DETA_TAU_TAU[deltaIEta];
  
      int iPhi = data->tauIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_TAU_TAU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

		{
	        iEta = data->tauIEta.at(jj);
    deltaIEta = abs(iEta - data->jetIEta.at(index));
      unsigned int deltaEta = LUT_DETA_TAU_TAU[deltaIEta];
  
      int iPhi = data->tauIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_TAU_TAU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

	    	    
	      pass = true;
	      break;

	  }

      if (pass) break;
    }
    if (pass) break;
  }

  return pass;

}


                        





bool
InvariantMass3_i374
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      const int idx2 = candidates.at(set.at(indicies.at(2)));
      {
        // MU2: ET >= 5 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 5)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

      }

      {
        // MU1p5: ET >= 4 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 4)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx2) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx2)) & 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx2) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not(fabs(data->muonChg.at(idx0)+data->muonChg.at(idx1)+data->muonChg.at(idx2)) == fabs(data->muonChg.at(idx0)))) continue;
          // 0.0 <= mass <= 12.0
      int ab_iEta = data->muonIEtaAtVtx.at(idx0);
    int ab_deltaIEta = abs(ab_iEta - data->muonIEtaAtVtx.at(idx1));
  
      int ab_iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int ab_deltaIPhi = abs(ab_iPhi - data->muonIPhiAtVtx.at(idx1));
  if (ab_deltaIPhi >= 288) ab_deltaIPhi = 2*288 - ab_deltaIPhi;
  
    const long long ab_coshDeltaEta = LUT_COSH_DETA_MU_MU[ab_deltaIEta];
    const long long ab_cosDeltaPhi = LUT_COS_DPHI_MU_MU[ab_deltaIPhi];
    const long long ab_pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long ab_pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long ab_mass2 = ab_pt0 * ab_pt1 * (ab_coshDeltaEta - ab_cosDeltaPhi);
      int ac_iEta = data->muonIEtaAtVtx.at(idx0);
    int ac_deltaIEta = abs(ac_iEta - data->muonIEtaAtVtx.at(idx2));
  
      int ac_iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int ac_deltaIPhi = abs(ac_iPhi - data->muonIPhiAtVtx.at(idx2));
  if (ac_deltaIPhi >= 288) ac_deltaIPhi = 2*288 - ac_deltaIPhi;
  
    const long long ac_coshDeltaEta = LUT_COSH_DETA_MU_MU[ac_deltaIEta];
    const long long ac_cosDeltaPhi = LUT_COS_DPHI_MU_MU[ac_deltaIPhi];
    const long long ac_pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long ac_pt1 = LUT_MU_ET[data->muonIEt.at(idx2)];
    const long long ac_mass2 = ac_pt0 * ac_pt1 * (ac_coshDeltaEta - ac_cosDeltaPhi);
      int bc_iEta = data->muonIEtaAtVtx.at(idx1);
    int bc_deltaIEta = abs(bc_iEta - data->muonIEtaAtVtx.at(idx2));
  
      int bc_iPhi = data->muonIPhiAtVtx.at(idx1);
  
  unsigned int bc_deltaIPhi = abs(bc_iPhi - data->muonIPhiAtVtx.at(idx2));
  if (bc_deltaIPhi >= 288) bc_deltaIPhi = 2*288 - bc_deltaIPhi;
  
    const long long bc_coshDeltaEta = LUT_COSH_DETA_MU_MU[bc_deltaIEta];
    const long long bc_cosDeltaPhi = LUT_COS_DPHI_MU_MU[bc_deltaIPhi];
    const long long bc_pt0 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long bc_pt1 = LUT_MU_ET[data->muonIEt.at(idx2)];
    const long long bc_mass2 = bc_pt0 * bc_pt1 * (bc_coshDeltaEta - bc_cosDeltaPhi);
    const long long mass3 = ab_mass2 + ac_mass2 + bc_mass2;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 72000000; // 72.0 * 10^6
      if (not ((minimum <= mass3) and (mass3 <= maximum))) continue;
    }

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


                            





bool
InvariantMass3_i411
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      const int idx2 = candidates.at(set.at(indicies.at(2)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

      }

      {
        // MU2p5: ET >= 6 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 6)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx2) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx2)) & 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx2) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not(fabs(data->muonChg.at(idx0)+data->muonChg.at(idx1)+data->muonChg.at(idx2)) == fabs(data->muonChg.at(idx0)))) continue;
          // 0.0 <= mass <= 12.0
      int ab_iEta = data->muonIEtaAtVtx.at(idx0);
    int ab_deltaIEta = abs(ab_iEta - data->muonIEtaAtVtx.at(idx1));
  
      int ab_iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int ab_deltaIPhi = abs(ab_iPhi - data->muonIPhiAtVtx.at(idx1));
  if (ab_deltaIPhi >= 288) ab_deltaIPhi = 2*288 - ab_deltaIPhi;
  
    const long long ab_coshDeltaEta = LUT_COSH_DETA_MU_MU[ab_deltaIEta];
    const long long ab_cosDeltaPhi = LUT_COS_DPHI_MU_MU[ab_deltaIPhi];
    const long long ab_pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long ab_pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long ab_mass2 = ab_pt0 * ab_pt1 * (ab_coshDeltaEta - ab_cosDeltaPhi);
      int ac_iEta = data->muonIEtaAtVtx.at(idx0);
    int ac_deltaIEta = abs(ac_iEta - data->muonIEtaAtVtx.at(idx2));
  
      int ac_iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int ac_deltaIPhi = abs(ac_iPhi - data->muonIPhiAtVtx.at(idx2));
  if (ac_deltaIPhi >= 288) ac_deltaIPhi = 2*288 - ac_deltaIPhi;
  
    const long long ac_coshDeltaEta = LUT_COSH_DETA_MU_MU[ac_deltaIEta];
    const long long ac_cosDeltaPhi = LUT_COS_DPHI_MU_MU[ac_deltaIPhi];
    const long long ac_pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long ac_pt1 = LUT_MU_ET[data->muonIEt.at(idx2)];
    const long long ac_mass2 = ac_pt0 * ac_pt1 * (ac_coshDeltaEta - ac_cosDeltaPhi);
      int bc_iEta = data->muonIEtaAtVtx.at(idx1);
    int bc_deltaIEta = abs(bc_iEta - data->muonIEtaAtVtx.at(idx2));
  
      int bc_iPhi = data->muonIPhiAtVtx.at(idx1);
  
  unsigned int bc_deltaIPhi = abs(bc_iPhi - data->muonIPhiAtVtx.at(idx2));
  if (bc_deltaIPhi >= 288) bc_deltaIPhi = 2*288 - bc_deltaIPhi;
  
    const long long bc_coshDeltaEta = LUT_COSH_DETA_MU_MU[bc_deltaIEta];
    const long long bc_cosDeltaPhi = LUT_COS_DPHI_MU_MU[bc_deltaIPhi];
    const long long bc_pt0 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long bc_pt1 = LUT_MU_ET[data->muonIEt.at(idx2)];
    const long long bc_mass2 = bc_pt0 * bc_pt1 * (bc_coshDeltaEta - bc_cosDeltaPhi);
    const long long mass3 = ab_mass2 + ac_mass2 + bc_mass2;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 72000000; // 72.0 * 10^6
      if (not ((minimum <= mass3) and (mass3 <= maximum))) continue;
    }

    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


                                      

  





bool
InvariantMassOvRm_i275
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
        
          {
        // TAU45: ET >= 90 at BX = 0
        if (not (data->tauIEt.at(ii) >= 90)) continue;

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(ii)) & 1)) continue;

      }

    reference.emplace_back(ii);
  }
  if (not reference.size()) return false;

  bool pass = false;
  size_t nobj0 = 0;

  //Loop over leg1
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;

      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(ii) >= 70)) continue;

      }


    size_t nobj1 = 0;
    //Loop over leg2, starting from index of leg1 + 1 to avoid double counting pairs.
    for (size_t jj = ii+1 ; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
            

      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(jj) >= 70)) continue;

      }

//Check for correlation requirements between leg1 and leg2
  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 450.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
  
      int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(ii)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(jj)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 10125000000; // 101250.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


        

                    // 0.00 <= DeltaR <= 0.20
      const long long cutDeltaR2Min = 0; // 0.0 * 10^6
      const long long cutDeltaR2Max = 41000; // 0.041 * 10^6
       
	//Loop over saved reference objects (leg3)
	for (size_t kk = 0; kk < reference.size(); kk++)
	  {
	    const int index = reference.at(kk);

	    //Check for correlation conditions between leg1-leg3 and leg2-leg3
	    	    	    	    		{
	        iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
      int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_JET[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

		{
	        iEta = data->jetIEta.at(jj);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
      int iPhi = data->jetIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_JET[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

	    	    
	      pass = true;
	      break;

	  }

      if (pass) break;
    }
    if (pass) break;
  }

  return pass;

}

	
          
                                  

  





bool
InvariantMassOvRm_i276
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
        
          {
        // TAU40: ET >= 80 at BX = 0
        if (not (data->tauIEt.at(ii) >= 80)) continue;

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(ii)) & 1)) continue;

      }

    reference.emplace_back(ii);
  }
  if (not reference.size()) return false;

  bool pass = false;
  size_t nobj0 = 0;

  //Loop over leg1
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;

      {
        // JET80: ET >= 160 at BX = 0
        if (not (data->jetIEt.at(ii) >= 160)) continue;

      }


    size_t nobj1 = 0;
    //Loop over leg2, starting from index of leg1 + 1 to avoid double counting pairs.
    for (size_t jj = ii+1 ; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
            

      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(jj) >= 60)) continue;

      }

//Check for correlation requirements between leg1 and leg2
  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 420.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
  
      int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(ii)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(jj)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 8820000000; // 88200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


        

                    // 0.00 <= DeltaR <= 0.20
      const long long cutDeltaR2Min = 0; // 0.0 * 10^6
      const long long cutDeltaR2Max = 41000; // 0.041 * 10^6
       
	//Loop over saved reference objects (leg3)
	for (size_t kk = 0; kk < reference.size(); kk++)
	  {
	    const int index = reference.at(kk);

	    //Check for correlation conditions between leg1-leg3 and leg2-leg3
	    	    	    	    		{
	        iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
      int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_JET[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

		{
	        iEta = data->jetIEta.at(jj);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
      int iPhi = data->jetIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_JET[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

	    	    
	      pass = true;
	      break;

	  }

      if (pass) break;
    }
    if (pass) break;
  }

  return pass;

}

	
          
                                  

  





bool
InvariantMassOvRm_i368
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
        
          {
        // TAU45: ET >= 90 at BX = 0
        if (not (data->tauIEt.at(ii) >= 90)) continue;

        const auto eta = data->tauIEta.at(ii);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

    reference.emplace_back(ii);
  }
  if (not reference.size()) return false;

  bool pass = false;
  size_t nobj0 = 0;

  //Loop over leg1
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;

      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(ii) >= 70)) continue;

      }


    size_t nobj1 = 0;
    //Loop over leg2, starting from index of leg1 + 1 to avoid double counting pairs.
    for (size_t jj = ii+1 ; jj < data->jetBx.size(); jj++)
    {
      if (not (data->jetBx.at(jj) == 0)) continue;
      nobj1++;
            

      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(jj) >= 70)) continue;

      }

//Check for correlation requirements between leg1 and leg2
  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 450.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->jetIEta.at(jj));
  
      int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(ii)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(jj)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 10125000000; // 101250.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


        

                    // 0.00 <= DeltaR <= 0.50
      const long long cutDeltaR2Min = 0; // 0.0 * 10^6
      const long long cutDeltaR2Max = 250000; // 0.25 * 10^6
       
	//Loop over saved reference objects (leg3)
	for (size_t kk = 0; kk < reference.size(); kk++)
	  {
	    const int index = reference.at(kk);

	    //Check for correlation conditions between leg1-leg3 and leg2-leg3
	    	    	    	    		{
	        iEta = data->jetIEta.at(ii);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
      int iPhi = data->jetIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_JET[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

		{
	        iEta = data->jetIEta.at(jj);
    deltaIEta = abs(iEta - data->tauIEta.at(index));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
      int iPhi = data->jetIPhi.at(jj);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(index));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_JET_JET[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
	      if ((cutDeltaR2Min <= deltaR2) and (deltaR2 <= cutDeltaR2Max)) continue;
		}

	    	    
	      pass = true;
	      break;

	  }

      if (pass) break;
    }
    if (pass) break;
  }

  return pass;

}

	
          
                            






bool
InvariantMass_i202
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
        if (nobj0 > 12) break;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // TAU28: ET >= 56 at BX = 0
        if (not (data->tauIEt.at(idx0) >= 56)) continue;

        const auto eta = data->tauIEta.at(idx0);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // TAU28: ET >= 56 at BX = 0
        if (not (data->tauIEt.at(idx1) >= 56)) continue;

        const auto eta = data->tauIEta.at(idx1);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 90.0
      iEta = data->tauIEta.at(idx0);
    deltaIEta = abs(iEta - data->tauIEta.at(idx1));
  
      int iPhi = data->tauIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_TAU_TAU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_TAU_TAU[deltaIPhi];
    const long long pt0 = LUT_TAU_ET[data->tauIEt.at(idx0)];
    const long long pt1 = LUT_TAU_ET[data->tauIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 0; // 0.0 * 10^5
      const long long maximum = 405000000; // 4050.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
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
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
        if (nobj0 > 12) break;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // TAU28: ET >= 56 at BX = 0
        if (not (data->tauIEt.at(idx0) >= 56)) continue;

        const auto eta = data->tauIEta.at(idx0);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // TAU28: ET >= 56 at BX = 0
        if (not (data->tauIEt.at(idx1) >= 56)) continue;

        const auto eta = data->tauIEta.at(idx1);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 80.0
      iEta = data->tauIEta.at(idx0);
    deltaIEta = abs(iEta - data->tauIEta.at(idx1));
  
      int iPhi = data->tauIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_TAU_TAU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_TAU_TAU[deltaIPhi];
    const long long pt0 = LUT_TAU_ET[data->tauIEt.at(idx0)];
    const long long pt1 = LUT_TAU_ET[data->tauIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 0; // 0.0 * 10^5
      const long long maximum = 320000000; // 3200.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i204
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
        if (nobj0 > 12) break;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // TAU30: ET >= 60 at BX = 0
        if (not (data->tauIEt.at(idx0) >= 60)) continue;

        const auto eta = data->tauIEta.at(idx0);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // TAU30: ET >= 60 at BX = 0
        if (not (data->tauIEt.at(idx1) >= 60)) continue;

        const auto eta = data->tauIEta.at(idx1);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 90.0
      iEta = data->tauIEta.at(idx0);
    deltaIEta = abs(iEta - data->tauIEta.at(idx1));
  
      int iPhi = data->tauIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_TAU_TAU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_TAU_TAU[deltaIPhi];
    const long long pt0 = LUT_TAU_ET[data->tauIEt.at(idx0)];
    const long long pt1 = LUT_TAU_ET[data->tauIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 0; // 0.0 * 10^5
      const long long maximum = 405000000; // 4050.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
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
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj0++;
        if (nobj0 > 12) break;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // TAU30: ET >= 60 at BX = 0
        if (not (data->tauIEt.at(idx0) >= 60)) continue;

        const auto eta = data->tauIEta.at(idx0);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // TAU30: ET >= 60 at BX = 0
        if (not (data->tauIEt.at(idx1) >= 60)) continue;

        const auto eta = data->tauIEta.at(idx1);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 80.0
      iEta = data->tauIEta.at(idx0);
    deltaIEta = abs(iEta - data->tauIEta.at(idx1));
  
      int iPhi = data->tauIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->tauIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_TAU_TAU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_TAU_TAU[deltaIPhi];
    const long long pt0 = LUT_TAU_ET[data->tauIEt.at(idx0)];
    const long long pt1 = LUT_TAU_ET[data->tauIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 0; // 0.0 * 10^5
      const long long maximum = 320000000; // 3200.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i248
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 150.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 1125000000; // 11250.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1500; // 1.5 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i249
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1500; // 1.5 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


          // 200.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 2000000000; // 20000.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i250
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1500; // 1.5 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


          // 250.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 3125000000; // 31250.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i251
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 300.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 4500000000; // 45000.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1500; // 1.5 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i252
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 330.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 5445000000; // 54450.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1500; // 1.5 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i253
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 360.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 6480000000; // 64800.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


          // 0.0 <= DeltaEta <= 1.5
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
      unsigned int deltaEta = LUT_DETA_JET_JET[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1500; // 1.5 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i255
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 60)) continue;

      }

      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 60)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i258
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 70)) continue;

      }

      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 70)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i260
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 80)) continue;

      }

      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 80)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i262
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 90)) continue;

      }

      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 90)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i264
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.6969999999999996 <= eta <= 2.6969999999999996
        const bool etaWindow0 = ((-62 <= eta) and (eta <= 61));

        if (not (etaWindow0)) continue;

      }

      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.6969999999999996 <= eta <= 2.6969999999999996
        const bool etaWindow0 = ((-62 <= eta) and (eta <= 61));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
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
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 120)) continue;

      }

      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 120)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
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
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 120)) continue;

      }

      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow1 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i267
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 120)) continue;

      }

      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.6969999999999996 <= eta <= 2.6969999999999996
        const bool etaWindow0 = ((-62 <= eta) and (eta <= 61));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i268
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow1 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.6969999999999996 <= eta <= 2.6969999999999996
        const bool etaWindow0 = ((-62 <= eta) and (eta <= 61));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i269
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow1 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow1 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
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
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -2.6969999999999996 <= eta <= 2.6969999999999996
        const bool etaWindow0 = ((-62 <= eta) and (eta <= 61));

        if (not (etaWindow0)) continue;

      }

      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.6969999999999996 <= eta <= 2.6969999999999996
        const bool etaWindow0 = ((-62 <= eta) and (eta <= 61));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i271
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 120)) continue;

      }

      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow1 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i272
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 120)) continue;

      }

      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.6969999999999996 <= eta <= 2.6969999999999996
        const bool etaWindow0 = ((-62 <= eta) and (eta <= 61));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i273
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow1 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -2.6969999999999996 <= eta <= 2.6969999999999996
        const bool etaWindow0 = ((-62 <= eta) and (eta <= 61));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i274
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx0);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow1 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx1);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow1 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 620.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 19220000000; // 192200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                            






bool
InvariantMass_i277
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // JET80: ET >= 160 at BX = 0
        if (not (data->jetIEt.at(idx0) >= 160)) continue;

      }

      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx1) >= 60)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 420.0 <= mass <= 151982.0
      iEta = data->jetIEta.at(idx0);
    deltaIEta = abs(iEta - data->jetIEta.at(idx1));
  
      int iPhi = data->jetIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->jetIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_JET_JET[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_JET_JET[deltaIPhi];
    const long long pt0 = LUT_JET_ET[data->jetIEt.at(idx0)];
    const long long pt1 = LUT_JET_ET[data->jetIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 8820000000; // 88200.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                    





bool
InvariantMass_i37
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 1.0 <= mass <= 151982.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
    const long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 500000; // 0.5 * 10^6
      const long long maximum = 11549264162000000; // 11549264162.0 * 10^6
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i373
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 7.0 <= mass <= 18.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
    const long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 24500000; // 24.5 * 10^6
      const long long maximum = 162000000; // 162.0 * 10^6
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


          if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i377
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx0)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.3000625 <= eta <= 2.3000624999999997
        const bool etaWindow0 = ((-211 <= eta) and (eta <= 211));

        if (not (etaWindow0)) continue;

      }

      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx1)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.3000625 <= eta <= 2.3000624999999997
        const bool etaWindow0 = ((-211 <= eta) and (eta <= 211));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 14.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
    const long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 98000000; // 98.0 * 10^6
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


          if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                            






bool
InvariantMass_i378
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG7p5: ET >= 15 at BX = 0
        if (not (data->egIEt.at(idx0) >= 15)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      {
        // EG7p5: ET >= 15 at BX = 0
        if (not (data->egIEt.at(idx1) >= 15)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 20.0
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_EG_EG[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_EG_EG[deltaIPhi];
    const long long pt0 = LUT_EG_ET[data->egIEt.at(idx0)];
    const long long pt1 = LUT_EG_ET[data->egIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 0; // 0.0 * 10^5
      const long long maximum = 20000000; // 200.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                    





bool
InvariantMass_i379
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 11)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx0)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.3000625 <= eta <= 2.3000624999999997
        const bool etaWindow0 = ((-211 <= eta) and (eta <= 211));

        if (not (etaWindow0)) continue;

      }

      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 11)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx1)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.3000625 <= eta <= 2.3000624999999997
        const bool etaWindow0 = ((-211 <= eta) and (eta <= 211));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 8.0 <= mass <= 14.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
    const long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 32000000; // 32.0 * 10^6
      const long long maximum = 98000000; // 98.0 * 10^6
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


          if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                            






bool
InvariantMass_i380
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    size_t nobj0 = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj0++;
             candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // EG3: ET >= 6 at BX = 0
        if (not (data->egIEt.at(idx0) >= 6)) continue;

        const auto eta = data->egIEta.at(idx0);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      {
        // EG3: ET >= 6 at BX = 0
        if (not (data->egIEt.at(idx1) >= 6)) continue;

        const auto eta = data->egIEta.at(idx1);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= mass <= 20.0
      iEta = data->egIEta.at(idx0);
    deltaIEta = abs(iEta - data->egIEta.at(idx1));
  
      int iPhi = data->egIPhi.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->egIPhi.at(idx1));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_EG_EG[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_EG_EG[deltaIPhi];
    const long long pt0 = LUT_EG_ET[data->egIEt.at(idx0)];
    const long long pt1 = LUT_EG_ET[data->egIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 0; // 0.0 * 10^5
      const long long maximum = 20000000; // 200.0 * 10^5
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


	
          
                    





bool
InvariantMass_i44
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU15: ET >= 31 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 31)) continue;

      }

      {
        // MU7: ET >= 15 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 15)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 1.0 <= mass <= 151982.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
    const long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 500000; // 0.5 * 10^6
      const long long maximum = 11549264162000000; // 11549264162.0 * 10^6
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i58
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 7.0 <= mass <= 151982.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
    const long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 24500000; // 24.5 * 10^6
      const long long maximum = 11549264162000000; // 11549264162.0 * 10^6
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i70
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 11)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx0)) & 1)) continue;

      }

      {
        // MU2p5: ET >= 6 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 6)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx1)) & 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 5.0 <= mass <= 17.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
    const long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 12500000; // 12.5 * 10^6
      const long long maximum = 144500000; // 144.5 * 10^6
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


          if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i71
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 11)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx0)) & 1)) continue;

      }

      {
        // MU2p5: ET >= 6 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 6)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx1)) & 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 5.0 <= mass <= 17.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
    const long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 12500000; // 12.5 * 10^6
      const long long maximum = 144500000; // 144.5 * 10^6
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


          if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
                    





bool
InvariantMass_i73
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 11)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

      }

      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 0.0 <= mass <= 9.0
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
  
    const long long coshDeltaEta = LUT_COSH_DETA_MU_MU[deltaIEta];
    const long long cosDeltaPhi = LUT_COS_DPHI_MU_MU[deltaIPhi];
    const long long pt0 = LUT_MU_ET[data->muonIEt.at(idx0)];
    const long long pt1 = LUT_MU_ET[data->muonIEt.at(idx1)];
    const long long mass2 = pt0 * pt1 * (coshDeltaEta - cosDeltaPhi);
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 40500000; // 40.5 * 10^6
      if (not ((minimum <= mass2) and (mass2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}


    
      





bool
MuonMuonCorrelation_i110
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx0)) & 1)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx1)) & 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.60
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 2561000; // 2.561 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i112
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx0)) & 1)) continue;

      }

      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx1)) & 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.60
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 2561000; // 2.561 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i344
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj0 = 0;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == -1)) continue;
    nobj0++;

        const int idx0 = ii;
      {
        // MU3-1: ET >= 7 at BX = -1
        if (not (data->muonIEt.at(idx0) >= 7)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -1.2016874999999998 <= eta <= 1.2016875
        const bool etaWindow0 = ((-110 <= eta) and (eta <= 110));

        const auto phi = data->muonIPhiAtVtx.at(idx0);
        // 0.5235987755982988 <= phi <= 2.6179938779914944
        const bool phiWindow0 = ((48 <= phi) and (phi <= 239));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

        if (not (phiWindow0)) continue;

      }


    size_t nobj1 = 0;
    for (size_t jj = 0; jj < data->muonBx.size(); jj++)
    {
      if (not (data->muonBx.at(jj) == 0)) continue;
      nobj1++;
      
            const int idx1 = jj;
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 7)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -1.2016874999999998 <= eta <= 1.2016875
        const bool etaWindow0 = ((-110 <= eta) and (eta <= 110));

        const auto phi = data->muonIPhiAtVtx.at(idx1);
        // 3.665191429188092 <= phi <= 5.759586531581287
        const bool phiWindow0 = ((336 <= phi) and (phi <= 527));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

        if (not (phiWindow0)) continue;

      }

  
        // 2.618 <= DeltaPhi <= 3.142
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    {
      const long long minimum = 2618; // 2.618 * 10^3
      const long long maximum = 3142; // 3.142 * 10^3
      if (not ((minimum <= deltaPhi) and (deltaPhi <= maximum))) continue;
    }


    
      pass = true;
      break;
    }
    if (pass) break;
  }

  return pass;
}


      





bool
MuonMuonCorrelation_i351
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1500; // 1.5 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


          if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i352
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 0.0 <= DeltaEta <= 1.6
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1600; // 1.6 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i353
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 7)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 7)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 1961000; // 1.961 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i372
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -1.4083124999999999 <= eta <= 1.4083124999999999
        const bool etaWindow0 = ((-129 <= eta) and (eta <= 129));

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -1.4083124999999999 <= eta <= 1.4083124999999999
        const bool etaWindow0 = ((-129 <= eta) and (eta <= 129));

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 0.0 <= DeltaEta <= 1.6
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1600; // 1.6 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i406
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.6
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1600; // 1.6 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i407
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.0 <= DeltaEta <= 1.5
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
    {
      const long long minimum = 0; // 0.0 * 10^3
      const long long maximum = 1500; // 1.5 * 10^3
      if (not ((minimum <= deltaEta) and (deltaEta <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i46
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 1961000; // 1.961 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i47
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -2.0064375 <= eta <= 2.0064374999999997
        const bool etaWindow0 = ((-184 <= eta) and (eta <= 184));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 1961000; // 1.961 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i50
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 1961000; // 1.961 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i51
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 1961000; // 1.961 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i52
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx0);
        // -1.4083124999999999 <= eta <= 1.4083124999999999
        const bool etaWindow0 = ((-129 <= eta) and (eta <= 129));

        if (not (etaWindow0)) continue;

      }

      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx1);
        // -1.4083124999999999 <= eta <= 1.4083124999999999
        const bool etaWindow0 = ((-129 <= eta) and (eta <= 129));

        if (not (etaWindow0)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 0.00 <= DeltaR <= 1.40
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 1961000; // 1.961 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i54
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU4: ET >= 9 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 9)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

      }

      {
        // MU4: ET >= 9 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 9)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 0.00 <= DeltaR <= 1.20
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 1441000; // 1.441 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      





bool
MuonMuonCorrelation_i56
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;
     candidates.emplace_back(ii);
  }

  bool pass = false;
  if (candidates.size() < 2) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 2);
  const auto& permutation = PermutationFactory::get(2);

  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      const int idx0 = candidates.at(set.at(indicies.at(0)));
      const int idx1 = candidates.at(set.at(indicies.at(1)));
      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx0) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx0)) & 1)) continue;

      }

      {
        // MU4p5: ET >= 10 at BX = 0
        if (not (data->muonIEt.at(idx1) >= 10)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx1)) & 1)) continue;

      }

  int iEta = -9999999; unsigned int deltaIEta = 9999999;
        if (data->muonChg.at(idx0) == 0) continue;  // charge valid bit not set
    if (data->muonChg.at(idx1) == 0) continue;  // charge valid bit not set
    // opposite-sign (os)
    if (not (data->muonChg.at(idx0) != data->muonChg.at(idx1))) continue;
          // 0.00 <= DeltaR <= 1.20
      iEta = data->muonIEtaAtVtx.at(idx0);
    deltaIEta = abs(iEta - data->muonIEtaAtVtx.at(idx1));
      unsigned int deltaEta = LUT_DETA_MU_MU[deltaIEta];
  
      int iPhi = data->muonIPhiAtVtx.at(idx0);
  
  unsigned int deltaIPhi = abs(iPhi - data->muonIPhiAtVtx.at(idx1));
  if (deltaIPhi >= 288) deltaIPhi = 2*288 - deltaIPhi;
        const unsigned int deltaPhi = LUT_DPHI_MU_MU[deltaIPhi];
  
    const long long deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
    {
      const long long minimum = 0; // 0.0 * 10^6
      const long long maximum = 1441000; // 1.441 * 10^6
      if (not ((minimum <= deltaR2) and (deltaR2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}



      

bool
MuonShower0_i375
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

    if (data->nMuonShowers >= 1) {

        if (data->muonShowerBx.at(0) == 0)
      {
        if (data->muonShowerOneNominal.at(0))
          {
            pass = true;
          }
        }
  }
  return pass;
}

      

bool
MuonShower1_i376
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

    if (data->nMuonShowers >= 1) {

        if (data->muonShowerBx.at(0) == 0)
      {
        if (data->muonShowerOneTight.at(0))
          {
            pass = true;
          }
        }
  }
  return pass;
}

      


bool
QuadJET_i285
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET95: ET >= 190 at BX = 0
        if (not (data->jetIEt.at(idx) >= 190)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET75: ET >= 150 at BX = 0
        if (not (data->jetIEt.at(idx) >= 150)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET65: ET >= 130 at BX = 0
        if (not (data->jetIEt.at(idx) >= 130)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // JET20: ET >= 40 at BX = 0
        if (not (data->jetIEt.at(idx) >= 40)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i288
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i289
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET80: ET >= 160 at BX = 0
        if (not (data->jetIEt.at(idx) >= 160)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i290
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET80: ET >= 160 at BX = 0
        if (not (data->jetIEt.at(idx) >= 160)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET50: ET >= 100 at BX = 0
        if (not (data->jetIEt.at(idx) >= 100)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // JET45: ET >= 90 at BX = 0
        if (not (data->jetIEt.at(idx) >= 90)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3055 <= eta <= 2.3055
        const bool etaWindow0 = ((-53 <= eta) and (eta <= 52));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i386
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET36: ET >= 72 at BX = 0
        if (not (data->jetIEt.at(idx) >= 72)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET36: ET >= 72 at BX = 0
        if (not (data->jetIEt.at(idx) >= 72)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET36: ET >= 72 at BX = 0
        if (not (data->jetIEt.at(idx) >= 72)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // JET36: ET >= 72 at BX = 0
        if (not (data->jetIEt.at(idx) >= 72)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i388
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET70: ET >= 140 at BX = 0
        if (not (data->jetIEt.at(idx) >= 140)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3924999999999996 <= eta <= 2.3924999999999996
        const bool etaWindow0 = ((-55 <= eta) and (eta <= 54));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET55: ET >= 110 at BX = 0
        if (not (data->jetIEt.at(idx) >= 110)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3924999999999996 <= eta <= 2.3924999999999996
        const bool etaWindow0 = ((-55 <= eta) and (eta <= 54));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3924999999999996 <= eta <= 2.3924999999999996
        const bool etaWindow0 = ((-55 <= eta) and (eta <= 54));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(idx) >= 70)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3924999999999996 <= eta <= 2.3924999999999996
        const bool etaWindow0 = ((-55 <= eta) and (eta <= 54));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadJET_i389
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET70: ET >= 140 at BX = 0
        if (not (data->jetIEt.at(idx) >= 140)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3924999999999996 <= eta <= 2.3924999999999996
        const bool etaWindow0 = ((-55 <= eta) and (eta <= 54));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET55: ET >= 110 at BX = 0
        if (not (data->jetIEt.at(idx) >= 110)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3924999999999996 <= eta <= 2.3924999999999996
        const bool etaWindow0 = ((-55 <= eta) and (eta <= 54));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3924999999999996 <= eta <= 2.3924999999999996
        const bool etaWindow0 = ((-55 <= eta) and (eta <= 54));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.3924999999999996 <= eta <= 2.3924999999999996
        const bool etaWindow0 = ((-55 <= eta) and (eta <= 54));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadMU_i75
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadMU_i76
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
QuadMU_i77
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 4) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 4);
  const auto& permutation = PermutationFactory::get(4);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(3)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i121
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG8: ET >= 16 at BX = 0
        if (not (data->egIEt.at(idx) >= 16)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG15: ET >= 30 at BX = 0
        if (not (data->egIEt.at(idx) >= 30)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i123
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG26: ET >= 52 at BX = 0
        if (not (data->egIEt.at(idx) >= 52)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // 2.5229999999999997 <= eta <= 5.0
        const bool etaWindow0 = ((58 <= eta) and (eta <= 114));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -5.0 <= eta <= -2.5229999999999997
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -59));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i129
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG34: ET >= 68 at BX = 0
        if (not (data->egIEt.at(idx) >= 68)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i130
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG36: ET >= 72 at BX = 0
        if (not (data->egIEt.at(idx) >= 72)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i131
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG38: ET >= 76 at BX = 0
        if (not (data->egIEt.at(idx) >= 76)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i132
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG40: ET >= 80 at BX = 0
        if (not (data->egIEt.at(idx) >= 80)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i133
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG42: ET >= 84 at BX = 0
        if (not (data->egIEt.at(idx) >= 84)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i134
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG45: ET >= 90 at BX = 0
        if (not (data->egIEt.at(idx) >= 90)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i135
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG50: ET >= 100 at BX = 0
        if (not (data->egIEt.at(idx) >= 100)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i136
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG60: ET >= 120 at BX = 0
        if (not (data->egIEt.at(idx) >= 120)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i137
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG26: ET >= 52 at BX = 0
        if (not (data->egIEt.at(idx) >= 52)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i138
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG26: ET >= 52 at BX = 0
        if (not (data->egIEt.at(idx) >= 52)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i139
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // 2.5229999999999997 <= eta <= 5.0
        const bool etaWindow0 = ((58 <= eta) and (eta <= 114));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i140
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -5.0 <= eta <= -2.5229999999999997
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -59));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i141
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i142
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i143
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i144
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG30: ET >= 60 at BX = 0
        if (not (data->egIEt.at(idx) >= 60)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i145
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG30: ET >= 60 at BX = 0
        if (not (data->egIEt.at(idx) >= 60)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i146
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG24: ET >= 48 at BX = 0
        if (not (data->egIEt.at(idx) >= 48)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i147
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG24: ET >= 48 at BX = 0
        if (not (data->egIEt.at(idx) >= 48)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i148
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG26: ET >= 52 at BX = 0
        if (not (data->egIEt.at(idx) >= 52)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i149
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG26: ET >= 52 at BX = 0
        if (not (data->egIEt.at(idx) >= 52)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i150
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG26: ET >= 52 at BX = 0
        if (not (data->egIEt.at(idx) >= 52)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i151
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // 2.5229999999999997 <= eta <= 5.0
        const bool etaWindow0 = ((58 <= eta) and (eta <= 114));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i152
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -5.0 <= eta <= -2.5229999999999997
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -59));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i153
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i154
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i155
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG28: ET >= 56 at BX = 0
        if (not (data->egIEt.at(idx) >= 56)) continue;

        const auto eta = data->egIEta.at(idx);
        // -1.5225 <= eta <= 1.5225
        const bool etaWindow0 = ((-35 <= eta) and (eta <= 34));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i156
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG30: ET >= 60 at BX = 0
        if (not (data->egIEt.at(idx) >= 60)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i157
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG30: ET >= 60 at BX = 0
        if (not (data->egIEt.at(idx) >= 60)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i158
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG32: ET >= 64 at BX = 0
        if (not (data->egIEt.at(idx) >= 64)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i159
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG32: ET >= 64 at BX = 0
        if (not (data->egIEt.at(idx) >= 64)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i160
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG34: ET >= 68 at BX = 0
        if (not (data->egIEt.at(idx) >= 68)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i182
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG24: ET >= 48 at BX = 0
        if (not (data->egIEt.at(idx) >= 48)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i184
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG26: ET >= 52 at BX = 0
        if (not (data->egIEt.at(idx) >= 52)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i185
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG30: ET >= 60 at BX = 0
        if (not (data->egIEt.at(idx) >= 60)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i78
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG23: ET >= 46 at BX = 0
        if (not (data->egIEt.at(idx) >= 46)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i79
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG20: ET >= 40 at BX = 0
        if (not (data->egIEt.at(idx) >= 40)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i80
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG10: ET >= 20 at BX = 0
        if (not (data->egIEt.at(idx) >= 20)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i81
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG20: ET >= 40 at BX = 0
        if (not (data->egIEt.at(idx) >= 40)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i82
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG23: ET >= 46 at BX = 0
        if (not (data->egIEt.at(idx) >= 46)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xc
        if (not ((0xc >> data->egIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleEG_i89
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG9: ET >= 18 at BX = 0
        if (not (data->egIEt.at(idx) >= 18)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      
bool
SingleETMHF_i100
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF40: ET >= 80 at BX = 0
      if (not (data->sumIEt.at(ii) >= 80)) continue;
      
  
    }



    pass = true;
    break;
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
        {
                    // ETMHF50: ET >= 100 at BX = 0
      if (not (data->sumIEt.at(ii) >= 100)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i118
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF60: ET >= 120 at BX = 0
      if (not (data->sumIEt.at(ii) >= 120)) continue;
      
  
    }



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
        {
                    // ETMHF70: ET >= 140 at BX = 0
      if (not (data->sumIEt.at(ii) >= 140)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i241
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF80: ET >= 160 at BX = 0
      if (not (data->sumIEt.at(ii) >= 160)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i242
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF90: ET >= 180 at BX = 0
      if (not (data->sumIEt.at(ii) >= 180)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i303
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF100: ET >= 200 at BX = 0
      if (not (data->sumIEt.at(ii) >= 200)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i304
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF110: ET >= 220 at BX = 0
      if (not (data->sumIEt.at(ii) >= 220)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i305
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF120: ET >= 240 at BX = 0
      if (not (data->sumIEt.at(ii) >= 240)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i306
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF130: ET >= 260 at BX = 0
      if (not (data->sumIEt.at(ii) >= 260)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i307
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF140: ET >= 280 at BX = 0
      if (not (data->sumIEt.at(ii) >= 280)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i308
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF150: ET >= 300 at BX = 0
      if (not (data->sumIEt.at(ii) >= 300)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETMHF_i420
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEtHF)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETMHF30: ET >= 60 at BX = 0
      if (not (data->sumIEt.at(ii) >= 60)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETM_i301
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETM120: ET >= 240 at BX = 0
      if (not (data->sumIEt.at(ii) >= 240)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETM_i302
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMissingEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETM150: ET >= 300 at BX = 0
      if (not (data->sumIEt.at(ii) >= 300)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETT_i298
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETT1200: ET >= 2400 at BX = 0
      if (not (data->sumIEt.at(ii) >= 2400)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETT_i299
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETT1600: ET >= 3200 at BX = 0
      if (not (data->sumIEt.at(ii) >= 3200)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleETT_i300
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalEt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // ETT2000: ET >= 4000 at BX = 0
      if (not (data->sumIEt.at(ii) >= 4000)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleEXT_i309
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i310
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i311
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i312
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i313
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i314
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i316
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i317
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i324
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i325
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i326
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i327
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i328
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i329
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i330
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i331
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i332
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i333
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i334
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i335
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i336
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i337
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i338
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i339
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i340
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i341
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i342
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i343
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i345
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i346
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i347
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i348
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i349
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleEXT_i350
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // Before the start of Run 3, the decision was always set FALSE;
  // This is set to be always TRUE in the emulator. This is important to do for Heavy Ion studies (as requested by the Heavy Ions team).
  //  bool pass = true;
  
  // However, this affects the rate quite significantly for pp collisions studies, and therefore has to be set to FALSE for pp studies.
  bool pass = false;
  
  return pass;
}

      
bool
SingleHTT_i102
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT240: ET >= 480 at BX = 0
      if (not (data->sumIEt.at(ii) >= 480)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i103
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT250: ET >= 500 at BX = 0
      if (not (data->sumIEt.at(ii) >= 500)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i115
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT60: ET >= 120 at BX = 0
      if (not (data->sumIEt.at(ii) >= 120)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i119
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT220: ET >= 440 at BX = 0
      if (not (data->sumIEt.at(ii) >= 440)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i120
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT260: ET >= 520 at BX = 0
      if (not (data->sumIEt.at(ii) >= 520)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i183
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT100: ET >= 200 at BX = 0
      if (not (data->sumIEt.at(ii) >= 200)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i187
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT280: ET >= 560 at BX = 0
      if (not (data->sumIEt.at(ii) >= 560)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i188
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT300: ET >= 600 at BX = 0
      if (not (data->sumIEt.at(ii) >= 600)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i189
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT320: ET >= 640 at BX = 0
      if (not (data->sumIEt.at(ii) >= 640)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i190
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT340: ET >= 680 at BX = 0
      if (not (data->sumIEt.at(ii) >= 680)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i291
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT120: ET >= 240 at BX = 0
      if (not (data->sumIEt.at(ii) >= 240)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i292
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT160: ET >= 320 at BX = 0
      if (not (data->sumIEt.at(ii) >= 320)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i293
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT200: ET >= 400 at BX = 0
      if (not (data->sumIEt.at(ii) >= 400)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i294
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT255: ET >= 510 at BX = 0
      if (not (data->sumIEt.at(ii) >= 510)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i295
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT360: ET >= 720 at BX = 0
      if (not (data->sumIEt.at(ii) >= 720)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i296
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT400: ET >= 800 at BX = 0
      if (not (data->sumIEt.at(ii) >= 800)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleHTT_i297
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kTotalHt)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // HTT450: ET >= 900 at BX = 0
      if (not (data->sumIEt.at(ii) >= 900)) continue;
      
  
    }



    pass = true;
    break;
  }

  return pass;
}

      


bool
SingleJET_i116
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i217
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(idx) >= 70)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i218
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i219
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET90: ET >= 180 at BX = 0
        if (not (data->jetIEt.at(idx) >= 180)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i220
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET120: ET >= 240 at BX = 0
        if (not (data->jetIEt.at(idx) >= 240)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i221
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET180: ET >= 360 at BX = 0
        if (not (data->jetIEt.at(idx) >= 360)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i222
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET200: ET >= 400 at BX = 0
        if (not (data->jetIEt.at(idx) >= 400)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i223
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(idx) >= 70)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i224
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET90: ET >= 180 at BX = 0
        if (not (data->jetIEt.at(idx) >= 180)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i225
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET120: ET >= 240 at BX = 0
        if (not (data->jetIEt.at(idx) >= 240)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i226
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET140: ET >= 280 at BX = 0
        if (not (data->jetIEt.at(idx) >= 280)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i227
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET160: ET >= 320 at BX = 0
        if (not (data->jetIEt.at(idx) >= 320)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i228
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET180: ET >= 360 at BX = 0
        if (not (data->jetIEt.at(idx) >= 360)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i229
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(idx) >= 70)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i230
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET35: ET >= 70 at BX = 0
        if (not (data->jetIEt.at(idx) >= 70)) continue;

        const auto eta = data->jetIEta.at(idx);
        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow0 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i231
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i232
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        const auto eta = data->jetIEta.at(idx);
        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow0 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i233
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET90: ET >= 180 at BX = 0
        if (not (data->jetIEt.at(idx) >= 180)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i234
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET90: ET >= 180 at BX = 0
        if (not (data->jetIEt.at(idx) >= 180)) continue;

        const auto eta = data->jetIEta.at(idx);
        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow0 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i235
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET120: ET >= 240 at BX = 0
        if (not (data->jetIEt.at(idx) >= 240)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET120: ET >= 240 at BX = 0
        if (not (data->jetIEt.at(idx) >= 240)) continue;

        const auto eta = data->jetIEta.at(idx);
        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow0 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET8: ET >= 16 at BX = 0
        if (not (data->jetIEt.at(idx) >= 16)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.9579999999999997 <= eta <= -1.392
        const bool etaWindow0 = ((-68 <= eta) and (eta <= -33));

        // 1.392 <= eta <= 2.9579999999999997
        const bool etaWindow1 = ((32 <= eta) and (eta <= 67));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET10: ET >= 20 at BX = 0
        if (not (data->jetIEt.at(idx) >= 20)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.9579999999999997 <= eta <= -1.392
        const bool etaWindow0 = ((-68 <= eta) and (eta <= -33));

        // 1.392 <= eta <= 2.9579999999999997
        const bool etaWindow1 = ((32 <= eta) and (eta <= 67));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i239
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET12: ET >= 24 at BX = 0
        if (not (data->jetIEt.at(idx) >= 24)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.9579999999999997 <= eta <= -1.392
        const bool etaWindow0 = ((-68 <= eta) and (eta <= -33));

        // 1.392 <= eta <= 2.9579999999999997
        const bool etaWindow1 = ((32 <= eta) and (eta <= 67));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i263
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET115: ET >= 230 at BX = 0
        if (not (data->jetIEt.at(idx) >= 230)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i286
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET20: ET >= 40 at BX = 0
        if (not (data->jetIEt.at(idx) >= 40)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -5.0 <= eta <= -3.0014999999999996
        const bool etaWindow0 = ((-115 <= eta) and (eta <= -70));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET20: ET >= 40 at BX = 0
        if (not (data->jetIEt.at(idx) >= 40)) continue;

        const auto eta = data->jetIEta.at(idx);
        // 3.0014999999999996 <= eta <= 5.0
        const bool etaWindow0 = ((69 <= eta) and (eta <= 114));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i319
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET20: ET >= 40 at BX = 0
        if (not (data->jetIEt.at(idx) >= 40)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i320
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET43: ET >= 86 at BX = 0
        if (not (data->jetIEt.at(idx) >= 86)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i321
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET46: ET >= 92 at BX = 0
        if (not (data->jetIEt.at(idx) >= 92)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i369
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET40: ET >= 80 at BX = 0
        if (not (data->jetIEt.at(idx) >= 80)) continue;

        // displaced jet bit : 0x1
        if (not ((0x1 & data->jetHwQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i370
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET60: ET >= 120 at BX = 0
        if (not (data->jetIEt.at(idx) >= 120)) continue;

        // displaced jet bit : 0x1
        if (not ((0x1 & data->jetHwQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i371
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET50: ET >= 100 at BX = 0
        if (not (data->jetIEt.at(idx) >= 100)) continue;

        // displaced jet bit : 0x1
        if (not ((0x1 & data->jetHwQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i409
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET55: ET >= 110 at BX = 0
        if (not (data->jetIEt.at(idx) >= 110)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i410
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET70: ET >= 140 at BX = 0
        if (not (data->jetIEt.at(idx) >= 140)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i412
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET70: ET >= 140 at BX = 0
        if (not (data->jetIEt.at(idx) >= 140)) continue;

        // displaced jet bit : 0x1
        if (not ((0x1 & data->jetHwQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i91
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET30: ET >= 60 at BX = 0
        if (not (data->jetIEt.at(idx) >= 60)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleJET_i99
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET100: ET >= 200 at BX = 0
        if (not (data->jetIEt.at(idx) >= 200)) continue;

        const auto eta = data->jetIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      
bool
SingleMBT0HFM_i323
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFM0)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // MBT0HFM1: Count >= 1 at BX = 0
      if (not (data->sumIEt.at(ii) >= 1)) continue;
          
  
    }



    pass = true;
    break;
  }

  return pass;
}

      
bool
SingleMBT0HFP_i322
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  bool pass = false;

  for (size_t ii = 0; ii < data->sumBx.size(); ii++)
  {
    if (not (data->sumType.at(ii) == L1Analysis::kMinBiasHFP0)) continue;
    if (not (data->sumBx.at(ii) == 0)) continue;
        {
                    // MBT0HFP1: Count >= 1 at BX = 0
      if (not (data->sumIEt.at(ii) >= 1)) continue;
          
  
    }



    pass = true;
    break;
  }

  return pass;
}

      


bool
SingleMU_i0
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i1
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i10
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU7: ET >= 15 at BX = 0
        if (not (data->muonIEt.at(idx) >= 15)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU7: ET >= 15 at BX = 0
        if (not (data->muonIEt.at(idx) >= 15)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU12: ET >= 25 at BX = 0
        if (not (data->muonIEt.at(idx) >= 25)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU12: ET >= 25 at BX = 0
        if (not (data->muonIEt.at(idx) >= 25)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // 0.7993125 <= eta <= 1.2451875
        const bool etaWindow0 = ((74 <= eta) and (eta <= 114));

        // -1.2451875 <= eta <= -0.7993125
        const bool etaWindow1 = ((-114 <= eta) and (eta <= -74));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU12: ET >= 25 at BX = 0
        if (not (data->muonIEt.at(idx) >= 25)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // 1.2451875 <= eta <= 2.45
        const bool etaWindow0 = ((115 <= eta) and (eta <= 225));

        // -2.45 <= eta <= -1.2451875
        const bool etaWindow1 = ((-225 <= eta) and (eta <= -115));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU15: ET >= 31 at BX = 0
        if (not (data->muonIEt.at(idx) >= 31)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU18: ET >= 37 at BX = 0
        if (not (data->muonIEt.at(idx) >= 37)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU20: ET >= 41 at BX = 0
        if (not (data->muonIEt.at(idx) >= 41)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU22: ET >= 45 at BX = 0
        if (not (data->muonIEt.at(idx) >= 45)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i2
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // 0.7993125 <= eta <= 1.2451875
        const bool etaWindow0 = ((74 <= eta) and (eta <= 114));

        // -1.2451875 <= eta <= -0.7993125
        const bool etaWindow1 = ((-114 <= eta) and (eta <= -74));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU22: ET >= 45 at BX = 0
        if (not (data->muonIEt.at(idx) >= 45)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i206
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU18: ET >= 37 at BX = 0
        if (not (data->muonIEt.at(idx) >= 37)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.1043125000000003 <= eta <= 2.1043125
        const bool etaWindow0 = ((-193 <= eta) and (eta <= 193));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i209
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU22: ET >= 45 at BX = 0
        if (not (data->muonIEt.at(idx) >= 45)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -2.1043125000000003 <= eta <= 2.1043125
        const bool etaWindow0 = ((-193 <= eta) and (eta <= 193));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU22: ET >= 45 at BX = 0
        if (not (data->muonIEt.at(idx) >= 45)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // 0.7993125 <= eta <= 1.2451875
        const bool etaWindow0 = ((74 <= eta) and (eta <= 114));

        // -1.2451875 <= eta <= -0.7993125
        const bool etaWindow1 = ((-114 <= eta) and (eta <= -74));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU22: ET >= 45 at BX = 0
        if (not (data->muonIEt.at(idx) >= 45)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // 1.2451875 <= eta <= 2.45
        const bool etaWindow0 = ((115 <= eta) and (eta <= 225));

        // -2.45 <= eta <= -1.2451875
        const bool etaWindow1 = ((-225 <= eta) and (eta <= -115));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU25: ET >= 51 at BX = 0
        if (not (data->muonIEt.at(idx) >= 51)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU6: ET >= 13 at BX = 0
        if (not (data->muonIEt.at(idx) >= 13)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU7: ET >= 15 at BX = 0
        if (not (data->muonIEt.at(idx) >= 15)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU8: ET >= 17 at BX = 0
        if (not (data->muonIEt.at(idx) >= 17)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU9: ET >= 19 at BX = 0
        if (not (data->muonIEt.at(idx) >= 19)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i278
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU8: ET >= 17 at BX = 0
        if (not (data->muonIEt.at(idx) >= 17)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i28
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU10: ET >= 21 at BX = 0
        if (not (data->muonIEt.at(idx) >= 21)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i29
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU12: ET >= 25 at BX = 0
        if (not (data->muonIEt.at(idx) >= 25)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i3
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // 1.2451875 <= eta <= 2.45
        const bool etaWindow0 = ((115 <= eta) and (eta <= 225));

        // -2.45 <= eta <= -1.2451875
        const bool etaWindow1 = ((-225 <= eta) and (eta <= -115));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i30
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU14: ET >= 29 at BX = 0
        if (not (data->muonIEt.at(idx) >= 29)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i31
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU16: ET >= 33 at BX = 0
        if (not (data->muonIEt.at(idx) >= 33)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i315
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.4083124999999999 <= eta <= 1.4083124999999999
        const bool etaWindow0 = ((-129 <= eta) and (eta <= 129));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i318
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.1038124999999999 <= eta <= 1.1038124999999999
        const bool etaWindow0 = ((-101 <= eta) and (eta <= 101));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i32
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU18: ET >= 37 at BX = 0
        if (not (data->muonIEt.at(idx) >= 37)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i4
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i414
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU22: ET >= 45 at BX = 0
        if (not (data->muonIEt.at(idx) >= 45)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i415
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU22: ET >= 45 at BX = 0
        if (not (data->muonIEt.at(idx) >= 45)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i5
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -0.7993125 <= eta <= 0.7993125
        const bool etaWindow0 = ((-73 <= eta) and (eta <= 73));

        if (not (etaWindow0)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // 0.7993125 <= eta <= 1.2451875
        const bool etaWindow0 = ((74 <= eta) and (eta <= 114));

        // -1.2451875 <= eta <= -0.7993125
        const bool etaWindow1 = ((-114 <= eta) and (eta <= -74));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // 1.2451875 <= eta <= 2.45
        const bool etaWindow0 = ((115 <= eta) and (eta <= 225));

        // -2.45 <= eta <= -1.2451875
        const bool etaWindow1 = ((-225 <= eta) and (eta <= -115));

        if (not (etaWindow0 or etaWindow1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i83
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU6: ET >= 13 at BX = 0
        if (not (data->muonIEt.at(idx) >= 13)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleMU_i98
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

        const auto eta = data->muonIEtaAtVtx.at(idx);
        // -1.5061874999999998 <= eta <= 1.5061875
        const bool etaWindow0 = ((-138 <= eta) and (eta <= 138));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i194
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU120: ET >= 240 at BX = 0
        if (not (data->tauIEt.at(idx) >= 240)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i195
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU130: ET >= 260 at BX = 0
        if (not (data->tauIEt.at(idx) >= 260)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i207
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU24: ET >= 48 at BX = 0
        if (not (data->tauIEt.at(idx) >= 48)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i208
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU26: ET >= 52 at BX = 0
        if (not (data->tauIEt.at(idx) >= 52)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i210
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU28: ET >= 56 at BX = 0
        if (not (data->tauIEt.at(idx) >= 56)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i211
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU30: ET >= 60 at BX = 0
        if (not (data->tauIEt.at(idx) >= 60)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i212
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU32: ET >= 64 at BX = 0
        if (not (data->tauIEt.at(idx) >= 64)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i213
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU34: ET >= 68 at BX = 0
        if (not (data->tauIEt.at(idx) >= 68)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i214
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU36: ET >= 72 at BX = 0
        if (not (data->tauIEt.at(idx) >= 72)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i215
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU40: ET >= 80 at BX = 0
        if (not (data->tauIEt.at(idx) >= 80)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i216
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU70: ET >= 140 at BX = 0
        if (not (data->tauIEt.at(idx) >= 140)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
SingleTAU_i387
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->tauBx.size(); ii++)
  {
    if (not (data->tauBx.at(ii) == 0)) continue;
    nobj++;
    if (nobj > 12) break;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 1) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 1);
  const auto& permutation = PermutationFactory::get(1);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // TAU52: ET >= 104 at BX = 0
        if (not (data->tauIEt.at(idx) >= 104)) continue;

        const auto eta = data->tauIEta.at(idx);
        // -2.1315 <= eta <= 2.1315
        const bool etaWindow0 = ((-49 <= eta) and (eta <= 48));

        // isolation : 0xe
        if (not ((0xe >> data->tauIso.at(idx)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

                    




bool
TransverseMass_i161
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
    bool pass = false;
  size_t nobj = 0;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;
        
              {
        // EG32: ET >= 64 at BX = 0
        if (not (data->egIEt.at(ii) >= 64)) continue;

        const auto eta = data->egIEta.at(ii);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        // isolation : 0xa
        if (not ((0xa >> data->egIso.at(ii)) & 1)) continue;

        if (not (etaWindow0)) continue;

      }


    for (size_t jj = 0; jj < data->sumBx.size(); jj++)
    {
      if (not (data->sumType.at(jj) == L1Analysis::kMissingEt)) continue;
      if (not (data->sumBx.at(jj) == 0)) continue;
          {
                    // ETM10: ET >= 20 at BX = 0
      if (not (data->sumIEt.at(jj) >= 20)) continue;
      
  
    }


              // 40.0 <= Mt <= 151982.0
      int iPhi = data->egIPhi.at(ii);
  
  unsigned int deltaIPhi = abs(iPhi - data->sumIPhi.at(jj));
  if (deltaIPhi >= 72) deltaIPhi = 2*72 - deltaIPhi;
  
    const long long cosDeltaPhi = LUT_COS_DPHI_EG_ETM[deltaIPhi];
    const long long pt0 = LUT_EG_ET[data->egIEt.at(ii)];
    const long long pt1 = LUT_ETM_ET[data->sumIEt.at(jj)];
    const long long mt2 = pt0 * pt1 * (1000 - cosDeltaPhi);
    {
      const long long minimum = 80000000; // 800.0 * 10^5
      const long long maximum = 1154926416200000; // 11549264162.0 * 10^5
      if (not ((minimum <= mt2) and (mt2 <= maximum))) continue;
    }


    
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}
    
      


bool
TripleEG_i174
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG16: ET >= 32 at BX = 0
        if (not (data->egIEt.at(idx) >= 32)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // EG8: ET >= 16 at BX = 0
        if (not (data->egIEt.at(idx) >= 16)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_i175
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG16: ET >= 32 at BX = 0
        if (not (data->egIEt.at(idx) >= 32)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG15: ET >= 30 at BX = 0
        if (not (data->egIEt.at(idx) >= 30)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // EG8: ET >= 16 at BX = 0
        if (not (data->egIEt.at(idx) >= 16)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_i176
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG18: ET >= 36 at BX = 0
        if (not (data->egIEt.at(idx) >= 36)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG17: ET >= 34 at BX = 0
        if (not (data->egIEt.at(idx) >= 34)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // EG8: ET >= 16 at BX = 0
        if (not (data->egIEt.at(idx) >= 16)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_i177
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG18: ET >= 36 at BX = 0
        if (not (data->egIEt.at(idx) >= 36)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG18: ET >= 36 at BX = 0
        if (not (data->egIEt.at(idx) >= 36)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // EG12: ET >= 24 at BX = 0
        if (not (data->egIEt.at(idx) >= 24)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleEG_i178
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->egBx.size(); ii++)
  {
    if (not (data->egBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // EG16: ET >= 32 at BX = 0
        if (not (data->egIEt.at(idx) >= 32)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // EG16: ET >= 32 at BX = 0
        if (not (data->egIEt.at(idx) >= 32)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // EG16: ET >= 32 at BX = 0
        if (not (data->egIEt.at(idx) >= 32)) continue;

        const auto eta = data->egIEta.at(idx);
        // -2.5229999999999997 <= eta <= 2.5229999999999997
        const bool etaWindow0 = ((-58 <= eta) and (eta <= 57));

        if (not (etaWindow0)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleJET_i279
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET95: ET >= 190 at BX = 0
        if (not (data->jetIEt.at(idx) >= 190)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET75: ET >= 150 at BX = 0
        if (not (data->jetIEt.at(idx) >= 150)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET65: ET >= 130 at BX = 0
        if (not (data->jetIEt.at(idx) >= 130)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleJET_i281
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET100: ET >= 200 at BX = 0
        if (not (data->jetIEt.at(idx) >= 200)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET80: ET >= 160 at BX = 0
        if (not (data->jetIEt.at(idx) >= 160)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET70: ET >= 140 at BX = 0
        if (not (data->jetIEt.at(idx) >= 140)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleJET_i283
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->jetBx.size(); ii++)
  {
    if (not (data->jetBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // JET105: ET >= 210 at BX = 0
        if (not (data->jetIEt.at(idx) >= 210)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // JET85: ET >= 170 at BX = 0
        if (not (data->jetIEt.at(idx) >= 170)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // JET75: ET >= 150 at BX = 0
        if (not (data->jetIEt.at(idx) >= 150)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i405
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU2: ET >= 5 at BX = 0
        if (not (data->muonIEt.at(idx) >= 5)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU1p5: ET >= 4 at BX = 0
        if (not (data->muonIEt.at(idx) >= 4)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i59
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i60
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i61
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i62
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i63
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i64
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i65
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU3p5: ET >= 8 at BX = 0
        if (not (data->muonIEt.at(idx) >= 8)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU2p5: ET >= 6 at BX = 0
        if (not (data->muonIEt.at(idx) >= 6)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i66
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i67
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i68
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i69
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU3p5: ET >= 8 at BX = 0
        if (not (data->muonIEt.at(idx) >= 8)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU2p5: ET >= 6 at BX = 0
        if (not (data->muonIEt.at(idx) >= 6)) continue;

        // quality : 0xfff0
        if (not ((0xfff0 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i72
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU4: ET >= 9 at BX = 0
        if (not (data->muonIEt.at(idx) >= 9)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU2p5: ET >= 6 at BX = 0
        if (not (data->muonIEt.at(idx) >= 6)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
      pass = true;
      break;
    }

    if (pass) break;
  }

  return pass;
}

      


bool
TripleMU_i74
(L1Analysis::L1AnalysisL1UpgradeDataFormat* data)
{
  // collect candidates of same BX offset
  size_t nobj = 0;
  std::vector<int> candidates;
  for (size_t ii = 0; ii < data->muonBx.size(); ii++)
  {
    if (not (data->muonBx.at(ii) == 0)) continue;
    nobj++;

    candidates.emplace_back(ii);
  }

  bool pass = false;

  // sufficient candidates found?
  if (candidates.size() < 3) return pass;

  const auto& combination = CombinationFactory::get(candidates.size(), 3);
  const auto& permutation = PermutationFactory::get(3);

  // match combinations
  for (const auto& set: combination)
  {
    for (const auto& indicies: permutation)
    {
      int idx = -1;
      idx = candidates.at(set.at(indicies.at(0)));
      {
        // MU5: ET >= 11 at BX = 0
        if (not (data->muonIEt.at(idx) >= 11)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(1)));
      {
        // MU3: ET >= 7 at BX = 0
        if (not (data->muonIEt.at(idx) >= 7)) continue;

        // quality : 0xf000
        if (not ((0xf000 >> data->muonQual.at(idx)) & 1)) continue;

      }

      idx = candidates.at(set.at(indicies.at(2)));
      {
        // MU0: ET >= 1 at BX = 0
        if (not (data->muonIEt.at(idx) >= 1)) continue;

        // quality : 0xff00
        if (not ((0xff00 >> data->muonQual.at(idx)) & 1)) continue;

      }

      
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
  return SingleEXT_i311(data) or ( not SingleEXT_i311(data));
}
bool
L1_BPTX_AND_Ref1_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i337(data);
}
bool
L1_BPTX_AND_Ref3_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i338(data);
}
bool
L1_BPTX_AND_Ref4_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i339(data);
}
bool
L1_BPTX_BeamGas_B1_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i342(data);
}
bool
L1_BPTX_BeamGas_B2_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i343(data);
}
bool
L1_BPTX_BeamGas_Ref1_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i340(data);
}
bool
L1_BPTX_BeamGas_Ref2_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i341(data);
}
bool
L1_BPTX_NotOR_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i333(data);
}
bool
L1_BPTX_OR_Ref3_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i334(data);
}
bool
L1_BPTX_OR_Ref4_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i335(data);
}
bool
L1_BPTX_RefAND_VME(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i336(data);
}
bool
L1_BptxMinus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i325(data);
}
bool
L1_BptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i314(data);
}
bool
L1_BptxPlus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i324(data);
}
bool
L1_BptxXOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return (SingleEXT_i324(data) and ( not SingleEXT_i325(data))) or (SingleEXT_i325(data) and ( not SingleEXT_i324(data)));
}
bool
L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i344(data);
}
bool
L1_DoubleEG10_er1p2_dR_Max0p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i400(data);
}
bool
L1_DoubleEG10p5_er1p2_dR_Max0p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i401(data);
}
bool
L1_DoubleEG11_er1p2_dR_Max0p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i402(data);
}
bool
L1_DoubleEG4_er1p2_dR_Max0p9(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i404(data);
}
bool
L1_DoubleEG4p5_er1p2_dR_Max0p9(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i403(data);
}
bool
L1_DoubleEG5_er1p2_dR_Max0p9(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i390(data);
}
bool
L1_DoubleEG5p5_er1p2_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i391(data);
}
bool
L1_DoubleEG6_er1p2_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i392(data);
}
bool
L1_DoubleEG6p5_er1p2_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i393(data);
}
bool
L1_DoubleEG7_er1p2_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i394(data);
}
bool
L1_DoubleEG7p5_er1p2_dR_Max0p7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i395(data);
}
bool
L1_DoubleEG8_er1p2_dR_Max0p7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i396(data);
}
bool
L1_DoubleEG8er2p5_HTT260er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i186(data) and SingleHTT_i120(data);
}
bool
L1_DoubleEG8er2p5_HTT280er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i186(data) and SingleHTT_i187(data);
}
bool
L1_DoubleEG8er2p5_HTT300er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i186(data) and SingleHTT_i188(data);
}
bool
L1_DoubleEG8er2p5_HTT320er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i186(data) and SingleHTT_i189(data);
}
bool
L1_DoubleEG8er2p5_HTT340er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i186(data) and SingleHTT_i190(data);
}
bool
L1_DoubleEG8p5_er1p2_dR_Max0p7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i397(data);
}
bool
L1_DoubleEG9_er1p2_dR_Max0p7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i398(data);
}
bool
L1_DoubleEG9p5_er1p2_dR_Max0p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i399(data);
}
bool
L1_DoubleEG_15_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i162(data);
}
bool
L1_DoubleEG_20_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i163(data);
}
bool
L1_DoubleEG_22_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i164(data);
}
bool
L1_DoubleEG_25_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i165(data);
}
bool
L1_DoubleEG_25_14_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i166(data);
}
bool
L1_DoubleEG_27_14_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i167(data);
}
bool
L1_DoubleEG_LooseIso16_LooseIso12_er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i354(data);
}
bool
L1_DoubleEG_LooseIso18_LooseIso12_er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i355(data);
}
bool
L1_DoubleEG_LooseIso20_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i168(data);
}
bool
L1_DoubleEG_LooseIso20_LooseIso12_er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i356(data);
}
bool
L1_DoubleEG_LooseIso22_10_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i169(data);
}
bool
L1_DoubleEG_LooseIso22_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i170(data);
}
bool
L1_DoubleEG_LooseIso22_LooseIso12_er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i357(data);
}
bool
L1_DoubleEG_LooseIso25_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i171(data);
}
bool
L1_DoubleEG_LooseIso25_LooseIso12_er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i358(data);
}
bool
L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTauOvRm_i381(data);
}
bool
L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTauOvRm_i408(data);
}
bool
L1_DoubleIsoTau28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i197(data);
}
bool
L1_DoubleIsoTau28er2p1_Mass_Max80(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i203(data);
}
bool
L1_DoubleIsoTau28er2p1_Mass_Max90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i202(data);
}
bool
L1_DoubleIsoTau30er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i198(data);
}
bool
L1_DoubleIsoTau30er2p1_Mass_Max80(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i205(data);
}
bool
L1_DoubleIsoTau30er2p1_Mass_Max90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i204(data);
}
bool
L1_DoubleIsoTau32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i199(data);
}
bool
L1_DoubleIsoTau34er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i200(data);
}
bool
L1_DoubleIsoTau35er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i367(data);
}
bool
L1_DoubleIsoTau36er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i201(data);
}
bool
L1_DoubleJet100er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i246(data);
}
bool
L1_DoubleJet100er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i243(data);
}
bool
L1_DoubleJet112er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i247(data);
}
bool
L1_DoubleJet120er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i244(data);
}
bool
L1_DoubleJet150er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i245(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i248(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i249(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i250(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i251(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i252(data);
}
bool
L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i253(data);
}
bool
L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMassOvRm_i275(data);
}
bool
L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMassOvRm_i368(data);
}
bool
L1_DoubleJet40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i117(data);
}
bool
L1_DoubleJet_100_30_DoubleJet30_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i256(data) and InvariantMass_i255(data);
}
bool
L1_DoubleJet_110_35_DoubleJet35_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i257(data) and InvariantMass_i258(data);
}
bool
L1_DoubleJet_115_40_DoubleJet40_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i259(data) and InvariantMass_i260(data);
}
bool
L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i263(data) and (InvariantMass_i264(data) or InvariantMass_i265(data) or InvariantMass_i266(data) or InvariantMass_i267(data) or InvariantMass_i268(data) or InvariantMass_i269(data));
}
bool
L1_DoubleJet_120_45_DoubleJet45_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i261(data) and InvariantMass_i262(data);
}
bool
L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i220(data) and (InvariantMass_i270(data) or InvariantMass_i265(data) or InvariantMass_i271(data) or InvariantMass_i272(data) or InvariantMass_i273(data) or InvariantMass_i274(data));
}
bool
L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i277(data) and DoubleMU_i35(data);
}
bool
L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMassOvRm_i276(data);
}
bool
L1_DoubleJet_80_30_Mass_Min420_Mu8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i277(data) and SingleMU_i278(data);
}
bool
L1_DoubleJet_90_30_DoubleJet30_Mass_Min620(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i254(data) and InvariantMass_i255(data);
}
bool
L1_DoubleLLPJet40(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleJET_i413(data);
}
bool
L1_DoubleLooseIsoEG22er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i172(data);
}
bool
L1_DoubleLooseIsoEG24er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleEG_i173(data);
}
bool
L1_DoubleMu0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i34(data);
}
bool
L1_DoubleMu0_Mass_Min1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i37(data);
}
bool
L1_DoubleMu0_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i33(data);
}
bool
L1_DoubleMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i35(data);
}
bool
L1_DoubleMu0_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i36(data);
}
bool
L1_DoubleMu0_Upt15_Upt7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i359(data) or DoubleMU_i360(data) or DoubleMU_i361(data) or DoubleMU_i362(data);
}
bool
L1_DoubleMu0_Upt5_Upt5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i416(data) or DoubleMU_i417(data) or DoubleMU_i418(data) or DoubleMU_i419(data);
}
bool
L1_DoubleMu0_Upt6_IP_Min1_Upt4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i363(data) or DoubleMU_i364(data) or DoubleMU_i365(data) or DoubleMU_i366(data);
}
bool
L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i110(data) and CaloMuonCorrelation_i111(data);
}
bool
L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i372(data);
}
bool
L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i52(data);
}
bool
L1_DoubleMu0er1p5_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i48(data);
}
bool
L1_DoubleMu0er1p5_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i49(data);
}
bool
L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i51(data);
}
bool
L1_DoubleMu0er1p5_SQ_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i50(data);
}
bool
L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i351(data);
}
bool
L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i352(data);
}
bool
L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i47(data);
}
bool
L1_DoubleMu0er2p0_SQ_dEta_Max1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i407(data);
}
bool
L1_DoubleMu0er2p0_SQ_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i406(data);
}
bool
L1_DoubleMu0er2p0_SQ_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i46(data);
}
bool
L1_DoubleMu18er2p1_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i45(data);
}
bool
L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i377(data) and InvariantMass_i378(data);
}
bool
L1_DoubleMu3_SQ_ETMHF30_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleETMHF_i420(data) and SingleHTT_i115(data);
}
bool
L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleETMHF_i420(data) and (SingleJET_i116(data) or DoubleJET_i117(data));
}
bool
L1_DoubleMu3_SQ_ETMHF40_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleETMHF_i100(data) and SingleHTT_i115(data);
}
bool
L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleETMHF_i100(data) and (SingleJET_i116(data) or DoubleJET_i117(data));
}
bool
L1_DoubleMu3_SQ_ETMHF50_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleETMHF_i101(data) and SingleHTT_i115(data);
}
bool
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleETMHF_i101(data) and SingleJET_i116(data);
}
bool
L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleETMHF_i101(data) and (SingleJET_i116(data) or DoubleJET_i117(data));
}
bool
L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleETMHF_i118(data) and SingleJET_i116(data);
}
bool
L1_DoubleMu3_SQ_HTT220er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleHTT_i119(data);
}
bool
L1_DoubleMu3_SQ_HTT240er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleHTT_i102(data);
}
bool
L1_DoubleMu3_SQ_HTT260er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i114(data) and SingleHTT_i120(data);
}
bool
L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i112(data) and CaloMuonCorrelation_i113(data);
}
bool
L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i353(data);
}
bool
L1_DoubleMu4_SQ_EG9er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i88(data) and SingleEG_i89(data);
}
bool
L1_DoubleMu4_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i53(data);
}
bool
L1_DoubleMu4_SQ_OS_dR_Max1p2(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i54(data);
}
bool
L1_DoubleMu4p5_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i55(data);
}
bool
L1_DoubleMu4p5_SQ_OS_dR_Max1p2(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonMuonCorrelation_i56(data);
}
bool
L1_DoubleMu4p5er2p0_SQ_OS(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i57(data);
}
bool
L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i373(data);
}
bool
L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i58(data);
}
bool
L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i379(data) and InvariantMass_i380(data);
}
bool
L1_DoubleMu5_SQ_EG9er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i90(data) and SingleEG_i89(data);
}
bool
L1_DoubleMu8_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i38(data);
}
bool
L1_DoubleMu9_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i39(data);
}
bool
L1_DoubleMu_12_5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i40(data);
}
bool
L1_DoubleMu_15_5_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i41(data);
}
bool
L1_DoubleMu_15_7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i42(data);
}
bool
L1_DoubleMu_15_7_Mass_Min1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass_i44(data);
}
bool
L1_DoubleMu_15_7_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleMU_i43(data);
}
bool
L1_DoubleTau70er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return DoubleTAU_i196(data);
}
bool
L1_ETM120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETM_i301(data);
}
bool
L1_ETM150(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETM_i302(data);
}
bool
L1_ETMHF100(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i303(data);
}
bool
L1_ETMHF100_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i303(data) and SingleHTT_i115(data);
}
bool
L1_ETMHF110(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i304(data);
}
bool
L1_ETMHF110_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i304(data) and SingleHTT_i115(data);
}
bool
L1_ETMHF110_HTT60er_NotSecondBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i304(data) and SingleHTT_i115(data) and ((SingleEXT_i309(data)) or ( not SingleEXT_i310(data)) or ( not SingleEXT_i311(data)) or ( not SingleEXT_i312(data)) or ( not SingleEXT_i313(data)));
}
bool
L1_ETMHF120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i305(data);
}
bool
L1_ETMHF120_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i305(data) and SingleHTT_i115(data);
}
bool
L1_ETMHF120_NotSecondBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i305(data) and ((SingleEXT_i309(data)) or ( not SingleEXT_i310(data)) or ( not SingleEXT_i311(data)) or ( not SingleEXT_i312(data)) or ( not SingleEXT_i313(data)));
}
bool
L1_ETMHF130(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i306(data);
}
bool
L1_ETMHF130_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i306(data) and SingleHTT_i115(data);
}
bool
L1_ETMHF140(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i307(data);
}
bool
L1_ETMHF150(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i308(data);
}
bool
L1_ETMHF70(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i240(data);
}
bool
L1_ETMHF70_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i240(data) and SingleHTT_i115(data);
}
bool
L1_ETMHF80(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i241(data);
}
bool
L1_ETMHF80_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i241(data) and SingleHTT_i115(data);
}
bool
L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i241(data) and CaloEsumCorrelation_i421(data);
}
bool
L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i241(data) and CaloEsumCorrelation_i422(data);
}
bool
L1_ETMHF90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i242(data);
}
bool
L1_ETMHF90_HTT60er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i242(data) and SingleHTT_i115(data);
}
bool
L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i242(data) and CaloEsumCorrelation_i382(data);
}
bool
L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i242(data) and CaloEsumCorrelation_i383(data);
}
bool
L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i242(data) and CaloEsumCorrelation_i384(data);
}
bool
L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETMHF_i242(data) and CaloEsumCorrelation_i385(data);
}
bool
L1_ETT1200(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETT_i298(data);
}
bool
L1_ETT1600(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETT_i299(data);
}
bool
L1_ETT2000(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleETT_i300(data);
}
bool
L1_FirstBunchAfterTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i329(data) and SingleEXT_i310(data) and ( not SingleEXT_i314(data)) and ( not SingleEXT_i317(data)) and ( not SingleEXT_i328(data));
}
bool
L1_FirstBunchBeforeTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_i309(data)) and ( not SingleEXT_i316(data)) and ( not SingleEXT_i314(data)) and SingleEXT_i312(data) and SingleEXT_i313(data);
}
bool
L1_FirstBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_i309(data)) and ( not SingleEXT_i316(data)) and SingleEXT_i311(data) and SingleEXT_i312(data) and SingleEXT_i313(data);
}
bool
L1_FirstCollisionInOrbit(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i332(data);
}
bool
L1_FirstCollisionInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i331(data);
}
bool
L1_HCAL_LaserMon_Trig(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i345(data);
}
bool
L1_HCAL_LaserMon_Veto(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i346(data);
}
bool
L1_HTT120_SingleLLPJet40(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i369(data) and SingleHTT_i291(data);
}
bool
L1_HTT120er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i291(data);
}
bool
L1_HTT160_SingleLLPJet50(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i371(data) and SingleHTT_i292(data);
}
bool
L1_HTT160er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i292(data);
}
bool
L1_HTT200_SingleLLPJet60(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i370(data) and SingleHTT_i293(data);
}
bool
L1_HTT200er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i293(data);
}
bool
L1_HTT240_SingleLLPJet70(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i412(data) and SingleHTT_i102(data);
}
bool
L1_HTT255er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i294(data);
}
bool
L1_HTT280er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i187(data);
}
bool
L1_HTT280er_QuadJet_70_55_40_35_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i187(data) and QuadJET_i388(data);
}
bool
L1_HTT320er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i189(data);
}
bool
L1_HTT320er_QuadJet_70_55_40_40_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i189(data) and QuadJET_i389(data);
}
bool
L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i189(data) and QuadJET_i289(data);
}
bool
L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i189(data) and QuadJET_i290(data);
}
bool
L1_HTT360er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i295(data);
}
bool
L1_HTT400er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i296(data);
}
bool
L1_HTT450er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleHTT_i297(data);
}
bool
L1_IsoEG32er2p5_Mt40(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TransverseMass_i161(data);
}
bool
L1_IsoTau52er2p1_QuadJet36er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadJET_i386(data) and SingleTAU_i387(data);
}
bool
L1_IsolatedBunch(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_i309(data)) and ( not SingleEXT_i316(data)) and SingleEXT_i311(data) and ( not SingleEXT_i317(data)) and ( not SingleEXT_i328(data));
}
bool
L1_LastBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i329(data) and SingleEXT_i310(data) and SingleEXT_i311(data) and ( not SingleEXT_i317(data)) and ( not SingleEXT_i328(data));
}
bool
L1_LastCollisionInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i330(data);
}
bool
L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i191(data);
}
bool
L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i193(data);
}
bool
L1_LooseIsoEG24er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i182(data) and SingleHTT_i183(data);
}
bool
L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i192(data);
}
bool
L1_LooseIsoEG26er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i184(data) and SingleHTT_i183(data);
}
bool
L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i179(data);
}
bool
L1_LooseIsoEG28er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i142(data) and SingleHTT_i183(data);
}
bool
L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i180(data);
}
bool
L1_LooseIsoEG30er2p1_HTT100er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i185(data) and SingleHTT_i183(data);
}
bool
L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloCaloCorrelation_i181(data);
}
bool
L1_MinimumBiasHF0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMBT0HFP_i322(data) and SingleMBT0HFM_i323(data);
}
bool
L1_MinimumBiasHF0_AND_BptxAND(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return (SingleMBT0HFP_i322(data) and SingleMBT0HFM_i323(data)) and SingleEXT_i311(data);
}
bool
L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i104(data) and CaloCaloCorrelation_i105(data);
}
bool
L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i108(data) and CaloCaloCorrelation_i109(data);
}
bool
L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i106(data) and CaloCaloCorrelation_i107(data);
}
bool
L1_Mu18er2p1_Tau24er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i206(data) and SingleTAU_i207(data);
}
bool
L1_Mu18er2p1_Tau26er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i206(data) and SingleTAU_i208(data);
}
bool
L1_Mu18er2p1_Tau26er2p1_Jet55(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i206(data) and SingleTAU_i208(data) and SingleJET_i409(data);
}
bool
L1_Mu18er2p1_Tau26er2p1_Jet70(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i206(data) and SingleTAU_i208(data) and SingleJET_i410(data);
}
bool
L1_Mu20_EG10er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i18(data) and SingleEG_i80(data);
}
bool
L1_Mu22er2p1_IsoTau28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i209(data) and SingleTAU_i210(data);
}
bool
L1_Mu22er2p1_IsoTau30er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i209(data) and SingleTAU_i211(data);
}
bool
L1_Mu22er2p1_IsoTau32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i209(data) and SingleTAU_i212(data);
}
bool
L1_Mu22er2p1_IsoTau34er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i209(data) and SingleTAU_i213(data);
}
bool
L1_Mu22er2p1_IsoTau36er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i209(data) and SingleTAU_i214(data);
}
bool
L1_Mu22er2p1_IsoTau40er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i209(data) and SingleTAU_i215(data);
}
bool
L1_Mu22er2p1_Tau70er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i209(data) and SingleTAU_i216(data);
}
bool
L1_Mu3_Jet120er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i97(data);
}
bool
L1_Mu3_Jet120er2p5_dR_Max0p8(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i96(data);
}
bool
L1_Mu3_Jet16er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i92(data);
}
bool
L1_Mu3_Jet30er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i9(data) and SingleJET_i91(data);
}
bool
L1_Mu3_Jet35er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i93(data);
}
bool
L1_Mu3_Jet60er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i94(data);
}
bool
L1_Mu3_Jet80er2p5_dR_Max0p4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return CaloMuonCorrelation_i95(data);
}
bool
L1_Mu3er1p5_Jet100er2p5_ETMHF30(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i98(data) and SingleJET_i99(data) and SingleETMHF_i420(data);
}
bool
L1_Mu3er1p5_Jet100er2p5_ETMHF40(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i98(data) and SingleJET_i99(data) and SingleETMHF_i100(data);
}
bool
L1_Mu3er1p5_Jet100er2p5_ETMHF50(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i98(data) and SingleJET_i99(data) and SingleETMHF_i101(data);
}
bool
L1_Mu5_EG23er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i10(data) and SingleEG_i78(data);
}
bool
L1_Mu5_LooseIsoEG20er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i10(data) and SingleEG_i81(data);
}
bool
L1_Mu6_DoubleEG10er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i83(data) and DoubleEG_i84(data);
}
bool
L1_Mu6_DoubleEG12er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i83(data) and DoubleEG_i85(data);
}
bool
L1_Mu6_DoubleEG15er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i83(data) and DoubleEG_i86(data);
}
bool
L1_Mu6_DoubleEG17er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i83(data) and DoubleEG_i87(data);
}
bool
L1_Mu6_HTT240er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i83(data) and SingleHTT_i102(data);
}
bool
L1_Mu6_HTT250er(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i83(data) and SingleHTT_i103(data);
}
bool
L1_Mu7_EG20er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i12(data) and SingleEG_i79(data);
}
bool
L1_Mu7_EG23er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i12(data) and SingleEG_i78(data);
}
bool
L1_Mu7_LooseIsoEG20er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i12(data) and SingleEG_i81(data);
}
bool
L1_Mu7_LooseIsoEG23er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i12(data) and SingleEG_i82(data);
}
bool
L1_NotBptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return not SingleEXT_i314(data);
}
bool
L1_QuadJet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadJET_i288(data);
}
bool
L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadJET_i285(data) and DoubleJET_i280(data) and (SingleJET_i286(data) or SingleJET_i287(data));
}
bool
L1_QuadMu0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadMU_i76(data);
}
bool
L1_QuadMu0_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadMU_i75(data);
}
bool
L1_QuadMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return QuadMU_i77(data);
}
bool
L1_SecondBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return ( not SingleEXT_i309(data)) and SingleEXT_i310(data) and SingleEXT_i311(data) and SingleEXT_i312(data) and SingleEXT_i313(data);
}
bool
L1_SecondLastBunchInTrain(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i329(data) and SingleEXT_i310(data) and SingleEXT_i311(data) and SingleEXT_i312(data) and ( not SingleEXT_i328(data));
}
bool
L1_SingleEG10er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i80(data);
}
bool
L1_SingleEG15er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i122(data);
}
bool
L1_SingleEG26er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i123(data);
}
bool
L1_SingleEG28_FWD2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i124(data) or SingleEG_i125(data);
}
bool
L1_SingleEG28er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i128(data);
}
bool
L1_SingleEG28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i127(data);
}
bool
L1_SingleEG28er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i126(data);
}
bool
L1_SingleEG34er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i129(data);
}
bool
L1_SingleEG36er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i130(data);
}
bool
L1_SingleEG38er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i131(data);
}
bool
L1_SingleEG40er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i132(data);
}
bool
L1_SingleEG42er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i133(data);
}
bool
L1_SingleEG45er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i134(data);
}
bool
L1_SingleEG50(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i135(data);
}
bool
L1_SingleEG60(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i136(data);
}
bool
L1_SingleEG8er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i121(data);
}
bool
L1_SingleIsoEG24er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i147(data);
}
bool
L1_SingleIsoEG24er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i146(data);
}
bool
L1_SingleIsoEG26er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i150(data);
}
bool
L1_SingleIsoEG26er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i149(data);
}
bool
L1_SingleIsoEG26er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i148(data);
}
bool
L1_SingleIsoEG28_FWD2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i151(data) or SingleEG_i152(data);
}
bool
L1_SingleIsoEG28er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i155(data);
}
bool
L1_SingleIsoEG28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i154(data);
}
bool
L1_SingleIsoEG28er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i153(data);
}
bool
L1_SingleIsoEG30er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i157(data);
}
bool
L1_SingleIsoEG30er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i156(data);
}
bool
L1_SingleIsoEG32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i159(data);
}
bool
L1_SingleIsoEG32er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i158(data);
}
bool
L1_SingleIsoEG34er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i160(data);
}
bool
L1_SingleIsoTau32er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i212(data);
}
bool
L1_SingleJet10erHE(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i238(data);
}
bool
L1_SingleJet120(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i220(data);
}
bool
L1_SingleJet120_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i235(data) or SingleJET_i236(data);
}
bool
L1_SingleJet120er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i225(data);
}
bool
L1_SingleJet12erHE(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i239(data);
}
bool
L1_SingleJet140er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i226(data);
}
bool
L1_SingleJet140er2p5_ETMHF70(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i226(data) and SingleETMHF_i240(data);
}
bool
L1_SingleJet140er2p5_ETMHF80(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i226(data) and SingleETMHF_i241(data);
}
bool
L1_SingleJet140er2p5_ETMHF90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i226(data) and SingleETMHF_i242(data);
}
bool
L1_SingleJet160er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i227(data);
}
bool
L1_SingleJet180(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i221(data);
}
bool
L1_SingleJet180er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i228(data);
}
bool
L1_SingleJet200(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i222(data);
}
bool
L1_SingleJet20er2p5_NotBptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i319(data) and ( not SingleEXT_i314(data));
}
bool
L1_SingleJet20er2p5_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i319(data) and ( not SingleEXT_i316(data)) and ( not SingleEXT_i314(data)) and ( not SingleEXT_i317(data));
}
bool
L1_SingleJet35(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i217(data);
}
bool
L1_SingleJet35_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i229(data) or SingleJET_i230(data);
}
bool
L1_SingleJet35er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i223(data);
}
bool
L1_SingleJet43er2p5_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i320(data) and ( not SingleEXT_i316(data)) and ( not SingleEXT_i314(data)) and ( not SingleEXT_i317(data));
}
bool
L1_SingleJet46er2p5_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i321(data) and ( not SingleEXT_i316(data)) and ( not SingleEXT_i314(data)) and ( not SingleEXT_i317(data));
}
bool
L1_SingleJet60(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i218(data);
}
bool
L1_SingleJet60_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i231(data) or SingleJET_i232(data);
}
bool
L1_SingleJet60er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i116(data);
}
bool
L1_SingleJet8erHE(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i237(data);
}
bool
L1_SingleJet90(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i219(data);
}
bool
L1_SingleJet90_FWD3p0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i233(data) or SingleJET_i234(data);
}
bool
L1_SingleJet90er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleJET_i224(data);
}
bool
L1_SingleLooseIsoEG26er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i138(data);
}
bool
L1_SingleLooseIsoEG26er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i137(data);
}
bool
L1_SingleLooseIsoEG28_FWD2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i139(data) or SingleEG_i140(data);
}
bool
L1_SingleLooseIsoEG28er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i143(data);
}
bool
L1_SingleLooseIsoEG28er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i142(data);
}
bool
L1_SingleLooseIsoEG28er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i141(data);
}
bool
L1_SingleLooseIsoEG30er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i145(data);
}
bool
L1_SingleLooseIsoEG30er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEG_i144(data);
}
bool
L1_SingleMu0_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i6(data);
}
bool
L1_SingleMu0_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i5(data);
}
bool
L1_SingleMu0_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i8(data);
}
bool
L1_SingleMu0_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i7(data);
}
bool
L1_SingleMu10er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i28(data);
}
bool
L1_SingleMu12_DQ_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i13(data);
}
bool
L1_SingleMu12_DQ_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i15(data);
}
bool
L1_SingleMu12_DQ_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i14(data);
}
bool
L1_SingleMu12er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i29(data);
}
bool
L1_SingleMu14er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i30(data);
}
bool
L1_SingleMu15_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i16(data);
}
bool
L1_SingleMu16er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i31(data);
}
bool
L1_SingleMu18(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i17(data);
}
bool
L1_SingleMu18er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i32(data);
}
bool
L1_SingleMu20(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i18(data);
}
bool
L1_SingleMu22(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i19(data);
}
bool
L1_SingleMu22_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i20(data);
}
bool
L1_SingleMu22_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i415(data);
}
bool
L1_SingleMu22_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i22(data);
}
bool
L1_SingleMu22_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i21(data);
}
bool
L1_SingleMu22_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i414(data);
}
bool
L1_SingleMu25(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i23(data);
}
bool
L1_SingleMu3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i9(data);
}
bool
L1_SingleMu5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i10(data);
}
bool
L1_SingleMu6er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i24(data);
}
bool
L1_SingleMu7(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i12(data);
}
bool
L1_SingleMu7_DQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i11(data);
}
bool
L1_SingleMu7er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i25(data);
}
bool
L1_SingleMu8er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i26(data);
}
bool
L1_SingleMu9er1p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i27(data);
}
bool
L1_SingleMuCosmics(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i0(data);
}
bool
L1_SingleMuCosmics_BMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i1(data);
}
bool
L1_SingleMuCosmics_EMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i3(data);
}
bool
L1_SingleMuCosmics_OMTF(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i2(data);
}
bool
L1_SingleMuOpen(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i4(data);
}
bool
L1_SingleMuOpen_NotBptxOR(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i4(data) and ( not SingleEXT_i314(data));
}
bool
L1_SingleMuOpen_er1p1_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i318(data) and ( not SingleEXT_i316(data)) and ( not SingleEXT_i314(data)) and ( not SingleEXT_i317(data));
}
bool
L1_SingleMuOpen_er1p4_NotBptxOR_3BX(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleMU_i315(data) and ( not SingleEXT_i316(data)) and ( not SingleEXT_i314(data)) and ( not SingleEXT_i317(data));
}
bool
L1_SingleMuShower_Nominal(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonShower0_i375(data);
}
bool
L1_SingleMuShower_Tight(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return MuonShower1_i376(data);
}
bool
L1_SingleTau120er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i194(data);
}
bool
L1_SingleTau130er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i195(data);
}
bool
L1_SingleTau70er2p1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleTAU_i216(data);
}
bool
L1_TOTEM_1(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i347(data);
}
bool
L1_TOTEM_2(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i348(data);
}
bool
L1_TOTEM_3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i349(data);
}
bool
L1_TOTEM_4(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i350(data);
}
bool
L1_TripleEG16er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i178(data);
}
bool
L1_TripleEG_16_12_8_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i174(data);
}
bool
L1_TripleEG_16_15_8_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i175(data);
}
bool
L1_TripleEG_18_17_8_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i176(data);
}
bool
L1_TripleEG_18_18_12_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleEG_i177(data);
}
bool
L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleJET_i281(data) and DoubleJET_i282(data);
}
bool
L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleJET_i283(data) and DoubleJET_i284(data);
}
bool
L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleJET_i279(data) and DoubleJET_i280(data);
}
bool
L1_TripleMu0(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i60(data);
}
bool
L1_TripleMu0_OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i59(data);
}
bool
L1_TripleMu0_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i61(data);
}
bool
L1_TripleMu3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i62(data);
}
bool
L1_TripleMu3_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i63(data);
}
bool
L1_TripleMu_2SQ_1p5SQ_0OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i405(data);
}
bool
L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass3_i374(data);
}
bool
L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return InvariantMass3_i411(data);
}
bool
L1_TripleMu_5SQ_3SQ_0OQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i64(data);
}
bool
L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i64(data) and InvariantMass_i73(data);
}
bool
L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i74(data) and InvariantMass_i73(data);
}
bool
L1_TripleMu_5_3_3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i66(data);
}
bool
L1_TripleMu_5_3_3_SQ(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i67(data);
}
bool
L1_TripleMu_5_3p5_2p5(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i65(data);
}
bool
L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i65(data) and InvariantMass_i71(data);
}
bool
L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i69(data) and InvariantMass_i70(data);
}
bool
L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i72(data) and InvariantMass_i71(data);
}
bool
L1_TripleMu_5_5_3(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return TripleMU_i68(data);
}
bool
L1_UnpairedBunchBptxMinus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i327(data);
}
bool
L1_UnpairedBunchBptxPlus(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i326(data);
}
bool
L1_ZeroBias(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i311(data);
}
bool
L1_ZeroBias_copy(L1Analysis::L1AnalysisL1UpgradeDataFormat* data, L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  return SingleEXT_i311(data);
}


std::string getNameFromId(const size_t index)
{
  static const std::map<size_t, std::string> id2name = {
    {458, "L1_AlwaysTrue"},
    {486, "L1_BPTX_AND_Ref1_VME"},
    {487, "L1_BPTX_AND_Ref3_VME"},
    {488, "L1_BPTX_AND_Ref4_VME"},
    {491, "L1_BPTX_BeamGas_B1_VME"},
    {492, "L1_BPTX_BeamGas_B2_VME"},
    {489, "L1_BPTX_BeamGas_Ref1_VME"},
    {490, "L1_BPTX_BeamGas_Ref2_VME"},
    {482, "L1_BPTX_NotOR_VME"},
    {483, "L1_BPTX_OR_Ref3_VME"},
    {484, "L1_BPTX_OR_Ref4_VME"},
    {485, "L1_BPTX_RefAND_VME"},
    {467, "L1_BptxMinus"},
    {464, "L1_BptxOR"},
    {466, "L1_BptxPlus"},
    {465, "L1_BptxXOR"},
    {494, "L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142"},
    {212, "L1_DoubleEG10_er1p2_dR_Max0p6"},
    {213, "L1_DoubleEG10p5_er1p2_dR_Max0p6"},
    {214, "L1_DoubleEG11_er1p2_dR_Max0p6"},
    {200, "L1_DoubleEG4_er1p2_dR_Max0p9"},
    {201, "L1_DoubleEG4p5_er1p2_dR_Max0p9"},
    {202, "L1_DoubleEG5_er1p2_dR_Max0p9"},
    {203, "L1_DoubleEG5p5_er1p2_dR_Max0p8"},
    {204, "L1_DoubleEG6_er1p2_dR_Max0p8"},
    {205, "L1_DoubleEG6p5_er1p2_dR_Max0p8"},
    {206, "L1_DoubleEG7_er1p2_dR_Max0p8"},
    {207, "L1_DoubleEG7p5_er1p2_dR_Max0p7"},
    {208, "L1_DoubleEG8_er1p2_dR_Max0p7"},
    {247, "L1_DoubleEG8er2p5_HTT260er"},
    {248, "L1_DoubleEG8er2p5_HTT280er"},
    {249, "L1_DoubleEG8er2p5_HTT300er"},
    {250, "L1_DoubleEG8er2p5_HTT320er"},
    {251, "L1_DoubleEG8er2p5_HTT340er"},
    {209, "L1_DoubleEG8p5_er1p2_dR_Max0p7"},
    {210, "L1_DoubleEG9_er1p2_dR_Max0p7"},
    {211, "L1_DoubleEG9p5_er1p2_dR_Max0p6"},
    {215, "L1_DoubleEG_15_10_er2p5"},
    {216, "L1_DoubleEG_20_10_er2p5"},
    {217, "L1_DoubleEG_22_10_er2p5"},
    {218, "L1_DoubleEG_25_12_er2p5"},
    {219, "L1_DoubleEG_25_14_er2p5"},
    {220, "L1_DoubleEG_27_14_er2p5"},
    {225, "L1_DoubleEG_LooseIso16_LooseIso12_er1p5"},
    {226, "L1_DoubleEG_LooseIso18_LooseIso12_er1p5"},
    {221, "L1_DoubleEG_LooseIso20_10_er2p5"},
    {227, "L1_DoubleEG_LooseIso20_LooseIso12_er1p5"},
    {222, "L1_DoubleEG_LooseIso22_10_er2p5"},
    {223, "L1_DoubleEG_LooseIso22_12_er2p5"},
    {228, "L1_DoubleEG_LooseIso22_LooseIso12_er1p5"},
    {224, "L1_DoubleEG_LooseIso25_12_er2p5"},
    {229, "L1_DoubleEG_LooseIso25_LooseIso12_er1p5"},
    {283, "L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5"},
    {284, "L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5"},
    {268, "L1_DoubleIsoTau28er2p1"},
    {275, "L1_DoubleIsoTau28er2p1_Mass_Max80"},
    {274, "L1_DoubleIsoTau28er2p1_Mass_Max90"},
    {269, "L1_DoubleIsoTau30er2p1"},
    {277, "L1_DoubleIsoTau30er2p1_Mass_Max80"},
    {276, "L1_DoubleIsoTau30er2p1_Mass_Max90"},
    {270, "L1_DoubleIsoTau32er2p1"},
    {271, "L1_DoubleIsoTau34er2p1"},
    {272, "L1_DoubleIsoTau35er2p1"},
    {273, "L1_DoubleIsoTau36er2p1"},
    {345, "L1_DoubleJet100er2p3_dEta_Max1p6"},
    {341, "L1_DoubleJet100er2p5"},
    {346, "L1_DoubleJet112er2p3_dEta_Max1p6"},
    {342, "L1_DoubleJet120er2p5"},
    {343, "L1_DoubleJet150er2p5"},
    {348, "L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5"},
    {349, "L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5"},
    {350, "L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5"},
    {351, "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5"},
    {352, "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5"},
    {353, "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5"},
    {363, "L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp"},
    {362, "L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5"},
    {340, "L1_DoubleJet40er2p5"},
    {356, "L1_DoubleJet_100_30_DoubleJet30_Mass_Min620"},
    {357, "L1_DoubleJet_110_35_DoubleJet35_Mass_Min620"},
    {358, "L1_DoubleJet_115_40_DoubleJet40_Mass_Min620"},
    {360, "L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28"},
    {359, "L1_DoubleJet_120_45_DoubleJet45_Mass_Min620"},
    {361, "L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28"},
    {366, "L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ"},
    {364, "L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp"},
    {365, "L1_DoubleJet_80_30_Mass_Min420_Mu8"},
    {355, "L1_DoubleJet_90_30_DoubleJet30_Mass_Min620"},
    {383, "L1_DoubleLLPJet40"},
    {230, "L1_DoubleLooseIsoEG22er2p1"},
    {231, "L1_DoubleLooseIsoEG24er2p1"},
    {36, "L1_DoubleMu0"},
    {39, "L1_DoubleMu0_Mass_Min1"},
    {35, "L1_DoubleMu0_OQ"},
    {37, "L1_DoubleMu0_SQ"},
    {38, "L1_DoubleMu0_SQ_OS"},
    {50, "L1_DoubleMu0_Upt15_Upt7"},
    {48, "L1_DoubleMu0_Upt5_Upt5"},
    {49, "L1_DoubleMu0_Upt6_IP_Min1_Upt4"},
    {138, "L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8"},
    {62, "L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6"},
    {61, "L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4"},
    {57, "L1_DoubleMu0er1p5_SQ"},
    {58, "L1_DoubleMu0er1p5_SQ_OS"},
    {60, "L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4"},
    {59, "L1_DoubleMu0er1p5_SQ_dR_Max1p4"},
    {56, "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5"},
    {55, "L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6"},
    {52, "L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4"},
    {54, "L1_DoubleMu0er2p0_SQ_dEta_Max1p5"},
    {53, "L1_DoubleMu0er2p0_SQ_dEta_Max1p6"},
    {51, "L1_DoubleMu0er2p0_SQ_dR_Max1p4"},
    {47, "L1_DoubleMu18er2p1_SQ"},
    {112, "L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20"},
    {141, "L1_DoubleMu3_SQ_ETMHF30_HTT60er"},
    {144, "L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5"},
    {142, "L1_DoubleMu3_SQ_ETMHF40_HTT60er"},
    {145, "L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5"},
    {143, "L1_DoubleMu3_SQ_ETMHF50_HTT60er"},
    {147, "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5"},
    {146, "L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5"},
    {148, "L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5"},
    {150, "L1_DoubleMu3_SQ_HTT220er"},
    {151, "L1_DoubleMu3_SQ_HTT240er"},
    {152, "L1_DoubleMu3_SQ_HTT260er"},
    {139, "L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8"},
    {63, "L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4"},
    {109, "L1_DoubleMu4_SQ_EG9er2p5"},
    {64, "L1_DoubleMu4_SQ_OS"},
    {65, "L1_DoubleMu4_SQ_OS_dR_Max1p2"},
    {66, "L1_DoubleMu4p5_SQ_OS"},
    {67, "L1_DoubleMu4p5_SQ_OS_dR_Max1p2"},
    {68, "L1_DoubleMu4p5er2p0_SQ_OS"},
    {70, "L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18"},
    {69, "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7"},
    {113, "L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20"},
    {110, "L1_DoubleMu5_SQ_EG9er2p5"},
    {40, "L1_DoubleMu8_SQ"},
    {41, "L1_DoubleMu9_SQ"},
    {42, "L1_DoubleMu_12_5"},
    {43, "L1_DoubleMu_15_5_SQ"},
    {44, "L1_DoubleMu_15_7"},
    {46, "L1_DoubleMu_15_7_Mass_Min1"},
    {45, "L1_DoubleMu_15_7_SQ"},
    {267, "L1_DoubleTau70er2p1"},
    {416, "L1_ETM120"},
    {417, "L1_ETM150"},
    {421, "L1_ETMHF100"},
    {430, "L1_ETMHF100_HTT60er"},
    {422, "L1_ETMHF110"},
    {431, "L1_ETMHF110_HTT60er"},
    {444, "L1_ETMHF110_HTT60er_NotSecondBunchInTrain"},
    {423, "L1_ETMHF120"},
    {432, "L1_ETMHF120_HTT60er"},
    {443, "L1_ETMHF120_NotSecondBunchInTrain"},
    {424, "L1_ETMHF130"},
    {433, "L1_ETMHF130_HTT60er"},
    {425, "L1_ETMHF140"},
    {426, "L1_ETMHF150"},
    {418, "L1_ETMHF70"},
    {427, "L1_ETMHF70_HTT60er"},
    {419, "L1_ETMHF80"},
    {428, "L1_ETMHF80_HTT60er"},
    {334, "L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1"},
    {335, "L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6"},
    {420, "L1_ETMHF90"},
    {429, "L1_ETMHF90_HTT60er"},
    {336, "L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1"},
    {337, "L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6"},
    {338, "L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1"},
    {339, "L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6"},
    {410, "L1_ETT1200"},
    {411, "L1_ETT1600"},
    {412, "L1_ETT2000"},
    {477, "L1_FirstBunchAfterTrain"},
    {472, "L1_FirstBunchBeforeTrain"},
    {473, "L1_FirstBunchInTrain"},
    {480, "L1_FirstCollisionInOrbit"},
    {479, "L1_FirstCollisionInTrain"},
    {500, "L1_HCAL_LaserMon_Trig"},
    {501, "L1_HCAL_LaserMon_Veto"},
    {384, "L1_HTT120_SingleLLPJet40"},
    {398, "L1_HTT120er"},
    {385, "L1_HTT160_SingleLLPJet50"},
    {399, "L1_HTT160er"},
    {386, "L1_HTT200_SingleLLPJet60"},
    {400, "L1_HTT200er"},
    {387, "L1_HTT240_SingleLLPJet70"},
    {401, "L1_HTT255er"},
    {402, "L1_HTT280er"},
    {388, "L1_HTT280er_QuadJet_70_55_40_35_er2p5"},
    {403, "L1_HTT320er"},
    {389, "L1_HTT320er_QuadJet_70_55_40_40_er2p5"},
    {390, "L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3"},
    {391, "L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3"},
    {404, "L1_HTT360er"},
    {405, "L1_HTT400er"},
    {406, "L1_HTT450er"},
    {197, "L1_IsoEG32er2p5_Mt40"},
    {298, "L1_IsoTau52er2p1_QuadJet36er2p5"},
    {471, "L1_IsolatedBunch"},
    {476, "L1_LastBunchInTrain"},
    {478, "L1_LastCollisionInTrain"},
    {257, "L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3"},
    {259, "L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3"},
    {241, "L1_LooseIsoEG24er2p1_HTT100er"},
    {258, "L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3"},
    {242, "L1_LooseIsoEG26er2p1_HTT100er"},
    {238, "L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3"},
    {243, "L1_LooseIsoEG28er2p1_HTT100er"},
    {239, "L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3"},
    {244, "L1_LooseIsoEG30er2p1_HTT100er"},
    {240, "L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3"},
    {462, "L1_MinimumBiasHF0"},
    {461, "L1_MinimumBiasHF0_AND_BptxAND"},
    {134, "L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6"},
    {136, "L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6"},
    {135, "L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6"},
    {279, "L1_Mu18er2p1_Tau24er2p1"},
    {280, "L1_Mu18er2p1_Tau26er2p1"},
    {281, "L1_Mu18er2p1_Tau26er2p1_Jet55"},
    {282, "L1_Mu18er2p1_Tau26er2p1_Jet70"},
    {99, "L1_Mu20_EG10er2p5"},
    {285, "L1_Mu22er2p1_IsoTau28er2p1"},
    {286, "L1_Mu22er2p1_IsoTau30er2p1"},
    {287, "L1_Mu22er2p1_IsoTau32er2p1"},
    {288, "L1_Mu22er2p1_IsoTau34er2p1"},
    {289, "L1_Mu22er2p1_IsoTau36er2p1"},
    {290, "L1_Mu22er2p1_IsoTau40er2p1"},
    {291, "L1_Mu22er2p1_Tau70er2p1"},
    {126, "L1_Mu3_Jet120er2p5_dR_Max0p4"},
    {125, "L1_Mu3_Jet120er2p5_dR_Max0p8"},
    {121, "L1_Mu3_Jet16er2p5_dR_Max0p4"},
    {119, "L1_Mu3_Jet30er2p5"},
    {122, "L1_Mu3_Jet35er2p5_dR_Max0p4"},
    {123, "L1_Mu3_Jet60er2p5_dR_Max0p4"},
    {124, "L1_Mu3_Jet80er2p5_dR_Max0p4"},
    {127, "L1_Mu3er1p5_Jet100er2p5_ETMHF30"},
    {128, "L1_Mu3er1p5_Jet100er2p5_ETMHF40"},
    {129, "L1_Mu3er1p5_Jet100er2p5_ETMHF50"},
    {96, "L1_Mu5_EG23er2p5"},
    {100, "L1_Mu5_LooseIsoEG20er2p5"},
    {104, "L1_Mu6_DoubleEG10er2p5"},
    {105, "L1_Mu6_DoubleEG12er2p5"},
    {106, "L1_Mu6_DoubleEG15er2p5"},
    {107, "L1_Mu6_DoubleEG17er2p5"},
    {131, "L1_Mu6_HTT240er"},
    {132, "L1_Mu6_HTT250er"},
    {97, "L1_Mu7_EG20er2p5"},
    {98, "L1_Mu7_EG23er2p5"},
    {101, "L1_Mu7_LooseIsoEG20er2p5"},
    {102, "L1_Mu7_LooseIsoEG23er2p5"},
    {463, "L1_NotBptxOR"},
    {382, "L1_QuadJet60er2p5"},
    {376, "L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0"},
    {90, "L1_QuadMu0"},
    {89, "L1_QuadMu0_OQ"},
    {91, "L1_QuadMu0_SQ"},
    {474, "L1_SecondBunchInTrain"},
    {475, "L1_SecondLastBunchInTrain"},
    {160, "L1_SingleEG10er2p5"},
    {161, "L1_SingleEG15er2p5"},
    {162, "L1_SingleEG26er2p5"},
    {163, "L1_SingleEG28_FWD2p5"},
    {166, "L1_SingleEG28er1p5"},
    {165, "L1_SingleEG28er2p1"},
    {164, "L1_SingleEG28er2p5"},
    {167, "L1_SingleEG34er2p5"},
    {168, "L1_SingleEG36er2p5"},
    {169, "L1_SingleEG38er2p5"},
    {170, "L1_SingleEG40er2p5"},
    {171, "L1_SingleEG42er2p5"},
    {172, "L1_SingleEG45er2p5"},
    {173, "L1_SingleEG50"},
    {174, "L1_SingleEG60"},
    {159, "L1_SingleEG8er2p5"},
    {184, "L1_SingleIsoEG24er1p5"},
    {183, "L1_SingleIsoEG24er2p1"},
    {187, "L1_SingleIsoEG26er1p5"},
    {186, "L1_SingleIsoEG26er2p1"},
    {185, "L1_SingleIsoEG26er2p5"},
    {188, "L1_SingleIsoEG28_FWD2p5"},
    {191, "L1_SingleIsoEG28er1p5"},
    {190, "L1_SingleIsoEG28er2p1"},
    {189, "L1_SingleIsoEG28er2p5"},
    {193, "L1_SingleIsoEG30er2p1"},
    {192, "L1_SingleIsoEG30er2p5"},
    {195, "L1_SingleIsoEG32er2p1"},
    {194, "L1_SingleIsoEG32er2p5"},
    {196, "L1_SingleIsoEG34er2p5"},
    {262, "L1_SingleIsoTau32er2p1"},
    {329, "L1_SingleJet10erHE"},
    {312, "L1_SingleJet120"},
    {327, "L1_SingleJet120_FWD3p0"},
    {319, "L1_SingleJet120er2p5"},
    {330, "L1_SingleJet12erHE"},
    {320, "L1_SingleJet140er2p5"},
    {331, "L1_SingleJet140er2p5_ETMHF70"},
    {332, "L1_SingleJet140er2p5_ETMHF80"},
    {333, "L1_SingleJet140er2p5_ETMHF90"},
    {321, "L1_SingleJet160er2p5"},
    {313, "L1_SingleJet180"},
    {322, "L1_SingleJet180er2p5"},
    {314, "L1_SingleJet200"},
    {451, "L1_SingleJet20er2p5_NotBptxOR"},
    {452, "L1_SingleJet20er2p5_NotBptxOR_3BX"},
    {309, "L1_SingleJet35"},
    {324, "L1_SingleJet35_FWD3p0"},
    {316, "L1_SingleJet35er2p5"},
    {453, "L1_SingleJet43er2p5_NotBptxOR_3BX"},
    {454, "L1_SingleJet46er2p5_NotBptxOR_3BX"},
    {310, "L1_SingleJet60"},
    {325, "L1_SingleJet60_FWD3p0"},
    {317, "L1_SingleJet60er2p5"},
    {328, "L1_SingleJet8erHE"},
    {311, "L1_SingleJet90"},
    {326, "L1_SingleJet90_FWD3p0"},
    {318, "L1_SingleJet90er2p5"},
    {176, "L1_SingleLooseIsoEG26er1p5"},
    {175, "L1_SingleLooseIsoEG26er2p5"},
    {177, "L1_SingleLooseIsoEG28_FWD2p5"},
    {180, "L1_SingleLooseIsoEG28er1p5"},
    {179, "L1_SingleLooseIsoEG28er2p1"},
    {178, "L1_SingleLooseIsoEG28er2p5"},
    {182, "L1_SingleLooseIsoEG30er1p5"},
    {181, "L1_SingleLooseIsoEG30er2p5"},
    {6, "L1_SingleMu0_BMTF"},
    {5, "L1_SingleMu0_DQ"},
    {8, "L1_SingleMu0_EMTF"},
    {7, "L1_SingleMu0_OMTF"},
    {30, "L1_SingleMu10er1p5"},
    {13, "L1_SingleMu12_DQ_BMTF"},
    {15, "L1_SingleMu12_DQ_EMTF"},
    {14, "L1_SingleMu12_DQ_OMTF"},
    {31, "L1_SingleMu12er1p5"},
    {32, "L1_SingleMu14er1p5"},
    {16, "L1_SingleMu15_DQ"},
    {33, "L1_SingleMu16er1p5"},
    {17, "L1_SingleMu18"},
    {34, "L1_SingleMu18er1p5"},
    {18, "L1_SingleMu20"},
    {21, "L1_SingleMu22"},
    {22, "L1_SingleMu22_BMTF"},
    {20, "L1_SingleMu22_DQ"},
    {24, "L1_SingleMu22_EMTF"},
    {23, "L1_SingleMu22_OMTF"},
    {19, "L1_SingleMu22_OQ"},
    {25, "L1_SingleMu25"},
    {9, "L1_SingleMu3"},
    {10, "L1_SingleMu5"},
    {26, "L1_SingleMu6er1p5"},
    {12, "L1_SingleMu7"},
    {11, "L1_SingleMu7_DQ"},
    {27, "L1_SingleMu7er1p5"},
    {28, "L1_SingleMu8er1p5"},
    {29, "L1_SingleMu9er1p5"},
    {0, "L1_SingleMuCosmics"},
    {1, "L1_SingleMuCosmics_BMTF"},
    {3, "L1_SingleMuCosmics_EMTF"},
    {2, "L1_SingleMuCosmics_OMTF"},
    {4, "L1_SingleMuOpen"},
    {446, "L1_SingleMuOpen_NotBptxOR"},
    {448, "L1_SingleMuOpen_er1p1_NotBptxOR_3BX"},
    {447, "L1_SingleMuOpen_er1p4_NotBptxOR_3BX"},
    {92, "L1_SingleMuShower_Nominal"},
    {93, "L1_SingleMuShower_Tight"},
    {264, "L1_SingleTau120er2p1"},
    {265, "L1_SingleTau130er2p1"},
    {263, "L1_SingleTau70er2p1"},
    {503, "L1_TOTEM_1"},
    {504, "L1_TOTEM_2"},
    {505, "L1_TOTEM_3"},
    {506, "L1_TOTEM_4"},
    {236, "L1_TripleEG16er2p5"},
    {232, "L1_TripleEG_16_12_8_er2p5"},
    {233, "L1_TripleEG_16_15_8_er2p5"},
    {234, "L1_TripleEG_18_17_8_er2p5"},
    {235, "L1_TripleEG_18_18_12_er2p5"},
    {373, "L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5"},
    {374, "L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5"},
    {372, "L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5"},
    {72, "L1_TripleMu0"},
    {71, "L1_TripleMu0_OQ"},
    {73, "L1_TripleMu0_SQ"},
    {75, "L1_TripleMu3"},
    {76, "L1_TripleMu3_SQ"},
    {74, "L1_TripleMu_2SQ_1p5SQ_0OQ"},
    {82, "L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12"},
    {83, "L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12"},
    {77, "L1_TripleMu_5SQ_3SQ_0OQ"},
    {87, "L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9"},
    {88, "L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9"},
    {79, "L1_TripleMu_5_3_3"},
    {80, "L1_TripleMu_5_3_3_SQ"},
    {78, "L1_TripleMu_5_3p5_2p5"},
    {85, "L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17"},
    {84, "L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17"},
    {86, "L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17"},
    {81, "L1_TripleMu_5_5_3"},
    {469, "L1_UnpairedBunchBptxMinus"},
    {468, "L1_UnpairedBunchBptxPlus"},
    {459, "L1_ZeroBias"},
    {460, "L1_ZeroBias_copy"}
  };

  const auto rc = id2name.find(index);
  if (rc == id2name.end())
  {
    std::ostringstream oss;
    oss << "no such algorithm index: " << index << ", in menu: L1Menu_Collisions2022_v1_4_0\n";
    throw std::runtime_error(oss.str());
  }
  return rc->second;
}


int getIdFromName(const std::string& name)
{
  static const std::map<std::string, int> name2id = {
  {"L1_AlwaysTrue", 458},
  {"L1_BPTX_AND_Ref1_VME", 486},
  {"L1_BPTX_AND_Ref3_VME", 487},
  {"L1_BPTX_AND_Ref4_VME", 488},
  {"L1_BPTX_BeamGas_B1_VME", 491},
  {"L1_BPTX_BeamGas_B2_VME", 492},
  {"L1_BPTX_BeamGas_Ref1_VME", 489},
  {"L1_BPTX_BeamGas_Ref2_VME", 490},
  {"L1_BPTX_NotOR_VME", 482},
  {"L1_BPTX_OR_Ref3_VME", 483},
  {"L1_BPTX_OR_Ref4_VME", 484},
  {"L1_BPTX_RefAND_VME", 485},
  {"L1_BptxMinus", 467},
  {"L1_BptxOR", 464},
  {"L1_BptxPlus", 466},
  {"L1_BptxXOR", 465},
  {"L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", 494},
  {"L1_DoubleEG10_er1p2_dR_Max0p6", 212},
  {"L1_DoubleEG10p5_er1p2_dR_Max0p6", 213},
  {"L1_DoubleEG11_er1p2_dR_Max0p6", 214},
  {"L1_DoubleEG4_er1p2_dR_Max0p9", 200},
  {"L1_DoubleEG4p5_er1p2_dR_Max0p9", 201},
  {"L1_DoubleEG5_er1p2_dR_Max0p9", 202},
  {"L1_DoubleEG5p5_er1p2_dR_Max0p8", 203},
  {"L1_DoubleEG6_er1p2_dR_Max0p8", 204},
  {"L1_DoubleEG6p5_er1p2_dR_Max0p8", 205},
  {"L1_DoubleEG7_er1p2_dR_Max0p8", 206},
  {"L1_DoubleEG7p5_er1p2_dR_Max0p7", 207},
  {"L1_DoubleEG8_er1p2_dR_Max0p7", 208},
  {"L1_DoubleEG8er2p5_HTT260er", 247},
  {"L1_DoubleEG8er2p5_HTT280er", 248},
  {"L1_DoubleEG8er2p5_HTT300er", 249},
  {"L1_DoubleEG8er2p5_HTT320er", 250},
  {"L1_DoubleEG8er2p5_HTT340er", 251},
  {"L1_DoubleEG8p5_er1p2_dR_Max0p7", 209},
  {"L1_DoubleEG9_er1p2_dR_Max0p7", 210},
  {"L1_DoubleEG9p5_er1p2_dR_Max0p6", 211},
  {"L1_DoubleEG_15_10_er2p5", 215},
  {"L1_DoubleEG_20_10_er2p5", 216},
  {"L1_DoubleEG_22_10_er2p5", 217},
  {"L1_DoubleEG_25_12_er2p5", 218},
  {"L1_DoubleEG_25_14_er2p5", 219},
  {"L1_DoubleEG_27_14_er2p5", 220},
  {"L1_DoubleEG_LooseIso16_LooseIso12_er1p5", 225},
  {"L1_DoubleEG_LooseIso18_LooseIso12_er1p5", 226},
  {"L1_DoubleEG_LooseIso20_10_er2p5", 221},
  {"L1_DoubleEG_LooseIso20_LooseIso12_er1p5", 227},
  {"L1_DoubleEG_LooseIso22_10_er2p5", 222},
  {"L1_DoubleEG_LooseIso22_12_er2p5", 223},
  {"L1_DoubleEG_LooseIso22_LooseIso12_er1p5", 228},
  {"L1_DoubleEG_LooseIso25_12_er2p5", 224},
  {"L1_DoubleEG_LooseIso25_LooseIso12_er1p5", 229},
  {"L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5", 283},
  {"L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5", 284},
  {"L1_DoubleIsoTau28er2p1", 268},
  {"L1_DoubleIsoTau28er2p1_Mass_Max80", 275},
  {"L1_DoubleIsoTau28er2p1_Mass_Max90", 274},
  {"L1_DoubleIsoTau30er2p1", 269},
  {"L1_DoubleIsoTau30er2p1_Mass_Max80", 277},
  {"L1_DoubleIsoTau30er2p1_Mass_Max90", 276},
  {"L1_DoubleIsoTau32er2p1", 270},
  {"L1_DoubleIsoTau34er2p1", 271},
  {"L1_DoubleIsoTau35er2p1", 272},
  {"L1_DoubleIsoTau36er2p1", 273},
  {"L1_DoubleJet100er2p3_dEta_Max1p6", 345},
  {"L1_DoubleJet100er2p5", 341},
  {"L1_DoubleJet112er2p3_dEta_Max1p6", 346},
  {"L1_DoubleJet120er2p5", 342},
  {"L1_DoubleJet150er2p5", 343},
  {"L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", 348},
  {"L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", 349},
  {"L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", 350},
  {"L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", 351},
  {"L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", 352},
  {"L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", 353},
  {"L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", 363},
  {"L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5", 362},
  {"L1_DoubleJet40er2p5", 340},
  {"L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", 356},
  {"L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", 357},
  {"L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", 358},
  {"L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", 360},
  {"L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", 359},
  {"L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", 361},
  {"L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", 366},
  {"L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", 364},
  {"L1_DoubleJet_80_30_Mass_Min420_Mu8", 365},
  {"L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", 355},
  {"L1_DoubleLLPJet40", 383},
  {"L1_DoubleLooseIsoEG22er2p1", 230},
  {"L1_DoubleLooseIsoEG24er2p1", 231},
  {"L1_DoubleMu0", 36},
  {"L1_DoubleMu0_Mass_Min1", 39},
  {"L1_DoubleMu0_OQ", 35},
  {"L1_DoubleMu0_SQ", 37},
  {"L1_DoubleMu0_SQ_OS", 38},
  {"L1_DoubleMu0_Upt15_Upt7", 50},
  {"L1_DoubleMu0_Upt5_Upt5", 48},
  {"L1_DoubleMu0_Upt6_IP_Min1_Upt4", 49},
  {"L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", 138},
  {"L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", 62},
  {"L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", 61},
  {"L1_DoubleMu0er1p5_SQ", 57},
  {"L1_DoubleMu0er1p5_SQ_OS", 58},
  {"L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", 60},
  {"L1_DoubleMu0er1p5_SQ_dR_Max1p4", 59},
  {"L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", 56},
  {"L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", 55},
  {"L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", 52},
  {"L1_DoubleMu0er2p0_SQ_dEta_Max1p5", 54},
  {"L1_DoubleMu0er2p0_SQ_dEta_Max1p6", 53},
  {"L1_DoubleMu0er2p0_SQ_dR_Max1p4", 51},
  {"L1_DoubleMu18er2p1_SQ", 47},
  {"L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20", 112},
  {"L1_DoubleMu3_SQ_ETMHF30_HTT60er", 141},
  {"L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5", 144},
  {"L1_DoubleMu3_SQ_ETMHF40_HTT60er", 142},
  {"L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5", 145},
  {"L1_DoubleMu3_SQ_ETMHF50_HTT60er", 143},
  {"L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", 147},
  {"L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", 146},
  {"L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", 148},
  {"L1_DoubleMu3_SQ_HTT220er", 150},
  {"L1_DoubleMu3_SQ_HTT240er", 151},
  {"L1_DoubleMu3_SQ_HTT260er", 152},
  {"L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", 139},
  {"L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4", 63},
  {"L1_DoubleMu4_SQ_EG9er2p5", 109},
  {"L1_DoubleMu4_SQ_OS", 64},
  {"L1_DoubleMu4_SQ_OS_dR_Max1p2", 65},
  {"L1_DoubleMu4p5_SQ_OS", 66},
  {"L1_DoubleMu4p5_SQ_OS_dR_Max1p2", 67},
  {"L1_DoubleMu4p5er2p0_SQ_OS", 68},
  {"L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18", 70},
  {"L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", 69},
  {"L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20", 113},
  {"L1_DoubleMu5_SQ_EG9er2p5", 110},
  {"L1_DoubleMu8_SQ", 40},
  {"L1_DoubleMu9_SQ", 41},
  {"L1_DoubleMu_12_5", 42},
  {"L1_DoubleMu_15_5_SQ", 43},
  {"L1_DoubleMu_15_7", 44},
  {"L1_DoubleMu_15_7_Mass_Min1", 46},
  {"L1_DoubleMu_15_7_SQ", 45},
  {"L1_DoubleTau70er2p1", 267},
  {"L1_ETM120", 416},
  {"L1_ETM150", 417},
  {"L1_ETMHF100", 421},
  {"L1_ETMHF100_HTT60er", 430},
  {"L1_ETMHF110", 422},
  {"L1_ETMHF110_HTT60er", 431},
  {"L1_ETMHF110_HTT60er_NotSecondBunchInTrain", 444},
  {"L1_ETMHF120", 423},
  {"L1_ETMHF120_HTT60er", 432},
  {"L1_ETMHF120_NotSecondBunchInTrain", 443},
  {"L1_ETMHF130", 424},
  {"L1_ETMHF130_HTT60er", 433},
  {"L1_ETMHF140", 425},
  {"L1_ETMHF150", 426},
  {"L1_ETMHF70", 418},
  {"L1_ETMHF70_HTT60er", 427},
  {"L1_ETMHF80", 419},
  {"L1_ETMHF80_HTT60er", 428},
  {"L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1", 334},
  {"L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6", 335},
  {"L1_ETMHF90", 420},
  {"L1_ETMHF90_HTT60er", 429},
  {"L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1", 336},
  {"L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6", 337},
  {"L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1", 338},
  {"L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6", 339},
  {"L1_ETT1200", 410},
  {"L1_ETT1600", 411},
  {"L1_ETT2000", 412},
  {"L1_FirstBunchAfterTrain", 477},
  {"L1_FirstBunchBeforeTrain", 472},
  {"L1_FirstBunchInTrain", 473},
  {"L1_FirstCollisionInOrbit", 480},
  {"L1_FirstCollisionInTrain", 479},
  {"L1_HCAL_LaserMon_Trig", 500},
  {"L1_HCAL_LaserMon_Veto", 501},
  {"L1_HTT120_SingleLLPJet40", 384},
  {"L1_HTT120er", 398},
  {"L1_HTT160_SingleLLPJet50", 385},
  {"L1_HTT160er", 399},
  {"L1_HTT200_SingleLLPJet60", 386},
  {"L1_HTT200er", 400},
  {"L1_HTT240_SingleLLPJet70", 387},
  {"L1_HTT255er", 401},
  {"L1_HTT280er", 402},
  {"L1_HTT280er_QuadJet_70_55_40_35_er2p5", 388},
  {"L1_HTT320er", 403},
  {"L1_HTT320er_QuadJet_70_55_40_40_er2p5", 389},
  {"L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", 390},
  {"L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", 391},
  {"L1_HTT360er", 404},
  {"L1_HTT400er", 405},
  {"L1_HTT450er", 406},
  {"L1_IsoEG32er2p5_Mt40", 197},
  {"L1_IsoTau52er2p1_QuadJet36er2p5", 298},
  {"L1_IsolatedBunch", 471},
  {"L1_LastBunchInTrain", 476},
  {"L1_LastCollisionInTrain", 478},
  {"L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", 257},
  {"L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", 259},
  {"L1_LooseIsoEG24er2p1_HTT100er", 241},
  {"L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", 258},
  {"L1_LooseIsoEG26er2p1_HTT100er", 242},
  {"L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", 238},
  {"L1_LooseIsoEG28er2p1_HTT100er", 243},
  {"L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", 239},
  {"L1_LooseIsoEG30er2p1_HTT100er", 244},
  {"L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", 240},
  {"L1_MinimumBiasHF0", 462},
  {"L1_MinimumBiasHF0_AND_BptxAND", 461},
  {"L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", 134},
  {"L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", 136},
  {"L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", 135},
  {"L1_Mu18er2p1_Tau24er2p1", 279},
  {"L1_Mu18er2p1_Tau26er2p1", 280},
  {"L1_Mu18er2p1_Tau26er2p1_Jet55", 281},
  {"L1_Mu18er2p1_Tau26er2p1_Jet70", 282},
  {"L1_Mu20_EG10er2p5", 99},
  {"L1_Mu22er2p1_IsoTau28er2p1", 285},
  {"L1_Mu22er2p1_IsoTau30er2p1", 286},
  {"L1_Mu22er2p1_IsoTau32er2p1", 287},
  {"L1_Mu22er2p1_IsoTau34er2p1", 288},
  {"L1_Mu22er2p1_IsoTau36er2p1", 289},
  {"L1_Mu22er2p1_IsoTau40er2p1", 290},
  {"L1_Mu22er2p1_Tau70er2p1", 291},
  {"L1_Mu3_Jet120er2p5_dR_Max0p4", 126},
  {"L1_Mu3_Jet120er2p5_dR_Max0p8", 125},
  {"L1_Mu3_Jet16er2p5_dR_Max0p4", 121},
  {"L1_Mu3_Jet30er2p5", 119},
  {"L1_Mu3_Jet35er2p5_dR_Max0p4", 122},
  {"L1_Mu3_Jet60er2p5_dR_Max0p4", 123},
  {"L1_Mu3_Jet80er2p5_dR_Max0p4", 124},
  {"L1_Mu3er1p5_Jet100er2p5_ETMHF30", 127},
  {"L1_Mu3er1p5_Jet100er2p5_ETMHF40", 128},
  {"L1_Mu3er1p5_Jet100er2p5_ETMHF50", 129},
  {"L1_Mu5_EG23er2p5", 96},
  {"L1_Mu5_LooseIsoEG20er2p5", 100},
  {"L1_Mu6_DoubleEG10er2p5", 104},
  {"L1_Mu6_DoubleEG12er2p5", 105},
  {"L1_Mu6_DoubleEG15er2p5", 106},
  {"L1_Mu6_DoubleEG17er2p5", 107},
  {"L1_Mu6_HTT240er", 131},
  {"L1_Mu6_HTT250er", 132},
  {"L1_Mu7_EG20er2p5", 97},
  {"L1_Mu7_EG23er2p5", 98},
  {"L1_Mu7_LooseIsoEG20er2p5", 101},
  {"L1_Mu7_LooseIsoEG23er2p5", 102},
  {"L1_NotBptxOR", 463},
  {"L1_QuadJet60er2p5", 382},
  {"L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", 376},
  {"L1_QuadMu0", 90},
  {"L1_QuadMu0_OQ", 89},
  {"L1_QuadMu0_SQ", 91},
  {"L1_SecondBunchInTrain", 474},
  {"L1_SecondLastBunchInTrain", 475},
  {"L1_SingleEG10er2p5", 160},
  {"L1_SingleEG15er2p5", 161},
  {"L1_SingleEG26er2p5", 162},
  {"L1_SingleEG28_FWD2p5", 163},
  {"L1_SingleEG28er1p5", 166},
  {"L1_SingleEG28er2p1", 165},
  {"L1_SingleEG28er2p5", 164},
  {"L1_SingleEG34er2p5", 167},
  {"L1_SingleEG36er2p5", 168},
  {"L1_SingleEG38er2p5", 169},
  {"L1_SingleEG40er2p5", 170},
  {"L1_SingleEG42er2p5", 171},
  {"L1_SingleEG45er2p5", 172},
  {"L1_SingleEG50", 173},
  {"L1_SingleEG60", 174},
  {"L1_SingleEG8er2p5", 159},
  {"L1_SingleIsoEG24er1p5", 184},
  {"L1_SingleIsoEG24er2p1", 183},
  {"L1_SingleIsoEG26er1p5", 187},
  {"L1_SingleIsoEG26er2p1", 186},
  {"L1_SingleIsoEG26er2p5", 185},
  {"L1_SingleIsoEG28_FWD2p5", 188},
  {"L1_SingleIsoEG28er1p5", 191},
  {"L1_SingleIsoEG28er2p1", 190},
  {"L1_SingleIsoEG28er2p5", 189},
  {"L1_SingleIsoEG30er2p1", 193},
  {"L1_SingleIsoEG30er2p5", 192},
  {"L1_SingleIsoEG32er2p1", 195},
  {"L1_SingleIsoEG32er2p5", 194},
  {"L1_SingleIsoEG34er2p5", 196},
  {"L1_SingleIsoTau32er2p1", 262},
  {"L1_SingleJet10erHE", 329},
  {"L1_SingleJet120", 312},
  {"L1_SingleJet120_FWD3p0", 327},
  {"L1_SingleJet120er2p5", 319},
  {"L1_SingleJet12erHE", 330},
  {"L1_SingleJet140er2p5", 320},
  {"L1_SingleJet140er2p5_ETMHF70", 331},
  {"L1_SingleJet140er2p5_ETMHF80", 332},
  {"L1_SingleJet140er2p5_ETMHF90", 333},
  {"L1_SingleJet160er2p5", 321},
  {"L1_SingleJet180", 313},
  {"L1_SingleJet180er2p5", 322},
  {"L1_SingleJet200", 314},
  {"L1_SingleJet20er2p5_NotBptxOR", 451},
  {"L1_SingleJet20er2p5_NotBptxOR_3BX", 452},
  {"L1_SingleJet35", 309},
  {"L1_SingleJet35_FWD3p0", 324},
  {"L1_SingleJet35er2p5", 316},
  {"L1_SingleJet43er2p5_NotBptxOR_3BX", 453},
  {"L1_SingleJet46er2p5_NotBptxOR_3BX", 454},
  {"L1_SingleJet60", 310},
  {"L1_SingleJet60_FWD3p0", 325},
  {"L1_SingleJet60er2p5", 317},
  {"L1_SingleJet8erHE", 328},
  {"L1_SingleJet90", 311},
  {"L1_SingleJet90_FWD3p0", 326},
  {"L1_SingleJet90er2p5", 318},
  {"L1_SingleLooseIsoEG26er1p5", 176},
  {"L1_SingleLooseIsoEG26er2p5", 175},
  {"L1_SingleLooseIsoEG28_FWD2p5", 177},
  {"L1_SingleLooseIsoEG28er1p5", 180},
  {"L1_SingleLooseIsoEG28er2p1", 179},
  {"L1_SingleLooseIsoEG28er2p5", 178},
  {"L1_SingleLooseIsoEG30er1p5", 182},
  {"L1_SingleLooseIsoEG30er2p5", 181},
  {"L1_SingleMu0_BMTF", 6},
  {"L1_SingleMu0_DQ", 5},
  {"L1_SingleMu0_EMTF", 8},
  {"L1_SingleMu0_OMTF", 7},
  {"L1_SingleMu10er1p5", 30},
  {"L1_SingleMu12_DQ_BMTF", 13},
  {"L1_SingleMu12_DQ_EMTF", 15},
  {"L1_SingleMu12_DQ_OMTF", 14},
  {"L1_SingleMu12er1p5", 31},
  {"L1_SingleMu14er1p5", 32},
  {"L1_SingleMu15_DQ", 16},
  {"L1_SingleMu16er1p5", 33},
  {"L1_SingleMu18", 17},
  {"L1_SingleMu18er1p5", 34},
  {"L1_SingleMu20", 18},
  {"L1_SingleMu22", 21},
  {"L1_SingleMu22_BMTF", 22},
  {"L1_SingleMu22_DQ", 20},
  {"L1_SingleMu22_EMTF", 24},
  {"L1_SingleMu22_OMTF", 23},
  {"L1_SingleMu22_OQ", 19},
  {"L1_SingleMu25", 25},
  {"L1_SingleMu3", 9},
  {"L1_SingleMu5", 10},
  {"L1_SingleMu6er1p5", 26},
  {"L1_SingleMu7", 12},
  {"L1_SingleMu7_DQ", 11},
  {"L1_SingleMu7er1p5", 27},
  {"L1_SingleMu8er1p5", 28},
  {"L1_SingleMu9er1p5", 29},
  {"L1_SingleMuCosmics", 0},
  {"L1_SingleMuCosmics_BMTF", 1},
  {"L1_SingleMuCosmics_EMTF", 3},
  {"L1_SingleMuCosmics_OMTF", 2},
  {"L1_SingleMuOpen", 4},
  {"L1_SingleMuOpen_NotBptxOR", 446},
  {"L1_SingleMuOpen_er1p1_NotBptxOR_3BX", 448},
  {"L1_SingleMuOpen_er1p4_NotBptxOR_3BX", 447},
  {"L1_SingleMuShower_Nominal", 92},
  {"L1_SingleMuShower_Tight", 93},
  {"L1_SingleTau120er2p1", 264},
  {"L1_SingleTau130er2p1", 265},
  {"L1_SingleTau70er2p1", 263},
  {"L1_TOTEM_1", 503},
  {"L1_TOTEM_2", 504},
  {"L1_TOTEM_3", 505},
  {"L1_TOTEM_4", 506},
  {"L1_TripleEG16er2p5", 236},
  {"L1_TripleEG_16_12_8_er2p5", 232},
  {"L1_TripleEG_16_15_8_er2p5", 233},
  {"L1_TripleEG_18_17_8_er2p5", 234},
  {"L1_TripleEG_18_18_12_er2p5", 235},
  {"L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", 373},
  {"L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", 374},
  {"L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", 372},
  {"L1_TripleMu0", 72},
  {"L1_TripleMu0_OQ", 71},
  {"L1_TripleMu0_SQ", 73},
  {"L1_TripleMu3", 75},
  {"L1_TripleMu3_SQ", 76},
  {"L1_TripleMu_2SQ_1p5SQ_0OQ", 74},
  {"L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12", 82},
  {"L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12", 83},
  {"L1_TripleMu_5SQ_3SQ_0OQ", 77},
  {"L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", 87},
  {"L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", 88},
  {"L1_TripleMu_5_3_3", 79},
  {"L1_TripleMu_5_3_3_SQ", 80},
  {"L1_TripleMu_5_3p5_2p5", 78},
  {"L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", 85},
  {"L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", 84},
  {"L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", 86},
  {"L1_TripleMu_5_5_3", 81},
  {"L1_UnpairedBunchBptxMinus", 469},
  {"L1_UnpairedBunchBptxPlus", 468},
  {"L1_ZeroBias", 459},
  {"L1_ZeroBias_copy", 460}
  };

  const auto rc = name2id.find(name);
  if (rc == name2id.end())
  {
    std::ostringstream oss;
    oss << "no such algorithm name: \"" << name << "\", in menu: L1Menu_Collisions2022_v1_4_0\n";
    throw std::runtime_error(oss.str());
  }
  return rc->second;
}


AlgorithmFunction getFuncFromId(const size_t index)
{
  static const std::map<size_t, AlgorithmFunction> id2func = {
    {458, &L1_AlwaysTrue},
    {486, &L1_BPTX_AND_Ref1_VME},
    {487, &L1_BPTX_AND_Ref3_VME},
    {488, &L1_BPTX_AND_Ref4_VME},
    {491, &L1_BPTX_BeamGas_B1_VME},
    {492, &L1_BPTX_BeamGas_B2_VME},
    {489, &L1_BPTX_BeamGas_Ref1_VME},
    {490, &L1_BPTX_BeamGas_Ref2_VME},
    {482, &L1_BPTX_NotOR_VME},
    {483, &L1_BPTX_OR_Ref3_VME},
    {484, &L1_BPTX_OR_Ref4_VME},
    {485, &L1_BPTX_RefAND_VME},
    {467, &L1_BptxMinus},
    {464, &L1_BptxOR},
    {466, &L1_BptxPlus},
    {465, &L1_BptxXOR},
    {494, &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142},
    {212, &L1_DoubleEG10_er1p2_dR_Max0p6},
    {213, &L1_DoubleEG10p5_er1p2_dR_Max0p6},
    {214, &L1_DoubleEG11_er1p2_dR_Max0p6},
    {200, &L1_DoubleEG4_er1p2_dR_Max0p9},
    {201, &L1_DoubleEG4p5_er1p2_dR_Max0p9},
    {202, &L1_DoubleEG5_er1p2_dR_Max0p9},
    {203, &L1_DoubleEG5p5_er1p2_dR_Max0p8},
    {204, &L1_DoubleEG6_er1p2_dR_Max0p8},
    {205, &L1_DoubleEG6p5_er1p2_dR_Max0p8},
    {206, &L1_DoubleEG7_er1p2_dR_Max0p8},
    {207, &L1_DoubleEG7p5_er1p2_dR_Max0p7},
    {208, &L1_DoubleEG8_er1p2_dR_Max0p7},
    {247, &L1_DoubleEG8er2p5_HTT260er},
    {248, &L1_DoubleEG8er2p5_HTT280er},
    {249, &L1_DoubleEG8er2p5_HTT300er},
    {250, &L1_DoubleEG8er2p5_HTT320er},
    {251, &L1_DoubleEG8er2p5_HTT340er},
    {209, &L1_DoubleEG8p5_er1p2_dR_Max0p7},
    {210, &L1_DoubleEG9_er1p2_dR_Max0p7},
    {211, &L1_DoubleEG9p5_er1p2_dR_Max0p6},
    {215, &L1_DoubleEG_15_10_er2p5},
    {216, &L1_DoubleEG_20_10_er2p5},
    {217, &L1_DoubleEG_22_10_er2p5},
    {218, &L1_DoubleEG_25_12_er2p5},
    {219, &L1_DoubleEG_25_14_er2p5},
    {220, &L1_DoubleEG_27_14_er2p5},
    {225, &L1_DoubleEG_LooseIso16_LooseIso12_er1p5},
    {226, &L1_DoubleEG_LooseIso18_LooseIso12_er1p5},
    {221, &L1_DoubleEG_LooseIso20_10_er2p5},
    {227, &L1_DoubleEG_LooseIso20_LooseIso12_er1p5},
    {222, &L1_DoubleEG_LooseIso22_10_er2p5},
    {223, &L1_DoubleEG_LooseIso22_12_er2p5},
    {228, &L1_DoubleEG_LooseIso22_LooseIso12_er1p5},
    {224, &L1_DoubleEG_LooseIso25_12_er2p5},
    {229, &L1_DoubleEG_LooseIso25_LooseIso12_er1p5},
    {283, &L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5},
    {284, &L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5},
    {268, &L1_DoubleIsoTau28er2p1},
    {275, &L1_DoubleIsoTau28er2p1_Mass_Max80},
    {274, &L1_DoubleIsoTau28er2p1_Mass_Max90},
    {269, &L1_DoubleIsoTau30er2p1},
    {277, &L1_DoubleIsoTau30er2p1_Mass_Max80},
    {276, &L1_DoubleIsoTau30er2p1_Mass_Max90},
    {270, &L1_DoubleIsoTau32er2p1},
    {271, &L1_DoubleIsoTau34er2p1},
    {272, &L1_DoubleIsoTau35er2p1},
    {273, &L1_DoubleIsoTau36er2p1},
    {345, &L1_DoubleJet100er2p3_dEta_Max1p6},
    {341, &L1_DoubleJet100er2p5},
    {346, &L1_DoubleJet112er2p3_dEta_Max1p6},
    {342, &L1_DoubleJet120er2p5},
    {343, &L1_DoubleJet150er2p5},
    {348, &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5},
    {349, &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5},
    {350, &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5},
    {351, &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5},
    {352, &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5},
    {353, &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5},
    {363, &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp},
    {362, &L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5},
    {340, &L1_DoubleJet40er2p5},
    {356, &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620},
    {357, &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620},
    {358, &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620},
    {360, &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28},
    {359, &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620},
    {361, &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28},
    {366, &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ},
    {364, &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp},
    {365, &L1_DoubleJet_80_30_Mass_Min420_Mu8},
    {355, &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620},
    {383, &L1_DoubleLLPJet40},
    {230, &L1_DoubleLooseIsoEG22er2p1},
    {231, &L1_DoubleLooseIsoEG24er2p1},
    {36, &L1_DoubleMu0},
    {39, &L1_DoubleMu0_Mass_Min1},
    {35, &L1_DoubleMu0_OQ},
    {37, &L1_DoubleMu0_SQ},
    {38, &L1_DoubleMu0_SQ_OS},
    {50, &L1_DoubleMu0_Upt15_Upt7},
    {48, &L1_DoubleMu0_Upt5_Upt5},
    {49, &L1_DoubleMu0_Upt6_IP_Min1_Upt4},
    {138, &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8},
    {62, &L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6},
    {61, &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4},
    {57, &L1_DoubleMu0er1p5_SQ},
    {58, &L1_DoubleMu0er1p5_SQ_OS},
    {60, &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4},
    {59, &L1_DoubleMu0er1p5_SQ_dR_Max1p4},
    {56, &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5},
    {55, &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6},
    {52, &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4},
    {54, &L1_DoubleMu0er2p0_SQ_dEta_Max1p5},
    {53, &L1_DoubleMu0er2p0_SQ_dEta_Max1p6},
    {51, &L1_DoubleMu0er2p0_SQ_dR_Max1p4},
    {47, &L1_DoubleMu18er2p1_SQ},
    {112, &L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20},
    {141, &L1_DoubleMu3_SQ_ETMHF30_HTT60er},
    {144, &L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5},
    {142, &L1_DoubleMu3_SQ_ETMHF40_HTT60er},
    {145, &L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5},
    {143, &L1_DoubleMu3_SQ_ETMHF50_HTT60er},
    {147, &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5},
    {146, &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5},
    {148, &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5},
    {150, &L1_DoubleMu3_SQ_HTT220er},
    {151, &L1_DoubleMu3_SQ_HTT240er},
    {152, &L1_DoubleMu3_SQ_HTT260er},
    {139, &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8},
    {63, &L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4},
    {109, &L1_DoubleMu4_SQ_EG9er2p5},
    {64, &L1_DoubleMu4_SQ_OS},
    {65, &L1_DoubleMu4_SQ_OS_dR_Max1p2},
    {66, &L1_DoubleMu4p5_SQ_OS},
    {67, &L1_DoubleMu4p5_SQ_OS_dR_Max1p2},
    {68, &L1_DoubleMu4p5er2p0_SQ_OS},
    {70, &L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18},
    {69, &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7},
    {113, &L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20},
    {110, &L1_DoubleMu5_SQ_EG9er2p5},
    {40, &L1_DoubleMu8_SQ},
    {41, &L1_DoubleMu9_SQ},
    {42, &L1_DoubleMu_12_5},
    {43, &L1_DoubleMu_15_5_SQ},
    {44, &L1_DoubleMu_15_7},
    {46, &L1_DoubleMu_15_7_Mass_Min1},
    {45, &L1_DoubleMu_15_7_SQ},
    {267, &L1_DoubleTau70er2p1},
    {416, &L1_ETM120},
    {417, &L1_ETM150},
    {421, &L1_ETMHF100},
    {430, &L1_ETMHF100_HTT60er},
    {422, &L1_ETMHF110},
    {431, &L1_ETMHF110_HTT60er},
    {444, &L1_ETMHF110_HTT60er_NotSecondBunchInTrain},
    {423, &L1_ETMHF120},
    {432, &L1_ETMHF120_HTT60er},
    {443, &L1_ETMHF120_NotSecondBunchInTrain},
    {424, &L1_ETMHF130},
    {433, &L1_ETMHF130_HTT60er},
    {425, &L1_ETMHF140},
    {426, &L1_ETMHF150},
    {418, &L1_ETMHF70},
    {427, &L1_ETMHF70_HTT60er},
    {419, &L1_ETMHF80},
    {428, &L1_ETMHF80_HTT60er},
    {334, &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1},
    {335, &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6},
    {420, &L1_ETMHF90},
    {429, &L1_ETMHF90_HTT60er},
    {336, &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1},
    {337, &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6},
    {338, &L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1},
    {339, &L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6},
    {410, &L1_ETT1200},
    {411, &L1_ETT1600},
    {412, &L1_ETT2000},
    {477, &L1_FirstBunchAfterTrain},
    {472, &L1_FirstBunchBeforeTrain},
    {473, &L1_FirstBunchInTrain},
    {480, &L1_FirstCollisionInOrbit},
    {479, &L1_FirstCollisionInTrain},
    {500, &L1_HCAL_LaserMon_Trig},
    {501, &L1_HCAL_LaserMon_Veto},
    {384, &L1_HTT120_SingleLLPJet40},
    {398, &L1_HTT120er},
    {385, &L1_HTT160_SingleLLPJet50},
    {399, &L1_HTT160er},
    {386, &L1_HTT200_SingleLLPJet60},
    {400, &L1_HTT200er},
    {387, &L1_HTT240_SingleLLPJet70},
    {401, &L1_HTT255er},
    {402, &L1_HTT280er},
    {388, &L1_HTT280er_QuadJet_70_55_40_35_er2p5},
    {403, &L1_HTT320er},
    {389, &L1_HTT320er_QuadJet_70_55_40_40_er2p5},
    {390, &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3},
    {391, &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3},
    {404, &L1_HTT360er},
    {405, &L1_HTT400er},
    {406, &L1_HTT450er},
    {197, &L1_IsoEG32er2p5_Mt40},
    {298, &L1_IsoTau52er2p1_QuadJet36er2p5},
    {471, &L1_IsolatedBunch},
    {476, &L1_LastBunchInTrain},
    {478, &L1_LastCollisionInTrain},
    {257, &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3},
    {259, &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3},
    {241, &L1_LooseIsoEG24er2p1_HTT100er},
    {258, &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3},
    {242, &L1_LooseIsoEG26er2p1_HTT100er},
    {238, &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3},
    {243, &L1_LooseIsoEG28er2p1_HTT100er},
    {239, &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3},
    {244, &L1_LooseIsoEG30er2p1_HTT100er},
    {240, &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3},
    {462, &L1_MinimumBiasHF0},
    {461, &L1_MinimumBiasHF0_AND_BptxAND},
    {134, &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6},
    {136, &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6},
    {135, &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6},
    {279, &L1_Mu18er2p1_Tau24er2p1},
    {280, &L1_Mu18er2p1_Tau26er2p1},
    {281, &L1_Mu18er2p1_Tau26er2p1_Jet55},
    {282, &L1_Mu18er2p1_Tau26er2p1_Jet70},
    {99, &L1_Mu20_EG10er2p5},
    {285, &L1_Mu22er2p1_IsoTau28er2p1},
    {286, &L1_Mu22er2p1_IsoTau30er2p1},
    {287, &L1_Mu22er2p1_IsoTau32er2p1},
    {288, &L1_Mu22er2p1_IsoTau34er2p1},
    {289, &L1_Mu22er2p1_IsoTau36er2p1},
    {290, &L1_Mu22er2p1_IsoTau40er2p1},
    {291, &L1_Mu22er2p1_Tau70er2p1},
    {126, &L1_Mu3_Jet120er2p5_dR_Max0p4},
    {125, &L1_Mu3_Jet120er2p5_dR_Max0p8},
    {121, &L1_Mu3_Jet16er2p5_dR_Max0p4},
    {119, &L1_Mu3_Jet30er2p5},
    {122, &L1_Mu3_Jet35er2p5_dR_Max0p4},
    {123, &L1_Mu3_Jet60er2p5_dR_Max0p4},
    {124, &L1_Mu3_Jet80er2p5_dR_Max0p4},
    {127, &L1_Mu3er1p5_Jet100er2p5_ETMHF30},
    {128, &L1_Mu3er1p5_Jet100er2p5_ETMHF40},
    {129, &L1_Mu3er1p5_Jet100er2p5_ETMHF50},
    {96, &L1_Mu5_EG23er2p5},
    {100, &L1_Mu5_LooseIsoEG20er2p5},
    {104, &L1_Mu6_DoubleEG10er2p5},
    {105, &L1_Mu6_DoubleEG12er2p5},
    {106, &L1_Mu6_DoubleEG15er2p5},
    {107, &L1_Mu6_DoubleEG17er2p5},
    {131, &L1_Mu6_HTT240er},
    {132, &L1_Mu6_HTT250er},
    {97, &L1_Mu7_EG20er2p5},
    {98, &L1_Mu7_EG23er2p5},
    {101, &L1_Mu7_LooseIsoEG20er2p5},
    {102, &L1_Mu7_LooseIsoEG23er2p5},
    {463, &L1_NotBptxOR},
    {382, &L1_QuadJet60er2p5},
    {376, &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0},
    {90, &L1_QuadMu0},
    {89, &L1_QuadMu0_OQ},
    {91, &L1_QuadMu0_SQ},
    {474, &L1_SecondBunchInTrain},
    {475, &L1_SecondLastBunchInTrain},
    {160, &L1_SingleEG10er2p5},
    {161, &L1_SingleEG15er2p5},
    {162, &L1_SingleEG26er2p5},
    {163, &L1_SingleEG28_FWD2p5},
    {166, &L1_SingleEG28er1p5},
    {165, &L1_SingleEG28er2p1},
    {164, &L1_SingleEG28er2p5},
    {167, &L1_SingleEG34er2p5},
    {168, &L1_SingleEG36er2p5},
    {169, &L1_SingleEG38er2p5},
    {170, &L1_SingleEG40er2p5},
    {171, &L1_SingleEG42er2p5},
    {172, &L1_SingleEG45er2p5},
    {173, &L1_SingleEG50},
    {174, &L1_SingleEG60},
    {159, &L1_SingleEG8er2p5},
    {184, &L1_SingleIsoEG24er1p5},
    {183, &L1_SingleIsoEG24er2p1},
    {187, &L1_SingleIsoEG26er1p5},
    {186, &L1_SingleIsoEG26er2p1},
    {185, &L1_SingleIsoEG26er2p5},
    {188, &L1_SingleIsoEG28_FWD2p5},
    {191, &L1_SingleIsoEG28er1p5},
    {190, &L1_SingleIsoEG28er2p1},
    {189, &L1_SingleIsoEG28er2p5},
    {193, &L1_SingleIsoEG30er2p1},
    {192, &L1_SingleIsoEG30er2p5},
    {195, &L1_SingleIsoEG32er2p1},
    {194, &L1_SingleIsoEG32er2p5},
    {196, &L1_SingleIsoEG34er2p5},
    {262, &L1_SingleIsoTau32er2p1},
    {329, &L1_SingleJet10erHE},
    {312, &L1_SingleJet120},
    {327, &L1_SingleJet120_FWD3p0},
    {319, &L1_SingleJet120er2p5},
    {330, &L1_SingleJet12erHE},
    {320, &L1_SingleJet140er2p5},
    {331, &L1_SingleJet140er2p5_ETMHF70},
    {332, &L1_SingleJet140er2p5_ETMHF80},
    {333, &L1_SingleJet140er2p5_ETMHF90},
    {321, &L1_SingleJet160er2p5},
    {313, &L1_SingleJet180},
    {322, &L1_SingleJet180er2p5},
    {314, &L1_SingleJet200},
    {451, &L1_SingleJet20er2p5_NotBptxOR},
    {452, &L1_SingleJet20er2p5_NotBptxOR_3BX},
    {309, &L1_SingleJet35},
    {324, &L1_SingleJet35_FWD3p0},
    {316, &L1_SingleJet35er2p5},
    {453, &L1_SingleJet43er2p5_NotBptxOR_3BX},
    {454, &L1_SingleJet46er2p5_NotBptxOR_3BX},
    {310, &L1_SingleJet60},
    {325, &L1_SingleJet60_FWD3p0},
    {317, &L1_SingleJet60er2p5},
    {328, &L1_SingleJet8erHE},
    {311, &L1_SingleJet90},
    {326, &L1_SingleJet90_FWD3p0},
    {318, &L1_SingleJet90er2p5},
    {176, &L1_SingleLooseIsoEG26er1p5},
    {175, &L1_SingleLooseIsoEG26er2p5},
    {177, &L1_SingleLooseIsoEG28_FWD2p5},
    {180, &L1_SingleLooseIsoEG28er1p5},
    {179, &L1_SingleLooseIsoEG28er2p1},
    {178, &L1_SingleLooseIsoEG28er2p5},
    {182, &L1_SingleLooseIsoEG30er1p5},
    {181, &L1_SingleLooseIsoEG30er2p5},
    {6, &L1_SingleMu0_BMTF},
    {5, &L1_SingleMu0_DQ},
    {8, &L1_SingleMu0_EMTF},
    {7, &L1_SingleMu0_OMTF},
    {30, &L1_SingleMu10er1p5},
    {13, &L1_SingleMu12_DQ_BMTF},
    {15, &L1_SingleMu12_DQ_EMTF},
    {14, &L1_SingleMu12_DQ_OMTF},
    {31, &L1_SingleMu12er1p5},
    {32, &L1_SingleMu14er1p5},
    {16, &L1_SingleMu15_DQ},
    {33, &L1_SingleMu16er1p5},
    {17, &L1_SingleMu18},
    {34, &L1_SingleMu18er1p5},
    {18, &L1_SingleMu20},
    {21, &L1_SingleMu22},
    {22, &L1_SingleMu22_BMTF},
    {20, &L1_SingleMu22_DQ},
    {24, &L1_SingleMu22_EMTF},
    {23, &L1_SingleMu22_OMTF},
    {19, &L1_SingleMu22_OQ},
    {25, &L1_SingleMu25},
    {9, &L1_SingleMu3},
    {10, &L1_SingleMu5},
    {26, &L1_SingleMu6er1p5},
    {12, &L1_SingleMu7},
    {11, &L1_SingleMu7_DQ},
    {27, &L1_SingleMu7er1p5},
    {28, &L1_SingleMu8er1p5},
    {29, &L1_SingleMu9er1p5},
    {0, &L1_SingleMuCosmics},
    {1, &L1_SingleMuCosmics_BMTF},
    {3, &L1_SingleMuCosmics_EMTF},
    {2, &L1_SingleMuCosmics_OMTF},
    {4, &L1_SingleMuOpen},
    {446, &L1_SingleMuOpen_NotBptxOR},
    {448, &L1_SingleMuOpen_er1p1_NotBptxOR_3BX},
    {447, &L1_SingleMuOpen_er1p4_NotBptxOR_3BX},
    {92, &L1_SingleMuShower_Nominal},
    {93, &L1_SingleMuShower_Tight},
    {264, &L1_SingleTau120er2p1},
    {265, &L1_SingleTau130er2p1},
    {263, &L1_SingleTau70er2p1},
    {503, &L1_TOTEM_1},
    {504, &L1_TOTEM_2},
    {505, &L1_TOTEM_3},
    {506, &L1_TOTEM_4},
    {236, &L1_TripleEG16er2p5},
    {232, &L1_TripleEG_16_12_8_er2p5},
    {233, &L1_TripleEG_16_15_8_er2p5},
    {234, &L1_TripleEG_18_17_8_er2p5},
    {235, &L1_TripleEG_18_18_12_er2p5},
    {373, &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5},
    {374, &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5},
    {372, &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5},
    {72, &L1_TripleMu0},
    {71, &L1_TripleMu0_OQ},
    {73, &L1_TripleMu0_SQ},
    {75, &L1_TripleMu3},
    {76, &L1_TripleMu3_SQ},
    {74, &L1_TripleMu_2SQ_1p5SQ_0OQ},
    {82, &L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12},
    {83, &L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12},
    {77, &L1_TripleMu_5SQ_3SQ_0OQ},
    {87, &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9},
    {88, &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9},
    {79, &L1_TripleMu_5_3_3},
    {80, &L1_TripleMu_5_3_3_SQ},
    {78, &L1_TripleMu_5_3p5_2p5},
    {85, &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17},
    {84, &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17},
    {86, &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17},
    {81, &L1_TripleMu_5_5_3},
    {469, &L1_UnpairedBunchBptxMinus},
    {468, &L1_UnpairedBunchBptxPlus},
    {459, &L1_ZeroBias},
    {460, &L1_ZeroBias_copy}
  };

  const auto rc = id2func.find(index);
  if (rc == id2func.end())
  {
    std::ostringstream oss;
    oss << "no such algorithm index: " << index << ", in menu: L1Menu_Collisions2022_v1_4_0\n";
    throw std::runtime_error(oss.str());
  }
  return rc->second;
}


AlgorithmFunction getFuncFromName(const std::string& name)
{
  static const std::map<std::string, AlgorithmFunction> name2func = {
    {"L1_AlwaysTrue", &L1_AlwaysTrue},
    {"L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME},
    {"L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME},
    {"L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME},
    {"L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME},
    {"L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME},
    {"L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME},
    {"L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME},
    {"L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME},
    {"L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME},
    {"L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME},
    {"L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME},
    {"L1_BptxMinus", &L1_BptxMinus},
    {"L1_BptxOR", &L1_BptxOR},
    {"L1_BptxPlus", &L1_BptxPlus},
    {"L1_BptxXOR", &L1_BptxXOR},
    {"L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142},
    {"L1_DoubleEG10_er1p2_dR_Max0p6", &L1_DoubleEG10_er1p2_dR_Max0p6},
    {"L1_DoubleEG10p5_er1p2_dR_Max0p6", &L1_DoubleEG10p5_er1p2_dR_Max0p6},
    {"L1_DoubleEG11_er1p2_dR_Max0p6", &L1_DoubleEG11_er1p2_dR_Max0p6},
    {"L1_DoubleEG4_er1p2_dR_Max0p9", &L1_DoubleEG4_er1p2_dR_Max0p9},
    {"L1_DoubleEG4p5_er1p2_dR_Max0p9", &L1_DoubleEG4p5_er1p2_dR_Max0p9},
    {"L1_DoubleEG5_er1p2_dR_Max0p9", &L1_DoubleEG5_er1p2_dR_Max0p9},
    {"L1_DoubleEG5p5_er1p2_dR_Max0p8", &L1_DoubleEG5p5_er1p2_dR_Max0p8},
    {"L1_DoubleEG6_er1p2_dR_Max0p8", &L1_DoubleEG6_er1p2_dR_Max0p8},
    {"L1_DoubleEG6p5_er1p2_dR_Max0p8", &L1_DoubleEG6p5_er1p2_dR_Max0p8},
    {"L1_DoubleEG7_er1p2_dR_Max0p8", &L1_DoubleEG7_er1p2_dR_Max0p8},
    {"L1_DoubleEG7p5_er1p2_dR_Max0p7", &L1_DoubleEG7p5_er1p2_dR_Max0p7},
    {"L1_DoubleEG8_er1p2_dR_Max0p7", &L1_DoubleEG8_er1p2_dR_Max0p7},
    {"L1_DoubleEG8er2p5_HTT260er", &L1_DoubleEG8er2p5_HTT260er},
    {"L1_DoubleEG8er2p5_HTT280er", &L1_DoubleEG8er2p5_HTT280er},
    {"L1_DoubleEG8er2p5_HTT300er", &L1_DoubleEG8er2p5_HTT300er},
    {"L1_DoubleEG8er2p5_HTT320er", &L1_DoubleEG8er2p5_HTT320er},
    {"L1_DoubleEG8er2p5_HTT340er", &L1_DoubleEG8er2p5_HTT340er},
    {"L1_DoubleEG8p5_er1p2_dR_Max0p7", &L1_DoubleEG8p5_er1p2_dR_Max0p7},
    {"L1_DoubleEG9_er1p2_dR_Max0p7", &L1_DoubleEG9_er1p2_dR_Max0p7},
    {"L1_DoubleEG9p5_er1p2_dR_Max0p6", &L1_DoubleEG9p5_er1p2_dR_Max0p6},
    {"L1_DoubleEG_15_10_er2p5", &L1_DoubleEG_15_10_er2p5},
    {"L1_DoubleEG_20_10_er2p5", &L1_DoubleEG_20_10_er2p5},
    {"L1_DoubleEG_22_10_er2p5", &L1_DoubleEG_22_10_er2p5},
    {"L1_DoubleEG_25_12_er2p5", &L1_DoubleEG_25_12_er2p5},
    {"L1_DoubleEG_25_14_er2p5", &L1_DoubleEG_25_14_er2p5},
    {"L1_DoubleEG_27_14_er2p5", &L1_DoubleEG_27_14_er2p5},
    {"L1_DoubleEG_LooseIso16_LooseIso12_er1p5", &L1_DoubleEG_LooseIso16_LooseIso12_er1p5},
    {"L1_DoubleEG_LooseIso18_LooseIso12_er1p5", &L1_DoubleEG_LooseIso18_LooseIso12_er1p5},
    {"L1_DoubleEG_LooseIso20_10_er2p5", &L1_DoubleEG_LooseIso20_10_er2p5},
    {"L1_DoubleEG_LooseIso20_LooseIso12_er1p5", &L1_DoubleEG_LooseIso20_LooseIso12_er1p5},
    {"L1_DoubleEG_LooseIso22_10_er2p5", &L1_DoubleEG_LooseIso22_10_er2p5},
    {"L1_DoubleEG_LooseIso22_12_er2p5", &L1_DoubleEG_LooseIso22_12_er2p5},
    {"L1_DoubleEG_LooseIso22_LooseIso12_er1p5", &L1_DoubleEG_LooseIso22_LooseIso12_er1p5},
    {"L1_DoubleEG_LooseIso25_12_er2p5", &L1_DoubleEG_LooseIso25_12_er2p5},
    {"L1_DoubleEG_LooseIso25_LooseIso12_er1p5", &L1_DoubleEG_LooseIso25_LooseIso12_er1p5},
    {"L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5", &L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5},
    {"L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5", &L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5},
    {"L1_DoubleIsoTau28er2p1", &L1_DoubleIsoTau28er2p1},
    {"L1_DoubleIsoTau28er2p1_Mass_Max80", &L1_DoubleIsoTau28er2p1_Mass_Max80},
    {"L1_DoubleIsoTau28er2p1_Mass_Max90", &L1_DoubleIsoTau28er2p1_Mass_Max90},
    {"L1_DoubleIsoTau30er2p1", &L1_DoubleIsoTau30er2p1},
    {"L1_DoubleIsoTau30er2p1_Mass_Max80", &L1_DoubleIsoTau30er2p1_Mass_Max80},
    {"L1_DoubleIsoTau30er2p1_Mass_Max90", &L1_DoubleIsoTau30er2p1_Mass_Max90},
    {"L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1},
    {"L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1},
    {"L1_DoubleIsoTau35er2p1", &L1_DoubleIsoTau35er2p1},
    {"L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1},
    {"L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6},
    {"L1_DoubleJet100er2p5", &L1_DoubleJet100er2p5},
    {"L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6},
    {"L1_DoubleJet120er2p5", &L1_DoubleJet120er2p5},
    {"L1_DoubleJet150er2p5", &L1_DoubleJet150er2p5},
    {"L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5},
    {"L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp},
    {"L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5", &L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5},
    {"L1_DoubleJet40er2p5", &L1_DoubleJet40er2p5},
    {"L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620},
    {"L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620},
    {"L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620},
    {"L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28},
    {"L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620},
    {"L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28},
    {"L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ},
    {"L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp},
    {"L1_DoubleJet_80_30_Mass_Min420_Mu8", &L1_DoubleJet_80_30_Mass_Min420_Mu8},
    {"L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620},
    {"L1_DoubleLLPJet40", &L1_DoubleLLPJet40},
    {"L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1},
    {"L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1},
    {"L1_DoubleMu0", &L1_DoubleMu0},
    {"L1_DoubleMu0_Mass_Min1", &L1_DoubleMu0_Mass_Min1},
    {"L1_DoubleMu0_OQ", &L1_DoubleMu0_OQ},
    {"L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ},
    {"L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS},
    {"L1_DoubleMu0_Upt15_Upt7", &L1_DoubleMu0_Upt15_Upt7},
    {"L1_DoubleMu0_Upt5_Upt5", &L1_DoubleMu0_Upt5_Upt5},
    {"L1_DoubleMu0_Upt6_IP_Min1_Upt4", &L1_DoubleMu0_Upt6_IP_Min1_Upt4},
    {"L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8},
    {"L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", &L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6},
    {"L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4},
    {"L1_DoubleMu0er1p5_SQ", &L1_DoubleMu0er1p5_SQ},
    {"L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS},
    {"L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4},
    {"L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4},
    {"L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5},
    {"L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6},
    {"L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4},
    {"L1_DoubleMu0er2p0_SQ_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_dEta_Max1p5},
    {"L1_DoubleMu0er2p0_SQ_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_dEta_Max1p6},
    {"L1_DoubleMu0er2p0_SQ_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_dR_Max1p4},
    {"L1_DoubleMu18er2p1_SQ", &L1_DoubleMu18er2p1_SQ},
    {"L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20", &L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20},
    {"L1_DoubleMu3_SQ_ETMHF30_HTT60er", &L1_DoubleMu3_SQ_ETMHF30_HTT60er},
    {"L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5},
    {"L1_DoubleMu3_SQ_ETMHF40_HTT60er", &L1_DoubleMu3_SQ_ETMHF40_HTT60er},
    {"L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5},
    {"L1_DoubleMu3_SQ_ETMHF50_HTT60er", &L1_DoubleMu3_SQ_ETMHF50_HTT60er},
    {"L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5},
    {"L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5},
    {"L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5},
    {"L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er},
    {"L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er},
    {"L1_DoubleMu3_SQ_HTT260er", &L1_DoubleMu3_SQ_HTT260er},
    {"L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8},
    {"L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4},
    {"L1_DoubleMu4_SQ_EG9er2p5", &L1_DoubleMu4_SQ_EG9er2p5},
    {"L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS},
    {"L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2},
    {"L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS},
    {"L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2},
    {"L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS},
    {"L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18},
    {"L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7},
    {"L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20", &L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20},
    {"L1_DoubleMu5_SQ_EG9er2p5", &L1_DoubleMu5_SQ_EG9er2p5},
    {"L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ},
    {"L1_DoubleMu9_SQ", &L1_DoubleMu9_SQ},
    {"L1_DoubleMu_12_5", &L1_DoubleMu_12_5},
    {"L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ},
    {"L1_DoubleMu_15_7", &L1_DoubleMu_15_7},
    {"L1_DoubleMu_15_7_Mass_Min1", &L1_DoubleMu_15_7_Mass_Min1},
    {"L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ},
    {"L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1},
    {"L1_ETM120", &L1_ETM120},
    {"L1_ETM150", &L1_ETM150},
    {"L1_ETMHF100", &L1_ETMHF100},
    {"L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er},
    {"L1_ETMHF110", &L1_ETMHF110},
    {"L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er},
    {"L1_ETMHF110_HTT60er_NotSecondBunchInTrain", &L1_ETMHF110_HTT60er_NotSecondBunchInTrain},
    {"L1_ETMHF120", &L1_ETMHF120},
    {"L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er},
    {"L1_ETMHF120_NotSecondBunchInTrain", &L1_ETMHF120_NotSecondBunchInTrain},
    {"L1_ETMHF130", &L1_ETMHF130},
    {"L1_ETMHF130_HTT60er", &L1_ETMHF130_HTT60er},
    {"L1_ETMHF140", &L1_ETMHF140},
    {"L1_ETMHF150", &L1_ETMHF150},
    {"L1_ETMHF70", &L1_ETMHF70},
    {"L1_ETMHF70_HTT60er", &L1_ETMHF70_HTT60er},
    {"L1_ETMHF80", &L1_ETMHF80},
    {"L1_ETMHF80_HTT60er", &L1_ETMHF80_HTT60er},
    {"L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1},
    {"L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6},
    {"L1_ETMHF90", &L1_ETMHF90},
    {"L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er},
    {"L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1},
    {"L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6},
    {"L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1", &L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1},
    {"L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6", &L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6},
    {"L1_ETT1200", &L1_ETT1200},
    {"L1_ETT1600", &L1_ETT1600},
    {"L1_ETT2000", &L1_ETT2000},
    {"L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain},
    {"L1_FirstBunchBeforeTrain", &L1_FirstBunchBeforeTrain},
    {"L1_FirstBunchInTrain", &L1_FirstBunchInTrain},
    {"L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit},
    {"L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain},
    {"L1_HCAL_LaserMon_Trig", &L1_HCAL_LaserMon_Trig},
    {"L1_HCAL_LaserMon_Veto", &L1_HCAL_LaserMon_Veto},
    {"L1_HTT120_SingleLLPJet40", &L1_HTT120_SingleLLPJet40},
    {"L1_HTT120er", &L1_HTT120er},
    {"L1_HTT160_SingleLLPJet50", &L1_HTT160_SingleLLPJet50},
    {"L1_HTT160er", &L1_HTT160er},
    {"L1_HTT200_SingleLLPJet60", &L1_HTT200_SingleLLPJet60},
    {"L1_HTT200er", &L1_HTT200er},
    {"L1_HTT240_SingleLLPJet70", &L1_HTT240_SingleLLPJet70},
    {"L1_HTT255er", &L1_HTT255er},
    {"L1_HTT280er", &L1_HTT280er},
    {"L1_HTT280er_QuadJet_70_55_40_35_er2p5", &L1_HTT280er_QuadJet_70_55_40_35_er2p5},
    {"L1_HTT320er", &L1_HTT320er},
    {"L1_HTT320er_QuadJet_70_55_40_40_er2p5", &L1_HTT320er_QuadJet_70_55_40_40_er2p5},
    {"L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3},
    {"L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3},
    {"L1_HTT360er", &L1_HTT360er},
    {"L1_HTT400er", &L1_HTT400er},
    {"L1_HTT450er", &L1_HTT450er},
    {"L1_IsoEG32er2p5_Mt40", &L1_IsoEG32er2p5_Mt40},
    {"L1_IsoTau52er2p1_QuadJet36er2p5", &L1_IsoTau52er2p1_QuadJet36er2p5},
    {"L1_IsolatedBunch", &L1_IsolatedBunch},
    {"L1_LastBunchInTrain", &L1_LastBunchInTrain},
    {"L1_LastCollisionInTrain", &L1_LastCollisionInTrain},
    {"L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3},
    {"L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3},
    {"L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er},
    {"L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3},
    {"L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er},
    {"L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3},
    {"L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er},
    {"L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3},
    {"L1_LooseIsoEG30er2p1_HTT100er", &L1_LooseIsoEG30er2p1_HTT100er},
    {"L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3},
    {"L1_MinimumBiasHF0", &L1_MinimumBiasHF0},
    {"L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND},
    {"L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6},
    {"L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6},
    {"L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6},
    {"L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1},
    {"L1_Mu18er2p1_Tau26er2p1", &L1_Mu18er2p1_Tau26er2p1},
    {"L1_Mu18er2p1_Tau26er2p1_Jet55", &L1_Mu18er2p1_Tau26er2p1_Jet55},
    {"L1_Mu18er2p1_Tau26er2p1_Jet70", &L1_Mu18er2p1_Tau26er2p1_Jet70},
    {"L1_Mu20_EG10er2p5", &L1_Mu20_EG10er2p5},
    {"L1_Mu22er2p1_IsoTau28er2p1", &L1_Mu22er2p1_IsoTau28er2p1},
    {"L1_Mu22er2p1_IsoTau30er2p1", &L1_Mu22er2p1_IsoTau30er2p1},
    {"L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1},
    {"L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1},
    {"L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1},
    {"L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1},
    {"L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1},
    {"L1_Mu3_Jet120er2p5_dR_Max0p4", &L1_Mu3_Jet120er2p5_dR_Max0p4},
    {"L1_Mu3_Jet120er2p5_dR_Max0p8", &L1_Mu3_Jet120er2p5_dR_Max0p8},
    {"L1_Mu3_Jet16er2p5_dR_Max0p4", &L1_Mu3_Jet16er2p5_dR_Max0p4},
    {"L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5},
    {"L1_Mu3_Jet35er2p5_dR_Max0p4", &L1_Mu3_Jet35er2p5_dR_Max0p4},
    {"L1_Mu3_Jet60er2p5_dR_Max0p4", &L1_Mu3_Jet60er2p5_dR_Max0p4},
    {"L1_Mu3_Jet80er2p5_dR_Max0p4", &L1_Mu3_Jet80er2p5_dR_Max0p4},
    {"L1_Mu3er1p5_Jet100er2p5_ETMHF30", &L1_Mu3er1p5_Jet100er2p5_ETMHF30},
    {"L1_Mu3er1p5_Jet100er2p5_ETMHF40", &L1_Mu3er1p5_Jet100er2p5_ETMHF40},
    {"L1_Mu3er1p5_Jet100er2p5_ETMHF50", &L1_Mu3er1p5_Jet100er2p5_ETMHF50},
    {"L1_Mu5_EG23er2p5", &L1_Mu5_EG23er2p5},
    {"L1_Mu5_LooseIsoEG20er2p5", &L1_Mu5_LooseIsoEG20er2p5},
    {"L1_Mu6_DoubleEG10er2p5", &L1_Mu6_DoubleEG10er2p5},
    {"L1_Mu6_DoubleEG12er2p5", &L1_Mu6_DoubleEG12er2p5},
    {"L1_Mu6_DoubleEG15er2p5", &L1_Mu6_DoubleEG15er2p5},
    {"L1_Mu6_DoubleEG17er2p5", &L1_Mu6_DoubleEG17er2p5},
    {"L1_Mu6_HTT240er", &L1_Mu6_HTT240er},
    {"L1_Mu6_HTT250er", &L1_Mu6_HTT250er},
    {"L1_Mu7_EG20er2p5", &L1_Mu7_EG20er2p5},
    {"L1_Mu7_EG23er2p5", &L1_Mu7_EG23er2p5},
    {"L1_Mu7_LooseIsoEG20er2p5", &L1_Mu7_LooseIsoEG20er2p5},
    {"L1_Mu7_LooseIsoEG23er2p5", &L1_Mu7_LooseIsoEG23er2p5},
    {"L1_NotBptxOR", &L1_NotBptxOR},
    {"L1_QuadJet60er2p5", &L1_QuadJet60er2p5},
    {"L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0},
    {"L1_QuadMu0", &L1_QuadMu0},
    {"L1_QuadMu0_OQ", &L1_QuadMu0_OQ},
    {"L1_QuadMu0_SQ", &L1_QuadMu0_SQ},
    {"L1_SecondBunchInTrain", &L1_SecondBunchInTrain},
    {"L1_SecondLastBunchInTrain", &L1_SecondLastBunchInTrain},
    {"L1_SingleEG10er2p5", &L1_SingleEG10er2p5},
    {"L1_SingleEG15er2p5", &L1_SingleEG15er2p5},
    {"L1_SingleEG26er2p5", &L1_SingleEG26er2p5},
    {"L1_SingleEG28_FWD2p5", &L1_SingleEG28_FWD2p5},
    {"L1_SingleEG28er1p5", &L1_SingleEG28er1p5},
    {"L1_SingleEG28er2p1", &L1_SingleEG28er2p1},
    {"L1_SingleEG28er2p5", &L1_SingleEG28er2p5},
    {"L1_SingleEG34er2p5", &L1_SingleEG34er2p5},
    {"L1_SingleEG36er2p5", &L1_SingleEG36er2p5},
    {"L1_SingleEG38er2p5", &L1_SingleEG38er2p5},
    {"L1_SingleEG40er2p5", &L1_SingleEG40er2p5},
    {"L1_SingleEG42er2p5", &L1_SingleEG42er2p5},
    {"L1_SingleEG45er2p5", &L1_SingleEG45er2p5},
    {"L1_SingleEG50", &L1_SingleEG50},
    {"L1_SingleEG60", &L1_SingleEG60},
    {"L1_SingleEG8er2p5", &L1_SingleEG8er2p5},
    {"L1_SingleIsoEG24er1p5", &L1_SingleIsoEG24er1p5},
    {"L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1},
    {"L1_SingleIsoEG26er1p5", &L1_SingleIsoEG26er1p5},
    {"L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1},
    {"L1_SingleIsoEG26er2p5", &L1_SingleIsoEG26er2p5},
    {"L1_SingleIsoEG28_FWD2p5", &L1_SingleIsoEG28_FWD2p5},
    {"L1_SingleIsoEG28er1p5", &L1_SingleIsoEG28er1p5},
    {"L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1},
    {"L1_SingleIsoEG28er2p5", &L1_SingleIsoEG28er2p5},
    {"L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1},
    {"L1_SingleIsoEG30er2p5", &L1_SingleIsoEG30er2p5},
    {"L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1},
    {"L1_SingleIsoEG32er2p5", &L1_SingleIsoEG32er2p5},
    {"L1_SingleIsoEG34er2p5", &L1_SingleIsoEG34er2p5},
    {"L1_SingleIsoTau32er2p1", &L1_SingleIsoTau32er2p1},
    {"L1_SingleJet10erHE", &L1_SingleJet10erHE},
    {"L1_SingleJet120", &L1_SingleJet120},
    {"L1_SingleJet120_FWD3p0", &L1_SingleJet120_FWD3p0},
    {"L1_SingleJet120er2p5", &L1_SingleJet120er2p5},
    {"L1_SingleJet12erHE", &L1_SingleJet12erHE},
    {"L1_SingleJet140er2p5", &L1_SingleJet140er2p5},
    {"L1_SingleJet140er2p5_ETMHF70", &L1_SingleJet140er2p5_ETMHF70},
    {"L1_SingleJet140er2p5_ETMHF80", &L1_SingleJet140er2p5_ETMHF80},
    {"L1_SingleJet140er2p5_ETMHF90", &L1_SingleJet140er2p5_ETMHF90},
    {"L1_SingleJet160er2p5", &L1_SingleJet160er2p5},
    {"L1_SingleJet180", &L1_SingleJet180},
    {"L1_SingleJet180er2p5", &L1_SingleJet180er2p5},
    {"L1_SingleJet200", &L1_SingleJet200},
    {"L1_SingleJet20er2p5_NotBptxOR", &L1_SingleJet20er2p5_NotBptxOR},
    {"L1_SingleJet20er2p5_NotBptxOR_3BX", &L1_SingleJet20er2p5_NotBptxOR_3BX},
    {"L1_SingleJet35", &L1_SingleJet35},
    {"L1_SingleJet35_FWD3p0", &L1_SingleJet35_FWD3p0},
    {"L1_SingleJet35er2p5", &L1_SingleJet35er2p5},
    {"L1_SingleJet43er2p5_NotBptxOR_3BX", &L1_SingleJet43er2p5_NotBptxOR_3BX},
    {"L1_SingleJet46er2p5_NotBptxOR_3BX", &L1_SingleJet46er2p5_NotBptxOR_3BX},
    {"L1_SingleJet60", &L1_SingleJet60},
    {"L1_SingleJet60_FWD3p0", &L1_SingleJet60_FWD3p0},
    {"L1_SingleJet60er2p5", &L1_SingleJet60er2p5},
    {"L1_SingleJet8erHE", &L1_SingleJet8erHE},
    {"L1_SingleJet90", &L1_SingleJet90},
    {"L1_SingleJet90_FWD3p0", &L1_SingleJet90_FWD3p0},
    {"L1_SingleJet90er2p5", &L1_SingleJet90er2p5},
    {"L1_SingleLooseIsoEG26er1p5", &L1_SingleLooseIsoEG26er1p5},
    {"L1_SingleLooseIsoEG26er2p5", &L1_SingleLooseIsoEG26er2p5},
    {"L1_SingleLooseIsoEG28_FWD2p5", &L1_SingleLooseIsoEG28_FWD2p5},
    {"L1_SingleLooseIsoEG28er1p5", &L1_SingleLooseIsoEG28er1p5},
    {"L1_SingleLooseIsoEG28er2p1", &L1_SingleLooseIsoEG28er2p1},
    {"L1_SingleLooseIsoEG28er2p5", &L1_SingleLooseIsoEG28er2p5},
    {"L1_SingleLooseIsoEG30er1p5", &L1_SingleLooseIsoEG30er1p5},
    {"L1_SingleLooseIsoEG30er2p5", &L1_SingleLooseIsoEG30er2p5},
    {"L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF},
    {"L1_SingleMu0_DQ", &L1_SingleMu0_DQ},
    {"L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF},
    {"L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF},
    {"L1_SingleMu10er1p5", &L1_SingleMu10er1p5},
    {"L1_SingleMu12_DQ_BMTF", &L1_SingleMu12_DQ_BMTF},
    {"L1_SingleMu12_DQ_EMTF", &L1_SingleMu12_DQ_EMTF},
    {"L1_SingleMu12_DQ_OMTF", &L1_SingleMu12_DQ_OMTF},
    {"L1_SingleMu12er1p5", &L1_SingleMu12er1p5},
    {"L1_SingleMu14er1p5", &L1_SingleMu14er1p5},
    {"L1_SingleMu15_DQ", &L1_SingleMu15_DQ},
    {"L1_SingleMu16er1p5", &L1_SingleMu16er1p5},
    {"L1_SingleMu18", &L1_SingleMu18},
    {"L1_SingleMu18er1p5", &L1_SingleMu18er1p5},
    {"L1_SingleMu20", &L1_SingleMu20},
    {"L1_SingleMu22", &L1_SingleMu22},
    {"L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF},
    {"L1_SingleMu22_DQ", &L1_SingleMu22_DQ},
    {"L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF},
    {"L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF},
    {"L1_SingleMu22_OQ", &L1_SingleMu22_OQ},
    {"L1_SingleMu25", &L1_SingleMu25},
    {"L1_SingleMu3", &L1_SingleMu3},
    {"L1_SingleMu5", &L1_SingleMu5},
    {"L1_SingleMu6er1p5", &L1_SingleMu6er1p5},
    {"L1_SingleMu7", &L1_SingleMu7},
    {"L1_SingleMu7_DQ", &L1_SingleMu7_DQ},
    {"L1_SingleMu7er1p5", &L1_SingleMu7er1p5},
    {"L1_SingleMu8er1p5", &L1_SingleMu8er1p5},
    {"L1_SingleMu9er1p5", &L1_SingleMu9er1p5},
    {"L1_SingleMuCosmics", &L1_SingleMuCosmics},
    {"L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF},
    {"L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF},
    {"L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF},
    {"L1_SingleMuOpen", &L1_SingleMuOpen},
    {"L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR},
    {"L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX},
    {"L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX},
    {"L1_SingleMuShower_Nominal", &L1_SingleMuShower_Nominal},
    {"L1_SingleMuShower_Tight", &L1_SingleMuShower_Tight},
    {"L1_SingleTau120er2p1", &L1_SingleTau120er2p1},
    {"L1_SingleTau130er2p1", &L1_SingleTau130er2p1},
    {"L1_SingleTau70er2p1", &L1_SingleTau70er2p1},
    {"L1_TOTEM_1", &L1_TOTEM_1},
    {"L1_TOTEM_2", &L1_TOTEM_2},
    {"L1_TOTEM_3", &L1_TOTEM_3},
    {"L1_TOTEM_4", &L1_TOTEM_4},
    {"L1_TripleEG16er2p5", &L1_TripleEG16er2p5},
    {"L1_TripleEG_16_12_8_er2p5", &L1_TripleEG_16_12_8_er2p5},
    {"L1_TripleEG_16_15_8_er2p5", &L1_TripleEG_16_15_8_er2p5},
    {"L1_TripleEG_18_17_8_er2p5", &L1_TripleEG_18_17_8_er2p5},
    {"L1_TripleEG_18_18_12_er2p5", &L1_TripleEG_18_18_12_er2p5},
    {"L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5},
    {"L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5},
    {"L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5},
    {"L1_TripleMu0", &L1_TripleMu0},
    {"L1_TripleMu0_OQ", &L1_TripleMu0_OQ},
    {"L1_TripleMu0_SQ", &L1_TripleMu0_SQ},
    {"L1_TripleMu3", &L1_TripleMu3},
    {"L1_TripleMu3_SQ", &L1_TripleMu3_SQ},
    {"L1_TripleMu_2SQ_1p5SQ_0OQ", &L1_TripleMu_2SQ_1p5SQ_0OQ},
    {"L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12", &L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12},
    {"L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12", &L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12},
    {"L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ},
    {"L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9},
    {"L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9},
    {"L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3},
    {"L1_TripleMu_5_3_3_SQ", &L1_TripleMu_5_3_3_SQ},
    {"L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5},
    {"L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17},
    {"L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17},
    {"L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17},
    {"L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3},
    {"L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus},
    {"L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus},
    {"L1_ZeroBias", &L1_ZeroBias},
    {"L1_ZeroBias_copy", &L1_ZeroBias_copy}
  };

  const auto rc = name2func.find(name);
  if (rc == name2func.end())
  {
    std::ostringstream oss;
    oss << "no such algorithm name: \"" << name << "\", in menu: L1Menu_Collisions2022_v1_4_0\n";
    throw std::runtime_error(oss.str());
  }
  return rc->second;
}


bool addFuncFromName(std::map<std::string, std::function<bool()>> &L1SeedFun,
                     L1Analysis::L1AnalysisL1UpgradeDataFormat* upgrade,
                     L1Analysis::L1AnalysisL1CaloTowerDataFormat* calo_tower)
{
  static const std::map<std::string, AlgorithmFunction> name2func = {
    {"L1_AlwaysTrue", &L1_AlwaysTrue},
    {"L1_BPTX_AND_Ref1_VME", &L1_BPTX_AND_Ref1_VME},
    {"L1_BPTX_AND_Ref3_VME", &L1_BPTX_AND_Ref3_VME},
    {"L1_BPTX_AND_Ref4_VME", &L1_BPTX_AND_Ref4_VME},
    {"L1_BPTX_BeamGas_B1_VME", &L1_BPTX_BeamGas_B1_VME},
    {"L1_BPTX_BeamGas_B2_VME", &L1_BPTX_BeamGas_B2_VME},
    {"L1_BPTX_BeamGas_Ref1_VME", &L1_BPTX_BeamGas_Ref1_VME},
    {"L1_BPTX_BeamGas_Ref2_VME", &L1_BPTX_BeamGas_Ref2_VME},
    {"L1_BPTX_NotOR_VME", &L1_BPTX_NotOR_VME},
    {"L1_BPTX_OR_Ref3_VME", &L1_BPTX_OR_Ref3_VME},
    {"L1_BPTX_OR_Ref4_VME", &L1_BPTX_OR_Ref4_VME},
    {"L1_BPTX_RefAND_VME", &L1_BPTX_RefAND_VME},
    {"L1_BptxMinus", &L1_BptxMinus},
    {"L1_BptxOR", &L1_BptxOR},
    {"L1_BptxPlus", &L1_BptxPlus},
    {"L1_BptxXOR", &L1_BptxXOR},
    {"L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142", &L1_CDC_SingleMu_3_er1p2_TOP120_DPHI2p618_3p142},
    {"L1_DoubleEG10_er1p2_dR_Max0p6", &L1_DoubleEG10_er1p2_dR_Max0p6},
    {"L1_DoubleEG10p5_er1p2_dR_Max0p6", &L1_DoubleEG10p5_er1p2_dR_Max0p6},
    {"L1_DoubleEG11_er1p2_dR_Max0p6", &L1_DoubleEG11_er1p2_dR_Max0p6},
    {"L1_DoubleEG4_er1p2_dR_Max0p9", &L1_DoubleEG4_er1p2_dR_Max0p9},
    {"L1_DoubleEG4p5_er1p2_dR_Max0p9", &L1_DoubleEG4p5_er1p2_dR_Max0p9},
    {"L1_DoubleEG5_er1p2_dR_Max0p9", &L1_DoubleEG5_er1p2_dR_Max0p9},
    {"L1_DoubleEG5p5_er1p2_dR_Max0p8", &L1_DoubleEG5p5_er1p2_dR_Max0p8},
    {"L1_DoubleEG6_er1p2_dR_Max0p8", &L1_DoubleEG6_er1p2_dR_Max0p8},
    {"L1_DoubleEG6p5_er1p2_dR_Max0p8", &L1_DoubleEG6p5_er1p2_dR_Max0p8},
    {"L1_DoubleEG7_er1p2_dR_Max0p8", &L1_DoubleEG7_er1p2_dR_Max0p8},
    {"L1_DoubleEG7p5_er1p2_dR_Max0p7", &L1_DoubleEG7p5_er1p2_dR_Max0p7},
    {"L1_DoubleEG8_er1p2_dR_Max0p7", &L1_DoubleEG8_er1p2_dR_Max0p7},
    {"L1_DoubleEG8er2p5_HTT260er", &L1_DoubleEG8er2p5_HTT260er},
    {"L1_DoubleEG8er2p5_HTT280er", &L1_DoubleEG8er2p5_HTT280er},
    {"L1_DoubleEG8er2p5_HTT300er", &L1_DoubleEG8er2p5_HTT300er},
    {"L1_DoubleEG8er2p5_HTT320er", &L1_DoubleEG8er2p5_HTT320er},
    {"L1_DoubleEG8er2p5_HTT340er", &L1_DoubleEG8er2p5_HTT340er},
    {"L1_DoubleEG8p5_er1p2_dR_Max0p7", &L1_DoubleEG8p5_er1p2_dR_Max0p7},
    {"L1_DoubleEG9_er1p2_dR_Max0p7", &L1_DoubleEG9_er1p2_dR_Max0p7},
    {"L1_DoubleEG9p5_er1p2_dR_Max0p6", &L1_DoubleEG9p5_er1p2_dR_Max0p6},
    {"L1_DoubleEG_15_10_er2p5", &L1_DoubleEG_15_10_er2p5},
    {"L1_DoubleEG_20_10_er2p5", &L1_DoubleEG_20_10_er2p5},
    {"L1_DoubleEG_22_10_er2p5", &L1_DoubleEG_22_10_er2p5},
    {"L1_DoubleEG_25_12_er2p5", &L1_DoubleEG_25_12_er2p5},
    {"L1_DoubleEG_25_14_er2p5", &L1_DoubleEG_25_14_er2p5},
    {"L1_DoubleEG_27_14_er2p5", &L1_DoubleEG_27_14_er2p5},
    {"L1_DoubleEG_LooseIso16_LooseIso12_er1p5", &L1_DoubleEG_LooseIso16_LooseIso12_er1p5},
    {"L1_DoubleEG_LooseIso18_LooseIso12_er1p5", &L1_DoubleEG_LooseIso18_LooseIso12_er1p5},
    {"L1_DoubleEG_LooseIso20_10_er2p5", &L1_DoubleEG_LooseIso20_10_er2p5},
    {"L1_DoubleEG_LooseIso20_LooseIso12_er1p5", &L1_DoubleEG_LooseIso20_LooseIso12_er1p5},
    {"L1_DoubleEG_LooseIso22_10_er2p5", &L1_DoubleEG_LooseIso22_10_er2p5},
    {"L1_DoubleEG_LooseIso22_12_er2p5", &L1_DoubleEG_LooseIso22_12_er2p5},
    {"L1_DoubleEG_LooseIso22_LooseIso12_er1p5", &L1_DoubleEG_LooseIso22_LooseIso12_er1p5},
    {"L1_DoubleEG_LooseIso25_12_er2p5", &L1_DoubleEG_LooseIso25_12_er2p5},
    {"L1_DoubleEG_LooseIso25_LooseIso12_er1p5", &L1_DoubleEG_LooseIso25_LooseIso12_er1p5},
    {"L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5", &L1_DoubleIsoTau26er2p1_Jet55_RmOvlp_dR0p5},
    {"L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5", &L1_DoubleIsoTau26er2p1_Jet70_RmOvlp_dR0p5},
    {"L1_DoubleIsoTau28er2p1", &L1_DoubleIsoTau28er2p1},
    {"L1_DoubleIsoTau28er2p1_Mass_Max80", &L1_DoubleIsoTau28er2p1_Mass_Max80},
    {"L1_DoubleIsoTau28er2p1_Mass_Max90", &L1_DoubleIsoTau28er2p1_Mass_Max90},
    {"L1_DoubleIsoTau30er2p1", &L1_DoubleIsoTau30er2p1},
    {"L1_DoubleIsoTau30er2p1_Mass_Max80", &L1_DoubleIsoTau30er2p1_Mass_Max80},
    {"L1_DoubleIsoTau30er2p1_Mass_Max90", &L1_DoubleIsoTau30er2p1_Mass_Max90},
    {"L1_DoubleIsoTau32er2p1", &L1_DoubleIsoTau32er2p1},
    {"L1_DoubleIsoTau34er2p1", &L1_DoubleIsoTau34er2p1},
    {"L1_DoubleIsoTau35er2p1", &L1_DoubleIsoTau35er2p1},
    {"L1_DoubleIsoTau36er2p1", &L1_DoubleIsoTau36er2p1},
    {"L1_DoubleJet100er2p3_dEta_Max1p6", &L1_DoubleJet100er2p3_dEta_Max1p6},
    {"L1_DoubleJet100er2p5", &L1_DoubleJet100er2p5},
    {"L1_DoubleJet112er2p3_dEta_Max1p6", &L1_DoubleJet112er2p3_dEta_Max1p6},
    {"L1_DoubleJet120er2p5", &L1_DoubleJet120er2p5},
    {"L1_DoubleJet150er2p5", &L1_DoubleJet150er2p5},
    {"L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min150_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min200_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5},
    {"L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", &L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5},
    {"L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp", &L1_DoubleJet35_Mass_Min450_IsoTau45_RmOvlp},
    {"L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5", &L1_DoubleJet35_Mass_Min450_IsoTau45er2p1_RmOvlp_dR0p5},
    {"L1_DoubleJet40er2p5", &L1_DoubleJet40er2p5},
    {"L1_DoubleJet_100_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_100_30_DoubleJet30_Mass_Min620},
    {"L1_DoubleJet_110_35_DoubleJet35_Mass_Min620", &L1_DoubleJet_110_35_DoubleJet35_Mass_Min620},
    {"L1_DoubleJet_115_40_DoubleJet40_Mass_Min620", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620},
    {"L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28", &L1_DoubleJet_115_40_DoubleJet40_Mass_Min620_Jet60TT28},
    {"L1_DoubleJet_120_45_DoubleJet45_Mass_Min620", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620},
    {"L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28", &L1_DoubleJet_120_45_DoubleJet45_Mass_Min620_Jet60TT28},
    {"L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ", &L1_DoubleJet_80_30_Mass_Min420_DoubleMu0_SQ},
    {"L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp", &L1_DoubleJet_80_30_Mass_Min420_IsoTau40_RmOvlp},
    {"L1_DoubleJet_80_30_Mass_Min420_Mu8", &L1_DoubleJet_80_30_Mass_Min420_Mu8},
    {"L1_DoubleJet_90_30_DoubleJet30_Mass_Min620", &L1_DoubleJet_90_30_DoubleJet30_Mass_Min620},
    {"L1_DoubleLLPJet40", &L1_DoubleLLPJet40},
    {"L1_DoubleLooseIsoEG22er2p1", &L1_DoubleLooseIsoEG22er2p1},
    {"L1_DoubleLooseIsoEG24er2p1", &L1_DoubleLooseIsoEG24er2p1},
    {"L1_DoubleMu0", &L1_DoubleMu0},
    {"L1_DoubleMu0_Mass_Min1", &L1_DoubleMu0_Mass_Min1},
    {"L1_DoubleMu0_OQ", &L1_DoubleMu0_OQ},
    {"L1_DoubleMu0_SQ", &L1_DoubleMu0_SQ},
    {"L1_DoubleMu0_SQ_OS", &L1_DoubleMu0_SQ_OS},
    {"L1_DoubleMu0_Upt15_Upt7", &L1_DoubleMu0_Upt15_Upt7},
    {"L1_DoubleMu0_Upt5_Upt5", &L1_DoubleMu0_Upt5_Upt5},
    {"L1_DoubleMu0_Upt6_IP_Min1_Upt4", &L1_DoubleMu0_Upt6_IP_Min1_Upt4},
    {"L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu0_dR_Max1p6_Jet90er2p5_dR_Max0p8},
    {"L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6", &L1_DoubleMu0er1p4_OQ_OS_dEta_Max1p6},
    {"L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p4_SQ_OS_dR_Max1p4},
    {"L1_DoubleMu0er1p5_SQ", &L1_DoubleMu0er1p5_SQ},
    {"L1_DoubleMu0er1p5_SQ_OS", &L1_DoubleMu0er1p5_SQ_OS},
    {"L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_OS_dR_Max1p4},
    {"L1_DoubleMu0er1p5_SQ_dR_Max1p4", &L1_DoubleMu0er1p5_SQ_dR_Max1p4},
    {"L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p5},
    {"L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_OS_dEta_Max1p6},
    {"L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_OS_dR_Max1p4},
    {"L1_DoubleMu0er2p0_SQ_dEta_Max1p5", &L1_DoubleMu0er2p0_SQ_dEta_Max1p5},
    {"L1_DoubleMu0er2p0_SQ_dEta_Max1p6", &L1_DoubleMu0er2p0_SQ_dEta_Max1p6},
    {"L1_DoubleMu0er2p0_SQ_dR_Max1p4", &L1_DoubleMu0er2p0_SQ_dR_Max1p4},
    {"L1_DoubleMu18er2p1_SQ", &L1_DoubleMu18er2p1_SQ},
    {"L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20", &L1_DoubleMu3_OS_er2p3_Mass_Max14_DoubleEG7p5_er2p1_Mass_Max20},
    {"L1_DoubleMu3_SQ_ETMHF30_HTT60er", &L1_DoubleMu3_SQ_ETMHF30_HTT60er},
    {"L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF30_Jet60er2p5_OR_DoubleJet40er2p5},
    {"L1_DoubleMu3_SQ_ETMHF40_HTT60er", &L1_DoubleMu3_SQ_ETMHF40_HTT60er},
    {"L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF40_Jet60er2p5_OR_DoubleJet40er2p5},
    {"L1_DoubleMu3_SQ_ETMHF50_HTT60er", &L1_DoubleMu3_SQ_ETMHF50_HTT60er},
    {"L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5},
    {"L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5", &L1_DoubleMu3_SQ_ETMHF50_Jet60er2p5_OR_DoubleJet40er2p5},
    {"L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5", &L1_DoubleMu3_SQ_ETMHF60_Jet60er2p5},
    {"L1_DoubleMu3_SQ_HTT220er", &L1_DoubleMu3_SQ_HTT220er},
    {"L1_DoubleMu3_SQ_HTT240er", &L1_DoubleMu3_SQ_HTT240er},
    {"L1_DoubleMu3_SQ_HTT260er", &L1_DoubleMu3_SQ_HTT260er},
    {"L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8", &L1_DoubleMu3_dR_Max1p6_Jet90er2p5_dR_Max0p8},
    {"L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4", &L1_DoubleMu3er2p0_SQ_OS_dR_Max1p4},
    {"L1_DoubleMu4_SQ_EG9er2p5", &L1_DoubleMu4_SQ_EG9er2p5},
    {"L1_DoubleMu4_SQ_OS", &L1_DoubleMu4_SQ_OS},
    {"L1_DoubleMu4_SQ_OS_dR_Max1p2", &L1_DoubleMu4_SQ_OS_dR_Max1p2},
    {"L1_DoubleMu4p5_SQ_OS", &L1_DoubleMu4p5_SQ_OS},
    {"L1_DoubleMu4p5_SQ_OS_dR_Max1p2", &L1_DoubleMu4p5_SQ_OS_dR_Max1p2},
    {"L1_DoubleMu4p5er2p0_SQ_OS", &L1_DoubleMu4p5er2p0_SQ_OS},
    {"L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_7to18},
    {"L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", &L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7},
    {"L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20", &L1_DoubleMu5_OS_er2p3_Mass_8to14_DoubleEG3er2p1_Mass_Max20},
    {"L1_DoubleMu5_SQ_EG9er2p5", &L1_DoubleMu5_SQ_EG9er2p5},
    {"L1_DoubleMu8_SQ", &L1_DoubleMu8_SQ},
    {"L1_DoubleMu9_SQ", &L1_DoubleMu9_SQ},
    {"L1_DoubleMu_12_5", &L1_DoubleMu_12_5},
    {"L1_DoubleMu_15_5_SQ", &L1_DoubleMu_15_5_SQ},
    {"L1_DoubleMu_15_7", &L1_DoubleMu_15_7},
    {"L1_DoubleMu_15_7_Mass_Min1", &L1_DoubleMu_15_7_Mass_Min1},
    {"L1_DoubleMu_15_7_SQ", &L1_DoubleMu_15_7_SQ},
    {"L1_DoubleTau70er2p1", &L1_DoubleTau70er2p1},
    {"L1_ETM120", &L1_ETM120},
    {"L1_ETM150", &L1_ETM150},
    {"L1_ETMHF100", &L1_ETMHF100},
    {"L1_ETMHF100_HTT60er", &L1_ETMHF100_HTT60er},
    {"L1_ETMHF110", &L1_ETMHF110},
    {"L1_ETMHF110_HTT60er", &L1_ETMHF110_HTT60er},
    {"L1_ETMHF110_HTT60er_NotSecondBunchInTrain", &L1_ETMHF110_HTT60er_NotSecondBunchInTrain},
    {"L1_ETMHF120", &L1_ETMHF120},
    {"L1_ETMHF120_HTT60er", &L1_ETMHF120_HTT60er},
    {"L1_ETMHF120_NotSecondBunchInTrain", &L1_ETMHF120_NotSecondBunchInTrain},
    {"L1_ETMHF130", &L1_ETMHF130},
    {"L1_ETMHF130_HTT60er", &L1_ETMHF130_HTT60er},
    {"L1_ETMHF140", &L1_ETMHF140},
    {"L1_ETMHF150", &L1_ETMHF150},
    {"L1_ETMHF70", &L1_ETMHF70},
    {"L1_ETMHF70_HTT60er", &L1_ETMHF70_HTT60er},
    {"L1_ETMHF80", &L1_ETMHF80},
    {"L1_ETMHF80_HTT60er", &L1_ETMHF80_HTT60er},
    {"L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p1},
    {"L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6", &L1_ETMHF80_SingleJet55er2p5_dPhi_Min2p6},
    {"L1_ETMHF90", &L1_ETMHF90},
    {"L1_ETMHF90_HTT60er", &L1_ETMHF90_HTT60er},
    {"L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p1},
    {"L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6", &L1_ETMHF90_SingleJet60er2p5_dPhi_Min2p6},
    {"L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1", &L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p1},
    {"L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6", &L1_ETMHF90_SingleJet80er2p5_dPhi_Min2p6},
    {"L1_ETT1200", &L1_ETT1200},
    {"L1_ETT1600", &L1_ETT1600},
    {"L1_ETT2000", &L1_ETT2000},
    {"L1_FirstBunchAfterTrain", &L1_FirstBunchAfterTrain},
    {"L1_FirstBunchBeforeTrain", &L1_FirstBunchBeforeTrain},
    {"L1_FirstBunchInTrain", &L1_FirstBunchInTrain},
    {"L1_FirstCollisionInOrbit", &L1_FirstCollisionInOrbit},
    {"L1_FirstCollisionInTrain", &L1_FirstCollisionInTrain},
    {"L1_HCAL_LaserMon_Trig", &L1_HCAL_LaserMon_Trig},
    {"L1_HCAL_LaserMon_Veto", &L1_HCAL_LaserMon_Veto},
    {"L1_HTT120_SingleLLPJet40", &L1_HTT120_SingleLLPJet40},
    {"L1_HTT120er", &L1_HTT120er},
    {"L1_HTT160_SingleLLPJet50", &L1_HTT160_SingleLLPJet50},
    {"L1_HTT160er", &L1_HTT160er},
    {"L1_HTT200_SingleLLPJet60", &L1_HTT200_SingleLLPJet60},
    {"L1_HTT200er", &L1_HTT200er},
    {"L1_HTT240_SingleLLPJet70", &L1_HTT240_SingleLLPJet70},
    {"L1_HTT255er", &L1_HTT255er},
    {"L1_HTT280er", &L1_HTT280er},
    {"L1_HTT280er_QuadJet_70_55_40_35_er2p5", &L1_HTT280er_QuadJet_70_55_40_35_er2p5},
    {"L1_HTT320er", &L1_HTT320er},
    {"L1_HTT320er_QuadJet_70_55_40_40_er2p5", &L1_HTT320er_QuadJet_70_55_40_40_er2p5},
    {"L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_45_40_er2p3},
    {"L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3", &L1_HTT320er_QuadJet_80_60_er2p1_50_45_er2p3},
    {"L1_HTT360er", &L1_HTT360er},
    {"L1_HTT400er", &L1_HTT400er},
    {"L1_HTT450er", &L1_HTT450er},
    {"L1_IsoEG32er2p5_Mt40", &L1_IsoEG32er2p5_Mt40},
    {"L1_IsoTau52er2p1_QuadJet36er2p5", &L1_IsoTau52er2p1_QuadJet36er2p5},
    {"L1_IsolatedBunch", &L1_IsolatedBunch},
    {"L1_LastBunchInTrain", &L1_LastBunchInTrain},
    {"L1_LastCollisionInTrain", &L1_LastCollisionInTrain},
    {"L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_IsoTau26er2p1_dR_Min0p3},
    {"L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3", &L1_LooseIsoEG22er2p1_Tau70er2p1_dR_Min0p3},
    {"L1_LooseIsoEG24er2p1_HTT100er", &L1_LooseIsoEG24er2p1_HTT100er},
    {"L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3", &L1_LooseIsoEG24er2p1_IsoTau27er2p1_dR_Min0p3},
    {"L1_LooseIsoEG26er2p1_HTT100er", &L1_LooseIsoEG26er2p1_HTT100er},
    {"L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG26er2p1_Jet34er2p5_dR_Min0p3},
    {"L1_LooseIsoEG28er2p1_HTT100er", &L1_LooseIsoEG28er2p1_HTT100er},
    {"L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG28er2p1_Jet34er2p5_dR_Min0p3},
    {"L1_LooseIsoEG30er2p1_HTT100er", &L1_LooseIsoEG30er2p1_HTT100er},
    {"L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3", &L1_LooseIsoEG30er2p1_Jet34er2p5_dR_Min0p3},
    {"L1_MinimumBiasHF0", &L1_MinimumBiasHF0},
    {"L1_MinimumBiasHF0_AND_BptxAND", &L1_MinimumBiasHF0_AND_BptxAND},
    {"L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6", &L1_Mu10er2p3_Jet32er2p3_dR_Max0p4_DoubleJet32er2p3_dEta_Max1p6},
    {"L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p1_dR_Max0p4_DoubleJet40er2p1_dEta_Max1p6},
    {"L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6", &L1_Mu12er2p3_Jet40er2p3_dR_Max0p4_DoubleJet40er2p3_dEta_Max1p6},
    {"L1_Mu18er2p1_Tau24er2p1", &L1_Mu18er2p1_Tau24er2p1},
    {"L1_Mu18er2p1_Tau26er2p1", &L1_Mu18er2p1_Tau26er2p1},
    {"L1_Mu18er2p1_Tau26er2p1_Jet55", &L1_Mu18er2p1_Tau26er2p1_Jet55},
    {"L1_Mu18er2p1_Tau26er2p1_Jet70", &L1_Mu18er2p1_Tau26er2p1_Jet70},
    {"L1_Mu20_EG10er2p5", &L1_Mu20_EG10er2p5},
    {"L1_Mu22er2p1_IsoTau28er2p1", &L1_Mu22er2p1_IsoTau28er2p1},
    {"L1_Mu22er2p1_IsoTau30er2p1", &L1_Mu22er2p1_IsoTau30er2p1},
    {"L1_Mu22er2p1_IsoTau32er2p1", &L1_Mu22er2p1_IsoTau32er2p1},
    {"L1_Mu22er2p1_IsoTau34er2p1", &L1_Mu22er2p1_IsoTau34er2p1},
    {"L1_Mu22er2p1_IsoTau36er2p1", &L1_Mu22er2p1_IsoTau36er2p1},
    {"L1_Mu22er2p1_IsoTau40er2p1", &L1_Mu22er2p1_IsoTau40er2p1},
    {"L1_Mu22er2p1_Tau70er2p1", &L1_Mu22er2p1_Tau70er2p1},
    {"L1_Mu3_Jet120er2p5_dR_Max0p4", &L1_Mu3_Jet120er2p5_dR_Max0p4},
    {"L1_Mu3_Jet120er2p5_dR_Max0p8", &L1_Mu3_Jet120er2p5_dR_Max0p8},
    {"L1_Mu3_Jet16er2p5_dR_Max0p4", &L1_Mu3_Jet16er2p5_dR_Max0p4},
    {"L1_Mu3_Jet30er2p5", &L1_Mu3_Jet30er2p5},
    {"L1_Mu3_Jet35er2p5_dR_Max0p4", &L1_Mu3_Jet35er2p5_dR_Max0p4},
    {"L1_Mu3_Jet60er2p5_dR_Max0p4", &L1_Mu3_Jet60er2p5_dR_Max0p4},
    {"L1_Mu3_Jet80er2p5_dR_Max0p4", &L1_Mu3_Jet80er2p5_dR_Max0p4},
    {"L1_Mu3er1p5_Jet100er2p5_ETMHF30", &L1_Mu3er1p5_Jet100er2p5_ETMHF30},
    {"L1_Mu3er1p5_Jet100er2p5_ETMHF40", &L1_Mu3er1p5_Jet100er2p5_ETMHF40},
    {"L1_Mu3er1p5_Jet100er2p5_ETMHF50", &L1_Mu3er1p5_Jet100er2p5_ETMHF50},
    {"L1_Mu5_EG23er2p5", &L1_Mu5_EG23er2p5},
    {"L1_Mu5_LooseIsoEG20er2p5", &L1_Mu5_LooseIsoEG20er2p5},
    {"L1_Mu6_DoubleEG10er2p5", &L1_Mu6_DoubleEG10er2p5},
    {"L1_Mu6_DoubleEG12er2p5", &L1_Mu6_DoubleEG12er2p5},
    {"L1_Mu6_DoubleEG15er2p5", &L1_Mu6_DoubleEG15er2p5},
    {"L1_Mu6_DoubleEG17er2p5", &L1_Mu6_DoubleEG17er2p5},
    {"L1_Mu6_HTT240er", &L1_Mu6_HTT240er},
    {"L1_Mu6_HTT250er", &L1_Mu6_HTT250er},
    {"L1_Mu7_EG20er2p5", &L1_Mu7_EG20er2p5},
    {"L1_Mu7_EG23er2p5", &L1_Mu7_EG23er2p5},
    {"L1_Mu7_LooseIsoEG20er2p5", &L1_Mu7_LooseIsoEG20er2p5},
    {"L1_Mu7_LooseIsoEG23er2p5", &L1_Mu7_LooseIsoEG23er2p5},
    {"L1_NotBptxOR", &L1_NotBptxOR},
    {"L1_QuadJet60er2p5", &L1_QuadJet60er2p5},
    {"L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0", &L1_QuadJet_95_75_65_20_DoubleJet_75_65_er2p5_Jet20_FWD3p0},
    {"L1_QuadMu0", &L1_QuadMu0},
    {"L1_QuadMu0_OQ", &L1_QuadMu0_OQ},
    {"L1_QuadMu0_SQ", &L1_QuadMu0_SQ},
    {"L1_SecondBunchInTrain", &L1_SecondBunchInTrain},
    {"L1_SecondLastBunchInTrain", &L1_SecondLastBunchInTrain},
    {"L1_SingleEG10er2p5", &L1_SingleEG10er2p5},
    {"L1_SingleEG15er2p5", &L1_SingleEG15er2p5},
    {"L1_SingleEG26er2p5", &L1_SingleEG26er2p5},
    {"L1_SingleEG28_FWD2p5", &L1_SingleEG28_FWD2p5},
    {"L1_SingleEG28er1p5", &L1_SingleEG28er1p5},
    {"L1_SingleEG28er2p1", &L1_SingleEG28er2p1},
    {"L1_SingleEG28er2p5", &L1_SingleEG28er2p5},
    {"L1_SingleEG34er2p5", &L1_SingleEG34er2p5},
    {"L1_SingleEG36er2p5", &L1_SingleEG36er2p5},
    {"L1_SingleEG38er2p5", &L1_SingleEG38er2p5},
    {"L1_SingleEG40er2p5", &L1_SingleEG40er2p5},
    {"L1_SingleEG42er2p5", &L1_SingleEG42er2p5},
    {"L1_SingleEG45er2p5", &L1_SingleEG45er2p5},
    {"L1_SingleEG50", &L1_SingleEG50},
    {"L1_SingleEG60", &L1_SingleEG60},
    {"L1_SingleEG8er2p5", &L1_SingleEG8er2p5},
    {"L1_SingleIsoEG24er1p5", &L1_SingleIsoEG24er1p5},
    {"L1_SingleIsoEG24er2p1", &L1_SingleIsoEG24er2p1},
    {"L1_SingleIsoEG26er1p5", &L1_SingleIsoEG26er1p5},
    {"L1_SingleIsoEG26er2p1", &L1_SingleIsoEG26er2p1},
    {"L1_SingleIsoEG26er2p5", &L1_SingleIsoEG26er2p5},
    {"L1_SingleIsoEG28_FWD2p5", &L1_SingleIsoEG28_FWD2p5},
    {"L1_SingleIsoEG28er1p5", &L1_SingleIsoEG28er1p5},
    {"L1_SingleIsoEG28er2p1", &L1_SingleIsoEG28er2p1},
    {"L1_SingleIsoEG28er2p5", &L1_SingleIsoEG28er2p5},
    {"L1_SingleIsoEG30er2p1", &L1_SingleIsoEG30er2p1},
    {"L1_SingleIsoEG30er2p5", &L1_SingleIsoEG30er2p5},
    {"L1_SingleIsoEG32er2p1", &L1_SingleIsoEG32er2p1},
    {"L1_SingleIsoEG32er2p5", &L1_SingleIsoEG32er2p5},
    {"L1_SingleIsoEG34er2p5", &L1_SingleIsoEG34er2p5},
    {"L1_SingleIsoTau32er2p1", &L1_SingleIsoTau32er2p1},
    {"L1_SingleJet10erHE", &L1_SingleJet10erHE},
    {"L1_SingleJet120", &L1_SingleJet120},
    {"L1_SingleJet120_FWD3p0", &L1_SingleJet120_FWD3p0},
    {"L1_SingleJet120er2p5", &L1_SingleJet120er2p5},
    {"L1_SingleJet12erHE", &L1_SingleJet12erHE},
    {"L1_SingleJet140er2p5", &L1_SingleJet140er2p5},
    {"L1_SingleJet140er2p5_ETMHF70", &L1_SingleJet140er2p5_ETMHF70},
    {"L1_SingleJet140er2p5_ETMHF80", &L1_SingleJet140er2p5_ETMHF80},
    {"L1_SingleJet140er2p5_ETMHF90", &L1_SingleJet140er2p5_ETMHF90},
    {"L1_SingleJet160er2p5", &L1_SingleJet160er2p5},
    {"L1_SingleJet180", &L1_SingleJet180},
    {"L1_SingleJet180er2p5", &L1_SingleJet180er2p5},
    {"L1_SingleJet200", &L1_SingleJet200},
    {"L1_SingleJet20er2p5_NotBptxOR", &L1_SingleJet20er2p5_NotBptxOR},
    {"L1_SingleJet20er2p5_NotBptxOR_3BX", &L1_SingleJet20er2p5_NotBptxOR_3BX},
    {"L1_SingleJet35", &L1_SingleJet35},
    {"L1_SingleJet35_FWD3p0", &L1_SingleJet35_FWD3p0},
    {"L1_SingleJet35er2p5", &L1_SingleJet35er2p5},
    {"L1_SingleJet43er2p5_NotBptxOR_3BX", &L1_SingleJet43er2p5_NotBptxOR_3BX},
    {"L1_SingleJet46er2p5_NotBptxOR_3BX", &L1_SingleJet46er2p5_NotBptxOR_3BX},
    {"L1_SingleJet60", &L1_SingleJet60},
    {"L1_SingleJet60_FWD3p0", &L1_SingleJet60_FWD3p0},
    {"L1_SingleJet60er2p5", &L1_SingleJet60er2p5},
    {"L1_SingleJet8erHE", &L1_SingleJet8erHE},
    {"L1_SingleJet90", &L1_SingleJet90},
    {"L1_SingleJet90_FWD3p0", &L1_SingleJet90_FWD3p0},
    {"L1_SingleJet90er2p5", &L1_SingleJet90er2p5},
    {"L1_SingleLooseIsoEG26er1p5", &L1_SingleLooseIsoEG26er1p5},
    {"L1_SingleLooseIsoEG26er2p5", &L1_SingleLooseIsoEG26er2p5},
    {"L1_SingleLooseIsoEG28_FWD2p5", &L1_SingleLooseIsoEG28_FWD2p5},
    {"L1_SingleLooseIsoEG28er1p5", &L1_SingleLooseIsoEG28er1p5},
    {"L1_SingleLooseIsoEG28er2p1", &L1_SingleLooseIsoEG28er2p1},
    {"L1_SingleLooseIsoEG28er2p5", &L1_SingleLooseIsoEG28er2p5},
    {"L1_SingleLooseIsoEG30er1p5", &L1_SingleLooseIsoEG30er1p5},
    {"L1_SingleLooseIsoEG30er2p5", &L1_SingleLooseIsoEG30er2p5},
    {"L1_SingleMu0_BMTF", &L1_SingleMu0_BMTF},
    {"L1_SingleMu0_DQ", &L1_SingleMu0_DQ},
    {"L1_SingleMu0_EMTF", &L1_SingleMu0_EMTF},
    {"L1_SingleMu0_OMTF", &L1_SingleMu0_OMTF},
    {"L1_SingleMu10er1p5", &L1_SingleMu10er1p5},
    {"L1_SingleMu12_DQ_BMTF", &L1_SingleMu12_DQ_BMTF},
    {"L1_SingleMu12_DQ_EMTF", &L1_SingleMu12_DQ_EMTF},
    {"L1_SingleMu12_DQ_OMTF", &L1_SingleMu12_DQ_OMTF},
    {"L1_SingleMu12er1p5", &L1_SingleMu12er1p5},
    {"L1_SingleMu14er1p5", &L1_SingleMu14er1p5},
    {"L1_SingleMu15_DQ", &L1_SingleMu15_DQ},
    {"L1_SingleMu16er1p5", &L1_SingleMu16er1p5},
    {"L1_SingleMu18", &L1_SingleMu18},
    {"L1_SingleMu18er1p5", &L1_SingleMu18er1p5},
    {"L1_SingleMu20", &L1_SingleMu20},
    {"L1_SingleMu22", &L1_SingleMu22},
    {"L1_SingleMu22_BMTF", &L1_SingleMu22_BMTF},
    {"L1_SingleMu22_DQ", &L1_SingleMu22_DQ},
    {"L1_SingleMu22_EMTF", &L1_SingleMu22_EMTF},
    {"L1_SingleMu22_OMTF", &L1_SingleMu22_OMTF},
    {"L1_SingleMu22_OQ", &L1_SingleMu22_OQ},
    {"L1_SingleMu25", &L1_SingleMu25},
    {"L1_SingleMu3", &L1_SingleMu3},
    {"L1_SingleMu5", &L1_SingleMu5},
    {"L1_SingleMu6er1p5", &L1_SingleMu6er1p5},
    {"L1_SingleMu7", &L1_SingleMu7},
    {"L1_SingleMu7_DQ", &L1_SingleMu7_DQ},
    {"L1_SingleMu7er1p5", &L1_SingleMu7er1p5},
    {"L1_SingleMu8er1p5", &L1_SingleMu8er1p5},
    {"L1_SingleMu9er1p5", &L1_SingleMu9er1p5},
    {"L1_SingleMuCosmics", &L1_SingleMuCosmics},
    {"L1_SingleMuCosmics_BMTF", &L1_SingleMuCosmics_BMTF},
    {"L1_SingleMuCosmics_EMTF", &L1_SingleMuCosmics_EMTF},
    {"L1_SingleMuCosmics_OMTF", &L1_SingleMuCosmics_OMTF},
    {"L1_SingleMuOpen", &L1_SingleMuOpen},
    {"L1_SingleMuOpen_NotBptxOR", &L1_SingleMuOpen_NotBptxOR},
    {"L1_SingleMuOpen_er1p1_NotBptxOR_3BX", &L1_SingleMuOpen_er1p1_NotBptxOR_3BX},
    {"L1_SingleMuOpen_er1p4_NotBptxOR_3BX", &L1_SingleMuOpen_er1p4_NotBptxOR_3BX},
    {"L1_SingleMuShower_Nominal", &L1_SingleMuShower_Nominal},
    {"L1_SingleMuShower_Tight", &L1_SingleMuShower_Tight},
    {"L1_SingleTau120er2p1", &L1_SingleTau120er2p1},
    {"L1_SingleTau130er2p1", &L1_SingleTau130er2p1},
    {"L1_SingleTau70er2p1", &L1_SingleTau70er2p1},
    {"L1_TOTEM_1", &L1_TOTEM_1},
    {"L1_TOTEM_2", &L1_TOTEM_2},
    {"L1_TOTEM_3", &L1_TOTEM_3},
    {"L1_TOTEM_4", &L1_TOTEM_4},
    {"L1_TripleEG16er2p5", &L1_TripleEG16er2p5},
    {"L1_TripleEG_16_12_8_er2p5", &L1_TripleEG_16_12_8_er2p5},
    {"L1_TripleEG_16_15_8_er2p5", &L1_TripleEG_16_15_8_er2p5},
    {"L1_TripleEG_18_17_8_er2p5", &L1_TripleEG_18_17_8_er2p5},
    {"L1_TripleEG_18_18_12_er2p5", &L1_TripleEG_18_18_12_er2p5},
    {"L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5", &L1_TripleJet_100_80_70_DoubleJet_80_70_er2p5},
    {"L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5", &L1_TripleJet_105_85_75_DoubleJet_85_75_er2p5},
    {"L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5", &L1_TripleJet_95_75_65_DoubleJet_75_65_er2p5},
    {"L1_TripleMu0", &L1_TripleMu0},
    {"L1_TripleMu0_OQ", &L1_TripleMu0_OQ},
    {"L1_TripleMu0_SQ", &L1_TripleMu0_SQ},
    {"L1_TripleMu3", &L1_TripleMu3},
    {"L1_TripleMu3_SQ", &L1_TripleMu3_SQ},
    {"L1_TripleMu_2SQ_1p5SQ_0OQ", &L1_TripleMu_2SQ_1p5SQ_0OQ},
    {"L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12", &L1_TripleMu_2SQ_1p5SQ_0OQ_Mass_Max12},
    {"L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12", &L1_TripleMu_3SQ_2p5SQ_0OQ_Mass_Max12},
    {"L1_TripleMu_5SQ_3SQ_0OQ", &L1_TripleMu_5SQ_3SQ_0OQ},
    {"L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0OQ_DoubleMu_5_3_SQ_OS_Mass_Max9},
    {"L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9", &L1_TripleMu_5SQ_3SQ_0_DoubleMu_5_3_SQ_OS_Mass_Max9},
    {"L1_TripleMu_5_3_3", &L1_TripleMu_5_3_3},
    {"L1_TripleMu_5_3_3_SQ", &L1_TripleMu_5_3_3_SQ},
    {"L1_TripleMu_5_3p5_2p5", &L1_TripleMu_5_3p5_2p5},
    {"L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_DoubleMu_5_2p5_OS_Mass_5to17},
    {"L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17", &L1_TripleMu_5_3p5_2p5_OQ_DoubleMu_5_2p5_OQ_OS_Mass_5to17},
    {"L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17", &L1_TripleMu_5_4_2p5_DoubleMu_5_2p5_OS_Mass_5to17},
    {"L1_TripleMu_5_5_3", &L1_TripleMu_5_5_3},
    {"L1_UnpairedBunchBptxMinus", &L1_UnpairedBunchBptxMinus},
    {"L1_UnpairedBunchBptxPlus", &L1_UnpairedBunchBptxPlus},
    {"L1_ZeroBias", &L1_ZeroBias},
    {"L1_ZeroBias_copy", &L1_ZeroBias_copy}
  };

  for (const auto pair : name2func)
  {
    L1SeedFun[pair.first] = std::bind(pair.second, upgrade, calo_tower);
  }

  return true;
}
// eof