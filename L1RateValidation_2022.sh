#!/bin/bash -ex

echo "Running automated Level-1 Trigger Rate validation script. Will compare rates menus using"
echo "reference and test GTs"
echo " "

###############################
starttime=$(date +%s.%N)
ARCH=slc7_amd64_gcc900
CMSREL=CMSSW_12_0_2
L1TTag=l1t-integration-v109.0
GT=120X_dataRun2_v2
Prescale=Prescale_2018_v2_1_0_Col_2.0.txt
nproc=`nproc`
sqlite1=$1 ##ref
sqlite2=$2
week=$3
year=$4
curdir=$PWD
username=$USER
pids=""
hasref=false
## ZeroBias Raw, Fill#6961 Run320065, LS1-225.6
filelist=('/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/1262EAC6-FD8D-E811-86BB-FA163E7FEBB7.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/92EB0BE7-FD8D-E811-BDEC-FA163E17A15F.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/5E19E964-FE8D-E811-9D38-FA163E7D9AFE.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/26A77465-FE8D-E811-B971-FA163E4200C7.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/2A857540-FF8D-E811-82C4-FA163E6011FE.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/185AD6BE-FD8D-E811-97D1-FA163E6AF4C5.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/54F0E73C-008E-E811-AAFA-FA163EBDCF4F.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/181C7B47-FF8D-E811-B626-02163E017F0F.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/C05C3D5B-FF8D-E811-8230-FA163EAE1B00.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/5E96232B-FE8D-E811-BE39-FA163EB3D34C.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/08AE8301-018E-E811-8793-FA163E17F3DA.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/201B654B-008E-E811-9422-02163E019EEA.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/94AD4962-FE8D-E811-94D6-02163E016CF8.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/9E6A1951-008E-E811-B29C-02163E01A04A.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/AE22FCF5-008E-E811-90C4-FA163E2C1E9E.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/10A0979D-018E-E811-B4B2-FA163E17A15F.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/008F1E66-008E-E811-8DCF-02163E010C36.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/F2F1445E-008E-E811-9E90-FA163E780F6E.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/C2C70E61-018E-E811-A6A4-FA163E06D84E.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/D43AB599-018E-E811-B24C-FA163E792AF4.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/EE166C90-038E-E811-AFBE-FA163E17F3DA.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/D47736E9-028E-E811-9913-FA163E31E513.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/12575FE9-028E-E811-820D-FA163EDA17B6.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/98625626-018E-E811-A48D-FA163E41C7B0.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/2E6A0D44-028E-E811-864D-FA163EABC162.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/C0A28D6C-028E-E811-BEE5-FA163E96DBF1.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/4450E207-048E-E811-9AC8-FA163EF250BA.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/2082E028-048E-E811-B3EF-FA163E254513.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/3C9C9177-058E-E811-8EE4-02163E017EB4.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/DA4EE943-048E-E811-968A-FA163E3B3CC6.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/680BB85F-068E-E811-B16B-FA163E3F1FDB.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/F0BC4E2C-058E-E811-A0D5-FA163E071249.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/6609EED8-058E-E811-8544-02163E010CFD.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/CCC49B3A-078E-E811-99C5-FA163EDC261A.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/588B001F-078E-E811-A25A-FA163EBD447D.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/80F1A88A-058E-E811-8F98-FA163E74976F.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/7CFA3727-078E-E811-895B-FA163E5CA922.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/44789C42-058E-E811-A789-02163E017F18.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/7ED5543C-078E-E811-B0EC-FA163E5193D5.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/6243C8AF-058E-E811-B3EB-FA163E1DC155.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/90EF3657-068E-E811-A473-02163E014D9D.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/0680DDF3-028E-E811-88C2-FA163E82CF22.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/EA95C67E-088E-E811-A8F5-FA163ED3058C.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/B0196D1B-088E-E811-A262-02163E017EEE.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/5E14504F-078E-E811-991A-FA163E8103B1.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/ACBBA428-078E-E811-8C81-FA163E733247.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/E6BB9037-098E-E811-AD51-FA163E6AE23B.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/4861F63D-098E-E811-A5D9-FA163E5B6746.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/EEABAB4C-098E-E811-BE1E-FA163E836BB3.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/04948A79-098E-E811-BEBC-FA163EA375DF.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/7AA6E404-0A8E-E811-9DA4-FA163E229782.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/0CB1B807-008E-E811-A9BD-FA163E8D7B5B.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/8095D777-098E-E811-A174-02163E01A0F8.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/086D114D-098E-E811-B1CD-02163E01800D.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/5A4D471B-0A8E-E811-A2FD-FA163EA70EE0.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/E6FFF613-0A8E-E811-A397-FA163E3F1FDB.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/94663FB1-0A8E-E811-A930-FA163E3BB583.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/88949DC2-0A8E-E811-A6C4-02163E019F56.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/725D1C32-0A8E-E811-A043-FA163E95ED85.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/56B25898-0B8E-E811-B8DB-FA163EE92492.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/52648233-0C8E-E811-A941-FA163E6FB005.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/A8057F26-0C8E-E811-A9E7-FA163E5DEA72.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/E6AE3F7A-0B8E-E811-9484-FA163E6F4C2E.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/6A1E2050-0C8E-E811-B6A4-FA163E3B8B97.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/54878731-0C8E-E811-9A40-FA163EE7E719.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/7A97134C-0C8E-E811-93E5-02163E01A03B.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/620308A3-0D8E-E811-9BAC-FA163E8F29FC.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/981AFAFE-088E-E811-BB3E-FA163EB7F77E.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/347275C2-0D8E-E811-9403-FA163EE7A04A.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/96B563DD-048E-E811-AAD9-02163E00B8DC.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/9C1FEEBF-0E8E-E811-AFC1-FA163E34F1D6.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/B6BF70D8-0E8E-E811-93B5-FA163E838299.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/2242FED4-0E8E-E811-9768-FA163ED24770.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/DAB74636-0F8E-E811-AC7A-FA163E5DEA72.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/C4350BE2-0E8E-E811-A05A-FA163EDFA1F1.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/28FAD4C2-0E8E-E811-AF2C-FA163E302745.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/1A1C15C1-0C8E-E811-A947-FA163E9A9A14.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/6C821B0E-108E-E811-B740-FA163E6CDF89.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/D093DFF2-0F8E-E811-A281-A4BF0114D040.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/20AA7316-108E-E811-98E3-FA163EAD1C70.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/4CB3C62E-108E-E811-9E32-FA163EBEB6D5.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/A8258709-108E-E811-B127-FA163EC648D9.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/107FD16F-118E-E811-995D-FA163E0EA436.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/38A7786E-118E-E811-B464-FA163E2B0A7A.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/12650447-108E-E811-A2DA-FA163E59D09E.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/1E7FBDF0-0E8E-E811-9B11-FA163E4F1151.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/123C7481-118E-E811-86C2-FA163E91FAAB.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/0ED5F44E-108E-E811-A33F-FA163E8E607A.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/A0DBBC09-138E-E811-A18B-FA163E27069D.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/E0D865C8-0E8E-E811-AB26-FA163ECA1847.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/52184C6F-118E-E811-A105-FA163E2C1E9E.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/7CB98C02-138E-E811-94EE-FA163E5F0A23.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/CCB36A39-138E-E811-BA6B-FA163E2D4503.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/48E2731C-138E-E811-81B6-FA163E3DDA6C.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/8AB241D2-138E-E811-9B24-FA163E1ED478.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/5A37C121-138E-E811-B8E1-FA163E792AF4.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/7063828D-128E-E811-A3E6-FA163E48F925.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/7A1BC243-148E-E811-8269-FA163E1190B1.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/1438B617-138E-E811-B190-FA163E6B9C65.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/B236E5B3-148E-E811-8F73-FA163E617795.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/BE2FB18C-158E-E811-BBCE-FA163E68CA49.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/080FE4BE-138E-E811-9D97-FA163E5043E5.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/10E03B5B-148E-E811-AEC9-FA163E75B351.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/80BEDB66-158E-E811-9CCC-02163E012FC0.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/F0BE7017-328E-E811-BE23-FA163E1C56E0.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/883CDB14-328E-E811-952A-FA163EB2DCAF.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/E8F0B714-328E-E811-96E6-FA163E4F129E.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/02DD7A19-328E-E811-8DA8-FA163ED07727.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/C41C181B-328E-E811-B164-02163E019F93.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/FE5C701B-328E-E811-A5A2-FA163E166081.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/F6A37338-328E-E811-93B9-FA163E335778.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/306E8117-328E-E811-861A-FA163E8CA0FA.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/E8034637-328E-E811-A399-FA163EAD2729.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/68869D46-328E-E811-A6F7-FA163E9D3E3D.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/50ED1F1B-328E-E811-A9DC-FA163EECEDCF.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/C2B6E519-328E-E811-9F9B-FA163E4C86A6.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/6A7C0118-328E-E811-8EBB-FA163E0B9045.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/54B3F115-328E-E811-A361-FA163E0CEA25.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/0A7D2017-328E-E811-AA9F-FA163ED1A85C.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/1ADA9822-328E-E811-8392-02163E017FBB.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/301EFD13-328E-E811-BDA5-02163E01A00B.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/20532837-328E-E811-B8D7-FA163EEB9FC0.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/C4589615-328E-E811-8E4C-FA163EF7D5EC.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/6C521621-328E-E811-9F7E-02163E01A169.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/DAA6A434-328E-E811-9505-FA163E1455BD.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/9C55BE18-328E-E811-8ADE-FA163E2B77C9.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/B2460E21-328E-E811-B80B-FA163E34934D.root'
'/store/data/Run2018C/ZeroBias/RAW/v1/000/320/065/00000/9A73591B-328E-E811-86A2-FA163E95D929.root'
)
#echo $filelist
#----------------------------------------------------------------------------#
#                                 Upload Ref                                 #
#----------------------------------------------------------------------------#
if [ -z "$WORKSPACE" ]
then
  WORKSPACE=${curdir}
fi

if [ -f ${WORKSPACE}/upload/$2 ]
then
  echo "dir is already existing"
  touch ${WORKSPACE}/upload/$2/.jenkins-upload
else
  mkdir -p ${WORKSPACE}/upload/$2
  touch ${WORKSPACE}/upload/$2/.jenkins-upload
fi
#----------------------------------------------------------------------------#
#                            Getting the reference                           #
#----------------------------------------------------------------------------#
#it may be gzipped infact ...
## prevent exit from failed wget
set +e 
#wget --no-check-certificate https://cmssdt.cern.ch/SDT/public/EcalLaserValidation/L1T_EcalLaserValidation/345982/L1TEcalValidation_2021_48_345982.tgz 
wget --no-check-certificate https://cmssdt.cern.ch/SDT/public/EcalLaserValidation/L1T_EcalLaserValidation/${sqlite1}/L1TEcalValidation_${year}_${week}_${sqlite1}.tgz 
if [ $? -ne 0 ]; then
  sqs="$sqlite1 $sqlite2"
else
  sqs=$sqlite2
  hasref=true
fi
#sqs="$sqlite1 $sqlite2"
echo "Running ECal validtion with ", $sqs

#----------------------------------------------------------------------------#
#                            Checkout L1 Emulator                            #
#----------------------------------------------------------------------------#
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=$ARCH
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
scramv1 project CMSSW $CMSREL
cd $CMSREL/src
eval `scramv1 runtime -sh`
git-cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline l1t-integration-CMSSW_12_0_2
git cms-merge-topic -u cms-l1t-offline:$L1TTag
#git cms-addpkg L1Trigger/L1TCommon
#git cms-addpkg L1Trigger/L1TMuon
#git clone https://github.com/cms-l1t-offline/L1Trigger-L1TMuon.git L1Trigger/L1TMuon/data
#git cms-addpkg L1Trigger/L1TCalorimeter
#git clone https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter.git L1Trigger/L1TCalorimeter/data

git cms-checkdeps -A -a

scram b -j ${nproc}

dur=$(echo "$(date +%s.%N) - $starttime" | bc)
printf "Execution time to L1T checkout: %.6f seconds" $dur
#----------------------------------------------------------------------------#
#                          Running the TP emulation                          #
#----------------------------------------------------------------------------#
echo running $GT
cmsDriver.py l1NtupleRAWEMU_2018 -s RAW2DIGI --era=Run2_2018  \
  --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMU \
  --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimEcalTP \
  --conditions=$GT -n 4000 --data --no_exec --no_output  \
  --filein=inputFiles \
  --python_filename=l1Ntuple_${GT}.py 
Nsq=`echo $sqs | awk -F ' ' '{print NF}'`
#Nfiles=$((wc -l <fileList_320065.txt))
Nfiles=${#filelist[@]}
NfpJ=`echo "${Nfiles} *${Nsq}/${nproc}" | bc`
NJ=`echo "${Nfiles}/${NfpJ}" | bc`
for sq in $sqs; do
  if [ ! -f EcalTPG_${sq}_moved_to_1.db ]; then
    wget http://cern.ch/ecaltrg/EcalLin/EcalTPG_${sq}_moved_to_1.db
  fi
  python ${curdir}/ModifyL1Ntuple.py --globalTag $GT --sqlite $sq
  for ((i = 0; i < $NJ; i++)); do
    let cnt1=$(($i*$NfpJ))
    args=`printf "inputFiles=%s " "${filelist[@]:$cnt1:$NfpJ}"`
    args+=`echo outputFile=L1Ntuple_${GT}_${sq}_${i}.root`
    timeout 1800 cmsRun l1Ntuple_${GT}_${sq}.py `echo $args` >& l1Ntuple_${GT}_${sq}_${i}.log  &
    pids="$pids $!"
  done
done
echo "Waiting for Ntuple production to finish......"
wait $pids
dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to L1Ntuple production: %.6f minutes" $dur

for sq in $sqs; do
  ls $PWD/L1Ntuple_${GT}_${sq}_*.root > L1Ntuple_${GT}_${sq}.list
  cp $PWD/L1Ntuple_${GT}_${sq}_*.log ${WORKSPACE}/upload/${2}/
done

################################


#----------------------------------------------------------------------------#
#                            Check out L1Menu code                           #
#----------------------------------------------------------------------------#
git clone --depth 1 https://github.com/cms-l1-dpg/L1MenuTools.git
cd L1MenuTools/rate-estimation/
cp $curdir/CompL1Rate.py  .
cp $curdir/menulib.cc .
cp $curdir/menulib.hh .
cp $curdir/$Prescale menu/
cp $curdir/run_lumi.csv menu/
cp $curdir/Selected_Seed.txt menu/
mkdir -p objs/include
make -j ${nproc}
#make comparePlots

:'
#----------------------------------------------------------------------------#
#                                 Lumi Table                                 #
#----------------------------------------------------------------------------#
export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda/bin:$PATH
pip install --user --upgrade brilws
cd menu
source GetLumi_setup.sh
./GetLumi.py
cd ..

dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to checkout and compile code: %.6f minutes" $dur
:'
#----------------------------------------------------------------------------#
#                                 Run L1Menu                                 #
#----------------------------------------------------------------------------#

for sq in $sqs; do
  ./testMenu2016 -u menu/run_lumi.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2544 --doPlotRate --doPlotEff --allPileUp --SelectCol 2E+34 >& L1Menu_${GT}_${sq}_emu.log &
#  ./testMenu2016 -u menu/run_lumi.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2544 --doPlotRate --doPlotEff >& L1Menu_${GT}_${sq}_emu.log &
  pids="$pids $!"
  ./testMenu2016 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu --allPileUp --SelectCol 2E+34 >& L1Seed_${GT}_${sq}_emu.log &
#  ./testMenu2016 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu >& L1Seed_${GT}_${sq}_emu.log &
  pids="$pids $!"
done
echo "Waiting for menu rate estimation to finish......"
wait $pids
dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time to L1Ntuple production: %.6f minutes" $dur
cp L1Menu_${GT}_*_emu.log ${WORKSPACE}/upload/${2}/
cp L1Seed_${GT}_*_emu.log ${WORKSPACE}/upload/${2}/

#----------------------------------------------------------------------------#
#                                Compare rate                                #
#----------------------------------------------------------------------------#

if $hasref; then
#  tar -xzvf $curdir/L1TEcalValidation_2021_48_345982.tgz -C results/
  tar -xzvf $curdir/L1TEcalValidation_${year}_${week}_${sqlite1}.tgz -C results/
fi

ls results/

python3 CompL1Rate.py --globalTag $GT --sqlite1 $sqlite1 --sqlite2 $sqlite2  | tee ${sqlite2}.log
#./comparePlots results/L1Seed*root

#----------------------------------------------------------------------------#
#                                 Upload Ref                                 #
#----------------------------------------------------------------------------#

#we need to make a tar gz of this one
cp results/L1Menu_${GT}_${sqlite2}_emu.csv ${WORKSPACE}/upload/${2}/
cp results/L1Seed_${GT}_${sqlite2}_emu.csv ${WORKSPACE}/upload/${2}/
cp results/L1Seed_${GT}_${sqlite2}_emu.root ${WORKSPACE}/upload/${2}/
cp ${sqlite2}.log ${WORKSPACE}/upload/${2}/
cp compRate.csv ${WORKSPACE}/upload/${2}/

#for i in nDiEGVsPt.gif nEGVsEta.gif nEGVsPt.gif nETMVsETM.gif nIsoEGVsPt.gif nJetVsPt.gif nETMVsETMHF.gif nQuadCenJetVsPt.gif nTauVsPt.gif; do
#  cp  results/comparePlots/Rate/${i}  ${WORKSPACE}/upload/${2}/
#done
cd ${WORKSPACE}/upload/${2}
tar -czvf ${WORKSPACE}/upload/L1TEcalValidation_${year}_${week}_${sqlite2}.tgz *

dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time of workflow: %.6f minutes" $dur
