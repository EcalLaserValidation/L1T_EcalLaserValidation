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
filelist=`cat fileList_320065.txt`

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
  --conditions=$GT -n -1 --data --no_exec --no_output  \
  --filein=inputFiles \
  --python_filename=l1Ntuple_${GT}.py 
Nsq=`echo $sqs | awk -F ' ' '{print NF}'`
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
    cmsRun l1Ntuple_${GT}_${sq}.py `echo $args` >& l1Ntuple_${GT}_${sq}_${i}.log  &
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
