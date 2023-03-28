#!/bin/bash -ex

echo "Running automated Level-1 Trigger Rate validation script. Will compare rates menus using"
echo "reference and test GTs"
echo " "

###############################
starttime=$(date +%s.%N)
ARCH=slc7_amd64_gcc11
CMSREL=CMSSW_13_0_0_pre4
L1TTag=l1t-integration-v147
GT=130X_dataRun3_Prompt_v1
Prescale=Prescale_2022_v1_4_0.csv
nproc=`nproc`
sqlite1=$1 ##ref
sqlite2=$2
week=$3
year=$4
curdir=$PWD
username=$USER
pids=""
hasref=false
## ZeroBias Raw, Fill#8128 Run357479, LS1-945
filelist=('/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/00c36f93-776b-4fd9-973a-2daf58e51b45.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/02caebd0-de05-4352-a7fa-8e37b75a6710.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/030fc479-ba14-4952-a146-12406a4b7c37.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/079b9639-a512-49b5-81c6-ea2ef7da2290.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/0a723b3a-434d-42e7-89f2-ed5e1ae96376.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/0c6909f5-9912-4b7a-b71e-dd4de1d27e61.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/0ec04dae-725e-4dfa-9387-afd7860ed335.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/0ee538d7-2c81-4359-b043-b17b4555354b.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/0f190e64-d0df-42ca-a5ed-2cd51a5ac05e.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/0f4bd30d-40d7-43e1-b1d7-aa132959d158.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/1510172f-58d9-48a3-99f7-3cdf4ee5591d.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/15cd6b6b-efd1-47ae-b874-7a21301296a7.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/15d89f1d-a4ea-4574-acce-340be8c351d6.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/16fd5126-9800-4f6e-9130-e676de6ae76d.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/18881777-2377-4b78-b3ca-e969c8a4fc2c.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/18b52df9-acba-488a-a847-c04710fad1ce.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/20ab59f8-0ce9-404a-9886-0232fee2920a.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/21d19d1a-fc61-4d05-a3af-3d6b6c017ac2.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/22cefcb2-5e80-4e06-8e29-cf6fbe6629d1.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/23dbd469-c9aa-46de-a44b-f3603f0f37e0.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/2585a350-ad16-48ba-a157-47f794bfbc6d.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/25df90c0-0351-46bd-b944-28ba019713f4.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/2983a1d3-678c-426b-92dd-1819fd8b8c72.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/2b8425aa-c777-482f-9f59-7f930355f8d0.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/2bbc03b8-c282-47eb-9117-e052f901b088.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/2c6bff89-1d07-4b92-b827-35bc07d59278.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/2ed56ab4-c930-4217-947f-3617f883669b.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/2ffd6c37-d5ad-4dab-b95c-46b5871303e9.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/3228cde5-be85-4a42-a953-456d62b9d1bd.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/33425f72-9952-4826-8625-17e57919842d.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/33b1c3b4-36f8-442e-a1e5-d0d04358e34a.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/349b2ad3-cef5-4f25-9008-d6097e351e43.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/351e2e1e-5063-46c7-b4c9-a2d5307e7ea0.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/3572034a-b249-4b5b-8c75-60b9b23b3ebc.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/3702eba4-1484-4255-9244-df4e11783363.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/382491f7-313b-47c0-9f08-12a05b3a218b.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/39ea4081-c238-4d94-a44c-fba73684ae05.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/3c60fd47-7dca-4015-bbdb-99c2ab18d16e.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/3dffea34-602b-4717-85da-24fe21daa9a9.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/3e385555-8632-40d4-a458-f9fe9409bab3.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/3e397fad-0c2c-4d99-b38c-e6294e26d552.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/4072daf2-8002-4776-8a7c-0e5c537b1f1e.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/42127cc2-409c-41b8-adac-7fc749d8af2f.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/441eeaf7-c326-4e55-8b01-7a324661906a.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/44db5ce2-1e34-4d6f-bfde-26a23571bccd.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/47ad5e41-acda-4ceb-92a9-c72ed47a7051.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/4a2751b1-cc39-4d72-8ac8-bd3b2b7a06ee.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/4e16eaf0-5619-476d-bbfe-e54b77e6fca7.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/50e0152a-f67f-460a-b19a-f678b946ddf7.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/5395779b-2210-4144-9c67-1bdb01396b95.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/58e40390-48e4-4982-a125-2050818563a0.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/59924cfa-16d2-483c-aeed-4a0d9e7b8aa3.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/5a4cbc5e-d0e7-40ce-be4f-8f5ff9f394c7.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/5a6b1bfb-7476-4f84-85eb-e9b6779f200a.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/60952d70-6611-4378-88d7-a84b2afdbd86.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/64b7daeb-9403-461c-ba57-6904e6fc4e38.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/6683f89f-f271-4dae-8b26-eafaaffb94df.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/6752d587-d8cd-4881-83ce-c868f4eeb870.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/675f7e5d-1218-415a-84a6-c5e6e6158740.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/676428eb-035c-4a63-9189-1ab8a66f6d83.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/67baea6f-3e95-41c4-86fc-3a9e18485542.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/6a49416a-b56c-4479-af51-51e78b6d1fb9.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/6aa4be47-81a1-44bd-a2ff-401881f2ef3d.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/6ebe5f94-dd4c-4619-9c9c-0410f605cc34.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/76cde6cc-500d-4151-a9cb-50dacf014e56.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/7ba78bfb-dffa-49b4-8184-38f76e66c216.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/7dadb556-0ba5-459f-ba5b-ac24b95aa4b6.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/7e1c269a-4dd0-4951-b118-7ff555e24eef.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/8347ad5d-e725-4987-917a-ac1aeeb7488f.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/852b96ab-dd7a-4705-94a6-848d79226894.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/85eeae32-6d28-4a11-8b16-c8b79cda0033.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/88469694-0d17-4095-96ea-eab345f6aefb.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/887bd37c-7809-4afa-a153-a95a07201888.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/898ceb0b-9ed0-491a-bc09-4c97899cf846.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/89b28bdb-1b74-41ec-ade9-6355ec7aaecb.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/89efb829-52f6-4098-88ea-6c1b29df89aa.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/8c1e818c-82f8-4f3e-b5d5-1ee66d228524.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/8c367426-8582-43f0-8a25-69d1cde0ae75.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/8d028cc7-fe2e-4819-a801-2e44b4f76895.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/8fbe8aac-da65-4c19-aadb-a01df03ed998.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/91cba411-0cac-4c78-b3cb-688c4434e720.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/94d549a1-6d6a-4f54-9dd6-6dd747e8d65d.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/95ab9a2d-2f17-43a5-9d82-422a7ffebb6f.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/96908f34-1924-4f07-bd17-5b2269765b35.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/9787bf5b-ebfe-4248-a826-b06bea69777f.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/97c83345-ea1f-47e6-82d8-66d204735f55.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/98553140-ecbd-40e6-9f8b-8b445dac8d08.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/9b707bf3-a68c-47a5-80ce-192b4aa19c0a.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/9ddcaa9d-2656-4302-9421-a43e7507357f.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/9e6949ae-ccd0-4d40-a82b-4cee04a601c4.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/a123cfd7-9f7d-425e-8727-765ecab93ef8.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/a2278a49-e2bf-4503-a99b-48d1c5e919bf.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/aa415e34-c778-48a2-9abf-518cd1c124e6.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/aae1ac71-f2a2-4c48-ae6f-9f195154baaf.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/aba20ccf-ee45-48cf-acde-db302ddd0866.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/ad9ad0b3-2071-40d9-a08e-63b17b54f91a.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/af165255-3886-44a8-9711-348df42cf63b.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/afa43839-4488-4e19-80b8-78938afcf9d9.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/aff477f3-6ef7-41bc-af85-2111dbfa0a87.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b1835866-29b9-4226-9dd3-92be539fd52a.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b3adb70f-ab55-4cb5-a2a0-e4962e55a22d.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b44b206d-acf0-4736-9acf-458629935708.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b49766f5-b751-4a4e-95f5-8de6ac0676b4.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b49efb56-2064-40c4-bfd5-6bc31bf08d73.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b5b46ccb-6957-4d41-80c4-567a94a4504b.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b5e4e3e8-63f5-4a57-a9e2-4900a313dde9.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b63ddc73-b2ab-4349-bcbf-efa3aec3740f.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b7369fd2-6a35-4e49-946a-e3e487e207cc.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/b9b5f70a-11b9-4172-a296-5894f2d2cd5c.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/bb4ea04f-d05b-40f0-973d-14d7e18bd3a8.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/bd200edb-b086-4ff4-9086-35ee867d62f3.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/bd9fe470-ebb8-4d3f-b04f-4c6c5c37f325.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c196c2be-4064-4a0f-9c32-c932dc83da06.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c25f764c-0986-4324-bc74-bd57813ecb91.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c2847113-68b8-4b4f-8479-71b46b0fe6bb.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c291569b-524b-43de-b51a-1e8a01d3c746.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c3a1e2d6-3b34-4e65-8402-ff6b0e83b5c3.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c8b096c9-060e-4e63-9390-ae34d72da07a.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c9228862-463e-4e13-abb6-34a462405625.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c95bb4c0-60ea-4459-817b-1fe35efb87cb.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c9629e48-b19a-42bc-b30b-973c2469f01a.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/c96f2b80-3528-4f56-9605-1ab5aa433ea2.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/ca1e92bc-429d-45d6-abd9-858401fc5e0d.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/ca9c1733-0e87-4240-b813-aeded75d99a0.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/cfbd47b5-3e07-41ac-b7c9-95bd264ff1ea.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/d14e9854-8b14-4ae1-81e0-00db4e60d9c6.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/d31700a0-b075-41c0-a09b-d2fbb4d6c32d.root'
'/store/data/Run2022C/ZeroBias/RAW/v1/000/357/479/00000/d9c41df5-3ed4-48a2-9b1d-0356d725eda8.root'
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
wget --no-check-certificate https://cmssdt.cern.ch/SDT/public/EcalLaserValidation/L1T_EcalLaserValidation/${sqlite1}/L1TEcalValidation_${sqlite1}.tgz 
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
git fetch cms-l1t-offline l1t-integration-CMSSW_13_0_0_pre4
git cms-merge-topic -u cms-l1t-offline:$L1TTag
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter.git L1Trigger/L1TCalorimeter/data

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
:'
cmsDriver.py l1NtupleRAWEMU_2018 -s RAW2DIGI --era=Run2_2018  \
  --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleRAWEMU \
  --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimEcalTP \
  --conditions=$GT -n 40 --data --no_exec --no_output  \
  --filein=inputFiles \
  --python_filename=l1Ntuple_${GT}.py
:'
cmsDriver.py l1Ntuple -s RAW2DIGI --python_filename=l1Ntuple_${GT}.py -n 8000 \
	     --no_output --no_exec --era=Run3 --data --conditions=$GT \
	     --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAW \
	     --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2022_v0_3 \
	     --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleAODRAWEMU \
	     --filein=inputFiles 

#for release CMSSW_13_0_0_pre4 and GT: 130X_dataRun3_Prompt_v1 I need to add this line SkipEvent = cms.untracked.vstring('ProductNotFound')
#replacing this one SkipEvent = cms.untracked.vstring()\March 24, 2023
var="SkipEvent = cms.untracked.vstring(\'ProductNotFound\')"
sed -i "s/SkipEvent = cms.untracked.vstring()/$var/g" l1Ntuple_${GT}.py

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
    echo "timeout 7200 cmsRun l1Ntuple_${GT}_${sq}.py `echo $args` >& l1Ntuple_${GT}_${sq}_${i}.log  &"
    timeout 7200 cmsRun l1Ntuple_${GT}_${sq}.py `echo $args` >& l1Ntuple_${GT}_${sq}_${i}.log  &
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
cp $curdir/Lumi_357479.csv menu/.
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
 echo " ./testMenu2016 -u menu/Lumi_357479.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2400 --doPlotRate --doPlotEff  --SelectCol 1.7E+34 >& L1Menu_${GT}_${sq}_emu.log &"
   ./testMenu2016 -u menu/Lumi_357479.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2400 --doPlotRate --doPlotEff  --SelectCol 1.7E+34 >& L1Menu_${GT}_${sq}_emu.log &
#  ./testMenu2016 -u menu/run_lumi.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2400 --doPlotRate --doPlotEff >& L1Menu_${GT}_${sq}_emu.log &
  pids="$pids $!"
echo "  ./testMenu2016 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu  --SelectCol 1.7E+34 >& L1Seed_${GT}_${sq}_emu.log &"
   ./testMenu2016 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu  --SelectCol 1.7E+34 >& L1Seed_${GT}_${sq}_emu.log &
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
echo " mkdir results"
  mkdir results
  tar -xzvf $curdir/L1TEcalValidation_${sqlite1}.tgz -C results/
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
tar -czvf ${WORKSPACE}/upload/L1TEcalValidation_${sqlite2}.tgz *

dur=$(echo "($(date +%s.%N) - $starttime)/60" | bc)
printf "Execution time of workflow: %.6f minutes" $dur
