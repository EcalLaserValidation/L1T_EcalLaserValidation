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
## ZeroBias Raw, Fill#8128 Run357479, LS1-945, 1.7E34 bx=2400
## ZeroBias Raw, Fill#8496 Run362760, LS1-728, 2E34, bx=2378
filelist=('/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/00247964-65bd-4b87-a54d-f42f4cbf1961.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0098b31e-e448-4927-81c3-210f4cc1edd0.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/00fc59b3-9e4b-4c71-b021-f84b55683987.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/01628edd-1ea2-41c3-8967-750a8931d873.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0188d0de-e836-4bca-881f-73f08b74ee25.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/02bd66ea-3d91-484e-8cda-6b4ccf46a522.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/02ea9d74-bcb1-42b3-926d-7da5e458e039.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/039202d5-1668-4229-884d-7634f5e7fe2c.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0722791e-bcfc-4404-8455-ec47a329c5e5.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0884a6c8-bf3d-4e50-9a01-e7520cb17ea0.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/09374b70-5cec-44d2-9678-b953fe7a2778.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0a20f8f7-4661-4663-aa73-6e33c861660a.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0a21566d-e1f4-4768-af0e-9706fd288789.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0c851fb9-104d-4086-9cfe-285d9193dd49.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0cbd6dda-7400-4407-8a12-1ac21306049b.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0d4c0c0f-5cb6-46b1-8139-92b8eff3b98e.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0d599627-1a85-400e-8ec4-82e896ede93a.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0e2d5183-232e-4358-acfb-ad26ac078b55.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0f284ef5-a4bf-42b9-8530-0a469045a0d4.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/0f297e79-e1fa-482c-a4fd-58f694e78f27.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/102f7d05-1523-4000-a179-4619bd183c68.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/10a5c99a-a6bd-4977-81f6-1ddf867b4d34.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/10aa2ca3-5f2e-4e7f-b2f9-530674b292fb.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/11f107de-343e-447b-b19a-76170dd7fd02.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/14148a10-94c0-4c7f-9976-a658431d6862.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/14927f06-e94e-4c32-9ccf-ea49467bbd22.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/156ac978-584d-4c0a-8982-2014ca5d0b39.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/16b08a5b-46f5-4ccd-8950-52dacbd8cbaa.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/17f74ecc-4a1c-42a8-9889-468d811b534c.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/181518ad-130b-4605-8545-06608a0e5661.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1907e0bf-0abf-4451-b13d-b0abe7746b0b.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1a07c522-16b2-4c1e-b218-a65fac9dd5ef.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1c6b1276-b033-435b-b8fc-7b6de2a28b20.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1cb86dc8-4378-480b-bccd-589a2ada79be.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1d79e7de-5044-42d2-9458-08f5216f9d2e.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1d80b1b2-115b-477a-94c6-ab1cd9741c3b.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1dab03e2-be53-4069-8b2b-d8597e167706.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1dda88ce-1ccb-42ec-b4da-53be8b15f734.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1f0dcf27-9799-4a40-8bf8-4cf937546416.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/1f9da000-6ae8-4b99-b819-1e9591e5e034.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/209c084d-ced8-4baa-9d36-962d21bd31ae.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/20c85896-81b7-4def-8735-bda94f8abeaa.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2186f3a0-635b-475b-9ba2-4cf82062a48e.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/21a393a2-3914-478f-9efc-401aa4773a5a.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/22719315-9342-4c3d-8712-be718c7f4f43.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/22b72030-304b-48eb-b0e3-4afa4330c12b.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/26d10446-43e7-473c-8835-c7b3944a6d33.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/27922fa9-0985-452b-9905-0bdf121a9376.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2a595349-2a0b-4ada-93d0-428f1a5629bc.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2b7bcf54-7a76-4181-8a32-a7d2f74d091b.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2b8874f6-ea9d-40b6-bb25-67a9170f9d5e.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2bb5353f-59fe-4ea9-aa35-28f29d7fc563.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2cb74bca-980e-43fe-9add-a40738edfbd3.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2cbc18a1-3cb1-4ef0-989b-b833b8ce5b80.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2d8a45a9-e885-4f80-b924-72dc78ea18ed.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2dbfeb37-58be-4e31-8b75-5eb33938e878.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/2eba975c-d599-4a50-9890-39c70b398dc2.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/319ce707-38ad-435b-91d4-5d1b4602fd09.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/31a8c8f5-eece-4cac-b8b2-772900b703c3.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/322d22ea-21cf-4260-921c-61589159faf3.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/348962b8-f973-4d2b-8146-f3811be62af1.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/34b0f2d9-e9bc-4a38-9dfc-d4c81b02627a.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/3566fc28-6b21-42ae-b725-8ec50934db95.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/36566a14-84bf-484d-a768-ccf56e4eaa89.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/367fea14-4066-4800-b880-6ab3a597dafe.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/373f6b9b-59f8-4e93-b984-1033cc0ee705.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/37bebac3-8cbe-4dd3-96c4-af340a6f37ed.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/3a7ff0eb-3a27-4c5a-b7c5-b9b5924fe9dd.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/3a85c382-9b34-4428-b58b-5158f73e01e6.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/3b2512b6-b29e-4750-af9b-d8ef806ba890.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/3f0ad04d-a5e1-47fb-8b44-87c605305301.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/3f0f9208-f974-424a-866d-d3ec1d67290f.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/3f315f9c-d5e1-4f63-9afc-b445d4f956cf.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/455cb47e-e9cf-4611-b3a7-9abce0880dd6.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/45684a29-d6e7-45d2-a764-f43c600690bd.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/456c40b1-4592-4851-a22c-c12ef19b131e.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/45843661-cb50-46a7-a6cd-c5721366eb1a.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/45b12c37-8005-4550-86ae-d80929e93ca3.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/46a88b32-995c-401c-a86b-4ecc929e9e9a.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/47514399-4a91-48c7-a06d-897c664b699d.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/49aecfdc-9db0-4f0b-ad74-152dee8b62f9.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/4a232d70-fd15-4f8a-99ab-e314b005dc69.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/4a483b32-b0dc-48d3-9963-a087ce41e280.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/4ba1f592-4945-4c32-9119-18712077f8f0.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/4bbf7e5c-d44b-4781-a481-87c448ab1c0f.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/4d144b85-a40b-40d3-b2b6-edfc05646e77.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/4e87c924-a973-4d01-ac09-cfffeedfe6b6.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/51c436da-3eee-47dc-9fb2-deab835eca51.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/52043393-38aa-42dd-be94-6ae88d004103.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/5210a4d2-af9e-459d-b462-f7f5f70a1007.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/53cd3f09-a369-4bd7-a6b1-dffb23c736ac.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/589ef122-8e0e-4906-bc08-ce53962af746.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/58ca5f93-4c5f-418c-b595-6d53bc614450.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/5a1d0498-56c3-44e7-b73b-4ab620d32e8e.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/5afcf3c7-d644-4853-9c4b-6a29972f059f.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/5b140b12-0ea9-4c52-82b9-455268cd4260.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/5cc82c3b-b22d-4351-8d8d-0a516fb50752.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/5d36e709-2509-400c-acbf-2a7cc86f6add.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/5f139a21-ef5e-4ca8-81d6-c02831a83063.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/610ea117-552b-45b3-bc4d-24be7bc0d676.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/6319ab0c-4ef5-4084-96a1-0fb9bdb97d74.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/648389bd-3bbc-4037-99ba-2b608d9e2795.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/65e4a825-a513-4609-9e9f-683c9f48a7ac.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/6686bd17-1cce-4e88-8296-8e06e7b10d64.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/671f662f-29eb-470a-a137-28bf6ed932f3.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/6969f8bb-3149-415b-85fc-eb8542b2d1e1.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/6e9125dc-37c6-410a-8e0c-8242ab31d09a.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/6fd638cf-10c8-43a1-8032-8fd0dce2796d.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/738e270d-150e-4ef2-9349-12c71dd09a05.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7564858b-f440-4451-9373-ae255ef39f35.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/75824daa-059f-4aa3-8850-a9be4fde3bc7.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/76bcfd36-9611-4aa6-94ad-e567b67e39b3.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7862b33a-4b1c-4dc6-a9d8-349f1f6868d8.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7a10cb02-092d-4950-b3c4-0b624cde25e8.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7a7d74b9-efa9-4144-a436-b6e0acb3b82b.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7ab65aa1-bd75-4a88-a3a3-4278ad52377e.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7afd215f-b22e-4afe-9c8c-288555194a29.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7b2fc94b-01c5-4aaa-a758-2574e32c8133.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7b92e99b-c442-42a7-9f34-9421072f1f73.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7b9a202e-2f86-4775-9f17-703cc4889f06.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7bb27814-6557-41ea-9d89-7f33b67883a9.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7c379ce7-a531-425d-90b9-892880ecfb4a.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7c3a319a-6205-46d0-b292-44b1cd9d6a27.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7ceca95c-b4b3-4f0f-a900-d8b76763e2ad.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7dc3836b-ec56-48c7-9427-8320f0e2fbad.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7e8dfb4e-b3b0-4334-ad88-8895f746c23e.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7f53a236-2195-45b7-9514-5d98cd5c22f8.root'
'/store/data/Run2022G/ZeroBias/RAW/v1/000/362/760/00000/7fc8e2bc-3205-4c2a-9fe1-e0abdf61509e.root'
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
cmsDriver.py l1Ntuple -s RAW2DIGI --python_filename=l1Ntuple_${GT}.py -n 100 \
	     --no_output --no_exec --era=Run3 --data --conditions=$GT \
	     --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAW \
	     --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2022_v0_6 \
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
NJ=1
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
  cp $PWD/L1Ntuple_${GT}_${sq}_*.log ${WORKSPACE}/upload/${2}/.
  #cp $PWD/L1Ntuple_${GT}_${sq}_*.root ${WORKSPACE}/upload/${2}/. 
  cp $PWD/L1Ntuple_${GT}_${sq}_*.py ${WORKSPACE}/upload/${2}/.
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
cp $curdir/Lumi_362760.csv menu/.
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
 echo " ./testMenu2016 -u menu/Lumi_362760.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2378 --doPlotRate --doPlotEff  --SelectCol 2E+34 >& L1Menu_${GT}_${sq}_emu.log &"
   ./testMenu2016 -u menu/Lumi_362760.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2378 --doPlotRate --doPlotEff  --SelectCol 2E+34 >& L1Menu_${GT}_${sq}_emu.log &
#  ./testMenu2016 -u menu/run_lumi.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2400 --doPlotRate --doPlotEff >& L1Menu_${GT}_${sq}_emu.log &
  pids="$pids $!"
echo "  ./testMenu2016 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu  --SelectCol 2E+34 >& L1Seed_${GT}_${sq}_emu.log &"
   ./testMenu2016 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu  --SelectCol 2EE+34 >& L1Seed_${GT}_${sq}_emu.log &
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
