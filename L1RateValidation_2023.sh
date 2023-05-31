#!/bin/bash -ex

echo "Running automated Level-1 Trigger Rate validation script. Will compare rates menus using"
echo "reference and test GTs"
echo " "

###############################
starttime=$(date +%s.%N)
ARCH=slc7_amd64_gcc11
CMSREL=CMSSW_13_0_1_pre4
L1TTag=l1t-integration-v156
GT=130X_dataRun3_Prompt_v3
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
filelist=('/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/00ea2980-5bbe-4e46-b5b8-12a6eff0be09.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/01ee253d-c826-4f2e-9670-a5632e72ae05.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/05aa244a-9aff-4d4e-a6d0-7e9434262adf.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/07cd29fe-59a0-4d91-ad97-459aa161f837.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/08f64a48-139e-44c0-ac45-fd7f7d2d1cd0.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/091e0230-0c16-4155-ad05-3e68bec128d7.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/097a6001-fd4e-4bf9-87ce-2882aa678323.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/09a17724-d2b8-4d5a-84a8-7888e1a205d1.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/09a17967-fd8a-451b-b6e9-660505fadb38.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/09ca0294-8100-4887-9181-620765f1d9c2.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/0af71738-204a-4699-a1a5-ca9281f89daa.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/12b54ccd-7368-4a5f-8839-1d7294b55334.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/16f81a41-4b08-45f5-9b8f-8ad17a9916bd.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/1f37bfd6-bb18-4009-a4cb-1f474dd37848.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/2503ae1d-78a4-400c-a102-4be31990e13d.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/2837e783-7baf-42c9-bf4c-ed931202fe46.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/296b66f2-24c3-4e3e-af59-3f52068d8c40.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/2971df46-5588-4234-b0d1-8c28fdca3697.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/2c4deb7a-dc34-4c78-9291-ac5d89d26d1e.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/2d403e84-0533-4bc9-88ce-8489c11cb621.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/2ed835c0-761a-4f1a-bfb3-13eadb0e6adf.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/2f40f3e2-b943-4197-9cb0-b014a525e107.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/2fab9d90-741d-4710-a627-c99b066357f1.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/3021bc6b-1666-49b0-a64d-2a5f12c75885.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/30cd235b-2f2e-4027-ae53-f5761879f928.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/327f3ec2-ef20-4752-87e5-868c9a56da20.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/33c916de-34dc-47b5-b1e5-1646c876b6ff.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/37a4232a-84fc-426e-880e-7998045f884d.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/3a6a6161-0015-461c-a465-4ff75b55698e.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/3adb2ebd-fea1-4bc7-828c-9b32153b36b6.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/3fe2ea7c-cd3f-4f1f-ad77-b04571e35a7e.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/4097201b-d20c-457e-94f9-054956fb3f06.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/41b8fa44-4569-4d65-9c54-8d13c2e1e85e.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/44a22c3f-065c-4fa2-84a5-6c8b3d075ebd.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/44b701fc-5a1c-4d71-ab6e-53461a0d5fde.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/46b45bef-1e7c-40cc-afe2-1b535a7ba924.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/4889b9a6-7874-4a73-9031-f1ad39397cdd.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/4b23d91c-acde-4a5a-b32c-ce258f0baf2c.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/4ea0326c-3bf0-4692-a1e6-e2612958d242.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/4f067b9e-1e2e-41cc-9858-3af9fcf50f76.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/4f2b9387-7eaa-43a6-869e-0b86ec890f95.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/505652e3-3e12-45b8-8ac3-ca3631c7bfa4.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/51f90dfb-650d-4851-b501-b2d07b29ce0b.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/54a3e9c1-9310-4b95-9c83-7574b9d51bda.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/55784c0d-a543-455a-bf25-68c44b99accd.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/574beacd-b663-4317-a8d1-ada2307d2177.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/578dd6f0-5050-4c2a-a66c-aedf1471b415.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/593c5071-e4d0-4519-9d46-11e9a2467b9b.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/5d293cec-1e82-4a47-8bc8-da4b77fe63a7.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/5d86f95b-9a01-49e3-bfe5-443a7c4e4054.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/5e02816c-c029-44ba-9d58-f43b05b45577.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/5ebe7e21-f179-4b3a-a0fb-f61052341dd0.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/6033cb52-b30d-4c8a-ab6b-a0cb1af92ece.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/607ed774-6669-496a-a7f0-976c0bc95259.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/608ff8dc-5a3f-4297-a948-aabc4c7c8d13.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/60bd87b3-e851-4fb1-a9a3-d6a7e36834c3.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/65d5e583-a027-46b4-a91d-82348bdd403e.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/682834e1-2ae4-4773-81cc-9b13cf043586.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/6901f49e-e856-491a-8a01-c4d28c35de3e.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/69f3aed8-71cd-4e89-be27-187ad1d1f975.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/6a95358e-8e5b-4a2e-afc7-42ab252bbbaa.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/7364feb4-fc5a-4ada-a119-9626c2f9bbd4.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/754411f5-2b87-44ea-9f1b-e07310a7fa0b.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/78088188-a459-40e5-91d7-1ee50804f8bb.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/7854d324-5af3-4a29-9f60-6a41af421235.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/7b37d169-45b0-4b90-8c48-c0ebbc055dbe.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/7b9dcc99-3514-4e89-ab60-b37109fdcf04.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/7d3374bc-6eeb-4709-ae2e-b44618bb7c97.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/7da1bfe5-da79-46a9-b07e-3f6ab0e6c94b.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/7ff96d2f-5ac2-437a-a30d-ea55c9841fdc.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/7ffdc7b5-8cc1-452b-a0b2-50c5a34a12d9.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/80938066-db6d-4367-b962-5d2f93185986.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/81bc24a5-3c37-494b-81fa-fc028deca1c4.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/823a7795-1d28-4cbe-bf37-db2479e8e204.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/830e52ed-7a4c-4943-939b-983cc5e771ac.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/86991493-673e-4e6c-a6df-df6d04383942.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/8806a5c5-0586-49a1-bd07-047330ee6841.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/8c96e9e7-94ab-4925-9ece-546126116378.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/8f733cc9-1cd2-457f-81f6-c02f0d4c4889.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/914963d8-706b-4416-9ef4-bd03ed7d2eb7.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/9174e3e6-f03f-4b1c-9f28-e4b49e1503db.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/9358fb72-ab78-4446-b6c3-63e0e0945b32.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/945a3f8a-48c4-4284-b174-f351f04171e0.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/968f8628-08f6-4d4e-9bc3-aa499885ea7c.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/975583da-794b-4856-9775-2611266ce024.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/978c97f7-c93c-42be-8baa-8072e391f73b.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/9938033c-ea6c-46d9-bd93-cbeb1040ddb2.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/99400733-ccdc-41c3-85f7-c31f600dad26.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/9a773c7e-d49c-4bed-968d-b77628d6d5b2.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/9b82dbfc-493f-47f4-91c3-57d63e98add2.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/9bb8d45c-b59e-4dfa-b05b-054570020489.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/9f9aa7d1-b4da-4219-8ca7-d1a3ccd39a2f.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/a196b180-8f92-4ed6-8e1b-f39636e907b1.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/a281390b-4799-43ea-8d79-a7fee9f5caaa.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/a3384e86-3967-4acf-8d4f-3b3c4e832142.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/a5f5f6fa-0efa-492a-b7a8-2ee90d726819.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/a9b41f24-45a7-41fe-a481-95808ae285b7.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/aa097a38-83b0-49fa-b4d5-156280329686.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/ad44e579-c548-4398-afe8-33f533a6ea7b.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/ae91192e-014e-472c-b715-92917e88626c.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/b3ab4b21-5c18-4904-8655-bf9d50dba14f.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/b7237a7e-81b6-4a41-84d7-e5a9b3529bd7.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/b7597c01-7d39-41c3-904b-e5aceb6a5dba.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/ba2fe34d-48d5-4b38-9767-7778b19134f4.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/ba52ce13-2a35-43e0-af22-c2df9fbfa94f.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/baa36af0-2bda-4e1d-9297-a0c83db219e1.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/baddf605-0b66-4c84-8abf-5c70e2013446.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/bb19df33-08f3-4fca-be42-801f9b599baf.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/bb69a755-20d2-446a-87f8-e3e040aef4ca.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/be89a3b0-32bb-4a73-b773-612649405bce.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/c00afbfe-1444-481a-8c2f-236c1d1787ab.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/c00e94ed-f91b-46ed-9c02-4517a068f961.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/c2b68a7a-eba3-437b-bb60-6cd71571dd02.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/c3db90e5-372e-454c-ad0f-b6f33d242e20.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/c4c98dbb-87f8-4090-8c0a-9fcb7c4bd20d.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/c5051fbb-2e7b-4d1a-ba7e-aebedfa18155.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/c64d5d94-8904-4fdc-8e01-78ff66908a89.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/ca3f9241-a04c-4ce8-ba9d-2b3181bf5d52.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/cb2666e8-e35a-45bd-9555-d3717d2eb014.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/cbe7ed5c-f944-4ab1-9f05-88a7a94f131f.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/cca91064-2a0c-474e-bf60-60dff8ffe6ec.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/cd1015bb-46c2-48da-af83-8bff119a8e8f.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/cee28f78-74ee-4cac-91c7-8a60f9cd0707.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/cee54c1f-097b-4dd9-a25a-20e22ba7dc42.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/d02b9a6e-a495-4795-aac2-509b624eac98.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/d0f0c32a-ce93-4337-b19b-a4573d31d547.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/d4b59d63-0e5e-4e5a-8b0d-7219fb209e15.root'
'/store/data/Run2023C/ZeroBias/RAW/v1/000/367/758/00000/d84fac77-86c2-4950-a1dc-0c291d197c54.root'
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
#cp l1Ntuple_130X_dataRun3_Prompt_v1.py $CMSREL/src/.
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
#now trying this with file: l1Ntuple_${GT}.py from the git repository 05 april 2023
cmsDriver.py l1Ntuple -s RAW2DIGI --python_filename=l1Ntuple_${GT}.py -n 4000 \
	     --no_output --no_exec --era=Run3 --data --conditions=$GT \
	     --customise=L1Trigger/Configuration/customiseReEmul.L1TReEmulFromRAWsimEcalTP  \
	     --customise=L1Trigger/Configuration/customiseSettings.L1TSettingsToCaloParams_2022_v0_6 \
	     --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleAODRAWEMU \
	     --filein=inputFiles 

#for release CMSSW_13_0_0_pre4 and GT: 130X_dataRun3_Prompt_v1 I need to add this line SkipEvent = cms.untracked.vstring('ProductNotFound')
#replacing this one SkipEvent = cms.untracked.vstring()\March 24, 2023
var="SkipEvent = cms.untracked.vstring(\'ProductNotFound\')"
sed -i "s/SkipEvent = cms.untracked.vstring()/$var/g" l1Ntuple_${GT}.py
#nevents=4000
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
  #cp l1Ntuple_${GT}.py l1Ntuple_${GT}_${sq}.py
  #sed -i "s/%ievents%/$nevents/g" l1Ntuple_${GT}_${sq}.py
  #sed -i "s/%iov%/${sq}/g" l1Ntuple_${GT}_${sq}.py
  cp l1Ntuple_${GT}_${sq}.py ${WORKSPACE}/upload/${2}/.
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
