#!/bin/bash -ex

echo "Running automated Level-1 Trigger Rate validation script. Will compare rates menus using"
echo "reference and test GTs"
echo " "

###############################
starttime=$(date +%s.%N)
printf "Start time" $starttime
ARCH=el8_amd64_gcc12
CMSREL=CMSSW_14_0_6
GT=140X_dataRun3_Prompt_v2
#Prescale=Prescale_2022_v1_4_0.csv
Prescale=Prescale_Collisions2023_v1_3_0.csv
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
## ZeroBias Raw, Fill#9605 Run380446, LS1-1048, 2E34, bx=2352
filelist=('/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/e748f661-688e-42fb-927a-7defdd5370eb.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/559633de-405f-4377-afb8-1f51c2ef8728.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/8438e081-f3b4-45be-b3ca-23f75bc7d611.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/4feb03c7-5cc8-4292-9adc-8f37bcfcf82b.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/66a210a6-e815-44d1-b81f-b01bd4ff6b9c.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/e6c4fb3e-8047-4209-8cf8-3e0ca811b004.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/0b79d9c9-d04e-44c9-8c48-2366277e662f.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/d17a7cf4-d63b-40d7-9acd-2960e56e0694.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/061e98d5-95b1-4df9-98d6-fe98474eb785.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/f2eb7387-8236-4520-8623-433bdcb2139f.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/7a33956c-154c-4ebe-87f9-23d740e3739f.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/8d3e9e1c-9a52-409a-9c0e-f04604ea681e.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/66023277-5b00-4b25-9130-2cc2ee1e2d25.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/f9476a38-1a88-4b38-9357-4a26fc450980.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/5827a1f3-f41a-443c-aa2b-f1370eb1019d.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/27241296-d72d-421b-b0f4-914c56e52377.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/3754c9da-48c4-4b2d-8276-82e66cdc5cf6.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/eeab22dd-138d-40d5-9aa0-e2799881f4f4.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/dbcd3e5b-1826-4eb8-a295-01b2c5c2c020.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/3b48f678-ec59-4253-b44a-2c5d01e69b95.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/8107352f-87e4-4f85-ba4e-82db8f2c6ab8.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/a422179d-f4d7-4dc6-81a6-767e9c9fa07e.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/c4a4b5e5-1dd0-4684-9a08-d918838e06ce.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/f8e2d478-ff2c-4bfc-9277-0eaf9ff8d007.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/ceed5b4d-4c3c-442d-8e56-e58164764e74.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/d18d33c7-15df-43a9-b5db-f9d9ac7312de.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/9d1e25dd-c2db-48a9-b297-973c9dc4e285.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/554d90c5-d80e-4214-a462-1525a7998b02.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/c4377f61-3188-4b1e-96a7-4c658ecc1ac3.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/3a1f8726-cb4c-4662-8846-1bbd20fdc742.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/e8fd6585-5d6f-41d6-9e03-d52bfd21a397.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/c587e421-41b6-4fd7-b235-0e9abee750a1.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/721c92c9-4298-4b59-8272-9c39fbf4ddea.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/0ac34e3f-c4a1-433d-a571-2ebb6a124e1c.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/85412681-7ba5-470d-9583-a2fb2cf29c0e.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/5a79e6a1-cba3-4539-b5ba-d285bbd09335.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/ec4976ad-ba98-43ec-8abb-3ae0bb936455.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/a4fdc96f-4f6e-47aa-99c9-911f9961ec23.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/2b9fbb6f-05c0-4963-ac42-8e2365a7e677.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/2d2bba85-e291-4278-8b33-8be850eeb1c1.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/3349fd69-7394-4f85-b60f-cdb50ac81c21.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/909aeb69-475c-46be-a3b8-de3d6a40fb5c.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/89061b20-89af-4950-bda9-3565c73758b0.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/f8da3f7f-6a50-4569-89b5-d314540bba50.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/578a1a71-4fe1-4e3d-a375-6a8d0c694066.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/49e9fed7-664d-4814-af79-37cbb071e5f2.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/c1f501f8-de5d-41b9-ac31-2d904c65b7d0.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/1a1d25da-96c9-44b9-a716-19f46eb8ed79.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/d0682c42-abcc-477d-a5d8-03d9eba86acb.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/5feecb48-0f2d-49cb-94e3-ec97212071be.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/6c234cec-ee96-4029-a3dd-e111ba0bbbde.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/b9601f47-3420-43bd-94cd-af6e1b5859ab.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/d6e09d32-c05e-4776-aaa5-6b913d82d76b.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/bc2620a9-b8fb-47f9-bdf8-c1b87f955fae.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/d5d7fa83-385c-4622-84f2-039fa158fcf7.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/549dc3cc-01f7-4ba4-8094-245f44520a57.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/37d8badb-1d23-4958-8752-ccb4d11904b4.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/4537e1a8-1ee1-4535-ad5a-3da8a411cf62.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/1e2885a6-4436-4338-b1e4-544bade4a63b.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/d7fe512b-6e6a-4670-85f8-b0ce6e4778a5.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/e7a4be02-ed46-4d41-86ae-d49b14bf6bf1.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/c1aca04d-b0bc-4ebe-9f50-ca42f1c5e1cf.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/62adb24c-de0d-4c3d-8a29-c7353d05f0bd.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/8c5bbf8b-a36b-45b6-8068-a55729bcfa14.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/400d119c-fcd7-4f5b-9d83-d445150dc730.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/bf0f94af-a5e9-4973-8901-687c198271d2.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/74455cbd-0aea-423c-bf3c-bf896322c65c.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/c1b232ef-69d0-4199-bd4e-78e06dd9a9ba.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/573ee606-4ba6-4417-87c7-246c2c230f6d.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/ae0690ec-1420-4ac8-a3bb-a1d5b9876444.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/b26cb766-c0eb-416a-a3ab-51bc4a10e877.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/9fa6e400-e108-46e7-8b3e-ba0cc55fdcd9.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/cd7ddf33-2a3c-46ae-8bf4-58771d2580ee.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/da8cb6d8-31e9-4dec-9aa9-0b78cd3097ab.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/37056fd8-89be-4aa2-be66-e6f7c8f135df.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/15679204-21d0-499e-b1fd-a90c4e176822.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/9dc8f1b0-0696-4752-9842-bfd8fa86789b.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/9b678242-10af-4ca0-b398-65bfc409e176.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/50b94785-13b7-4e4f-b67d-cba07cf489b5.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/0696229a-8b68-4ee7-b2d8-12bbdc7f2661.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/79ec2c86-279d-436f-a0d9-b44b9c30a9da.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/a514b2e2-51d2-4605-97e0-f0978870bb40.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/baed664b-799f-4895-acbb-d27772b97fb4.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/b5f4f34a-493a-4f6e-8908-512e8f5695fc.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/5ab3faad-49fd-450a-8622-53b54c2fc678.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/9a899e3b-1c0e-4b87-b9a3-92d4117196c9.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/e9f8f8b5-e431-4b3e-85bf-579a8c2a4a13.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/eb0ac7a3-dcc6-464b-b6d1-68d80510e03c.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/316af7c6-87a0-4bf0-a916-99518dbf30a1.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/8691d395-3caa-41d7-8c85-d302a4d289c2.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/8e7d3bd1-1cc3-40e8-b819-d9338d85f392.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/c0fd4460-b2e3-43ba-aebd-72375c4d85ab.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/7ceba859-dc64-4ac1-b1fb-7865f45769be.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/b5dade44-d32b-40d7-8be8-b89a7a7a147d.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/29610c04-840e-428a-8eec-07ffd83b15b7.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/81ab5980-566e-4a36-b763-4de78bb78fed.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/d8c7813d-7b99-4cb2-bedc-41f9dd5467f5.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/b3e52f84-582a-4dcd-b15b-876dc1cf63cc.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/3552ded9-1fa4-468b-a30b-49c4ec049bf3.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/d48660ab-7460-40ca-b409-73bf59ae4368.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/45eeeafb-1f81-490c-8c5f-72c364230d90.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/a6d89bac-fbc1-4c2d-a3ef-ddbb11cd578f.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/d4699c6b-273f-4f70-9d9c-20c3c5f4469f.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/eae5297c-0d01-420a-b75d-0345848fb4a9.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/74fe6877-a584-422e-a19c-65143ca9c292.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/7fc535eb-7c54-4848-9ed1-b806aa08bfbd.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/662d6692-1c88-440a-9f25-bfd78cef3ead.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/b5e1a9b0-3bc2-4670-9408-ab1af7cc6378.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/173503a1-daf1-4aac-b193-0fd2e5a9e267.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/6e425902-df3e-47f8-a0d5-1931a612f307.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/1ea36ed0-7e60-4aea-ac22-e044ce32baf0.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/80a8eae2-1853-487b-aa47-6f01ecfa212a.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/17d9f711-4fd3-465e-8afc-d03424f36843.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/e38d33ee-aa6a-4132-85e5-dd6562825ef2.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/edaa4c29-693c-46f2-918c-76590e9011f0.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/78854ed7-e5a9-4ce5-94cf-90e266df7067.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/92723561-7de9-41c5-8f35-6f14b046f508.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/ee4e1f78-4307-45ab-be36-b8f8d4d670b9.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/416c5d86-3594-4410-9655-eb85a1a3240c.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/741dd9e5-39ee-4269-92e9-c794c486eef5.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/69bfe3a5-284d-411f-978a-f61ad38348cf.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/e085f0b5-c7f9-4dcf-8c88-0dbb4768806e.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/34de7e36-92b8-4442-9051-ed44bde22711.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/7ddd3607-44e5-4fd9-8d12-1b861c1c0ace.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/7b2adffd-e2de-4481-b585-e1cb62cdc92d.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/9969f4ab-5c25-447a-90df-c68e8d738423.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/7a036e25-beb4-4503-ba17-3eae87729c47.root'
'/store/data/Run2024D/ZeroBias/RAW/v1/000/380/446/00000/c7bf9e6b-a2c4-49c2-9652-f6e0a4ca2150.root'
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
git cms-addpkg L1Trigger/L1TNtuples/python
git clone https://github.com/cms-l1t-offline/L1Trigger-L1TCalorimeter.git L1Trigger/L1TCalorimeter/data

sed 's/+l1UpgradeTfMuonTree/#+l1UpgradeTfMuonTree/g' L1Trigger/L1TNtuples/python/L1NtupleRAW_cff.py
sed 's/+l1TauRecoTree/#+l1TauRecoTree/g' L1Trigger/L1TNtuples/python/L1NtupleAOD_cff.py

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
	     --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleAODRAWEMU \
	     --filein=inputFiles 

#for release CMSSW_13_0_0_pre4 and GT: 130X_dataRun3_Prompt_v1 I need to add this line SkipEvent = cms.untracked.vstring('ProductNotFound')
#replacing this one SkipEvent = cms.untracked.vstring()\March 24, 2023
#var="SkipEvent = cms.untracked.vstring(\'ProductNotFound\')"
#sed -i "s/SkipEvent = cms.untracked.vstring()/$var/g" l1Ntuple_${GT}.py
var="TryToContinue = cms.untracked.vstring(\'ProductNotFound\')"
sed -i "s/TryToContinue = cms.untracked.vstring()/$var/g" l1Ntuple_${GT}.py
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
  python3.9 ${curdir}/ModifyL1Ntuple.py --globalTag $GT --sqlite $sq
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
 echo " ./testMenu2016 -u menu/Lumi_380446.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2352 --doPlotRate --doPlotEff  --SelectCol 2E+34 >& L1Menu_${GT}_${sq}_emu.log &"
   ./testMenu2016 -u menu/Lumi_380446.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2352 --doPlotRate --doPlotEff  --SelectCol 2E+34 >& L1Menu_${GT}_${sq}_emu.log &
#  ./testMenu2016 -u menu/run_lumi.csv -m menu/$Prescale -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Menu_${GT}_${sq}_emu -b 2400 --doPlotRate --doPlotEff >& L1Menu_${GT}_${sq}_emu.log &
  pids="$pids $!"
echo "  ./testMenu2016 -b 2352 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu  --SelectCol 2E+34 >& L1Seed_${GT}_${sq}_emu.log &"
   ./testMenu2016 -b 2352 --doPlotRate -m menu/Selected_Seed.txt -l ${CMSSW_BASE}/src/L1Ntuple_${GT}_${sq}.list -o L1Seed_${GT}_${sq}_emu  --SelectCol 2E+34 >& L1Seed_${GT}_${sq}_emu.log &
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

python3.9 CompL1Rate.py --globalTag $GT --sqlite1 $sqlite1 --sqlite2 $sqlite2  | tee ${sqlite2}.log
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
