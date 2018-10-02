#!/bin/bash

#Run gcloud init before executing script
#TODO: add error checks

cd tf

PROJECT_ID=${1:-treesapp}-$RANDOM
echo "--------------------------------" 
echo "YOUR PROJECT ID IS: "$PROJECT_ID
echo "--------------------------------"

gcloud projects create $PROJECT_ID

gcloud config set project $PROJECT_ID

gcloud iam service-accounts create terraform --display-name "terraform"

echo -n "On GCP Console, go to Billing tab and enable billing for $PROJECT_ID, press [ENTER] when complete"

read done

gcloud services enable compute.googleapis.com

gcloud iam service-accounts keys create credentials.json --iam-account=terraform@${PROJECT_ID}.iam.gserviceaccount.com

gcloud projects add-iam-policy-binding ${PROJECT_ID} --member serviceAccount:terraform@${PROJECT_ID}.iam.gserviceaccount.com --role roles/editor

terraform init

ARG="project_id="${PROJECT_ID}

terraform apply -auto-approve -var $ARG
 
