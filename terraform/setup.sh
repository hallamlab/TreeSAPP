#!/bin/bash

#Run gcloud init before executing script
#TODO: add error checks

PROJECT_ID=${1:-treesapp}-$RANDOM
echo "--------------------------------" 
echo "YOUR PROJECT ID IS: "$PROJECT_ID
echo "--------------------------------"

gcloud projects create $PROJECT_ID

gcloud config set project $PROJECT_ID

gcloud iam service-accounts create terraform --display-name "terraform"

echo -n "On GCP Console, go to Billing tab and enable billing for $PROJECT_ID, press [ENTER] when complete"

read done
echo "Creating virtual machine...this could take up to 10 minutes.."

gcloud services enable compute.googleapis.com

gcloud iam service-accounts keys create ./tf/credentials.json --iam-account=terraform@${PROJECT_ID}.iam.gserviceaccount.com

cd tf

gcloud projects add-iam-policy-binding ${PROJECT_ID} --member serviceAccount:terraform@${PROJECT_ID}.iam.gserviceaccount.com --role roles/editor

terraform init

ARG="project_id="${PROJECT_ID}

terraform apply -auto-approve -var $ARG
 
