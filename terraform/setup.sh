#!/bin/bash

# TERRAFORM AND GSDK INSTALLED CHECK
if ! [ -x "$(command -v terraform)" ]; then
   echo "ERROR: terraform is not installed or added to PATH..." >&2
   exit 1
fi

if ! [ -x "$(command -v gcloud)" ]; then
    echo "ERROR: Cloud SDK is not installed or added to PATH..." >&2
    exit 1
fi

# VERSION CHECKS
TERRAFORM_VER="$(terraform --version | head -n1 | cut -d"v" -f2)"
REQUIRED_TF_VER="0.10.0"

if [ $(printf '%s\n' "$REQUIRED_TF_VER" "$TERRAFORM_VER" | sort -V | head -n1) = $REQUIRED_TF_VER ]; then
    echo "Terraform version $TERRAFORM_VER..."
else
    echo "ERROR: Terraform version should be at least 0.10.0, please update!"
    exit 1
fi

# GCLOUD INIT CHECK
if [ -z $(gcloud config list account --format "value(core.account)" ) ]; then
    echo "ERROR: GCP account not found, use gcloud init to sign in.."
    exit 1
fi

PROJECT_ID=${1:-treesapp}-$RANDOM
echo "--------------------------------" 
echo "YOUR PROJECT ID IS: "$PROJECT_ID
echo "--------------------------------"

PROJECT_ID="treesapp-2948"
USER_NAME="$(whoami)"
 
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
NAME="username="${USER_NAME}

terraform apply -auto-approve -var $NAME -var $ARG
 
