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
    echo "Using Terraform version $TERRAFORM_VER..."
else
    echo "ERROR: Terraform version should be at least 0.10.0, please update!"
    exit 1
fi

# GCLOUD INIT CHECK
ACCOUNT=$(gcloud config list account --format "value(core.account)")

if [ -z $ACCOUNT ]; then
    echo "ERROR: GCP account not found, use gcloud init to sign in.."
    exit 1
else
    echo "Found GCP account: $ACCOUNT"
fi

PROJECT_ID=${1:-treesapp-$RANDOM}
PROJECT_NAME=treesapp
echo "--------------------------------" 
echo "YOUR PROJECT ID IS: "$PROJECT_ID
echo "--------------------------------"

USER_NAME="$(whoami)"
echo "Will install TreeSAPP and dependencies under $USER_NAME user on instance..."

# ------------------------------------------------------------------
# Create new project, service account, instance
if [ -z "$1" ]; then    
    gcloud projects create $PROJECT_ID
fi

gcloud config set project $PROJECT_ID

echo "-------------------------------------------------------------------------------" 

# Enable billing
echo -n "In GCP Console, go to Billing tab and enable billing for $PROJECT_ID, press [ENTER] when complete"
read done

echo "-------------------------------------------------------------------------------" 

echo "Adding permissions, this may take a few minutes..."

SERVICE_ACCOUNT=$(gcloud compute project-info describe --format='value(defaultServiceAccount)')

gcloud iam service-accounts keys create ./tf/credentials.json --iam-account=$SERVICE_ACCOUNT

gcloud services enable compute.googleapis.com

gcloud services enable iam.googleapis.com

# Check if instance with project name already exists
if [ -z "$(gcloud compute instances list --format='table(name)' | grep $PROJECT_NAME)" ]; then
    INSTANCE_NAME=$PROJECT_NAME
elif [ -z "$(gcloud compute instances list --format='table(name)' | grep $PROJECT_ID)" ]; then
    INSTANCE_NAME=$PROJECT_ID
else
    INSTANCE_NAME=treesapp-$RANDOM
    echo "Instance name will be $INSTANCE_NAME"
fi

echo "-------------------------------------------------------------------------------" 
echo "YOUR INSTANCE NAME IS: $INSTANCE_NAME"
echo "-------------------------------------------------------------------------------" 

cd tf

# Delete old terraform state files
rm -rf terraform*

terraform init

ARG="project_id="${PROJECT_ID}
NAME="username="${USER_NAME}
INSTANCE_NAME="instance_name="${INSTANCE_NAME}

terraform apply -auto-approve -var $NAME -var $ARG -var $INSTANCE_NAME

