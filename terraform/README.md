# Running TreeSAPP on Google Cloud Platform VM

## Setup GCP Account
If you do not already have an account, create one with your Google account [here](https://cloud.google.com/compute/) and click the blue "Try free" button.

## Installations

### Terraform
[Terraform](https://www.terraform.io/intro/index.html) will be used to provision a new project and instance on Google Cloud Platform. With Terraform, a state file will be generated that is saved to your local machine so that your resources can be easily recreated and redeployed. 

After [installing Terraform](https://www.terraform.io/intro/getting-started/install.html), verify that it is installed properly by running:

```
$ terraform
Usage: terraform [--version] [--help] <command> [args]
# ...
```

If you get an error, make sure your PATH variable is set correctly. See [here](https://stackoverflow.com/questions/14637979/how-to-permanently-set-path-on-linux-unix) for setting PATH on Linux and Mac and [here](https://stackoverflow.com/questions/1618280/where-can-i-set-path-to-make-exe-on-windows) for setting PATH on Windows.

### Cloud SDK
Cloud SDK is the command-line interface for Google Cloud Platform. You can install it [here](https://cloud.google.com/sdk/). 

## Running the Instance

Run the following commands to download the TreeSAPP repository from GitHub. 

```
git clone git@github.com:hallamlab/TreeSAPP.git
```
Next, move to the terraform directory by running:

```
cd TreeSAPP/terraform
```
Sign into your GCP account by running gcloud init (if you haven't already). 

Run the setup script and enable billing when prompted. 

```
./scripts/setup.sh
```

If you want to create an instance in an existing account, pass the project name as a parameter to the setup script like so:
```
./scripts/setup.sh [YOUR-PROJECT-NAME]
```

Your project ID and instance ID will be displayed in the script output:

```
-----------------------------------------------
YOUR PROJECT ID IS: [YOUR PROJECT NAME]-xxxxxx
-----------------------------------------------
YOUR INSTANCE ID IS: [INSTANCE-ID]
-----------------------------------------------

```

##### Enabling Billing
Go to your GCP console and in the upper left corner, select the project that you have just created. In the navigation menu, select the Billing tab and click 'Link a billing account'.  Return back to your terminal window and press [ENTER].

The script will now generate a key for your local machine to access your instance. This may take a few minutes.

Verify that you have the following output:

```
Apply complete! Resources: 1 added, 0 changed, 0 destroyed.

Outputs:

instance_id = https://www.googleapis.com/compute/v1/projects/[YOUR-PROJECT-ID-HERE]/zones/us-west1-b/instances/[YOUR-PROJECT-ID-HERE]
```

Now, you can ssh into your instance using:
```
gcloud compute ssh [YOUR-INSTANCE-ID-HERE]
```
The setup script will install all the dependencies TreeSAPP relies on and will take around 5-10 minutes to complete. 

Run the following to start TreeSAPP:
```
source ~/.bashrc
treesapp [INPUTS]
```


