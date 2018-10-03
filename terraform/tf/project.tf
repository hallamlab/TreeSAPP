variable "project_id" {
  default = "treesapp"
}

variable "username" {
  default = "treesapp"
}

provider "google" {
  region      = "us-west-1b"
  credentials = "${file("credentials.json")}"
  project     = "${var.project_id}"
}
