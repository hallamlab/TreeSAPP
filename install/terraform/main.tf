provider "google" {
  credentials = "${file("credentials.json")}"
  region      = "us-west1-b"
}
