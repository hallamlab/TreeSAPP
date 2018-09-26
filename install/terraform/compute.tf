resource "google_compute_instance" "default" {
  zone         = "us-west1-b"
  name         = "tftreesapp"
  machine_type = "n1-highcpu-8"

  boot_disk {
    initialize_params {
      image = "ubuntu-1604-xenial-v20180912"
      type  = "pd-standard"
      size  = 40
    }
  }

  metadata_startup_script = "${file("./scripts/user_data.sh")}"

  network_interface {
    network       = "default"
    access_config = {}
  }
}

output "instance_id" {
  value = "${google_compute_instance.default.self_link}"
}
