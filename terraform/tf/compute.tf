data "template_file" "init" {
  template = "${file("./scripts/user_data.sh")}"

  vars {
    username = "${var.username}"
  }
}

resource "google_compute_instance" "tftreesapp" {
  project      = "${var.project_id}"
  zone         = "us-west1-b"
  name         = "${var.instance_name}"
  machine_type = "n1-highcpu-8"

  boot_disk {
    initialize_params {
      image = "ubuntu-1604-xenial-v20180912"
      type  = "pd-standard"
      size  = 40
    }
  }

  metadata_startup_script = "${data.template_file.init.rendered}"

  network_interface {
    network       = "default"
    access_config = {}
  }
}

output "instance_id" {
  value = "${google_compute_instance.tftreesapp.self_link}"
}
