#!/usr/bin/perl
use strict;
use LWP::UserAgent;
use HTTP::Request::Common;
use Config::Simple;


use Data::Dumper;

my $upload_url = "https://itol.embl.de/batch_uploader.cgi";

my $configFile = $ARGV[0];

unless (-e $configFile) {
		print STDERR "Usage: iTOL_uploader.pl config_file\n\nPlease provide a config file with upload parameteres.\n";
		exit;
}

my ($zipFile, $APIkey, $projectName, $treeName, $treeDescription);

my $cfg = new Config::File::Simple($configFile);

$zipFile = $cfg->read('zipFile') if ($cfg->read('zipFile'));
$treeName = $cfg->read('treeName') if ($cfg->read('treeName'));
$treeDescription = $cfg->read('treeDescription') if ($cfg->read('treeDescription'));
$APIkey = $cfg->read('APIkey') if ($cfg->read('APIkey'));
$projectName = $cfg->read('projectName') if ($cfg->read('projectName'));

print "\niTOL batch uploader\n=====================\n";

unless ( -e $zipFile) {
  print STDERR "Please provide a zip file with the tree and optionally other annotation files\n";
 exit;
}

#prepare the  POST data
my %post_content;
$post_content{'zipFile'} = [ $zipFile ];
if ($treeName) {$post_content{'treeName'} = $treeName;}
if ($treeDescription) {$post_content{'treeDescription'} = $treeDescription;}
if ($APIkey) {
		$post_content{'APIkey'} = $APIkey;
		if (not $projectName)  {die "projectName is required when using APIkey";}
		$post_content{'projectName'} = $projectName;
}

#submit the data
my $ua  = LWP::UserAgent->new();
$ua->agent("iTOLbatchUploader4.0");
my $req = POST $upload_url, Content_Type => 'form-data', Content => [ %post_content ];
my $response = $ua->request($req);

if ($response->is_success()) {
		my @res = split(/\n/, $response->content);
		#check for an upload error
		if ($res[$#res] =~ /^ERR/) {
				print "Upload failed. iTOL returned the following error message:\n\n$res[$#res]\n\n";
				print "Full output:\n";
				print join("\n", @res);
				exit;
		}
		#upload without warnings, ID on first line
		if ($res[0] =~ /^SUCCESS: (\S+)/) {
				print "Upload successful. Your tree is accessible using the following iTOL tree ID:\n\n$1\n\n";
				exit;
		}
		if ($res[$#res] =~ /^SUCCESS: (\S+)/) {
            print "Upload successful. Warnings occured during upload:\n\n";
            for (0..$#res-1) {
                print $res[$_] . "\n";
            }

            print "\n\nYour tree is accessible using the following iTOL tree ID:\n\n$1\n\n";
            exit;
		} else {
            print "This shouldn't happen. iTOL did not return an error, but there is no tree ID. Please email the full dump below to ivica\@letunic.com:\n\n===DEBUG DUMP==\n";
            print join("\n", @res);
            print  "\n===END DUMP==\n";
		}
} else {
		print "iTOL returned a web server error. Full message follows:\n\n";
		print $response->as_string;
}
