#! /usr/bin/perl -w

#####################################################
##
## mltreemap_imagemaker.pl
##
#####################################################

use strict;
use warnings;

use lib 'lib';
use TREEMAP_ml_svg_visualizer;

print "Start MLTreeMap imagemaker v2.061\n";
my $mltreemap_imagemaker                             = mltreemap_imagemaker->new();
my $user_options                                     = $mltreemap_imagemaker->read_user_input(\@ARGV);
my $input_files                                      = $mltreemap_imagemaker->get_input_files($user_options);
my ($concatenated_input_files,$text_of_denominator,$used_colors)  = $mltreemap_imagemaker->concatenate_files($user_options,$input_files);
$mltreemap_imagemaker->run_the_imagemaker($user_options,$concatenated_input_files,$text_of_denominator,$used_colors);
print "Done.\n";


#####################################################
#####################################################

# package mltreemap_imagemaker

#####################################################
#####################################################

package mltreemap_imagemaker;

#####################################################
# new
#####################################################

sub new {
    my $package = shift;
    my $reference = {};
    bless ($reference, $package);
    return ($reference);    
}

#####################################################
# read_user_input
#####################################################

sub read_user_input {    
    my $package = shift;
    my $argv = shift;
    
    my $error_message = "\nERROR: incorrect input.\nexample input:\n./mltreemap_imagemaker.pl -i example_input/\n\n";
    $error_message .= "mandatory input option:\n-i: either a) a path to a directory with input files or b) to a specific input file.\n";
    $error_message .= "note that in case a) all input files of the same type (i.e. p, g, n, h, r, m, d) will be concatenated.\n\n";
    $error_message .= "optional input options:\n-o: output directory (default: output/)\n-d: use different colors for different datasets ";
    $error_message .= "(0: don't use this mode (default), 1: use this mode, 11: use this mode but concatenate percent values over datasets).";
    $error_message .= " See the user guide for a more detailed description.\n";
    $error_message .= "-t: display trees inclusive text and labels (1, default) or without text (0)\n";
    $error_message .= "-b: parameter for the size of the placement bubbles (default: 1)\n-r display 16s and 18s rRNA hits in different trees (2, default).";
    $error_message .= " or one tree of life (1)\n-h: high quality placement bubbles (1, default) or low quality bubbles (0). High quality bubbles work perfectly";
    $error_message .= "with inkscape but can cause trouble with Adobe Illustrator.\n";
        
    my %mandatory_options = (-i => 0);
    my %option_defaults = (-o => "output/", -d => 0, -b => 1, -r => 2, -h => 1, -t => 1);
    
    my %user_options = %option_defaults;
    my $current_option = "";
    my $arg_count = 0;
    
    foreach my $argument (@$argv) {      
        if ( $arg_count % 2 ) {
            #is input
            $user_options{$current_option} = $argument;
        } else {
            #is option denominator
            die "$error_message" unless (exists ($option_defaults{$argument}) || exists ($mandatory_options{$argument}));
        }
        $current_option = $argument;        
        $arg_count++;   
    }
    
    #second, do some (hard coded) checks on certain input variables:
    
    die "$error_message" unless ($user_options{"-b"} =~ /\A\d+\.*\d*\Z/);
    die "$error_message" unless (($user_options{"-d"} eq 0) || ($user_options{"-d"} eq 1) || ($user_options{"-d"} eq 11));
    die "$error_message" unless (($user_options{"-r"} eq 1) || ($user_options{"-r"} eq 2));
    die "$error_message" unless (($user_options{"-t"} eq 1) || ($user_options{"-t"} eq 0));
    
    #add a "/" to the output directory, if not done so by the user.
    
    $user_options{-o} = "$user_options{-o}"."/" unless ($user_options{-o} =~ /\/\Z/);
    #check if the mandatory options have been used:
    
    die "$error_message" unless $arg_count;
    foreach my $mandatory_option (sort {$a cmp $b} keys %mandatory_options) {
        die "$error_message" unless exists ($user_options{$mandatory_option});    
    }
        
    #create the output directories
        
    if (-e $user_options{-o}) {
        print "WARNING: Your output directory \"$user_options{-o}\" allready exists!\n";
        print "Overwrite [1], quit [2] or change directory [3]?\n";
        
        my $alarm_timer = 30;
        $SIG{ALRM} = sub { die "No input for more than $alarm_timer sek. Exit MLTreeMap.\n"};
        alarm $alarm_timer;
        my $answer = <STDIN>;
        $answer = "2" unless $answer; #quit if something goes wrong with <STDIN>.
        chomp $answer;
        alarm 0;
        while (($answer ne "1") && ($answer ne "2") && ($answer ne "3")) {
            print "invalid input. Please chose 1, 2 or 3.\n";
            alarm $alarm_timer;
            $answer = <STDIN>;
            chomp $answer;
            alarm 0;
        }
        if ($answer == 1) {
            print "Do you really want to overwrite the old output directory?\n";
            print "All data in it will be lost!\n";
            print "Yes [y] or no [n]?\n";         
            alarm $alarm_timer;
            my $answer2 = <STDIN>;
            chomp $answer2;
            alarm 0;
            while (($answer2 ne "y") && ($answer2 ne "n")) {
                print "invalid input. Please chose y or n.\n";
                alarm $alarm_timer;
                $answer2 = <STDIN>;
                chomp $answer2;
                alarm 0;
            }
            if ($answer2 eq "y") {
                my $rm_command = "rm -r $user_options{-o}";
                system ($rm_command);
            } else {
                die "Exit MLTreeMap\n";    
            }
        } elsif ($answer == 2) {
            die "Exit MLTreeMap\n";;
        } elsif ($answer == 3) {
            print "Please enter the new directory.\n";
            alarm 600;
            my $dir = <STDIN>;
            chomp $dir;
            alarm 0;
            $user_options{-o} = $dir;
        }
    }
    my $mkdir_command = "mkdir $user_options{-o}";     
    system ($mkdir_command) unless (-e $user_options{-o});
    
    #done.
    
    return (\%user_options);    
}

#####################################################
# get_input_files
#####################################################

sub get_input_files {
    my $package = shift;
    my $user_options = shift;
    
    my %input_files = ();
    my $input = $$user_options{-i};
    my $rRNA_display_option = $$user_options{-r};
    my $input2 = "/"."$input";    
    if ($input2 =~ /.*\/((.)_.*RAxML_.+.txt)\Z/) {
        if (-e $input) {
            my $filename = $1;
            my $denominator = $2;
            $denominator = "A" if (($denominator eq "b" || $denominator eq "a") && ($rRNA_display_option == 1));
            $input_files{$denominator}{$input} = $filename;    
        } else {
            die "ERROR: you entered a file as input ($input), but it does not exist!\n";    
        }
    } else {
        $input .= "/" unless ($input =~ /\/\Z/);
        opendir (PATH, "$input") or die "Error, your input directory ($input) does not exist!\n";
        my @files = readdir PATH;
        closedir (PATH);
        
        my %files2 = ();
        
        foreach my $filename (@files) {
            if ($filename =~ /\A(.)_.*RAxML_.+.txt\Z/) {
                my $denominator = $1;
                my $filename_long = "$input$filename";
                $denominator = "A" if (($denominator eq "b" || $denominator eq "a") && ($rRNA_display_option == 1));
                $input_files{$denominator}{$filename_long} = $filename; 
            }    
        }       
    }
    return (\%input_files);    
}

#####################################################
#####################################################

sub concatenate_files {
    
    my $package = shift;
    my $user_options = shift;
    my $input_files = shift;
    
    my $output_dir = $$user_options{-o};
    my $rRNA_display_option = $$user_options{-r};
    
    my %concatenated_input_files = ();
    my %text_of_denominator = ();
    
    my @colors = "";
    my %used_colors = ();
    
    open (IN, "tree_data/available_dataset_colors.txt") or die "ERROR: tree_data/available_dataset_colors.txt does not exist!\n";
    <IN>;
    my $count = 0;
    while (<IN>) {
        chomp $_;
        my ($red, $green, $blue) = split / /;
        next unless (defined $blue);
        my $color_tag = "rgb($red, $green, $blue)";
        $colors[$count] = $color_tag;
        $count++;
    }
        
    my $nr_of_colors = @colors + 0;
           
    foreach my $denominator (sort {$a cmp $b} keys %$input_files) {
        
        my $output_filename = "$denominator"."_concatenated_RAxML_outputs.txt";
        my $output_filename_long = "$output_dir"."$output_filename";
        open (OUTPUT, "> $output_filename_long") or die "ERROR: Can't create $output_filename_long\n";
        
        if ($$user_options{-d} > 0) {
            open (OUT2, "> $output_dir"."$denominator"."_color_legend.txt") or die "Error, can't create the color legend file\n";
        }
        
        $concatenated_input_files{$denominator}{$output_filename_long} = $output_filename;
        
        my %percentages_of_texts = ();
        my $nr_of_files = 0;

        foreach my $filename_long (sort {$a cmp $b} keys %{$$input_files{$denominator}}) {           
            $nr_of_files++;
            die "Error: this directory contains more datasets than can be displayed with different colors!\n" if ( ($nr_of_files > $nr_of_colors) && ($$user_options{-d} > 0));
            my $attachment = "colorcode_$colors[0]";
            $used_colors{"$colors[0]"} = 1;
            if ($$user_options{-d} > 0) {
                $attachment = "colorcode_$colors[$nr_of_files-1]";
                $used_colors{"$colors[$nr_of_files-1]"} = 1;
                print OUT2 "$colors[$nr_of_files-1]\t$filename_long\n";
            }
            open (FILE, "$filename_long") or die "can't open $filename_long";    
            
            while (<FILE>) {
                chomp $_;
            
                my $weight = "";
                my $text = "";
            
                if (/\APlacement weight (.+)%:(.+)/) {
                    $weight = $1;
                    $text = $2;
                    $text .= $attachment;
                } elsif (/\A# .+ analysis, (.+):/) {
                    $text_of_denominator{$denominator} = $1;
                    $text_of_denominator{$denominator} = "16s rRNA & 18s rRNA" if (($1 =~ /rRNA/) && ($rRNA_display_option == 1));
                    next;
                } elsif (/\A# Phylogenetic analysis(.*)/) {
                    $text_of_denominator{$denominator} = "MLTreeMap tree of life";
                    $text_of_denominator{$denominator} = "GEBA tree of life" if ($1 =~ /GEBA/);
                    $text_of_denominator{$denominator} = "Fungi" if ($1 =~ /fungi/);
                    next;
                } else {
                    next;    
                }
            
                if (exists $percentages_of_texts{$text}) {
                    $percentages_of_texts{$text} += $weight;
                } else {
                    $percentages_of_texts{$text} = $weight;    
                }     
            }
            close FILE;   
        }
        if ($$user_options{-d} > 0) {
            close OUT2;
        }

        my $check = 0;
        foreach my $text (sort {$a cmp $b} keys %percentages_of_texts) {
            my $weight = $percentages_of_texts{$text};
            my $relative_weight = $weight / $nr_of_files;
            $relative_weight = (int($relative_weight * 10000 + 0.5)) / 10000;
            print OUTPUT "Placement weight $relative_weight"."\%:$text\n";
            $check += $relative_weight;
        }    
        print "$denominator files: sum of percentages = $check\n";
        close OUTPUT;
    
    } 
    return(\%concatenated_input_files,\%text_of_denominator,\%used_colors);
}

#####################################################
#####################################################

sub run_the_imagemaker {
    my $package = shift;
    my $user_options = shift;
    my $concatenated_input_files = shift;
    my $text_of_denominator = shift;
    my $used_colors = shift;
    my $param_scale_bubble = $$user_options{-b};
    my $output_dir = $$user_options{-o};
    my $bubble_type = $$user_options{-h};
    my $color_mode = $$user_options{-d};
    my $text_mode = $$user_options{-t};
    
    foreach my $denominator (sort {$a cmp $b} keys %$concatenated_input_files) {
        foreach my $input_filename_long (sort {$a cmp $b} keys %{$$concatenated_input_files{$denominator}}) {
            my $input_filename = $$concatenated_input_files{$denominator}{$input_filename_long};
            
            my $message = TREEMAP_ml_svg_visualizer::run_visualisation ($input_filename,$input_filename_long,$output_dir,$denominator,$param_scale_bubble,$text_of_denominator,$bubble_type,$used_colors,$color_mode,$text_mode);
            print "$message\n";
        }    
    }
}

#####################################################
#####################################################
1;