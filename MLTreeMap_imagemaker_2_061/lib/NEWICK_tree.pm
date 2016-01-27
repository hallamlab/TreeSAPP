package NEWICK_tree;

use strict;
use warnings;

#########################################################################################
## constructor.
##
## this just makes the tree object and sets all values to 'empty'
#########################################################################################

sub new {

    my ($that) = @_;

    my $class = ref ($that) || $that;          
    my $self = {};                             
    bless $self, $class;

    $self->{node_id_counter} = -1;
    $self->{parent_of_node} = {};
    $self->{children_of_node} = {};
    $self->{terminal_children_of_node} = {};
    $self->{branch_length_of_node} = {};
    
    return $self;
}

####################################################################################################
## subroutine: read_tree_topology_from_file ()
##
## This reads the topology of the species-tree from a flat-file.
####################################################################################################

sub read_tree_topology_from_file {

    my ($self, $filename) = @_;

    my $input_concatenated = "";

    open (FH, "$filename") or die "cannot open inputfile '$filename'\n";

    while (<FH>) {
        chomp;
        s/\s+//g;
        $input_concatenated .= $_;
    }
    close FH;

    $self->parse_tree_topology_from_string ($input_concatenated);
}

####################################################################################################
## subroutine: parse_tree_topology_from_string ()
##
##
####################################################################################################

sub parse_tree_topology_from_string {

    my ($self, $tree_string) = @_;

    my @tree_input_symbols = split //, $tree_string;

    $self->species_tree_parse_node (\@tree_input_symbols);

    $self->connect_terminal_children_to_node (-1);

    $self->assign_cumulative_branchlength_to_node (-1, 0);

    $self->compute_average_cumulative_branch_length ();
}

##################################################################################################################
## subroutine: species_tree_parse_node ()
##
## This is a recursive routine to parse the species-tree topology out of an array of input symbols. The input
## symbols are brackets, letters, and numbers. The routine was initially written by Daniel, and is quite clever.
##################################################################################################################

sub species_tree_parse_node {

    my ($self, $inputarray) = @_;
    my $next_item;
    my $this_node = $self->{node_id_counter}--;
    my $nr = 0;
    

    my @children_of_this_node = ();

    ## first symbol must be open bracket.

    unless (($next_item = shift @$inputarray) eq "(") { die "#1 error parsing species tree\n"; }

    while (1) {  ## now loop over all subnodes of this node ...

        ## read next symbol: must be either open bracket or species_name.

        unless ($next_item = shift @$inputarray) { die "#2 error parsing species tree\n"; }

        if ($next_item eq "(") { ## open bracket ? call yourself recursively

            unshift @$inputarray, $next_item;
            my $sub_node = $self->species_tree_parse_node ($inputarray);
            $self->{parent_of_node}->{$sub_node} = $this_node;
            push @children_of_this_node, $sub_node;

        } elsif ($next_item =~ /\w+/) { ## word character ? this is a terminal node ('leaf' - species)

            my $node_summary = $next_item;
            while (($next_item = shift @$inputarray) =~ /[\w\.\:]+/) {
                $node_summary .= $next_item;
            }

	    my ($species, $branch_length) = (undef, undef);

	    if ($node_summary =~ /\A(\w+)\:([\d\.]+\z)/) {            ## are there branch-lengths ?
		($species, $branch_length) = ($1, $2);
	    } else {
		$species = $node_summary;
	    }
	    die "#4a error parsing species tree\n" unless defined $species;

            push @children_of_this_node, $species;
            $self->{parent_of_node}->{$species} = $this_node;
            $self->{branch_length_of_node}->{$species} = $branch_length if defined $branch_length;
	        $self->{terminal_nodes}->{$species} = 1;
            unless ($next_item) { die "#3 error parsing species tree\n"; }
            unshift @$inputarray, $next_item;
        
        }
        else { ## something else ? illegal.

            die "#4 error parsing species tree\n";
        }

        ## ok. one subnode has been parsed. next symbol must be either a comma
        ## (in which case we loop to the next subnode) or a closing bracket
        ## (in which case we are done).

        unless ($next_item = shift @$inputarray) { die "#5 error parsing species tree\n"; }

        next if ($next_item eq ",");
        next if ($next_item eq ":");

        if ($next_item eq ")") {

            ## code to get the branch_length here ...

            my $branch_summary = "";
            while ((@$inputarray) and (($next_item = shift @$inputarray) =~ /[\d\.\:]+/)) {
                $branch_summary .= $next_item;
            }
	    my ($bootstrap_support, $branch_length) = $branch_summary =~ /\A(\d*)\:([\d\.]+\z)/;
            $self->{branch_length_of_node}->{$this_node} = $branch_length if defined $branch_length;
	    $self->{bootstrap_support_of_node}->{$this_node} = $bootstrap_support if defined $bootstrap_support;
            unshift @$inputarray, $next_item;
            last;
        }

        ## anything else ? illegal.

        die "#6 error parsing species tree\n";
    }
    
    foreach my $child (@children_of_this_node) {
        $self->{children_of_node}->{$this_node}{$child} = 1;
            }
    
    return ($this_node);
}


####################################################################################################
## subroutine: connect_terminal_children_to_node ()
##
####################################################################################################

sub connect_terminal_children_to_node {

    my ($self, $node) = @_;

    if (($node =~ /\A[\-\d]+\z/) and ($node < 0)) {
	
	foreach my $child (keys %{$self->{children_of_node}->{$node}}) {
	    
	    if (($child =~ /\A[\-\d]+\z/) and ($child < 0)) {
		
		$self->connect_terminal_children_to_node ($child);
		foreach my $subchild (keys %{$self->{terminal_children_of_node}->{$child}}) {
		    $self->{terminal_children_of_node}->{$node}{$subchild} = 1;
		}
		
	    } else {
		
		$self->{terminal_children_of_node}->{$node}{$child} = 1;
		$self->{terminal_children_of_node}->{$child}{$child} = 1;
	    }
	}

    } else {       ## if the current node is a terminal node, attach it to itself.

	$self->{terminal_children_of_node}->{$node}{$node} = 1;
	$self->{terminal_children_of_node}->{$node}{$node} = 1;
    }
}

####################################################################################################
## subroutine: assign_cumulative_branchlength_to_node ()
##
####################################################################################################

sub assign_cumulative_branchlength_to_node {

    my ($self, $node, $previous_cumulative_length) = @_;

    my $branch_length_this_node = 0;
    $branch_length_this_node = $self->{branch_length_of_node}->{$node} if exists $self->{branch_length_of_node}->{$node};

    my $new_cumulative_length = $branch_length_this_node + $previous_cumulative_length;
    
    $self->{cumulative_branch_length_of_node}{$node} = $new_cumulative_length;

    if (($node =~ /\A[\-\d]+\z/) and ($node < 0)) {
	
	foreach my $child (keys %{$self->{children_of_node}->{$node}}) {

	    $self->assign_cumulative_branchlength_to_node ($child, $new_cumulative_length);
	}
    } 
}

####################################################################################################
## subroutine: compute_average_cumulative_branch_length ()
##
####################################################################################################
	    
sub compute_average_cumulative_branch_length {

    my ($self) = @_;

    my $valid_node_counter = 0;
    my $branch_length_sum = 0;

    foreach my $node (keys %{$self->{cumulative_branch_length_of_node}}) {
	
	next if (($node =~ /\A[\-\d]+\z/) and ($node < 0));       ## internal nodes not considered here.
	
	$valid_node_counter += 1;
	$branch_length_sum += $self->{cumulative_branch_length_of_node}->{$node};
    }

    $self->{average_cumulative_branch_length} = $branch_length_sum / $valid_node_counter;
}

####################################################################################################
## subroutine: connect_terminal_children_to_node ()
##
####################################################################################################

sub print_as_verbose_text {

    my ($self) = @_;

    foreach my $node (sort {$b <=> $a} keys %{$self->{children_of_node}}) {
    
	my $branch_length = "undefined";
        $branch_length = $self->{branch_length_of_node}->{$node} if exists $self->{branch_length_of_node}->{$node};
	my @children = sort {$a <=> $b} keys %{$self->{children_of_node}->{$node}};
	my $nr_children = scalar @children;
	my @terminal_children = sort {$a <=> $b} keys %{$self->{terminal_children_of_node}->{$node}};
	my $nr_terminal_children = scalar @terminal_children;
	
	print "node $node has branch-length $branch_length, $nr_children children [";
	print join ".", @children;
	print "], and $nr_terminal_children terminal children: [";
	print join ".", @terminal_children;
	print "]\n";
    }

    foreach my $node (sort {$a <=> $b} keys %{$self->{parent_of_node}}) {
	next unless $node > 0;
	my $branch_length = "undefined";
	$branch_length = $self->{branch_length_of_node}->{$node} if exists $self->{branch_length_of_node}->{$node};
	print "node $node has branch-length $branch_length\n";
    }
}

#######################################################################################################
## subroutine: get_species_to_species_tree_distance ()
##
## given two species, the cumulative branch-length needed to traverse from one to the other is 
## computed and provided by this routine. The routine will return -1 in case of errors.
#######################################################################################################

sub get_species_to_species_tree_distance {

    my ($self, $species1, $species2) = @_;

    return -1 unless defined $species1;
    return -1 unless defined $species2;

    return 0 if $species1 eq $species2;

    next unless exists $self->{parent_of_node}{$species1};
    next unless exists $self->{parent_of_node}{$species2};

    my %parent_traversion_count = ();

    my $node1 = $species1; $parent_traversion_count{$node1} = 1;
    while ((defined $node1) and (not ($node1 eq '-1'))) { $node1 = $self->{parent_of_node}{$node1}; $parent_traversion_count{$node1} += 1; }

    my $node2 = $species2; $parent_traversion_count{$node2} = 1;
    while ((defined $node2) and (not ($node2 eq '-1'))) { $node2 = $self->{parent_of_node}{$node2}; $parent_traversion_count{$node2} += 1; }

    my $cumulative_tree_distance = 0;

    $node1 = $species1; 
    while ($parent_traversion_count{$node1} < 2) { 
	$cumulative_tree_distance += $self->{branch_length_of_node}{$node1}; 
	$node1 = $self->{parent_of_node}{$node1};
    }

    $node2 = $species2; 
    while ($parent_traversion_count{$node2} < 2) { 
	$cumulative_tree_distance += $self->{branch_length_of_node}{$node2}; 
	$node2 = $self->{parent_of_node}{$node2};
    }

    return $cumulative_tree_distance;
}

#######################################################################################################
## subroutine: optimize_species_display_order ()
##
## this makes the tree optically nicer by swapping subtrees such that the internal nodes are maximally 
## 'staggered'. This only affects the displaying of the tree, not its topology or biological meaning. 
## If necessary, manual intervention into the display order is also possible, in the code below 
## (note: this will fuck up if the underlying tree is changed !!)
#######################################################################################################

sub optimize_species_display_order {

    my ($self) = @_;

    $self->{optimized_display_order_counter} = 1;

    $self->optimize_node_recursively (-1);
}

####################################################################################################
## subroutine: optimize_node_recursive ()
##
####################################################################################################

sub optimize_node_recursively {

    my ($self, $node) = @_;

    unless (($node =~ /\A[\-\d]+\z/) and ($node < 0)) {

	$self->{species_display_order}->{$node} = $self->{optimized_display_order_counter};
	$self->{optimized_display_order_counter} += 1;
	return;
    }

    ## ok, this node is an internal one - sort its children by 'terminal_size', and recurse down into them.

    my @children_this_node = keys %{$self->{children_of_node}->{$node}};
    my %nr_terminal_children_per_child = ();

    foreach my $child (@children_this_node) {
	my $nr_children = scalar keys %{$self->{terminal_children_of_node}->{$child}};
	$nr_terminal_children_per_child{$child} = $nr_children;
    }

    foreach my $child (sort { my $sort_value_a = $nr_terminal_children_per_child{$a};
			      my $sort_value_b = $nr_terminal_children_per_child{$b};
			      if ($sort_value_a == $sort_value_b) { $sort_value_a = $a; $sort_value_b = $b; }
			      return $sort_value_a <=> $sort_value_b;
			  } @children_this_node) {

	$self->optimize_node_recursively ($child);
    }
}

####################################################################################################
## subroutine: find_last_common_ancestor ()
##
####################################################################################################

sub find_last_common_ancestor {
    
    my ($self, $node1, $node2) = @_;
    return $node1 if $node1 == $node2;
    return $node1 if $self->is_in_parental_line ($node1, $node2);
    return $node2 if $self->is_in_parental_line ($node2, $node1);
    while (1) {
	return -1 unless exists $self->{parent_of_node}{$node1};
	$node1 = $self->{parent_of_node}{$node1};
	return $node1 if $self->is_in_parental_line ($node1, $node2);
	return -1 if $node1 == -1;
    }
}

####################################################################################################
## subroutine: is_in_parental_line ()
##
####################################################################################################

sub is_in_parental_line {

    my ($self, $putative_parent_node, $query_node) = @_;
    return 1 if $putative_parent_node == $query_node;
    while (1) {
	$query_node = $self->{parent_of_node}{$query_node};
	return 1 if $putative_parent_node == $query_node;
	last if $query_node == -1;
    }
    return 0;
}

    
1;
