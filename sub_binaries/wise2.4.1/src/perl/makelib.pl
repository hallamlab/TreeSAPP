#!/usr/local/bin/perl


%dontcopy = ( 
	     'psw' => 1,
	     'pswdb' => 1,
	     'dyc' => 1,
	     'genewisedb' => 1,
	     'estwisedb' => 1,
	     'genewise' => 1,
	     'estwise' => 1,
	     'oldpostwise' => 1,
	     'dba' => 1,
	     'hmmalign' => 1,
	     'hmmsearch' => 1,
	     'hmmconvert' => 1,
	     'hmmpfam' => 1,
	     'hmmbuild' => 1,
	     );

@dirs = ( '../base' , '../dynlibsrc' , '../models' , '../HMMer2' );



foreach $dir ( @dirs ) {
    opendir(DIR,$dir) || die "Cannot open $dir";
    @files = readdir(DIR);
    closedir(DIR);

    shift @files;
    shift @files;

    foreach ( @files ) {
	/^(\S+).o$/ || next;
	$name = $1;
	if( exists $dontcopy{$name} ) {
	    next;
	}
	
	$dest = "$dir";
	$dest =~ s/^.*\///g;
	$dest .= "_";
	$dest .= $_;
	push(@objs,$dest);
	&copy_file_if_needed("$dir/$_","Wise2/libs/$dest");
    }

}

open(MA,">Wise2/libs/Makefile");


$files = join("\\\n\t",@objs);

print MA "\nOBJS = $files\n\nlibwise2.a : \$(OBJS)\n\tar -ru libwise2.a \$(OBJS)\n";
    

sub copy_file_if_needed {
    my $file = shift;
    my $dest  = shift;

    if(!(-e $dest) ||  (-M $file) < (-M $dest) ) {
	print STDERR "Updating $dest\n";
	system("cp $file $dest");
    } else {
	print STDERR "Skipping $dest\n";
    }

}


