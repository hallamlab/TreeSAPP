use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	=> 'Wise2',
    'VERSION'	=> '0.1',
    'LIBS'	=> ['-lm'],   # e.g., '-lm' 
    'DEFINE'	=> '-DPOSIX -DNOERROR',     # e.g., '-DHAVE_SOMETHING' 
    'INC'	=> '-I./libs',     # e.g., '-I/usr/include/other'
    'MYEXTLIB'  => 'libs/libwise2$(LIB_EXT)' 
);


