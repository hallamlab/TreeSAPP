#!/usr/local/bin/perl


while(<>) {
    s/module_/!!MODULE/g;
    s/object_/!!OBJECT/g;
    s/accessing_/!!ACCESS/g;
    s/\_/\\\_/g;
    s/\&/\\\&/g;
    s/\$/\\\$/g;
    s/\'/\\\'/g;
    s/\#/\\\#/g;
    s/!!OBJECT/object_/g;
    s/!!ACCESS/accessing_/g;
    s/!!MODULE/module_/g;
    s/module_sw\\_wrap/module_sw_wrap/g;
    print;
}
