jdkWin = /cygdrive/c/Program Files/Java/jdk1.5.0_04
jdkLinux = /usr/java/jdk1.5.0_04

.PHONY: javah genewiseJNIWin genewiseJNI

javah:
	javah -jni -classpath ../../../java/src -o genewiseJNI.h edu.stanford.hiv.algorithm.sequence.align.GenewiseAlign

genewiseJNIWin:
	gcc -mno-cygwin -I"${jdkWin}/include" -I"${jdkWin}/include/win32" -Wl,--add-stdcall-alias -shared -o genewiseJNI.dll genewiseJNI.c

genewiseJNI:
	gcc -I"${jdkLinux}/include" -I"${jdkLinux}/include/linux" -fpic -shared -o genewiseJNI.so genewiseJNI.c
