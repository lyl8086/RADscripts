#!/bin/bash

# Functions
Path()                 
{
    oldIFS=$IFS
    IFS=:
    for d in $PATH
    do
        if [ $d == $1 ]; 
        then 
            return 0
            break
        fi
    done
    IFS=oldIFS
    return 1
}

# Test requirements.

# perl modules
p=(Parallel::ForkManager Array::Shuffle threads Storable Getopt::Long)
for i in ${p[*]}; do
    perl -le "use $i" >&/dev/null || \
    (echo "Please install Perl module $i." && \
    echo "Run (sudo) cpan install $i" && \
    exit 1) || exit 1
done

# stacks.
s=(ustacks cstacks sstacks)
for i in ${s[*]}; do
    if [ ! `which $i 2>/dev/null` ]; then
        echo "$i is not installed!"
        echo "Please download and install it, see http://catchenlab.life.illinois.edu/stacks/"
        exit 1
    fi
done
cstacks --version 2>/tmp/RADpipeline.log
v=`cut -d ' ' -f2 /tmp/RADpipeline.log | tr -d "a-zA-Z"`
rm /tmp/RADpipeline.log
if [ `awk -v v=$v 'BEGIN {if (v<1.45) {print 1} else {print 0}}'` -ne 0 ];
then
    echo "You need to update stacks to 1.45 or later."
    exit 1
fi
sstacks --version 2>/tmp/RADpipeline.log
v=`cut -d ' ' -f2 /tmp/RADpipeline.log | tr -d "a-zA-Z"`
rm /tmp/RADpipeline.log
if [ `awk -v v=$v 'BEGIN {if (v<1.45) {print 1} else {print 0}}'` -ne 0 ];
then
    echo "You need to update stacks to 1.45 or later."
    exit 1
fi
# cap3
if [ ! `which cap3 2>/dev/null` ];
then
echo "cap3 is not installed!"
echo "Please download and install it, see http://seq.cs.iastate.edu/cap3.html"
exit 1
fi
if [ ! `which gnuplot 2>/dev/null` ];
then
echo "please install gnuplot!"
exit 1
fi
# install.
exepath="$HOME/bin"

if [ ! -d "$exepath" ];
then
    echo "Make dir $exepath!"
    mkdir $exepath
fi
cp -a RADassembler bin/ $exepath
chmod +x $exepath/RADassembler

Path "$exepath"
if [ $? != 0 ];
then
    cp ~/.bashrc ~/.bashrc.bak
    echo -e "export PATH=$exepath:\$PATH" >> ~/.bashrc
    echo "Please run 'source ~/.bashrc' first!"
fi

echo "Successfully Installed!"
echo "Run 'RADassembler' for a test"
