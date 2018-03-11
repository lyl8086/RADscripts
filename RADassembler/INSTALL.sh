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
perl -le "use Parallel::ForkManager" >& /dev/null
if [ $? != 0 ];
then
	echo "Please install Perl module Parallel::ForkManager."
    echo "Run (sudo) cpan install Parallel::ForkManager"
	exit 1
fi
perl -le "use Array::Shuffle" >& /dev/null
if [ $? != 0 ];
then
	echo "Please install Perl module Array::Shuffle."
    echo "Run (sudo) cpan install Array::Shuffle"
	exit 1
fi
# stacks.
if [ ! `which ustacks 2>/dev/null` ];
then
	echo "ustacks is not installed!"
	echo "Please download and install it, see http://catchenlab.life.illinois.edu/stacks/"
	exit 1
fi
if [ ! `which cstacks 2>/dev/null` ];
then
	echo "cstacks is not installed!"
	echo "Please download and install it, see http://catchenlab.life.illinois.edu/stacks/"
	exit 1
else
	cstacks --version 2>/tmp/RADpipeline.log
	v=`cut -d ' ' -f2 /tmp/RADpipeline.log | tr -d "a-zA-Z"`
    rm /tmp/RADpipeline.log
	if [ `awk -v v=$v 'BEGIN {if (v<1.45) {print 1} else {print 0}}'` -ne 0 ];
	then
		echo "You need to update stacks to 1.45 or later."
		exit 1
	fi
fi
if [ ! `which sstacks 2>/dev/null` ];
then
	echo "sstacks is not installed!"
	echo "Please download and install it, see http://catchenlab.life.illinois.edu/stacks/"
	exit 1
else
	sstacks --version 2>/tmp/RADpipeline.log
	v=`cut -d ' ' -f2 /tmp/RADpipeline.log | tr -d "a-zA-Z"`
    rm /tmp/RADpipeline.log
	if [ `awk -v v=$v 'BEGIN {if (v<1.45) {print 1} else {print 0}}'` -ne 0 ];
	then
		echo "You need to update stacks to 1.45 or later."
		exit 1
	fi
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
