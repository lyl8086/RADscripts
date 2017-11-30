#/bin/bash

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
# bc
if [ ! `which bc 2>/dev/null` ];
then
	echo "Please install bc"
	exit 1
fi

# perl Parallel::ForkManage
perldoc Parallel::ForkManager >& /dev/null
if [ $? != 0 ];
then
	echo "Please install Perl module Parallel::ForkManager."
    echo "Run cpan install Parallel::ForkManager"
	exit 1
fi

# stacks.
if [ ! `which ustacks 2>/dev/null` ];
then
	echo "ustacks does not exist!"
	echo "Please download it at http://catchenlab.life.illinois.edu/stacks/"
	exit 1
fi
if [ ! `which cstacks 2>/dev/null` ];
then
	echo "cstacks does not exist!"
	echo "Please download it at http://catchenlab.life.illinois.edu/stacks/"
	exit 1
else
	cstacks --version 2>/tmp/RADpipeline.log
	v=`cut -d ' ' -f2 /tmp/RADpipeline.log`
	if [ `echo "$v<1.45"|bc` -ne 0 ];
	then
		echo "You need to update stacks to 1.45 or later."
		exit 1
	fi
fi
if [ ! `which sstacks 2>/dev/null` ];
then
	echo "sstacks does not exist!"
	echo "Please download it at http://catchenlab.life.illinois.edu/stacks/"
	exit 1
else
	sstacks --version 2>/tmp/RADpipeline.log
	v=`cut -d ' ' -f2 /tmp/RADpipeline.log`
	if [ `echo "$v<1.45"|bc` -ne 0 ];
	then
		echo "You need to update stacks to 1.45 or later."
		exit 1
	fi
fi

# cap3
if [ ! `which cap3 2>/dev/null` ];
then
echo "cap3 does not exist!"
echo "Please download it at http://seq.cs.iastate.edu/cap3.html"
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
