#/bin/bash

exepath="$HOME/bin/"

if [ ! -d "$exepath" ];
then
	echo "Make dir $exepath!"
	mkdir $exepath
fi

cp -a RADassembly bin/ $exepath

chmod +x $exepath/RADassembly
cp ~/.bashrc ~/.bashrc.bak
echo -e "export PATH=$exepath:\$PATH" >> ~/.bashrc

echo "Please run source ~/.bashrc first!"




