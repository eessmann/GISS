# !/bin/sh

if test -d test; then
  :
else
    exit 1
fi

cd test
for file in *.xmgr; do
    grace -param ../divfig.par $file -nosafe -batch ../makeref
    mv -f /tmp/makeref.xmgr ../$1/$file
done
cd ..
