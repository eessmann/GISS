#! /bin/bash

for file in faq/*.html; do
awk '{
    if ($1 == "<PRE>") {
	inpre = 1;
	npre = 0;
    }
    else if ($1 == "</PRE>")
	inpre = 0;
    if (inpre) {
	if (NF > 0 || npre > 1)
	    print $0;
	npre++;
    }
    else
	print $0;
}' < $file | sed 's/<TT> /<TT>/g' > /tmp/`basename $file`
mv -f /tmp/`basename $file` $file
done
