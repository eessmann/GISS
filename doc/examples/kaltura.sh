# insert the Kaltura.org javascript for HTML5 video
# http://www.kaltura.org/project/HTML5_Video_Media_JavaScript_Library

if grep -q "</video>" $1.html; then
    sed 's/<\/HEAD>/<script type="text\/javascript" src="http:\/\/html5.kaltura.org\/js"><\/script>\n<\/HEAD>/' < $1.html > tmp
    mv -f tmp $1.html
fi
