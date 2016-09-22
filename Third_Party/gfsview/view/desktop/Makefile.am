## Process this file with automake to produce Makefile.in

mimedir = $(datadir)/mime/packages
mime_DATA = gfsview.xml

desktopdir = $(datadir)/applications

DESKTOPS = gfsview2D.desktop gfsview3D.desktop gfsview.desktop

desktop_DATA = $(DESKTOPS)

iconsdir = $(datadir)/icons/hicolor/48x48/mimetypes

ICONS = application-gfsview.png \
	application-gfsview-2D.png \
	application-gfsview-3D.png

icons_DATA = $(ICONS)

EXTRA_DIST = gfsview.xml $(DESKTOPS) $(ICONS)
