%define	alphatag %(date +%Y%m%d)
%define current %(gerris2D -V 2>&1 | head -1 | cut -d' ' -f6)

Summary: The Gerris Flow Solver
Name: gerris
%if "%{current}" == ""
Version: 1.1.2
%else
Version: %{current}
%endif
Release: 3.%{alphatag}cvs%{?dist}
License: GPLv2
Group: Applications/Engineering
URL: http://gfs.sourceforge.net
Packager: Ivan Adam Vari <i.vari@niwa.co.nz>
Source0: %{name}-stable.tar.gz
Provides: %{name}
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
Requires: gts-snapshot-devel >= 0.7.6 pkgconfig gcc sed gawk m4
BuildRequires: gts-snapshot-devel >= 0.7.6 glib2-devel geomview
BuildRequires: glibc-devel automake libtool gtkglext-devel
BuildRequires: startup-notification-devel
%if 0%{?fedora_version}
BuildRequires: libXt-devel netpbm-devel
%elif 0%{?suse_version}
BuildRequires: xorg-x11-devel libnetpbm
%endif

%description
Gerris is an Open Source Free Software library for the solution of the 
partial differential equations describing fluid flow. The source code 
is available free of charge under the Free Software GPL license.
Gerris is supported by NIWA (National Institute of Water and Atmospheric research) 
and by the Marsden Fund of the Royal Society of New Zealand.
The code is written entirely in C and uses both the GLib Library and 
the GTS Library for geometrical functions and object-oriented programming. 


%prep
%setup -q -n %{name}-stable


%build
RPM_OPT_FLAGS="$RPM_OPT_FLAGS -fPIC -DPIC"
if [ -x ./configure ]; then
    CFLAGS="$RPM_OPT_FLAGS" ./configure --prefix=%{_prefix} \
	    --libdir=%{_prefix}/%_lib
else
    CFLAGS="$RPM_OPT_FLAGS" sh autogen.sh --prefix=%{_prefix} \
	    --libdir=%{_prefix}/%_lib
fi

%{__make} %{?_smp_mflags}


%install
rm -rf $RPM_BUILD_ROOT
mkdir $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT


%clean
rm -rf $RPM_BUILD_ROOT


%post
/sbin/ldconfig
touch --no-create %{_datadir}/icons/hicolor
if [ -x %{_bindir}/gtk-update-icon-cache ]; then
 %{_bindir}/gtk-update-icon-cache --quiet %{_datadir}/icons/hicolor || :
elif [ -x /opt/gnome/bin/gtk-update-icon-cache ]; then
 /opt/gnome/bin/gtk-update-icon-cache -t --quiet %{_datadir}/icons/hicolor || :
fi
update-mime-database %{_datadir}/mime &> /dev/null || :



%postun
/sbin/ldconfig
touch --no-create %{_datadir}/icons/hicolor
if [ -x %{_bindir}/gtk-update-icon-cache ]; then
 %{_bindir}/gtk-update-icon-cache --quiet %{_datadir}/icons/hicolor || :
elif [ -x /opt/gnome/bin/gtk-update-icon-cache ]; then
 /opt/gnome/bin/gtk-update-icon-cache -t --quiet %{_datadir}/icons/hicolor || :
fi
update-mime-database %{_datadir}/mime &> /dev/null || :


%files
%defattr(-,root,root,-)
%doc NEWS README TODO COPYING
%{_bindir}/*
%{_includedir}/*.h
%{_includedir}/gerris/*.h
%{_libdir}/*.so.*
%{_libdir}/*.so
%{_libdir}/*.a
%{_libdir}/*.la
%{_libdir}/gerris/*
%{_libdir}/pkgconfig/*.pc
%{_datadir}/mime/packages/*.xml
%{_datadir}/icons/hicolor/48x48/mimetypes/*.png


%changelog
* Mon Nov 12 2007 Ivan Adam Vari <i.vari@niwa.co.nz>
- Fixed package (install) dependencies

* Mon Oct 1 2007 Ivan Adam Vari <i.vari@niwa.co.nz>
- Removed unnecessary version specifications for some 
  build requirements
- Added SLEx/SuSE compatibilty
- Added 64bit compatibility
- Updated %post, %postun scriptlets

* Tue May 1 2007 Ivan Adam Vari <i.vari@niwa.co.nz>
- Initial rpm release based on Fedora/Redhat Linux
