%define alphatag %(date +%Y%m%d)
%define gerris_version %(pkg-config gerris2D --modversion --silence-errors)
%define current %(gfsview-batch2D -V 2>&1 | head -1 | cut -d ' ' -f6)

Summary: GfsView visualises Gerris simulation files (development snapshot)
Name: gfsview-snapshot
%if "%{current}" == ""
Version: 0.5.0
%else
Version: %{current}
%endif
Release: 10.%{alphatag}cvs%{?dist}
License: GPLv2
%if 0%{?suse_version}
Group: Productivity/Scientific/Other
%else
Group: Applications/Engineering
%endif
URL: http://gfs.sourceforge.net
Packager: Matthieu Castellazzi <m.castellazzi@niwa.co.nz>
Source0: %{name}.tar.gz
BuildRoot:  %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
%if 0%{?suse_version}
BuildRequires: update-desktop-files
%if 0%{?suse_version} >= 1110
BuildRequires: ftgl-devel
%endif
%else
BuildRequires: ftgl-devel
%endif
%if 0%{?centos_version} || 0%{?fedora} 
BuildRequires: mesa-libOSMesa-devel
%endif
%if 0%{?fedora_version}
%if %{?fedora_version} > 15
BuildRequires: atlas
%endif
%endif
# For all distros
Requires: zenity, gzip
BuildRequires: gerris-snapshot >= %{gerris_version}
BuildRequires: gtk2-devel, gtkglext-devel, startup-notification-devel
BuildRequires: gerris-snapshot-devel >= %{gerris_version}
BuildRequires: openmpi openmpi-devel

%description
GfsView displays the results of 2D and 3D Gerris simulations.

A brief summary of its main features:

    * Scalar and vector cross-sections.
    * Isosurfaces.
    * Streamlines.
    * User-defined functions.
    * Fast adaptive display (using the multiresolution data
      representation of Gerris).
    * Scriptable.
    * Offline image generation.
    * Quality PostScript, PDF and bitmap outputs.

See http://gfs.sf.net for more information and documentation.


%prep
%setup -q -n %{name}


%build

# if we have centos or rhel, set mpi-selector
%if 0%{?rhel_version} || 0%{?centos_version}
%if %{?centos_version} < 6
mpi-selector --set $(mpi-selector --list)
source /etc/profile.d/mpi-selector.sh
%endif
%endif

RPM_OPT_FLAGS="$RPM_OPT_FLAGS -fPIC -DPIC"
if [ -x ./configure ]; then
    CFLAGS="$RPM_OPT_FLAGS" ./configure \
	--prefix=%{_prefix} \
	--libdir=%{_prefix}/%_lib \
	--mandir=%{_mandir} \
	--disable-static
else
    CFLAGS="$RPM_OPT_FLAGS" sh autogen.sh \
	--prefix=%{_prefix} \
	--libdir=%{_prefix}/%_lib \
	--mandir=%{_mandir} \
	--disable-static
fi

%{__make}


%install

# if we have centos or rhel, set mpi-selector
%if 0%{?rhel_version} || 0%{?centos_version}
%if %{?centos_version} < 6
mpi-selector --set $(mpi-selector --list)
source /etc/profile.d/mpi-selector.sh
%endif
%endif

rm -rf $RPM_BUILD_ROOT
mkdir $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT

# Comply shared library policy
find $RPM_BUILD_ROOT -name *.la -exec rm -f {} \;

# Comply static build policy
find $RPM_BUILD_ROOT -name *.a -exec rm -f {} \;

%if 0%{?suse_version}
%suse_update_desktop_file -n $RPM_BUILD_ROOT/%{_datadir}/applications/gfsview2D.desktop
%suse_update_desktop_file -n $RPM_BUILD_ROOT/%{_datadir}/applications/gfsview3D.desktop
%suse_update_desktop_file -n $RPM_BUILD_ROOT/%{_datadir}/applications/gfsview.desktop
%endif

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
%defattr(-,root,root)
%doc COPYING NEWS README TODO
%{_bindir}/*
%{_libdir}/libgfsgl?D.so*
%dir %{_libdir}/gerris
%{_libdir}/gerris/libgfsview*.so
%{_datadir}/applications/*.desktop
%dir %{_datadir}/gfsview
%dir %{_datadir}/gfsview/pixmaps
%{_datadir}/gfsview/pixmaps/*.png
%dir %{_datadir}/gfsview/fonts
%{_datadir}/gfsview/fonts/*.ttf
%dir %{_datadir}/icons/hicolor
%dir %{_datadir}/icons/hicolor/48x48
%dir %{_datadir}/icons/hicolor/48x48/mimetypes
%{_datadir}/icons/hicolor/48x48/mimetypes/*.png
%{_datadir}/mime/packages/*.xml

%changelog
* Tue Feb 16 2010 Matthieu Castellazzi <m.castellazzi@niwa.co.nz> - 10
- add gerris-snapshot-devel dependency

* Thu Nov 05 2009 Matthieu Castellazzi <m.castellazzi@niwa.co.nz> - 9
- rename the package to gfsview-snapshot

* Tue Nov 03 2009 Matthieu Castellazzi <m.castellazzi@niwa.co.nz> - 8
- first attempt to use buildservice for gerris
- add --silence-errors options to define %current

* Wed Jul 16 2008 Ivan Adam Vari <i.vari@niwa.co.nz> - 7
- Fixed wrong version specification
  Some other changes found in debian packages

* Thu Jul 3 2008 Ivan Adam Vari <i.vari@niwa.co.nz> - 6
- Fixed typo in %files section (attr)

* Thu May 15 2008 Ivan Adam Vari <i.vari@niwa.co.nz> - 5
- Added fedora 8 support for x86 (32bit only)
- Removed libtool config files to comply with shared
  library policy
- Removed static build bits to comply with shared
  library policy
- Fixed dependencies

* Mon Jan 7 2008 Ivan Adam Vari <i.vari@niwa.co.nz> - 4
- Removed %{?_smp_mflags} from make due to intermittent
  build errors on some SMP systems

* Mon Nov 12 2007 Stephane Popinet <popinet@users.sf.net> - 3
- Fixed package (install) dependencies

* Mon Oct 1 2007 Ivan Adam Vari <i.vari@niwa.co.nz> - 2
- Removed unnecessary version specifications for some
  build requirements
- Added SLEx/SuSE compatibilty
- Added 64bit compatibility
- Updated %post, %postun scriptlets

* Tue May 1 2007 Ivan Adam Vari <i.vari@niwa.co.nz> - 1
- Initial rpm release based on Fedora/Redhat Linux.
