%define	alphatag %(date +%Y%m%d)
%define current %(pkg-config gerris2D --modversion)
%define gts_version %(pkg-config gts --modversion)

Summary: The Gerris Flow Solver (development snapshot)
Name: gerris
%if "%{current}" == ""
Version: 1.3.2
%else
Version: %{current}
%endif
Release: 9.%{alphatag}cvs%{?dist}
License: GPLv2
# SuSE should have this macro set. If doubt specify in ~/.rpmmacros
%if 0%{?suse_version}
Group: Productivity/Scientific/Other
%endif
# For Fedora you must specify fedora_version in your ~/.rpmmacros file
%if 0%{?fedora_version}
Group: Applications/Engineering
%endif
URL: http://gfs.sourceforge.net
Packager: Ivan Adam Vari <i.vari@niwa.co.nz>
Source0: %{name}-stable.tar.gz
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root-%(%{__id_u} -n)
%if 0%{?fedora_version}
Requires: proj gsl netcdf
BuildRequires: netcdf-devel proj-devel gcc-gfortran
%endif
%if 0%{?suse_version}
Requires: libproj0 gsl libnetcdf-4
BuildRequires: libnetcdf-devel libproj-devel gcc-fortran
%endif
# For both distros
Requires: gts-snapshot-devel >= %{gts_version} pkgconfig gcc sed gawk m4
BuildRequires: glibc-devel automake libtool gsl-devel gts-snapshot-devel >= %{gts_version}


%description
Gerris is an Open Source Free Software library for the solution of the 
partial differential equations describing fluid flow.

Gerris is supported by NIWA (National Institute of Water and Atmospheric
research) and by the Marsden Fund of the Royal Society of New Zealand.

A brief summary of its main (current) features:

    * Quadtree-based (Octree in 3D) spatial discretisation with
      automatic and dynamic refinement.
    * Multigrid Poisson solver.
    * Second-order Godunov type advection scheme.
    * Solves the time-dependent incompressible Euler, Stokes ans Navier-Stokes
      equations.
    * Support for complex solid boundaries (automatic locally-refined
      mesh generation).
      
See http://gfs.sf.net for more information and documentation.


%prep
%setup -q -n %{name}-stable


%build
RPM_OPT_FLAGS="$RPM_OPT_FLAGS -fPIC -DPIC"
%if 0%{?suse_version}
if [ -x ./configure ]; then
    CFLAGS="$RPM_OPT_FLAGS" ./configure \
	--prefix=%{_prefix} \
	--libdir=%{_prefix}/%_lib \
	--mandir=%{_mandir} \
	--disable-mpi \
	--disable-static
else
    CFLAGS="$RPM_OPT_FLAGS" sh autogen.sh \
	--prefix=%{_prefix} \
	--libdir=%{_prefix}/%_lib \
	--mandir=%{_mandir} \
	--disable-mpi \
	--disable-static
fi
%endif
%if 0%{?fedora_version}
if [ -x ./configure ]; then
    CFLAGS="$RPM_OPT_FLAGS" \
    CPPFLAGS="-I%{_includedir}/netcdf-3" ./configure \
	--prefix=%{_prefix} \
	--libdir=%{_prefix}/%_lib \
	--mandir=%{_mandir} \
	--disable-mpi \
	--disable-static
else
    CFLAGS="$RPM_OPT_FLAGS" \
    CPPFLAGS="-I%{_includedir}/netcdf-3" sh autogen.sh \
	--prefix=%{_prefix} \
	--libdir=%{_prefix}/%_lib \
	--mandir=%{_mandir} \
	--disable-mpi \
	--disable-static
fi
%endif

%{__make}


%install
rm -rf $RPM_BUILD_ROOT
mkdir $RPM_BUILD_ROOT
make install DESTDIR=$RPM_BUILD_ROOT

# Comply shared library policy
find $RPM_BUILD_ROOT -name *.la -exec rm -f {} \;

# Comply static build policy
find $RPM_BUILD_ROOT -name *.a -exec rm -f {} \;


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
%doc NEWS README TODO COPYING
%{_bindir}/*
%{_includedir}/*.h
%dir %{_includedir}/gerris
%{_includedir}/gerris/*.h
%{_libdir}/*.so.*
%{_libdir}/*.so
%dir %{_libdir}/gerris
%{_libdir}/gerris/*
%{_libdir}/pkgconfig/*.pc
%dir %{_datadir}/gerris
%{_datadir}/gerris/gfs.lang
%{_datadir}/gerris/gerris.dic
%{_datadir}/mime/packages/*.xml
%{_datadir}/icons/hicolor/48x48/mimetypes/*.png
%{_mandir}/man1/*.gz


%changelog
* Thu Jan 29 2009 Ivan Adam Vari <i.vari@niwa.co.nz> - 9
- Version change (1.3.1 -> 1.3.2) related minor fixes
- Added fortran dependency

* Wed Jul 16 2008 Ivan Adam Vari <i.vari@niwa.co.nz> - 8
- Version change (1.2.0 -> 1.3.0) related minor fixes
- Added macro for gts version specification
  Some other changes found in debian packages

* Thu Jul 3 2008 Ivan Adam Vari <i.vari@niwa.co.nz> - 7
- Fixed typo in %files section (attr)
- Added new file gfs.lang to %files section
- Disabled MPI according to debian build rules

* Thu May 15 2008 Ivan Adam Vari <i.vari@niwa.co.nz> - 6
- Added fedora 8 support for x86 (32bit only)
- Removed libtool config files to comply with shared
  library policy
- Removed static build bits to comply with shared
  library policy
- Fixed dependencies

* Mon May 12 2008 Ivan Adam Vari <i.vari@niwa.co.nz> - 5
- Added new package dependencies, minor fixes

* Mon Jan 7 2008 Ivan Adam Vari <i.vari@niwa.co.nz> - 4
- Removed %{?_smp_mflags} from make due to intermittent
  build errors on some SMP systems

* Mon Nov 12 2007 Ivan Adam Vari <i.vari@niwa.co.nz> - 3
- Fixed package (install) dependencies

* Mon Oct 1 2007 Ivan Adam Vari <i.vari@niwa.co.nz> - 2
- Removed unnecessary version specifications for some
  build requirements
- Added SLEx/SuSE compatibilty
- Added 64bit compatibility
- Updated %post, %postun scriptlets

* Tue May 1 2007 Ivan Adam Vari <i.vari@niwa.co.nz> - 1
- Initial rpm release based on Fedora/Redhat Linux
