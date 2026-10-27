# autoconf file inspiration from 
# https://github.com/USCbiostats/software-dev/wiki/Setting-up-optional-OpenMP-support
# EDIT: who are we kidding we didn't take inspiration from the above, we just shamelessly stole and minimally
#       modified the RcppArmadillo configure.ac file that others have perfected.
#       https://github.com/RcppCore/RcppArmadillo/blob/master/configure.ac
#

AC_PREREQ([2.69])
AC_INIT(fido, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

AC_ARG_WITH([mkl],
    [AS_HELP_STRING([--with-mkl=auto|no|PATH],
        [detect MKL (default), disable it, or require the MKL installation at PATH])],
    [], [with_mkl=auto])
case "${with_mkl}" in
    yes|'') AC_MSG_ERROR([--with-mkl requires auto, no, or an installation path]) ;;
esac

## Set R_HOME, respecting an environment variable if one is set
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi

## Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

## We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP
AC_PROG_CXX

## Is R already configured to compile things using OpenMP without
## any extra hand-holding?
openmp_already_works="no"

## default to not even thinking about OpenMP as Armadillo wants a pragma
## variant available if and only if C++11 is used with g++ 5.4 or newer
can_use_openmp="no"

## Ensure TMPDIR is set.
AC_MSG_CHECKING([whether we have a suitable tempdir])
TMPDIR=$("${R_HOME}/bin/R" --vanilla --slave -e "cat(dirname(tempdir()))")
AC_MSG_RESULT([${TMPDIR}])

## Check if R is configured to compile programs using OpenMP out-of-the-box.
AC_MSG_CHECKING([whether R CMD SHLIB can already compile programs using OpenMP])

## Create private directory in TMPDIR.
BUILDDIR="${TMPDIR}/rcpparmadillo-$$-$RANDOM"
mkdir -p "${BUILDDIR}"

owd=$(pwd)
cd "${BUILDDIR}"

# MPN: These updates were from a (very, very useful) help question posted to r-package-devel
# Credit goes to Ivan Krylov.

cat > test-omp.cpp <<EOF
#include <omp.h>
extern "C" void configtest(int * arg) {
  *arg = omp_get_num_threads();
}
EOF
# Without the following you're relying on the GNU/Linux-like behaviour
# w.r.t. undefined symbols (see WRE 1.2.1.1):
cat > Makevars <<EOF
PKG_CXXFLAGS = \$(SHLIB_OPENMP_CXXFLAGS)
PKG_LIBS = \$(SHLIB_OPENMP_CXXFLAGS)
EOF

# compile and link the shared object
"${R_HOME}/bin/R" CMD SHLIB test-omp.cpp
# load it and use a function from the inside
"${R_HOME}/bin/Rscript" -e '
 dyn.load(paste0("test-omp", .Platform$dynlib.ext))
 .C("configtest", arg = integer(1))$arg
'
if test $? -eq 0; then
 # everything seems to have succeeded
 openmp_already_works="yes"
fi

## Go back home.
cd "${owd}"
rm -rf "${BUILDDIR}"

## additional reporting/diagnostics for mac don't trust the above approach for mac
RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')
AC_MSG_CHECKING([for mac OS])
if test x"${RSysinfoName}" = x"Darwin"; then
   AC_MSG_RESULT([found, not trusting "already works"])
   openmp_already_works="no"
else
   AC_MSG_RESULT([not found])
fi # end mac check

## If the above checks failed, then perform other heuristics
## based on the compiler version, etc.
if test x"${openmp_already_works}" = x"no"; then

    ## Check the C++ compiler using the CXX value set

    ## If it is g++, we have GXX set so let's examine it
    if test "${GXX}" = yes; then
        AC_MSG_CHECKING([whether g++ version is sufficient])
        gxx_version=$(${CXX} -v 2>&1 | awk '/^.*g.. version/ {print $3}')
        case ${gxx_version} in
            1.*|2.*|3.*|4.0.*|4.1.*|4.2.*|4.3.*|4.4.*|4.5.*|4.6.*|4.7.0|4.7.1)
                 AC_MSG_RESULT([no])
                 AC_MSG_WARN([Only g++ version 4.7.2 or greater can be used with RcppArmadillo.])
                 AC_MSG_ERROR([Please use a different compiler.])
            ;;
            4.7.*|4.8.*|4.9.*|5.0*|5.1*|5.2*|5.3*)
                 AC_MSG_RESULT([yes, but without OpenMP as version ${gxx_version} (Armadillo constraint)])
                 ## we know this one is bad
                 can_use_openmp="no"
            ;;
            5.4*|5.5*|5.6*|5.7*|5.8*|5.9*|6.*|7.*|8.*|9.*|10.*|11.*|12.*)
                 AC_MSG_RESULT([yes, with OpenMP as version ${gxx_version}])
                 ## we know this one is good, yay
                 can_use_openmp="yes"
            ;;
            *)
                 AC_MSG_RESULT([almost])
                 AC_MSG_WARN([Compiler self-identifies as being compliant with GNUC extensions but is not g++.])
                 ## we know nothing, so no
                 can_use_openmp="no"
            ;;
        esac
    fi

    ## Check for Apple LLVM

    AC_MSG_CHECKING([for macOS])
    RSysinfoName=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(Sys.info()[["sysname"]])')

    if test x"${RSysinfoName}" = x"Darwin"; then
        AC_MSG_RESULT([found])
        AC_MSG_CHECKING([for macOS Apple compiler])

        apple_compiler=$($CXX --version 2>&1 | grep -i -c -e 'apple llvm')

        if test x"${apple_compiler}" = x"1"; then
            AC_MSG_RESULT([found])
            AC_MSG_WARN([OpenMP unavailable and turned off.])
            can_use_openmp="no"
        else
            AC_MSG_RESULT([not found])
            AC_MSG_CHECKING([for clang compiler])
            clang_compiler=$($CXX --version 2>&1 | grep -i -c -e 'clang ')

            if test x"${clang_compiler}" = x"1"; then
                AC_MSG_RESULT([found])
                AC_MSG_CHECKING([for OpenMP compatible version of clang])
                clang_version=$(${CXX} -v 2>&1 | awk '/^.*clang version/ {print $3}')

                case ${clang_version} in
                    4.*|5.*|6.*|7.*|8.*|9.*|10.*|11.*)
                        AC_MSG_RESULT([found and suitable])
                        can_use_openmp="yes"
                    ;;
                    *)
                        AC_MSG_RESULT([not found])
                        AC_MSG_WARN([OpenMP unavailable and turned off.])
                        can_use_openmp="no"
                    ;;
                esac
            else
                AC_MSG_RESULT([not found])
                AC_MSG_WARN([unsupported macOS build detected; if anything breaks, you keep the pieces.])
            fi
        fi
    else
        AC_MSG_RESULT([no])
    fi

fi # if test x"${openmp_already_works}" = x"no"

## Check for suitable LAPACK_LIBS
AC_MSG_CHECKING([LAPACK_LIBS])

## external LAPACK has the required function
lapack=$(${R_HOME}/bin/R CMD config LAPACK_LIBS)
hasRlapack=$(echo ${lapack} | grep lRlapack)

## in what follows below we substitute both side of the define/undef
## while this may seem a little unusual we do it to fully reproduce the
## previous bash-based implementation

if test x"${hasRlapack}" = x""; then
    ## We are using a full Lapack and can use zgbsv -- so #undef remains
    AC_MSG_RESULT([system LAPACK found])
    arma_lapack="#undef ARMA_CRIPPLED_LAPACK"
else
    ## We are using R's subset of Lapack and CANNOT use zgbsv etc, so we mark it
    AC_MSG_RESULT([R-supplied partial LAPACK found])
    AC_MSG_WARN([Some complex-valued LAPACK functions may not be available])
    arma_lapack="#define ARMA_CRIPPLED_LAPACK 1"
fi

## Default the OpenMP flag to the empty string.
## If and only if OpenMP is found, expand to $(SHLIB_OPENMP_CXXFLAGS)
openmp_flag=""

## Set the fallback, by default it is nope
arma_have_openmp="#define ARMA_DONT_USE_OPENMP 1"

if test x"${openmp_already_works}" = x"yes"; then
    arma_have_openmp="#define ARMA_USE_OPENMP"
    ## Added by JDS 
    openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
fi

if test x"${can_use_openmp}" = x"yes"; then
    AC_MSG_CHECKING([for OpenMP])
    ## if R has -fopenmp we should be good
    allldflags=$(${R_HOME}/bin/R CMD config --ldflags)
    hasOpenMP=$(echo ${allldflags} | grep -- -fopenmp)
    if test x"${hasOpenMP}" = x""; then
        AC_MSG_RESULT([missing])
        arma_have_openmp="#define ARMA_DONT_USE_OPENMP"
    else
        AC_MSG_RESULT([found and suitable])
        arma_have_openmp="#define ARMA_USE_OPENMP"
        openmp_flag='$(SHLIB_OPENMP_CXXFLAGS)'
    fi
fi


mkl_cppflags=""
mkl_libs=""
mkl_reason="not found via MKLROOT or pkg-config"
mkl_source=""
mkl_root=""
mkl_arch=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(R.version$arch)')

if test "${with_mkl}" = no; then
    mkl_reason="disabled by --with-mkl=no"
elif test "${RSysinfoName}" = Darwin && { test "${mkl_arch}" = aarch64 || test "${mkl_arch}" = arm64; }; then
    mkl_reason="unavailable on Apple Silicon"
elif test "${RSysinfoName}" != Linux || test "${mkl_arch}" != x86_64; then
    mkl_reason="MKL is currently tested only on Linux x86_64"
else
    if test "${with_mkl}" != auto; then
        mkl_root="${with_mkl}"
        mkl_source="--with-mkl"
    elif test -n "${MKLROOT}"; then
        mkl_root="${MKLROOT}"
        mkl_source="MKLROOT"
    elif command -v pkg-config >/dev/null 2>&1 &&
         pkg-config --exists mkl-dynamic-lp64-gomp; then
        mkl_root=$(pkg-config --variable=prefix mkl-dynamic-lp64-gomp)
        mkl_source="pkg-config"
    fi

    if test -n "${mkl_root}"; then
        mkl_cr=$(printf '\r')
        case "${mkl_root}" in
            *'$'*|*'#'*|*'`'*|*'\'*|*'"'*|*"'"*|*','*|*"${as_nl}"*|*"${mkl_cr}"*)
                mkl_reason="unsupported characters in MKL path"
                mkl_root=""
                ;;
        esac
        if test -n "${mkl_root}"; then
            if test -d "${mkl_root}"; then
                mkl_root=$(cd "${mkl_root}" && pwd -P)
            else
                mkl_reason="${mkl_source} directory does not exist"
                mkl_root=""
            fi
        fi
        case "${mkl_root}" in
            *'$'*|*'#'*|*'`'*|*'\'*|*'"'*|*"'"*|*','*|*"${as_nl}"*|*"${mkl_cr}"*)
                mkl_reason="unsupported characters in resolved MKL path"
                mkl_root=""
                ;;
        esac
        if test -n "${mkl_root}"; then
            mkl_candidate_libs=""
            mkl_libdir=""
            for mkl_directory in "${mkl_root}/lib" "${mkl_root}/lib/intel64" "${mkl_root}/lib/intel64_lin"; do
                mkl_candidate_libs=""
                for mkl_library in mkl_intel_lp64 mkl_gnu_thread mkl_core; do
                    mkl_library_path=""
                    for mkl_suffix in .so .so.2; do
                        if test -f "${mkl_directory}/lib${mkl_library}${mkl_suffix}"; then
                            mkl_library_path="${mkl_directory}/lib${mkl_library}${mkl_suffix}"
                            break
                        fi
                    done
                    if test -z "${mkl_library_path}"; then
                        mkl_candidate_libs=""
                        break
                    fi
                    mkl_candidate_libs="${mkl_candidate_libs} \"${mkl_library_path}\""
                done
                if test -n "${mkl_candidate_libs}"; then
                    mkl_libdir="${mkl_directory}"
                    break
                fi
            done
            if test ! -f "${mkl_root}/include/mkl.h"; then
                mkl_reason="${mkl_source} is missing include/mkl.h"
            elif test -z "${mkl_libdir}"; then
                mkl_reason="${mkl_source} is missing LP64 GNU shared libraries"
            elif test -z "${openmp_flag}"; then
                mkl_reason="OpenMP is unavailable"
            else
                mkl_candidate_cppflags="-I\"${mkl_root}/include\""
                mkl_candidate_libs="-Wl,--no-as-needed${mkl_candidate_libs} -Wl,-rpath,\"${mkl_libdir}\" -lgomp -lpthread -lm -ldl"
                mkl_eigen_include=$("${R_HOME}/bin/Rscript" --vanilla -e 'cat(system.file("include", package="RcppEigen"))')
                mkl_probe_dir=$(mktemp -d "${TMPDIR}/fido-mkl-XXXXXX")
                if test -z "${mkl_probe_dir}"; then
                    mkl_reason="could not create MKL probe directory"
                else
                    cat > "${mkl_probe_dir}/fido-mkl-probe.cpp" <<'EOF'
#include <mkl.h>
#include <Eigen/Dense>
#include <omp.h>
#if !defined(_OPENMP) || !defined(__GNUC__) || defined(__clang__)
#error MKL GNU threading requires GNU C++ with OpenMP
#endif
static_assert(sizeof(MKL_INT) == 4, "MKL requires LP64");
extern "C" void fido_mkl_config_probe(int* status) {
    const int previous_limit = mkl_set_num_threads_local(1);
    const int local_limit = mkl_get_max_threads();
    mkl_set_num_threads_local(previous_limit);
    if (local_limit != 1) {
        *status = -1;
        return;
    }
    Eigen::Matrix2d precision = Eigen::Matrix2d::Zero();
    precision(0, 0) = 4.0;
    precision(1, 1) = 9.0;
    *status = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', 2, precision.data(), 2);
    if (*status != 0) return;
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(Eigen::Matrix2d::Identity());
    if (precision(0, 0) != 2.0 || precision(1, 1) != 3.0 ||
        omp_get_max_threads() < 1 || solver.info() != Eigen::Success ||
        (solver.eigenvalues().array() - 1.0).abs().maxCoeff() > 1e-12) {
        *status = -1;
    }
}
EOF
                    cat > "${mkl_probe_dir}/Makevars" <<EOF
PKG_CPPFLAGS = ${mkl_candidate_cppflags} -I"${mkl_eigen_include}" -DFIDO_USE_MKL -DEIGEN_USE_MKL_ALL
PKG_CXXFLAGS = ${openmp_flag}
PKG_LIBS = ${openmp_flag} ${mkl_candidate_libs}
EOF
                    if (cd "${mkl_probe_dir}" &&
                        "${R_HOME}/bin/R" CMD SHLIB fido-mkl-probe.cpp > probe.log 2>&1 &&
                        "${R_HOME}/bin/Rscript" --vanilla -e '
                            probe <- dyn.load(paste0("fido-mkl-probe", .Platform$dynlib.ext))
                            result <- .C(getNativeSymbolInfo("fido_mkl_config_probe", probe), status=integer(1))
                            stopifnot(identical(result$status, 0L))
                        ' >> probe.log 2>&1); then
                        mkl_cppflags="${mkl_candidate_cppflags}"
                        mkl_libs="${mkl_candidate_libs}"
                    else
                        mkl_reason="compile, link, or R load/run probe failed; see config.log"
                    fi
                    cat "${mkl_probe_dir}/probe.log" >&AS_MESSAGE_LOG_FD
                    rm -rf "${mkl_probe_dir}"
                fi
            fi
        fi
    fi
fi

if test -n "${mkl_libs}"; then
    AC_MSG_NOTICE([MKL enabled: LP64, GNU threading (${mkl_source})])
elif test "${with_mkl}" != auto && test "${with_mkl}" != no; then
    AC_MSG_ERROR([MKL requested but unavailable: ${mkl_reason}])
else
    AC_MSG_NOTICE([MKL disabled: ${mkl_reason}])
fi

AC_SUBST([MKL_CPPFLAGS], ["${mkl_cppflags}"])
AC_SUBST([MKL_LIBS], ["${mkl_libs}"])

## now use all these
AC_SUBST([ARMA_LAPACK],["${arma_lapack}"])
AC_SUBST([ARMA_HAVE_OPENMP], ["${arma_have_openmp}"])
AC_SUBST([OPENMP_FLAG], ["${openmp_flag}"])
AC_CONFIG_FILES([inst/include/fidoGenerated.h src/Makevars])
AC_OUTPUT
