## Usage integrate_GLPK.sh [GLPK.tar.gz]
## Integrates the latest GLPK package
## theussl, 2008-07-08

## latest glpk tarball (Austrian mirror out of date)
#URL="http://gd.tuwien.ac.at/gnu/gnusrc/glpk/"
URL="ftp://ftp.gnu.org/gnu/glpk/"
latest="glpk-4.34.tar.gz"
## where to put source files and headers
DESTINATION=../src

## --------------------------------------------
## Usage

usage() {
cat << "USAGE_END"
Usage: integrate_GLPK.sh [-g ]
       integrate_GLPK.sh [-i GLPK_source.tar.gz]
       integrate_GLPK.sh [-c]
Get, integrate or clean GLPK sources

Options:
  -g, --get           get latest Version of GLPK
  -i, --integrate     integrate given GLPK sources
  -c, --clean         clean the R package's src directory

USAGE_END
        exit $1
}

## --------------------------------------------
## Read command line arguments

for x in "$@" ; do
    case "${x}" in
        -i|--integrate)        # integrate sources
             integrate=true
             ;;
        -c|--clean)            # clean sources
             clean=true
             ;;
        -g|--get)            # clean sources
             get=true
	     sources=$latest # sources gets latest
             ;;
        -*)  # invalid option
             echo "$0: Invalid switch \"${x}\"!  Run $0 without parameters for help."
             exit 1
             ;;
        *)   # this should be the tarball of glpk sources
             if [[ ! -z "${sources}" ]] ; then
                 echo "$0: Only one source file allowed: \"${sources}\"!  Run $0 without parameters for help."
                 exit 1
             fi
             sources="${x}"
             ;;
    esac
done

## --------------------------------------------
## input validation

if [[ ! ( ${integrate} || ${clean} || ${get}) ]] ; then
    echo "$0: No option given; nothing to do!"
    usage 1
    exit 1
fi

if [[ ( ${integrate} && ${clean} ) || ( ${get} && ${clean} )]] ; then
    echo "$0: --clean can only be used alone!  Run $0 without parameters for help."
    exit 1
fi

if [[ -z "${sources}" && $integrate ]] ; then
    echo "$0: No source file to integrate given!"
    usage 1
    exit 1
fi


## --------------------------------------------
## integrate GLPK sources to package

if [[ $get ]] ; then
    if [[ ! -s "${sources}" ]] ; then
	wget $URL/$sources
    else
	echo "$sources already available."
    fi
fi

if [[ $integrate ]] ; then
    
    if [[ ! -s "${sources}" ]] ; then
	echo "$0: Selected source file \"$sources\" is not available or zero!"
	usage 1
	exit 1
    fi
    GLPK=`basename $sources .tar.gz`
    SOURCEDIR=${GLPK}

    tar xzf $sources

    cp $SOURCEDIR/include/glp* $DESTINATION
    cp $SOURCEDIR/src/glp*     $DESTINATION
    if [[ -d $SOURCEDIR ]] ; then
	rm -rf $SOURCEDIR
    fi
    echo "NOTE: watch out for abort() statements in upstream code!"
    echo "NOTE: replacement = error(\"Execution aborted.\");"
fi


if [[ $clean ]] ; then
    rm $DESTINATION/glp*
fi

echo "done."

exit 0
