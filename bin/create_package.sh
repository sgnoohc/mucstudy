tar --exclude='*/.git/*' \
    -cvzf ${MY_MUCOLL_BASEDIR}/condor/package.tar.gz  \
    bin \
    condor \
    mucoll_software \
    steering \
    study \
    setup.sh
