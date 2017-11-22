FROM ubuntu:16.04

LABEL description="Docker image for Lacer/LoFreq workflow" \
      maintainer="Andreas Wilm <wilma@gis.a-star.edu.sg>"

ENV DEBIAN_FRONTEND noninteractive

# sources from biocontainers
RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list

ENV LACER_VERSION ea0e5e6
ENV BUILD_PACKAGES "build-essential zlib1g-dev wget ca-certificates git"

RUN apt-get update && \
    apt-get install -y --no-install-recommends $BUILD_PACKAGES && \
    cd /tmp && \
    echo "INFO Downloading and unpacking htslib" && \
    wget --quiet https://github.com/samtools/htslib/archive/1.1.tar.gz && \
    tar xzf 1.1.tar.gz && \
    cd htslib-1.1/ && \
    echo "INFO Compiling htslib" && \
    make >/dev/null && \
    cd .. && \
    rm 1.1.tar.gz && \
    echo "INFO Downloading and unpacking samtools" && \
    wget --quiet https://github.com/samtools/samtools/archive/1.1.tar.gz && \
    tar xzf 1.1.tar.gz && \
    cd samtools-1.1/ && \
    echo "INFO Compiling samtools" && \
    sed -i -e 's,-D_CURSES_LIB=1,-D_CURSES_LIB=0,' -e 's,$(LIBCURSES),,' Makefile && \
    make HTSDIR=../htslib-1.1/ >/dev/null && \
    cd .. && \
    rm 1.1.tar.gz && \
    echo "INFO Installing Lacer" && \
    git clone https://github.com/swainechen/lacer.git && \
    cd lacer/ && \
    git checkout $LACER_VERSION && \
    cd lacepr/ && \
    gcc -I../../samtools-1.1/ -I../../htslib-1.1/ lacepr.c -L../../samtools-1.1/ -L../../htslib-1.1/ \
        -lbam -l:libhts.a -lz -lpthread -lm -o lacepr && \
    cp lacepr ../lacer.pl /usr/local/bin && \
    echo "INFO Installing Lacer dependencies" && \
    apt-get install -y --no-install-recommends \
            libterm-progressbar-perl pdl libmemory-usage-perl \
            libgetopt-simple-perl libdata-dumper-simple-perl \
            libbio-samtools-perl && \
    echo yes | cpan PDL::Parallel::threads >/dev/null && \
    rm -rf $HOME/.cpan/ && \
    apt-get remove -y $BUILD_PACKAGES && \
    apt-get autoremove -y && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


CMD ["/bin/bash"]
