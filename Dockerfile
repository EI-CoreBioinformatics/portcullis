FROM maplesond/cppbuild:latest
LABEL maintainer="d.mapleson@gmail.com"

ARG VERSION
ENV VERSION ${VERSION:-1.2.0}
RUN echo ${VERSION}

COPY . /portcullis-src
WORKDIR /portcullis-src

RUN ./update_version.sh $VERSION
RUN ./autogen.sh && ./configure --prefix=/portcullis
RUN make clean && make -j4
RUN make -j4 check
RUN make install
RUN rm -rf /portcullis-src

ENV PATH=$PATH:/portcullis/bin
ENV PYTHONPATH=$PYTHONPATH:/portcullis/lib/python3.6/site-packages
WORKDIR /data
