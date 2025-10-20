FROM python:3.9-slim-bookworm AS iproverbuilder

RUN apt-get update \
    && apt-get install -y -q --no-install-recommends git ca-certificates python3 \
        opam bubblewrap libgmp-dev libz-dev \
    && update-ca-certificates
 
#  RUN apt-get install -y -q --no-install-recommends opam bubblewrap libgmp-dev libz-dev
 # Can disable sandboxing as we are running within a container here
# RUN opam --version

WORKDIR /usr/src
RUN git clone -b v3.9.2 --single-branch https://gitlab.com/korovin/iprover.git
WORKDIR /usr/src/iprover

RUN opam init -y --disable-sandboxing && opam update -y \
  && opam switch create 4.14.1+flambda --package=ocaml-variants.4.14.1+options,ocaml-option-flambda \
  && opam switch 4.14.1+flambda \
  && opam install -y --confirm-level=unsafe-yes ocamlfind ocamlgraph zarith yojson z3

ARG IPROVER_MAKE_TARGET="make"

RUN eval $(opam env) \
  && ./configure \
  && $IPROVER_MAKE_TARGET


ENV IPROVER_HOME="/usr/src/iprover"

WORKDIR /usr/src
RUN git clone -b pathways https://github.com/AlecGower/LGEMPlus.git && \
  mv LGEMPlus/src lgemplus && \
  rm -rf LGEMPlus

COPY py-requirements.txt /usr/src/py-requirements.txt
RUN pip install --no-cache-dir -r /usr/src/py-requirements.txt && \
    rm /usr/src/py-requirements.txt

LABEL maintainer.name="Alexander Gower"
LABEL maintainer.email="gower@chalmers.se"

CMD ["/bin/bash"]
