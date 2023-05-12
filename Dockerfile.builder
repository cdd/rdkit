FROM public.ecr.aws/docker/library/ubuntu:22.04 as build
ARG BUILD_PG
ARG BUILD_PY
ARG PACKAGE_NAME
ARG TEST

COPY install-deps install-deps
RUN ./install-deps ${BUILD_PG} ${BUILD_PY}

WORKDIR /rdkit
COPY . .
RUN --mount=type=cache,target=/rdkit/build ./build-deb ${BUILD_PG} ${BUILD_PY} ${PACKAGE_NAME} ${TEST}

FROM scratch as artifact
COPY --from=build /rdkit/*.deb .
