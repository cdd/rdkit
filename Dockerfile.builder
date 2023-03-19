FROM public.ecr.aws/docker/library/ubuntu:22.04 as build
ARG BUILD_PG
ARG BUILD_PY
ARG PACKAGE_NAME
COPY . .
RUN ./build-deb ${BUILD_PG} ${BUILD_PY} ${PACKAGE_NAME}

FROM scratch as artifact
COPY --from=build build/*.deb build/
