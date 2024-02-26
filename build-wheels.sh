#!/bin/bash
set -e -u -x

function repair_wheel {
    wheel="$1"
    if ! auditwheel show "$wheel"; then
        echo "Skipping non-platform wheel $wheel"
    else
        auditwheel repair "$wheel" -w /io/wheelhouse/
    fi
}

dnf install libarchive-devel -y
cd /io

# Compile wheels
for PYBIN in /opt/python/cp*/bin; do
    "${PYBIN}/python3" setup.py bdist_wheel
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    repair_wheel "$whl"
done
