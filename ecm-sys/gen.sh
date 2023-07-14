#!/bin/bash

bindgen wrapper.h -o src/bindings.rs --whitelist-function '^ecm_.*' --whitelist-var '^ECM_.*' -- -DUSE_ZLIB
