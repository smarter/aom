##
## Copyright (c) 2017, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
cmake_minimum_required(VERSION 3.5)

set(REQUIRED_ARGS "AOM_ROOT" "AOM_CONFIG_DIR" "AOM_TARGET_SYSTEM" "AOM_SYM_FILE"
    "CONFIG_AV1_DECODER" "CONFIG_AV1_ENCODER")

foreach (arg ${REQUIRED_ARGS})
  if ("${${arg}}" STREQUAL "")
    message(FATAL_ERROR "${arg} must not be empty.")
  endif ()
endforeach ()

include("${AOM_ROOT}/build/cmake/exports_sources.cmake")

if ("${AOM_TARGET_SYSTEM}" STREQUAL "Darwin")
  set(symbol_prefix "_")
elseif ("${AOM_TARGET_SYSTEM}" MATCHES "Windows\|MSYS" AND AOM_MSVC)
  set(symbol_prefix "_")
  file(WRITE "${AOM_SYM_FILE}"
       "LIBRARY libaom INITINSTANCE TERMINSTANCE\n"
       "DATA MULTIPLE NONSHARED\n"
       "EXPORTS\n")
else ()
  set(symbol_suffix ";")
endif ()

set(aom_sym_file "${AOM_SYM_FILE}")

if ("${AOM_TARGET_SYSTEM}" STREQUAL "Darwin")
  file(REMOVE "${aom_sym_file}")
elseif ("${AOM_TARGET_SYSTEM}" MATCHES "Windows\|MSYS")
  file(WRITE "${aom_sym_file}"
       "LIBRARY libaom INITINSTANCE TERMINSTANCE\n"
       "DATA MULTIPLE NONSHARED\n"
       "EXPORTS\n")
else ()
  file(WRITE "${aom_sym_file}" "{ global:\n")
endif ()

foreach (export_file ${AOM_EXPORTS_SOURCES})
  file(STRINGS "${export_file}" exported_file_data)
  set(exported_symbols "${exported_symbols} ${exported_file_data};")
  string(STRIP "${exported_symbols}" exported_symbols)
endforeach ()

foreach (exported_symbol ${exported_symbols})
  string(STRIP "${exported_symbol}" exported_symbol)
  string(REGEX REPLACE "text \|data " "" "exported_symbol" "${exported_symbol}")
  set(exported_symbol "${symbol_prefix}${exported_symbol}${symbol_suffix}")
  file(APPEND "${aom_sym_file}" "${exported_symbol}\n")
endforeach ()

if ("${aom_sym_file}" MATCHES "ver$")
  file(APPEND "${aom_sym_file}" " };")
endif ()