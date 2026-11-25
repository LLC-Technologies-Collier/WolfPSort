# Copyright 2026 Google LLC
# Licensed under the GNU General Public License v3 (GPLv3).
cmake_minimum_required(VERSION 3.16)
project(wolfPredict C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_executable(wolfPredict wolfPredict.c)
if(NOT WIN32)
    target_link_libraries(wolfPredict PRIVATE m)
endif()
