# -*- coding: utf-8 -*-
#
# Copyright 2012-2019 CNRS-UM LIRMM, CNRS-AIST JRL
#

from conans import python_requires
import os
import shutil
import subprocess
import sys

base = python_requires("Eigen3ToPython/1.0.0@gergondet/stable")


class SpaceVecAlgConan(base.Eigen3ToPythonConan):
    name = "SpaceVecAlg"
    version = "1.1.0"
    description = "Spatial Vector Algebra with the Eigen library"
    topics = ("robotics", "algebra", "eigen", "python")
    url = "https://github.com/jrl-umi3218/SpaceVecAlg"
    homepage = "https://github.com/jrl-umi3218/SpaceVecAlg"
    author = "Pierre Gergondet <pierre.gergondet@gmail.com>"
    license = "BSD-2-Clause"
    exports = ["LICENSE"]
    exports_sources = ["CMakeLists.txt", "conan/CMakeLists.txt", "binding/*", "cmake/*", "doc/*", "src/*", "tests/*"]
    generators = "cmake"
    settings = "os", "arch"
    options = { "python_version": ["2.7", "3.3", "3.4", "3.5", "3.6", "3.7"] }
    default_options = { "python_version": base.get_python_version() }

    requires = (
        "Eigen3ToPython/1.0.0@gergondet/stable"
    )

    def package_info(self):
        self.env_info.PYTHONPATH.append(self._extra_python_path())
