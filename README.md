[![Build Status](https://travis-ci.org/jorisv/SpaceVecAlg.svg?branch=master)](https://travis-ci.org/jorisv/SpaceVecAlg)

SpaceVecAlg aim to implement spatial vector algebra with the Eigen3 linear algebra library.

All this implementation is based on appendix A of Roy Featherstone Rigid Body Dynamics Algorithms book.

Pulling git subtree
------

To update sync cmake or .travis directory with their upstream git repository:

	git fetch git://github.com/jrl-umi3218/jrl-cmakemodules.git master
	git subtree pull --prefix cmake git://github.com/jrl-umi3218/jrl-cmakemodules.git master --squash

	git fetch git://github.com/jrl-umi3218/jrl-travis.git master
	git subtree pull --prefix .travis git://github.com/jrl-umi3218/jrl-travis.git master --squash
