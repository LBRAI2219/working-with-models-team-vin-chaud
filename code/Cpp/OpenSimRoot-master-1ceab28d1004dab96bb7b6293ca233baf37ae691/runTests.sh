#!/bin/bash
cd OpenSimRoot/tests/engine
./testEngine.sh "../../StaticBuild/OpenSimRoot"
rexe=$?
echo finished testing engine with error status $rexe
cd ../modules
./testModules.sh "../../StaticBuild/OpenSimRoot"
rexm=$?
echo finished testing engine with error status $rexm
cd ../../..
rex=$(($rexe+$rexm))
echo exiting with error status $rex
exit $rex


