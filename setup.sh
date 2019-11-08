set -xe
WD=$(pwd)
./copy_readmes_in_packages.sh
cd $WD/thermocepstrum
pip install .
cd $WD/thermocepstrum_gui
pip install .
cd $WD
